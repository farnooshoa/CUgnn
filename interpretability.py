"""Interpretability analysis — which copper genes and edges drive stage prediction.

Retrains one GAT on all 349 stage-labeled tumor samples (no CV split)
then extracts:
  1. Signed saliency  — which genes push the model toward late stage
  2. Attention edges  — which gene-gene connections the model focuses on
  3. Comparison table — saliency rank vs attention rank vs log2FC vs biology

Outputs:
  outputs/final_comparison/interp_saliency.csv
  outputs/final_comparison/interp_attention.csv
  outputs/final_comparison/interp_comparison_table.csv
  outputs/final_comparison/interp_report.md
"""
from __future__ import annotations
import sys
import random
from pathlib import Path
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from torch_geometric.loader import DataLoader

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from src.preprocessing import load_lihc_dataset
from src.graph_building import build_functional_graph
from src.gnn_models.dataset import build_graph_dataset
from src.gnn_models.models import GATGraphClassifier
from src.gnn_models.train import TrainConfig, _class_weights

OUT = ROOT / "outputs" / "final_comparison"
OUT.mkdir(parents=True, exist_ok=True)

SEED = 42


def set_seeds(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)


def classify_stage(s):
    if not isinstance(s, str): return None
    s = s.upper().replace("STAGE ", "").strip()
    if s in {"I", "II", "IIA", "IIB"}: return 0
    if s.startswith("III") or s.startswith("IV"): return 1
    return None


# ── known biology for comparison ───────────────────────────────────────────
# What we expect from copper biology + literature
# positive = higher expression → late stage
# negative = higher expression → early stage
KNOWN_BIOLOGY = {
    "LOX":    ("ECM crosslinking, invasion, metastasis", "late"),
    "LOXL2":  ("ECM remodeling, TGF-beta pathway, invasion", "late"),
    "LOXL1":  ("ECM crosslinking", "late"),
    "LOXL3":  ("ECM remodeling", "late"),
    "LOXL4":  ("ECM remodeling", "late"),
    "SPARC":  ("ECM remodeling, tumor progression", "late"),
    "ATP7B":  ("Cu efflux, cisplatin resistance, HCC progression", "late"),
    "ATP7A":  ("Cu efflux pump", "late"),
    "ATOX1":  ("Cu chaperone, nuclear signaling, proliferation", "late"),
    "CP":     ("Ferroxidase, serum marker elevated in HCC", "late"),
    "SLC31A1":("Cu importer CTR1, elevated in proliferating cells", "late"),
    "SOD1":   ("Antioxidant, protective in early stage", "early"),
    "PRNP":   ("Prion protein, Cu binding, tumor suppressor role", "early"),
    "MT3":    ("Metallothionein, protective", "early"),
    "DBH":    ("Dopamine beta-hydroxylase, low in HCC", "early"),
}


def main():
    print("[interp] loading dataset ...")
    set_seeds(SEED)

    ds    = load_lihc_dataset(require_real=True)
    genes = ds.expression.index.tolist()
    graph = build_functional_graph(genes, ds.copper_genes)
    md    = ds.metadata.loc[ds.expression.columns].copy()

    # stage subset
    stage_bin = md["stage"].map(classify_stage)
    keep      = (md["sample_type"] == "Tumor") & stage_bin.notna()
    md_stage  = md[keep].copy()
    y         = stage_bin[keep].astype(int).to_numpy()

    print(f"  n={len(y)}  early={int((y==0).sum())}  late={int((y==1).sum())}")

    # build dataset — z-score on all samples (full train)
    bundle = build_graph_dataset(ds, graph, zscore_train_mask=None)

    # overwrite labels with stage
    sample_ids = ds.expression.columns.tolist()
    data_list  = []
    for sid, label in zip(md_stage.index, y):
        idx = sample_ids.index(sid)
        d   = bundle.data_list[idx]
        d.y = torch.tensor([int(label)], dtype=torch.long)
        data_list.append(d)

    # ── train full model ───────────────────────────────────────────────────
    print("\n[interp] training GAT on all stage patients ...")
    cfg     = TrainConfig(model="gat", epochs=120)
    class_w = _class_weights(y, cfg.device)
    model   = GATGraphClassifier(
        bundle.in_dim, cfg.hidden, n_classes=2,
        n_heads=4, n_layers=cfg.n_layers, dropout=cfg.dropout,
    ).to(cfg.device)
    opt     = torch.optim.Adam(model.parameters(),
                               lr=cfg.lr, weight_decay=cfg.weight_decay)
    loss_fn = nn.CrossEntropyLoss(weight=class_w)
    dl      = DataLoader(data_list, batch_size=cfg.batch_size, shuffle=True)

    for ep in range(cfg.epochs):
        model.train()
        total_loss = 0.0
        for batch in dl:
            batch = batch.to(cfg.device)
            opt.zero_grad()
            loss = loss_fn(
                model(batch.x, batch.edge_index, batch.batch,
                      edge_weight=getattr(batch, "edge_weight", None)),
                batch.y)
            loss.backward()
            opt.step()
            total_loss += loss.item()
        if (ep + 1) % 20 == 0:
            print(f"  epoch {ep+1:3d}/{cfg.epochs}  loss={total_loss/len(dl):.4f}")

    # ── 1. SIGNED SALIENCY ─────────────────────────────────────────────────
    print("\n[interp] extracting signed saliency ...")
    model.eval()
    single_dl = DataLoader(data_list, batch_size=1, shuffle=False)

    # signed: positive = gene pushes toward late stage
    # negative = gene pushes toward early stage
    gene_saliency = np.zeros(len(bundle.gene_order), dtype=np.float64)
    n_late_graphs = 0

    for batch in single_dl:
        batch   = batch.to(cfg.device)
        batch.x = batch.x.detach().requires_grad_(True)
        logits  = model(batch.x, batch.edge_index, batch.batch,
                        edge_weight=getattr(batch, "edge_weight", None))
        # gradient of late-stage logit w.r.t. node features
        late_logit = logits[:, 1].sum()
        late_logit.backward()
        grad = batch.x.grad.detach().cpu().numpy()
        # signed mean gradient across feature dimensions per node
        gene_saliency += grad.mean(axis=1)
        n_late_graphs += 1

    gene_saliency /= n_late_graphs

    sal_df = pd.DataFrame({
        "gene":      bundle.gene_order,
        "saliency":  gene_saliency,
        "direction": ["→ late" if s > 0 else "→ early"
                      for s in gene_saliency],
        "abs_saliency": np.abs(gene_saliency),
    }).sort_values("abs_saliency", ascending=False).reset_index(drop=True)
    sal_df["saliency_rank"] = sal_df.index + 1
    sal_df.to_csv(OUT / "interp_saliency.csv", index=False)

    print("\n  Top 15 genes by saliency:")
    print(f"  {'rank':<5} {'gene':<12} {'saliency':>10} {'direction':<12}")
    print("  " + "-"*42)
    for _, r in sal_df.head(15).iterrows():
        print(f"  {int(r['saliency_rank']):<5} {r['gene']:<12} "
              f"{r['saliency']:>10.4f} {r['direction']:<12}")

    # ── 2. ATTENTION EDGES ─────────────────────────────────────────────────
    print("\n[interp] extracting attention weights ...")
    model.eval()
    agg_attention: dict[tuple[int,int], list[float]] = {}

    with torch.no_grad():
        for batch in single_dl:
            batch = batch.to(cfg.device)
            _     = model(batch.x, batch.edge_index, batch.batch,
                          edge_weight=getattr(batch, "edge_weight", None),
                          return_attention=True)
            attn_list = model.get_last_attention()
            if not attn_list:
                continue
            ei, alpha = attn_list[-1]
            ei    = ei.cpu().numpy()
            alpha = (alpha.mean(dim=1).cpu().numpy()
                     if alpha.ndim > 1 else alpha.cpu().numpy())
            for (s, t), a in zip(ei.T, alpha):
                key = (int(s), int(t))
                agg_attention.setdefault(key, []).append(float(a))

    n_nodes = len(bundle.gene_order)
    attn_rows = []
    for (s, t), vals in agg_attention.items():
        if s < n_nodes and t < n_nodes:
            attn_rows.append({
                "source":        bundle.gene_order[s],
                "target":        bundle.gene_order[t],
                "attention_mean": float(np.mean(vals)),
                "attention_std":  float(np.std(vals)),
            })

    attn_df = (pd.DataFrame(attn_rows)
               .sort_values("attention_mean", ascending=False)
               .reset_index(drop=True))
    attn_df["attention_rank"] = attn_df.index + 1

    # deduplicate undirected pairs — keep highest attention direction
    seen = set()
    dedup_rows = []
    for _, r in attn_df.iterrows():
        pair = tuple(sorted([r["source"], r["target"]]))
        if pair not in seen:
            seen.add(pair)
            dedup_rows.append(r)
    attn_dedup = pd.DataFrame(dedup_rows).reset_index(drop=True)
    attn_dedup.to_csv(OUT / "interp_attention.csv", index=False)

    print("\n  Top 15 attention edges:")
    print(f"  {'source':<12} {'target':<12} {'attention':>10}")
    print("  " + "-"*36)
    for _, r in attn_dedup.head(15).iterrows():
        print(f"  {r['source']:<12} {r['target']:<12} "
              f"{r['attention_mean']:>10.4f}")

    # ── 3. COMPARISON TABLE ────────────────────────────────────────────────
    print("\n[interp] building comparison table ...")

    # differential expression — tumor vs normal from whole dataset
    tumor_mask  = ds.tumor_mask
    normal_mask = ~tumor_mask
    expr        = ds.expression.to_numpy()
    log2fc      = expr[:, tumor_mask].mean(axis=1) - \
                  expr[:, normal_mask].mean(axis=1)
    de_df = pd.DataFrame({
        "gene":    ds.expression.index,
        "log2FC":  log2fc,
    }).set_index("gene")

    # build comparison for top 20 by saliency
    top_genes = sal_df.head(20)["gene"].tolist()
    comp_rows = []
    sal_rank_map = sal_df.set_index("gene")["saliency_rank"].to_dict()
    sal_dir_map  = sal_df.set_index("gene")["direction"].to_dict()
    sal_val_map  = sal_df.set_index("gene")["saliency"].to_dict()

    for gene in top_genes:
        fc  = de_df.loc[gene, "log2FC"] if gene in de_df.index else float("nan")
        bio = KNOWN_BIOLOGY.get(gene, ("—", "unknown"))
        comp_rows.append({
            "gene":           gene,
            "saliency_rank":  int(sal_rank_map.get(gene, 99)),
            "saliency_value": round(sal_val_map.get(gene, 0), 4),
            "model_direction": sal_dir_map.get(gene, "—"),
            "log2FC_tumor_vs_normal": round(fc, 3),
            "FC_direction":   "up in tumor" if fc > 0 else "down in tumor",
            "biology":        bio[0],
            "expected_stage": bio[1],
        })

    comp_df = pd.DataFrame(comp_rows)
    comp_df.to_csv(OUT / "interp_comparison_table.csv", index=False)

    # count agreements
    agrees = 0
    total  = 0
    for _, r in comp_df.iterrows():
        if r["expected_stage"] == "unknown":
            continue
        model_late = r["model_direction"] == "→ late"
        bio_late   = r["expected_stage"] == "late"
        total += 1
        if model_late == bio_late:
            agrees += 1

    print(f"\n  Model-biology agreement: {agrees}/{total} genes "
          f"({100*agrees/max(total,1):.0f}%)")

    # ── write report ───────────────────────────────────────────────────────
    report_lines = [
        "# Interpretability Analysis — Stage Classification",
        "",
        "## Method",
        "",
        "- **Signed saliency**: gradient of the late-stage logit w.r.t. node",
        "  features, averaged across all 349 stage-labeled tumor samples.",
        "  Positive = gene pushes model toward predicting late stage.",
        "  Negative = gene pushes model toward predicting early stage.",
        "- **Attention edges**: mean attention weight from the final GAT layer,",
        "  averaged across all samples. High attention = model focuses on this",
        "  gene-gene connection when making predictions.",
        "",
        "## Top 20 genes by saliency — comparison with biology",
        "",
        "| gene | saliency rank | model direction | log2FC (T/N) | expected stage | biology |",
        "|---|---:|---|---:|---|---|",
    ]
    for _, r in comp_df.iterrows():
        agree = ""
        if r["expected_stage"] != "unknown":
            model_late = r["model_direction"] == "→ late"
            bio_late   = r["expected_stage"] == "late"
            agree = " ✓" if model_late == bio_late else " ✗"
        report_lines.append(
            f"| **{r['gene']}** | {r['saliency_rank']} | "
            f"{r['model_direction']}{agree} | {r['log2FC_tumor_vs_normal']:+.3f} | "
            f"{r['expected_stage']} | {r['biology'][:50]} |"
        )

    report_lines += [
        "",
        f"**Model-biology agreement: {agrees}/{total} genes "
        f"({100*agrees/max(total,1):.0f}%)**",
        "",
        "## Top 15 attention edges",
        "",
        "| source | target | attention mean | in curated edges? |",
        "|---|---|---:|---|",
    ]

    # check if edge is in curated list
    from src.graph_building.build_graph import CURATED_CU_EDGES
    curated_pairs = {tuple(sorted([s, t])) for s, t, _ in CURATED_CU_EDGES}

    for _, r in attn_dedup.head(15).iterrows():
        pair    = tuple(sorted([r["source"], r["target"]]))
        in_cur  = "✓ curated" if pair in curated_pairs else "self-loop / compartment"
        report_lines.append(
            f"| {r['source']} | {r['target']} | "
            f"{r['attention_mean']:.4f} | {in_cur} |"
        )

    report_lines += [
        "",
        "## Interpretation",
        "",
        "The saliency map shows which copper genes the GAT relies on most",
        "when predicting AJCC stage. Agreement with known biology validates",
        "that the model is learning real copper biology, not statistical noise.",
        "",
        "The attention edges show which gene-gene connections the model",
        "focuses on. If high-attention edges match curated copper biology",
        "(e.g. ATOX1→ATP7B, LOX→SPARC, SLC31A1→ATOX1), this strongly",
        "supports the biological relevance of the graph structure.",
        "",
        "## Files produced",
        "- `interp_saliency.csv` — all genes ranked by saliency",
        "- `interp_attention.csv` — all attention edges ranked",
        "- `interp_comparison_table.csv` — top 20 genes with biology comparison",
        "- `interp_report.md` — this document",
    ]

    (OUT / "interp_report.md").write_text(
        "\n".join(report_lines), encoding="utf-8")

    print(f"\n[interp] wrote interp_saliency.csv")
    print(f"[interp] wrote interp_attention.csv")
    print(f"[interp] wrote interp_comparison_table.csv")
    print(f"[interp] wrote interp_report.md")
    print("\nDone.")


if __name__ == "__main__":
    main()
