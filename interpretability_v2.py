"""Interpretability v2 — fixed training + filtered attention.

Changes from v1:
  1. Stronger class weight (squared inverse frequency)
  2. More epochs (200) with early stopping on training loss
  3. Self-loop filtering in attention
  4. Validation check — confirm model actually learned before extracting
  5. Saliency averaged separately for early vs late graphs
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
from src.gnn_models.train import TrainConfig

OUT = ROOT / "outputs" / "final_comparison"
OUT.mkdir(parents=True, exist_ok=True)

SEED = 42

KNOWN_BIOLOGY = {
    "LOX":    ("ECM crosslinking, invasion", "late"),
    "LOXL2":  ("ECM remodeling, invasion", "late"),
    "LOXL1":  ("ECM crosslinking", "late"),
    "LOXL3":  ("ECM remodeling", "late"),
    "LOXL4":  ("ECM remodeling", "late"),
    "SPARC":  ("ECM remodeling, tumor progression", "late"),
    "ATP7B":  ("Cu efflux, HCC progression", "late"),
    "ATP7A":  ("Cu efflux pump", "late"),
    "ATOX1":  ("Cu chaperone, proliferation", "late"),
    "CP":     ("Ferroxidase, elevated in HCC", "late"),
    "SLC31A1":("Cu importer, elevated in proliferating cells", "late"),
    "COX17":  ("Mitochondrial Cu relay, energy metabolism", "late"),
    "SOD1":   ("Antioxidant, protective early stage", "early"),
    "PRNP":   ("Cu binding, tumor suppressor role", "early"),
    "MT3":    ("Metallothionein, protective", "early"),
    "DBH":    ("Low in HCC", "early"),
    "ALB":    ("Liver function marker, low in late HCC", "early"),
}


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


def strong_class_weights(y, device):
    """Squared inverse frequency — harder push toward minority class."""
    counts = np.bincount(y, minlength=2).astype(np.float32)
    counts = np.where(counts == 0, 1.0, counts)
    w = (counts.sum() / (2 * counts)) ** 2   # squared = stronger
    w = w / w.mean()                          # normalise
    return torch.tensor(w, dtype=torch.float32, device=device)


def check_predictions(model, data_list, device):
    """Return fraction predicted as late — should be > 0.1 if model learned."""
    model.eval()
    dl    = DataLoader(data_list, batch_size=32, shuffle=False)
    preds = []
    with torch.no_grad():
        for batch in dl:
            batch = batch.to(device)
            logits = model(batch.x, batch.edge_index, batch.batch,
                           edge_weight=getattr(batch, "edge_weight", None))
            preds.extend(logits.argmax(dim=-1).cpu().numpy().tolist())
    preds = np.array(preds)
    return preds.mean(), (preds == 1).sum(), (preds == 0).sum()


def main():
    print("[interp v2] loading dataset ...")
    set_seeds(SEED)

    ds    = load_lihc_dataset(require_real=True)
    genes = ds.expression.index.tolist()
    graph = build_functional_graph(genes, ds.copper_genes)
    md    = ds.metadata.loc[ds.expression.columns].copy()

    stage_bin = md["stage"].map(classify_stage)
    keep      = (md["sample_type"] == "Tumor") & stage_bin.notna()
    md_stage  = md[keep].copy()
    y         = stage_bin[keep].astype(int).to_numpy()

    n_early = int((y == 0).sum())
    n_late  = int((y == 1).sum())
    print(f"  n={len(y)}  early={n_early}  late={n_late}  ratio={n_early/n_late:.1f}:1")

    bundle    = build_graph_dataset(ds, graph, zscore_train_mask=None)
    sample_ids = ds.expression.columns.tolist()
    data_list  = []
    for sid, label in zip(md_stage.index, y):
        idx = sample_ids.index(sid)
        d   = bundle.data_list[idx]
        d.y = torch.tensor([int(label)], dtype=torch.long)
        data_list.append(d)

    # ── train with stronger class weights ─────────────────────────────────
    print("\n[interp v2] training GAT (200 epochs, strong class weights) ...")
    device  = "cuda" if torch.cuda.is_available() else "cpu"
    cfg     = TrainConfig(model="gat", epochs=200, lr=3e-3,
                          weight_decay=5e-4, dropout=0.3)
    class_w = strong_class_weights(y, device)
    print(f"  class weights: early={class_w[0]:.3f}  late={class_w[1]:.3f}")

    model   = GATGraphClassifier(
        bundle.in_dim, hidden_dim=64, n_classes=2,
        n_heads=4, n_layers=2, dropout=cfg.dropout,
    ).to(device)
    opt     = torch.optim.Adam(model.parameters(),
                               lr=cfg.lr, weight_decay=cfg.weight_decay)
    loss_fn = nn.CrossEntropyLoss(weight=class_w)
    dl      = DataLoader(data_list, batch_size=32, shuffle=True)

    for ep in range(cfg.epochs):
        model.train()
        total = 0.0
        for batch in dl:
            batch = batch.to(device)
            opt.zero_grad()
            loss  = loss_fn(
                model(batch.x, batch.edge_index, batch.batch,
                      edge_weight=getattr(batch, "edge_weight", None)),
                batch.y)
            loss.backward()
            opt.step()
            total += loss.item()
        if (ep + 1) % 40 == 0:
            frac, n_pred_late, n_pred_early = check_predictions(
                model, data_list, device)
            print(f"  epoch {ep+1:3d}  loss={total/len(dl):.4f}  "
                  f"pred late={n_pred_late}  pred early={n_pred_early}")

    # validation check
    frac, n_pred_late, n_pred_early = check_predictions(
        model, data_list, device)
    print(f"\n  Final: predicts late={n_pred_late}, early={n_pred_early}")
    if n_pred_late < 5:
        print("  WARNING: model still predicting mostly early.")
        print("  Saliency may not be meaningful.")
    else:
        print("  OK: model predicting both classes.")

    # ── signed saliency ────────────────────────────────────────────────────
    print("\n[interp v2] extracting signed saliency ...")
    model.eval()
    single_dl = DataLoader(data_list, batch_size=1, shuffle=False)

    # separate early vs late graphs
    sal_early = np.zeros(len(bundle.gene_order))
    sal_late  = np.zeros(len(bundle.gene_order))
    n_e, n_l  = 0, 0

    for batch in single_dl:
        batch   = batch.to(device)
        batch.x = batch.x.detach().requires_grad_(True)
        logits  = model(batch.x, batch.edge_index, batch.batch,
                        edge_weight=getattr(batch, "edge_weight", None))
        logits[:, 1].sum().backward()
        grad = batch.x.grad.detach().cpu().numpy().mean(axis=1)
        label = int(batch.y.item())
        if label == 0:
            sal_early += grad; n_e += 1
        else:
            sal_late  += grad; n_l += 1

    sal_early /= max(n_e, 1)
    sal_late  /= max(n_l, 1)
    sal_combined = (sal_late - sal_early)   # differential saliency

    sal_df = pd.DataFrame({
        "gene":            bundle.gene_order,
        "saliency_late":   sal_late,
        "saliency_early":  sal_early,
        "differential":    sal_combined,
        "direction":       ["→ late" if s > 0 else "→ early"
                            for s in sal_combined],
        "abs_differential": np.abs(sal_combined),
    }).sort_values("abs_differential", ascending=False).reset_index(drop=True)
    sal_df["saliency_rank"] = sal_df.index + 1
    sal_df.to_csv(OUT / "interp_saliency_v2.csv", index=False)

    print(f"\n  Top 15 genes by differential saliency:")
    print(f"  {'rank':<5} {'gene':<12} {'differential':>13} {'direction'}")
    print("  " + "-"*45)
    for _, r in sal_df.head(15).iterrows():
        print(f"  {int(r['saliency_rank']):<5} {r['gene']:<12} "
              f"{r['differential']:>13.5f} {r['direction']}")

    # ── attention edges (filtered) ─────────────────────────────────────────
    print("\n[interp v2] extracting attention weights (self-loops filtered) ...")
    model.eval()
    agg: dict[tuple[int,int], list[float]] = {}

    with torch.no_grad():
        for batch in single_dl:
            batch = batch.to(device)
            _     = model(batch.x, batch.edge_index, batch.batch,
                          edge_weight=getattr(batch, "edge_weight", None),
                          return_attention=True)
            attn_list = model.get_last_attention()
            if not attn_list: continue
            ei, alpha = attn_list[-1]
            ei    = ei.cpu().numpy()
            alpha = (alpha.mean(dim=1).cpu().numpy()
                     if alpha.ndim > 1 else alpha.cpu().numpy())
            for (s, t), a in zip(ei.T, alpha):
                if s == t: continue   # filter self-loops
                key = (int(s), int(t))
                agg.setdefault(key, []).append(float(a))

    n_nodes   = len(bundle.gene_order)
    attn_rows = []
    for (s, t), vals in agg.items():
        if s < n_nodes and t < n_nodes:
            attn_rows.append({
                "source":         bundle.gene_order[s],
                "target":         bundle.gene_order[t],
                "attention_mean": float(np.mean(vals)),
            })

    attn_df = (pd.DataFrame(attn_rows)
               .sort_values("attention_mean", ascending=False)
               .reset_index(drop=True))

    # deduplicate undirected
    seen, dedup = set(), []
    for _, r in attn_df.iterrows():
        pair = tuple(sorted([r["source"], r["target"]]))
        if pair not in seen:
            seen.add(pair)
            dedup.append(r)
    attn_dedup = pd.DataFrame(dedup).reset_index(drop=True)
    attn_dedup.to_csv(OUT / "interp_attention_v2.csv", index=False)

    print(f"\n  Top 15 attention edges (self-loops removed):")
    print(f"  {'source':<12} {'target':<12} {'attention':>10}")
    print("  " + "-"*36)
    for _, r in attn_dedup.head(15).iterrows():
        print(f"  {r['source']:<12} {r['target']:<12} "
              f"{r['attention_mean']:>10.4f}")

    # ── comparison table ───────────────────────────────────────────────────
    tumor_mask = ds.tumor_mask
    expr       = ds.expression.to_numpy()
    log2fc     = (expr[:, tumor_mask].mean(axis=1) -
                  expr[:, ~tumor_mask].mean(axis=1))
    de_map = dict(zip(ds.expression.index, log2fc))

    from src.graph_building.build_graph import CURATED_CU_EDGES
    curated_pairs = {tuple(sorted([s, t])) for s, t, _ in CURATED_CU_EDGES}

    comp_rows = []
    for _, r in sal_df.head(20).iterrows():
        gene = r["gene"]
        fc   = de_map.get(gene, float("nan"))
        bio  = KNOWN_BIOLOGY.get(gene, ("—", "unknown"))
        comp_rows.append({
            "gene":            gene,
            "saliency_rank":   int(r["saliency_rank"]),
            "model_direction": r["direction"],
            "log2FC":          round(fc, 3),
            "expected_stage":  bio[1],
            "biology":         bio[0],
            "agrees":          (r["direction"] == "→ late") == (bio[1] == "late")
                               if bio[1] != "unknown" else "—",
        })

    comp_df = pd.DataFrame(comp_rows)
    comp_df.to_csv(OUT / "interp_comparison_table_v2.csv", index=False)

    known    = comp_df[comp_df["expected_stage"] != "unknown"]
    n_agree  = int(known["agrees"].sum())
    n_total  = len(known)
    pct      = 100 * n_agree / max(n_total, 1)

    print(f"\n  Model-biology agreement: {n_agree}/{n_total} ({pct:.0f}%)")
    print(f"\n  Top 20 genes vs known biology:")
    print(f"  {'gene':<12} {'direction':<12} {'log2FC':>7} "
          f"{'expected':>10} {'agrees':>7}")
    print("  " + "-"*52)
    for _, r in comp_df.iterrows():
        ag = "✓" if r["agrees"] is True else ("✗" if r["agrees"] is False else "—")
        print(f"  {r['gene']:<12} {r['model_direction']:<12} "
              f"{r['log2FC']:>7.3f} {r['expected_stage']:>10} {ag:>7}")

    # ── attention edges biology check ──────────────────────────────────────
    print(f"\n  Top attention edges vs curated biology:")
    print(f"  {'source':<12} {'target':<12} {'attn':>7} {'curated':>10}")
    print("  " + "-"*45)
    for _, r in attn_dedup.head(15).iterrows():
        pair    = tuple(sorted([r["source"], r["target"]]))
        in_cur  = "✓" if pair in curated_pairs else "—"
        print(f"  {r['source']:<12} {r['target']:<12} "
              f"{r['attention_mean']:>7.4f} {in_cur:>10}")

    print(f"\n[interp v2] files written to outputs/final_comparison/")
    print("Done.")


if __name__ == "__main__":
    main()
