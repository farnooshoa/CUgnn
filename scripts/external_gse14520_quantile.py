"""Follow-up: quantile-normalise GSE14520 per-gene to TCGA distribution.

Cross-platform RNA-seq vs microarray comparison is the known hard case.
Per-gene quantile normalisation (rank-map each external sample to the TCGA
distribution for that gene) is the lightest-weight fix that usually halves
the gap.

If the 0.60 external AUC is driven by scale mismatch (not biology), this
script should bump it noticeably. If it doesn't, then the problem is
missing Cu-gene coverage (11/54 absent on the Affy platform), which cannot
be patched by normalisation.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from torch_geometric.loader import DataLoader
from torch_geometric.data import Data
from sklearn.metrics import (
    accuracy_score, balanced_accuracy_score, f1_score, roc_auc_score,
)

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.preprocessing import load_lihc_dataset
from src.graph_building import build_functional_graph
from src.gnn_models.dataset import build_graph_dataset, graph_to_edge_tensors
from src.gnn_models.models import GATGraphClassifier
from src.gnn_models.train import _class_weights, TrainConfig

OUT = ROOT / "outputs" / "final_comparison"


def quantile_normalize_to_reference(external: pd.DataFrame, reference: pd.DataFrame) -> pd.DataFrame:
    """Per-gene: map each external value to the matching-rank value in reference."""
    out = external.copy()
    for gene in external.index:
        if gene not in reference.index:
            continue
        ref_sorted = np.sort(reference.loc[gene].dropna().values)
        if len(ref_sorted) == 0:
            continue
        row = external.loc[gene]
        if row.isna().all():
            continue
        ranks = row.rank(pct=True, method="average")
        idx = np.clip((ranks * (len(ref_sorted) - 1)).round().fillna(0).astype(int),
                      0, len(ref_sorted) - 1)
        mapped = ref_sorted[idx.values]
        mapped[row.isna().values] = np.nan
        out.loc[gene] = mapped
    return out


def train_tcga_gat(ds, graph, cfg):
    bundle = build_graph_dataset(ds, graph, zscore_train_mask=np.ones(ds.n_samples, dtype=bool))
    y = np.array([int(d.y.item()) for d in bundle.data_list])
    class_w = _class_weights(y, cfg.device)
    model = GATGraphClassifier(bundle.in_dim, cfg.hidden, n_classes=2,
                                n_heads=4, n_layers=cfg.n_layers,
                                dropout=cfg.dropout).to(cfg.device)
    opt = torch.optim.Adam(model.parameters(), lr=cfg.lr, weight_decay=cfg.weight_decay)
    loss_fn = torch.nn.CrossEntropyLoss(weight=class_w)
    dl = DataLoader(bundle.data_list, batch_size=cfg.batch_size, shuffle=True)
    for _ in range(cfg.epochs):
        model.train()
        for batch in dl:
            batch = batch.to(cfg.device)
            opt.zero_grad()
            logits = model(batch.x, batch.edge_index, batch.batch,
                           edge_weight=getattr(batch, "edge_weight", None))
            loss = loss_fn(logits, batch.y)
            loss.backward(); opt.step()
    return model, bundle


def score_external(model, ds, graph, gse, gse_md, gene_order):
    cfg = TrainConfig(model="gat")
    tcga_mu = ds.expression.mean(axis=1)
    tcga_sd = ds.expression.std(axis=1).replace(0, 1)
    gse = gse.reindex(gene_order)
    for g in gse.index:
        if gse.loc[g].isna().all():
            gse.loc[g] = tcga_mu[g]
    gse = gse.fillna(gse.mean(axis=1).fillna(tcga_mu))
    z = gse.sub(tcga_mu, axis=0).div(tcga_sd, axis=0)

    cop = ds.copper_genes.set_index("gene_symbol")
    cats = ["transporter", "enzyme", "other_or_unknown"]
    cat_feats = np.zeros((len(gene_order), len(cats)), dtype=np.float32)
    for i, g in enumerate(gene_order):
        c = cop.loc[g, "functional_category"] if g in cop.index else "other_or_unknown"
        cat_feats[i, cats.index(c) if c in cats else 2] = 1.0
    edge_index, edge_weight = graph_to_edge_tensors(graph, gene_order)
    labels = (gse_md.loc[gse.columns, "sample_type"].str.lower() == "tumor").astype(int).to_numpy()

    data_list = []
    for i, sid in enumerate(gse.columns):
        x = np.concatenate([
            gse[sid].to_numpy().astype(np.float32).reshape(-1, 1),
            z[sid].to_numpy().astype(np.float32).reshape(-1, 1),
            cat_feats,
        ], axis=1)
        data_list.append(Data(x=torch.tensor(x, dtype=torch.float32),
                              edge_index=edge_index, edge_weight=edge_weight,
                              y=torch.tensor([int(labels[i])], dtype=torch.long)))
    dl = DataLoader(data_list, batch_size=cfg.batch_size, shuffle=False)
    model.eval()
    ys, ps, preds = [], [], []
    with torch.no_grad():
        for batch in dl:
            batch = batch.to(cfg.device)
            logits = model(batch.x, batch.edge_index, batch.batch,
                           edge_weight=getattr(batch, "edge_weight", None))
            prob = torch.softmax(logits, dim=-1)[:, 1].cpu().numpy()
            pr = logits.argmax(dim=-1).cpu().numpy()
            ys.extend(batch.y.cpu().numpy().tolist())
            ps.extend(prob.tolist())
            preds.extend(pr.tolist())
    ys = np.array(ys); ps = np.array(ps); preds = np.array(preds)
    return {"accuracy": accuracy_score(ys, preds),
            "balanced_accuracy": balanced_accuracy_score(ys, preds),
            "f1": f1_score(ys, preds, zero_division=0),
            "roc_auc": roc_auc_score(ys, ps) if len(set(ys)) > 1 else float("nan"),
            "n_tumor": int((ys == 1).sum()), "n_normal": int((ys == 0).sum())}


def main():
    ds = load_lihc_dataset(require_real=True)
    graph = build_functional_graph(ds.expression.index.tolist(), ds.copper_genes)
    gse_expr = pd.read_csv(ROOT / "data" / "gse14520_expression.tsv", sep="\t", index_col=0)
    gse_md = pd.read_csv(ROOT / "data" / "gse14520_metadata.tsv", sep="\t", index_col=0)

    cfg = TrainConfig(model="gat", epochs=80)
    print("[qn] training GAT on TCGA ...")
    model, _ = train_tcga_gat(ds, graph, cfg)

    # Non-normalised baseline
    print("[qn] scoring external (no quantile norm) ...")
    no_qn = score_external(model, ds, graph, gse_expr.copy(), gse_md, ds.expression.index.tolist())
    print("  ", no_qn)

    # With quantile normalisation to TCGA
    print("[qn] quantile-normalising GSE14520 to TCGA per-gene distribution ...")
    gse_qn = quantile_normalize_to_reference(gse_expr, ds.expression)
    print("[qn] scoring external (after quantile norm) ...")
    with_qn = score_external(model, ds, graph, gse_qn, gse_md, ds.expression.index.tolist())
    print("  ", with_qn)

    rows = [
        {"cohort": "TCGA-LIHC 5-fold CV", "normalization": "TCGA-internal",
         "roc_auc": 0.993, "balanced_accuracy": 0.913},
        {"cohort": "GSE14520 external", "normalization": "z-score only",
         **{k: no_qn[k] for k in ("roc_auc", "balanced_accuracy")}},
        {"cohort": "GSE14520 external", "normalization": "quantile + z-score",
         **{k: with_qn[k] for k in ("roc_auc", "balanced_accuracy")}},
    ]
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "external_gse14520_quantile_metrics.csv", index=False)
    print("\n" + df.round(4).to_string(index=False))

    with open(OUT / "external_gse14520.md", "a") as f:
        f.write(f"""

---

## Follow-up — quantile normalisation

Hypothesis: the internal → external AUC drop is driven by **platform scale mismatch** (Affymetrix RMA vs Illumina FPKM-UQ), not biology. Per-gene quantile-mapping of GSE14520 to TCGA distributions should close most of the gap if this hypothesis is right.

| normalisation | ROC-AUC | balanced acc |
|---|---:|---:|
| z-score only (TCGA µ/σ) | {no_qn['roc_auc']:.3f} | {no_qn['balanced_accuracy']:.3f} |
| **quantile + z-score** | **{with_qn['roc_auc']:.3f}** | **{with_qn['balanced_accuracy']:.3f}** |
| Δ | **{with_qn['roc_auc'] - no_qn['roc_auc']:+.3f}** | {with_qn['balanced_accuracy'] - no_qn['balanced_accuracy']:+.3f} |

### Verdict
- If Δ AUC > 0.10 → platform mismatch was the main issue; reporting with
  quantile-normalised inference is the honest comparison.
- If Δ AUC < 0.05 → the drop is **biological**, not technical: the 11 Cu
  genes missing from GSE14520 (notably mitochondrial MT-CO1, MT-CO2 and
  enzymes SCO1, ENOX1/2, HEPHL1, LOXL3/4, MEMO1) remove most of the
  features the TCGA-trained model relies on. Remediation requires a
  platform with full Cu-proteome coverage (e.g. RNA-seq cohort GSE36376
  or CPTAC LIHC).
""")
    print(f"[qn] appended to {OUT/'external_gse14520.md'}")


if __name__ == "__main__":
    main()
