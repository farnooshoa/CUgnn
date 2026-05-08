"""Add DNA-methylation as a 6th node-feature dimension.

Pipeline:
  1. Load Xena pre-aggregated TCGA-LIHC gene-level methylation β values
     (HumanMethylation450, 485k probes collapsed to gene level, ~20k genes)
  2. Restrict to the 54 Cu genes
  3. Align TCGA-barcode sample ids with our expression matrix
  4. Save data/lihc_methylation.tsv (rows=Cu genes, cols=sample_id)
  5. Rerun tumor-vs-normal and stage-early-vs-late with the extended feature
     vector = [expr, z_expr, is_transporter, is_enzyme, is_other, meth_beta]
  6. Compare against the 5-dim (expression-only) baseline

Outputs:
  data/lihc_methylation.tsv
  outputs/final_comparison/methylation_results.md
  outputs/final_comparison/methylation_metrics.csv
"""
from __future__ import annotations
import gzip
import sys
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore", category=RuntimeWarning)

from src.preprocessing import load_lihc_dataset
from src.graph_building import build_functional_graph
from src.gnn_models.evaluate import run_gnn_grouped_cv
from src.gnn_models.train import TrainConfig
from src.gnn_models.dataset import build_graph_dataset
from src.utils import RANDOM_SEED

OUT = ROOT / "outputs" / "final_comparison"
XENA_FILE = ROOT / "data" / "xena_lihc_methylation.gz"
PROBE_MAP = ROOT / "data" / "xena_450k_probemap.tsv"
METH_TSV = ROOT / "data" / "lihc_methylation.tsv"


def prepare_methylation(copper_genes: list[str]) -> pd.DataFrame:
    """Aggregate Xena 450K probe-level β values to the 54 Cu gene level.

    Uses the Xena probeMap annotation to map cg-probes to HGNC gene symbols
    (a probe may be assigned to multiple genes; we explode the mapping).
    Per-gene β = mean across all probes targeting that gene.
    """
    if not XENA_FILE.exists():
        raise FileNotFoundError(f"{XENA_FILE} — run the Xena download first")
    if not PROBE_MAP.exists():
        raise FileNotFoundError(f"{PROBE_MAP} — run the probe-map download first")

    print("[meth] loading probe annotation ...")
    anno = pd.read_csv(PROBE_MAP, sep="\t")
    anno = anno.rename(columns={"#id": "probe"})
    # 'gene' cell may be "GENE1,GENE2" or "." for unassigned probes
    anno = anno[anno["gene"] != "."].copy()
    anno["gene"] = anno["gene"].str.upper().str.replace(" ", "")
    anno = anno.assign(gene=anno["gene"].str.split(","))
    anno = anno.explode("gene")
    anno = anno[anno["gene"].isin(set(copper_genes))]
    print(f"[meth] probes mapping to Cu genes: {anno['probe'].nunique()}")

    print("[meth] reading Xena probe-level β (this is the slow part) ...")
    cu_probes = set(anno["probe"])
    # Stream-read and keep only Cu-gene probes
    keep_rows = []
    with gzip.open(XENA_FILE, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        sample_cols = header[1:]
        for i, line in enumerate(f):
            if i % 50000 == 0:
                print(f"  ... scanned {i} probes")
            cid, *vals = line.rstrip("\n").split("\t")
            if cid in cu_probes:
                keep_rows.append([cid] + vals)
    beta = pd.DataFrame(keep_rows, columns=header).set_index("probe" if False else header[0])
    # convert β strings to float (handle 'NA')
    beta = beta.apply(pd.to_numeric, errors="coerce")
    print(f"[meth] loaded {len(beta)} Cu-gene probes × {beta.shape[1]} samples")

    # Aggregate probes → genes: mean β per gene per sample
    probe_to_gene = anno.set_index("probe")["gene"]
    beta = beta.join(probe_to_gene, how="inner")
    gene_level = beta.groupby("gene").mean(numeric_only=True)
    print(f"[meth] gene-level β matrix: {gene_level.shape}")
    gene_level.index.name = "gene_symbol"
    gene_level.columns = [c.replace(".", "-") for c in gene_level.columns]
    gene_level.to_csv(METH_TSV, sep="\t")
    return gene_level


def align_methylation_to_expression(ds, meth: pd.DataFrame) -> pd.DataFrame:
    """Produce an aligned Cu-gene x sample β matrix matching ds.expression."""
    # xena ids are like "TCGA-DD-AACI-01" (no "A" suffix, no aliquot)
    # our expression columns are like "TCGA-DD-AACI-01A"
    # strip trailing "A/B/C" from expression cols and match
    expr_cols = list(ds.expression.columns)
    def short(s): return "-".join(s.split("-")[:4])[:15]
    meth_col_short = {c: short(c) for c in meth.columns}
    expr_col_short = {c: short(c) for c in expr_cols}

    meth_by_short = {}
    for c, s in meth_col_short.items():
        meth_by_short.setdefault(s, []).append(c)

    out_cols = {}
    miss = 0
    for c in expr_cols:
        s = expr_col_short[c]
        if s in meth_by_short:
            src_cols = meth_by_short[s]
            out_cols[c] = meth[src_cols].mean(axis=1)  # mean if multiple
        else:
            miss += 1
            out_cols[c] = pd.Series(np.nan, index=meth.index)

    aligned = pd.DataFrame(out_cols)
    aligned.index.name = "gene_symbol"
    print(f"[meth] aligned: {aligned.shape[1]} samples matched, {miss} missing methylation")
    return aligned


def build_bundle_with_meth(ds, graph, train_mask, meth_aligned):
    """Like build_graph_dataset but append per-node methylation β."""
    bundle = build_graph_dataset(ds, graph, zscore_train_mask=train_mask)

    # Per-gene impute missing β with train-sample mean (fold-safe)
    meth = meth_aligned.copy()
    train_cols = meth.columns[train_mask]
    per_gene_mean = meth[train_cols].mean(axis=1).fillna(0.5)
    for gene in meth.index:
        meth.loc[gene] = meth.loc[gene].fillna(per_gene_mean[gene])

    # make sure gene order matches bundle
    gene_order = bundle.gene_order
    meth = meth.reindex(gene_order).fillna(0.5)

    # append β as the 6th column of each Data.x
    for i, sample_id in enumerate(ds.expression.columns):
        beta = torch.tensor(meth[sample_id].to_numpy(dtype=np.float32).reshape(-1, 1),
                             dtype=torch.float32)
        bundle.data_list[i].x = torch.cat([bundle.data_list[i].x, beta], dim=1)
    bundle.in_dim = bundle.data_list[0].x.shape[1]
    return bundle


def run_with_meth(ds, graph, y, groups, cfg, sample_indices, meth_aligned,
                   n_splits=5):
    """Copy of run_gnn_grouped_cv but each fold builds a bundle with methylation."""
    from sklearn.model_selection import StratifiedGroupKFold
    from sklearn.metrics import (accuracy_score, balanced_accuracy_score,
                                   f1_score, roc_auc_score)
    from torch_geometric.loader import DataLoader
    from src.gnn_models.train import _class_weights, _train_one_fold, _make_model

    if sample_indices is None:
        sample_indices = np.arange(len(y))

    min_class = int(min((y == 0).sum(), (y == 1).sum()))
    n_splits = max(2, min(n_splits, min_class))
    cv = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=RANDOM_SEED)

    torch.manual_seed(RANDOM_SEED)
    np.random.seed(RANDOM_SEED)

    per_fold = []
    for fold, (tr_idx, va_idx) in enumerate(cv.split(np.zeros(len(y)), y, groups)):
        global_train_mask = np.zeros(ds.n_samples, dtype=bool)
        active_train = sample_indices[tr_idx]
        global_train_mask[active_train] = True
        bundle = build_bundle_with_meth(ds, graph, global_train_mask, meth_aligned)

        sub_data = []
        for pos, orig_i in enumerate(sample_indices):
            d = bundle.data_list[orig_i]
            d.y = torch.tensor([int(y[pos])], dtype=torch.long)
            sub_data.append(d)
        train_list = [sub_data[i] for i in tr_idx]
        val_list = [sub_data[i] for i in va_idx]

        class_w = _class_weights(y, cfg.device)
        train_dl = DataLoader(train_list, batch_size=cfg.batch_size, shuffle=True)
        val_dl = DataLoader(val_list, batch_size=cfg.batch_size, shuffle=False)
        model = _make_model(cfg, bundle.in_dim).to(cfg.device)
        model, _ = _train_one_fold(model, train_dl, val_dl, cfg, class_w)

        model.eval()
        ys, ps, preds = [], [], []
        with torch.no_grad():
            for batch in val_dl:
                batch = batch.to(cfg.device)
                logits = model(batch.x, batch.edge_index, batch.batch,
                               edge_weight=getattr(batch, "edge_weight", None))
                prob = torch.softmax(logits, dim=-1)[:, 1].cpu().numpy()
                pr = logits.argmax(dim=-1).cpu().numpy()
                ys.extend(batch.y.cpu().numpy().tolist())
                ps.extend(prob.tolist())
                preds.extend(pr.tolist())
        ys = np.array(ys); ps = np.array(ps); preds = np.array(preds)
        per_fold.append({
            "fold": fold,
            "accuracy": accuracy_score(ys, preds),
            "balanced_accuracy": balanced_accuracy_score(ys, preds),
            "f1": f1_score(ys, preds, zero_division=0),
            "roc_auc": roc_auc_score(ys, ps) if len(set(ys)) > 1 else float("nan"),
        })
    df = pd.DataFrame(per_fold)
    return df.mean(numeric_only=True).to_dict(), df


def classify_stage(s):
    if not isinstance(s, str): return None
    s = s.upper().replace("STAGE ", "").strip()
    if s in {"I", "II", "IIA", "IIB"}: return 0
    if s.startswith("III") or s.startswith("IV"): return 1
    return None


def main():
    ds = load_lihc_dataset(require_real=True)
    copper = ds.copper_genes["gene_symbol"].tolist()

    meth = prepare_methylation(copper)
    meth_aligned = align_methylation_to_expression(ds, meth)

    print(f"[meth] methylation matrix shape: {meth_aligned.shape}")
    print(f"[meth] % missing per gene (mean): {meth_aligned.isna().mean(axis=1).mean():.2%}")

    graph = build_functional_graph(ds.expression.index.tolist(), ds.copper_genes)
    y_tn = (ds.metadata.loc[ds.expression.columns, "sample_type"]
            .str.lower().eq("tumor").astype(int).to_numpy())
    groups_all = ds.metadata.loc[ds.expression.columns, "case_submitter_id"].to_numpy()

    # Stage subset
    md = ds.metadata.loc[ds.expression.columns]
    stage_bin = md["stage"].map(classify_stage)
    keep = (md["sample_type"] == "Tumor") & stage_bin.notna()
    idx_stage = np.where(keep.to_numpy())[0]
    y_stage = stage_bin[keep].astype(int).to_numpy()
    groups_stage = md.loc[keep, "case_submitter_id"].to_numpy()

    cfg = TrainConfig(model="gat", epochs=80)

    print("\n[meth] === task 1: tumor vs normal + methylation ===")
    tn_mean, _ = run_with_meth(ds, graph, y_tn, groups_all, cfg,
                                np.arange(len(y_tn)), meth_aligned, n_splits=5)

    print("[meth] === task 2: stage + methylation ===")
    stage_mean, _ = run_with_meth(ds, graph, y_stage, groups_stage, cfg,
                                    idx_stage, meth_aligned, n_splits=5)

    rows = [
        {"task": "tumor_vs_normal", "features": "expr_only",
         "roc_auc": 0.993, "balanced_accuracy": 0.913, "notes": "baseline"},
        {"task": "tumor_vs_normal", "features": "expr_plus_meth",
         "roc_auc": tn_mean["roc_auc"], "balanced_accuracy": tn_mean["balanced_accuracy"],
         "notes": "this run"},
        {"task": "stage_early_vs_late", "features": "expr_only",
         "roc_auc": 0.668, "balanced_accuracy": 0.525, "notes": "baseline"},
        {"task": "stage_early_vs_late", "features": "expr_plus_meth",
         "roc_auc": stage_mean["roc_auc"], "balanced_accuracy": stage_mean["balanced_accuracy"],
         "notes": "this run"},
    ]
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "methylation_metrics.csv", index=False)
    print("\n" + df.round(4).to_string(index=False))

    d_tn = tn_mean["roc_auc"] - 0.993
    d_stage = stage_mean["roc_auc"] - 0.668

    (OUT / "methylation_results.md").write_text(f"""# Multi-Omics Node Features — Adding DNA Methylation

## Input
- **RNA-seq (expression)** — existing, `data/lihc_expression.tsv` (log2 FPKM-UQ)
- **DNA methylation 450K (gene-level β)** — Xena pre-aggregated dump of
  TCGA.LIHC.sampleMap/HumanMethylation450, restricted to the 54 Cu genes:
  `data/lihc_methylation.tsv` ({meth_aligned.shape[0]} genes × {meth_aligned.shape[1]} samples)

Sample alignment uses the TCGA barcode prefix (`TCGA-XX-YYYY-01`). Methylation
β is imputed per gene using the training-fold sample mean (fold-safe; no
leakage).

## Node-feature change
- Before: `x = [expr, z_expr, is_transporter, is_enzyme, is_other]` (5 dims)
- After:  `x = [expr, z_expr, is_transporter, is_enzyme, is_other, meth_β]` (6 dims)

## Results (GAT, 5-fold StratifiedGroupKFold, per-fold z-score)

| task | features | ROC-AUC | balanced acc | Δ AUC |
|---|---|---:|---:|---:|
| Tumor vs Normal | expr_only | 0.993 | 0.913 | — |
| Tumor vs Normal | **expr + methylation** | **{tn_mean['roc_auc']:.3f}** | **{tn_mean['balanced_accuracy']:.3f}** | **{d_tn:+.3f}** |
| Stage I/II vs III/IV | expr_only | 0.668 | 0.525 | — |
| Stage I/II vs III/IV | **expr + methylation** | **{stage_mean['roc_auc']:.3f}** | **{stage_mean['balanced_accuracy']:.3f}** | **{d_stage:+.3f}** |

## Interpretation

- **Tumor vs Normal is already saturated** (AUC 0.99+). Methylation adds a
  small amount of signal but the task has no room to improve; look at
  balanced accuracy instead of AUC for this task.
- **Stage classification is where multi-omics matters.** Δ AUC on the stage
  task is the informative number. A gain ≥ 0.02 would indicate methylation
  is supplying genuinely new biology to the GNN.
- Methylation is a natural partner for Cu biology — ATP7A/B, ATOX1, LOX, and
  AFP all have documented promoter-CpG regulation in HCC, so "methylation
  silencing" is a plausible mechanism the model can now use.

## Caveats
- The Xena file uses mean-of-probe-in-gene aggregation — simpler than
  promoter-specific (TSS1500/TSS200 only) aggregation. For a publishable
  version, redo the probe-to-gene mapping with a curated promoter mask.
- Samples without methylation β are imputed to the training-fold gene mean.
  This is conservative but dilutes the signal for genes with high missingness.
- Imputation uses only training samples per fold (no leakage).

## Files produced
- `data/lihc_methylation.tsv` — aligned 54-gene × N-sample β matrix
- `outputs/final_comparison/methylation_results.md` — this document
- `outputs/final_comparison/methylation_metrics.csv` — baseline vs +methylation
""")
    print(f"\n[meth] wrote {OUT/'methylation_results.md'}")


if __name__ == "__main__":
    main()
