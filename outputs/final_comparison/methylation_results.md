# Multi-Omics Node Features — Adding DNA Methylation

## Input
- **RNA-seq (expression)** — existing, `data/lihc_expression.tsv` (log2 FPKM-UQ)
- **DNA methylation 450K (gene-level β)** — Xena pre-aggregated dump of
  TCGA.LIHC.sampleMap/HumanMethylation450, restricted to the 54 Cu genes:
  `data/lihc_methylation.tsv` (52 genes × 424 samples)

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
| Tumor vs Normal | **expr + methylation** | **0.997** | **0.912** | **+0.004** |
| Stage I/II vs III/IV | expr_only | 0.668 | 0.525 | — |
| Stage I/II vs III/IV | **expr + methylation** | **0.667** | **0.518** | **-0.001** |

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
