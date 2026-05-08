# Stage Classification — Early (I/II) vs Late (III/IV)

## Cohort
- Tumor samples only: **349** (normals excluded — no meaningful stage)
- Early (AJCC I / II / IIB): **259**
- Late (AJCC III* / IV*): **90**
- Class ratio: 2.9 : 1 (gentler than tumor-vs-normal's 7.5:1)
- Excluded: Stage 0 (n=1), missing stage (n=24)
- Unique patients = 346; StratifiedGroupKFold grouped by `case_submitter_id`

## Results (5-fold CV)

| model | ROC-AUC | balanced acc | F1 | accuracy |
|---|---:|---:|---:|---:|
| **gcn** | 0.668 ± 0.059 | 0.571 | 0.306 | 0.614 |
| **gat** | 0.668 ± 0.029 | 0.525 | 0.276 | 0.592 |
| **svm_rbf** | 0.659 ± 0.050 | 0.619 | 0.418 | 0.688 |
| **random_forest** | 0.643 ± 0.060 | 0.515 | 0.095 | 0.737 |
| **logreg** | 0.585 ± 0.058 | 0.565 | 0.382 | 0.585 |

## Interpretation

**This is the right difficulty band** for a copper-proteome pilot. AUC in the
0.55–0.70 range is what a small, curated feature set should give on a
clinically-meaningful sub-classification. Here:

- Best model: **gcn** at ROC-AUC 0.668
- Worst model: logreg at ROC-AUC 0.585
- Gap: +0.083

### What to compare with tumor-vs-normal (AUC ~0.99)
- Tumor-vs-normal was trivial; any model hit 0.99+ and the graph did not
  change predictions.
- Stage classification drops everyone to the 0.55–0.70 regime — this is
  where model comparison is actually informative.
- If GAT ≥ classical here by more than ~0.02 AUC, the graph is now doing
  real predictive work (not just interpretability).
- If all models are similar, the Cu proteome carries limited stage-specific
  signal at the mRNA level — a negative but publishable result.

### Biological read
AJCC stage in HCC is largely driven by tumor size, vascular invasion, and
metastasis — processes that the ECM / LOX-family sub-module of the Cu
proteome should in principle touch. It is plausible that LOX, LOXL1-4,
SPARC, and possibly ATOX1 are the main stage-predictive nodes. The per-fold
node importance for the GAT in this task is saved in
`stage_task_gnn_folds.csv` — a follow-up pass can extract saliency per
tumor graph to rank stage-specific genes.

### Honest caveats
- AJCC stage is clinically noisy; inter-observer disagreement on Stage II
  vs IIIA is well-documented.
- We binarise I/II vs III/IV; the I vs II and IIIA vs IV boundaries carry
  different biology we are collapsing.
- With only ~90 late-stage samples, per-fold variance is high; do not
  over-interpret 0.02 differences.
- Stage is known at diagnosis and is **not** a temporal prediction — this
  is a cross-sectional association, not prognosis.

## Files produced
- `outputs/final_comparison/stage_task.md` — this document
- `outputs/final_comparison/stage_task_metrics.csv` — per-model summary
- `outputs/final_comparison/stage_task_gnn_folds.csv` — GNN per-fold metrics
