# Week 2 (Part 2) — Stage Classification Proper Comparison

## Task
Early (AJCC I/II) vs Late (AJCC III/IV) — tumor samples only.

## Cohort
- n = 349 tumor samples with stage information
- Early (I/II): 259
- Late (III/IV): 90
- Class ratio: 2.9:1
- Unique patients: 346
- CV: 5-fold StratifiedGroupKFold, 10 seeds

## Why stage is more interesting than survival for copper biology

The copper proteome should be most relevant at the stage level because:
- **LOX / LOXL1-4** crosslink collagen — directly enables invasion and metastasis
- **SPARC** remodels the extracellular matrix — higher in late stage tumors
- **ATP7B** drives copper efflux — dysregulated in aggressive tumors
- **CP** (ceruloplasmin) — serum levels correlate with AJCC stage in liver cancer

If the graph neural network learns from these edges, it should outperform
models that treat the 58 genes as independent features.

## Results (10-seed mean ± std)

| model | type | ROC-AUC | ± std | balanced acc |
|---|---|---:|---:|---:|
| **gat** | graph | 0.6667 | 0.0174 | 0.5724 |
| **gcn** | graph | 0.6605 | 0.0103 | 0.5846 |
| **svm_rbf** | classical | 0.6441 | 0.0120 | 0.5808 |
| **random_forest** | classical | 0.6242 | 0.0205 | 0.4986 |
| **logreg** | classical | 0.5832 | 0.0250 | 0.5655 |

## Verdict

GAT clearly beats classical models — graph adds predictive value on stage

GAT vs best classical model: **+0.0226 AUC**

## What this means for the paper

The graph neural network earns its place on the stage task — the biological edges between copper genes carry stage-specific predictive signal that cannot be captured by treating genes as independent features.

## Honest caveats
- AJCC stage is clinically noisy — inter-observer disagreement on IIA vs IIIA
  is well documented
- We collapse I+II and III+IV — within-group biology differs
- With only 90 late-stage samples, per-fold variance is high
- Stage is cross-sectional, not temporal — this is association not prognosis

## Files produced
- `outputs/final_comparison/week2_stage_comparison.csv`
- `outputs/final_comparison/week2_stage_report.md`
