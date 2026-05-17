# Multi-Seed Stability Report

- torch           : 2.12.0+cpu
- seeds           : [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
- model           : GAT, 5-fold StratifiedGroupKFold

## Results

| task | AUC mean | AUC std | AUC min | AUC max |
|---|---:|---:|---:|---:|
| survival_3yr | **0.6799** | 0.0117 | 0.6623 | 0.6957 |

## Verdict
- std < 0.02 → **STABLE** — report mean ± std
- std 0.02–0.04 → **MODERATE** — report mean ± std, note limitation
- std > 0.04 → **UNSTABLE** — do not quote a single number

## Per-seed detail

        task  seed  roc_auc  balanced_accuracy
survival_3yr     1 0.662295           0.588972
survival_3yr     2 0.692398           0.572515
survival_3yr     3 0.680925           0.578739
survival_3yr     4 0.673768           0.591688
survival_3yr     5 0.693901           0.540351
survival_3yr     6 0.695684           0.620510
survival_3yr     7 0.675021           0.513701
survival_3yr     8 0.664439           0.608647
survival_3yr     9 0.677527           0.554219
survival_3yr    10 0.683124           0.538972