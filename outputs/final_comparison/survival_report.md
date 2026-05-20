# Week 2 — Proper Survival Analysis

## What changed vs the original pipeline

| Issue | Original | Fixed |
|---|---|---|
| Metric | ROC-AUC (binary classification) | **C-index** (survival ranking) |
| Censored patients | Dropped | **Included in Cox** |
| Baseline model | None | **Cox regression** |
| Stability | Single seed | **10 seeds, mean ± std** |

## Results

| model | C-index mean | C-index std | note |
|---|---:|---:|---|
| Random baseline | 0.500 | — | theoretical floor |
| Cox regression | 0.6385 | 0.0095 | uses all patients incl. censored |
| GAT | 0.5862 | 0.0196 | confirmed-outcome patients only |

GAT ROC-AUC (binary): **0.6155 ± 0.0301**

GAT vs Cox C-index gap: **-0.0523**

## Verdict

Cox beats GAT — simpler model wins on this task

## What the C-index means

A C-index of 0.5 is random guessing.
A C-index of 0.7 is considered good for a clinical survival model.
A C-index of 1.0 is perfect ranking.

Unlike ROC-AUC, the C-index:
- Does not require throwing away censored patients
- Measures ranking (who dies sooner) not just binary classification
- Is the standard metric used in clinical survival research

## Honest caveats

- Cox uses more patients than GAT (censored included vs excluded)
  so the comparison is not perfectly apples-to-apples
- GAT uses only 198 confirmed-outcome patients; Cox uses all available
- A truly fair comparison would implement a survival-aware GNN loss
  (e.g. DeepSurv or Cox-PH loss) — that is the next step
- With n=198 and 10 seeds, std around 0.01-0.02 is expected

## Files produced

- `outputs/final_comparison/week2_survival_comparison.csv`
- `outputs/final_comparison/week2_survival_report.md`
