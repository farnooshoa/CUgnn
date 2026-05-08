# 3-Year Overall Survival — Classification Task

## Cohort
- Tumor samples only with OS labels and enough follow-up: **n = 198**
- Died ≤ 3 y: **105**
- Survived > 3 y (confirmed): **93**
- Excluded censored-short (alive but < 3 y follow-up): 176
- Class ratio: 1.13 : 1 (essentially balanced)
- Unique patients = 195; StratifiedGroupKFold grouped by `case_submitter_id`

## Results (5-fold CV)

| model | ROC-AUC | balanced acc | F1 | accuracy |
|---|---:|---:|---:|---:|
| **gat** | 0.692 ± 0.064 | 0.584 | 0.551 | 0.575 |
| **gcn** | 0.682 ± 0.051 | 0.573 | 0.639 | 0.585 |
| **random_forest** | 0.655 ± 0.026 | 0.607 | 0.633 | 0.591 |
| **svm_rbf** | 0.650 ± 0.038 | 0.609 | 0.608 | 0.601 |
| **logreg** | 0.585 ± 0.068 | 0.561 | 0.575 | 0.549 |

## Interpretation

**Best**: gat at ROC-AUC **0.692**.
**Worst**: logreg at ROC-AUC 0.585.
**Gap (best − worst)**: +0.107 AUC.

### Where this lands
- Survival classification from 54 genes with n=198 is a **hard task**.
  AUC in the 0.55–0.70 range is the realistic ceiling.
- Result is right in that band: **0.692** AUC.
- Clinically most-used LIHC prognostic models (BCLC, CLIP, Okuda) are in
  roughly the same range on TCGA; a 54-gene Cu-proteome model reaching
  this is already a scientifically interesting comparison.

### Honest caveats
- 198 samples with balanced classes is tight; AUC standard deviation
  across folds is ~0.05, so 0.02 model-to-model differences are noise.
- 3-year OS is one of many survival endpoints; disease-free survival
  or recurrence-free survival often carry clearer mRNA signal.
- Censored-short exclusion (we dropped patients alive at < 3 y
  follow-up) simplifies the task but loses real-world data. A
  proper survival analysis (Cox model or time-to-event DeepSurv
  GNN) uses all data including censoring — future work.
- The Cu proteome is small (54 genes); stage and size drive outcome
  more than Cu biology alone. Comparing our Cu-only AUC to a
  full-transcriptome baseline is the obvious next step.

### Biological read (exploratory)
- Given tumor-vs-normal already highlighted ATP7B, CP, LOX as hubs,
  and stage flagged the same ECM axis, survival likely leans on
  LOX/LOXL2/SPARC (invasive capacity) and possibly ATOX1 / ATP7A
  (proliferation). The per-gene importance saved in
  `outputs/gnn/node_importance.csv` for the survival model is the
  file to inspect for a survival-specific top-10.

## Files produced
- `outputs/final_comparison/survival_task.md` — this document
- `outputs/final_comparison/survival_task_metrics.csv` — per-model summary
