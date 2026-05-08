# Model Comparison — TCGA-LIHC Cu Proteome (Real Data)

## 1. Predictive performance (5-fold CV, tumor vs normal)

| family | model | accuracy | balanced accuracy | F1 | ROC-AUC |
|---|---|---:|---:|---:|---:|
| Classical (tabular 54-gene vector) | logreg | 0.983 | 0.973 | 0.991 | 0.997 |
| Classical | random forest | 0.976 | 0.917 | 0.987 | 0.998 |
| Classical | **SVM (RBF)** | **0.986** | **0.983** | **0.992** | 0.998 |
| Graph-based (54-node per-patient graph) | GCN | 0.868 | 0.899 | 0.913 | 0.984 |
| Graph-based | **GAT** | **0.934** | **0.945** | **0.961** | **0.995** |

On raw prediction alone, **SVM narrowly wins** (AUC 0.998, balanced accuracy 0.983). GAT is the best graph model and is 0.003 AUC behind SVM — effectively a tie on a cohort of this size, but not a win.

This is the expected outcome for a small, well-separated tumor-vs-normal task with 54 features. The takeaway is **not** that the GNN is better; it is that the GNN reaches parity while offering structural information the tabular models cannot express.

## 2. Interpretability

| model | top feature / gene | structure | comment |
|---|---|---|---|
| logreg | signed coefficient per gene | flat | easy to read; no biology beyond the training data |
| random forest | Gini importance per gene | flat | same comment |
| SVM-RBF | permutation importance per gene | flat | no interpretable attribution for individual samples |
| GCN | gradient saliency per gene | fixed graph | adds graph context but no edge-level explanation |
| **GAT** | saliency + **per-edge attention** | fixed graph | attention edges recover textbook Cu-handling pairs (CCS↔SLC31A1, SOD3↔ATOX1, ATP7B↔COMMD1, MT-CO2↔SCO2) |

GAT attention is the only output in this pipeline that surfaces **known Cu biology** (not just expression statistics). For a preliminary / collaboration context this is the most valuable signal.

## 3. Practical usefulness for preliminary biological data generation

- **SVM** is the right choice if the only question is "can we classify tumor vs normal on the Cu proteome?" (The answer is clearly yes.)
- **GAT** is the right choice if the question is "which copper-axis interactions does this classification hinge on?" The per-sample attention and saliency outputs enable discussion of specific genes and edges with a biologist.
- **Classical + GNN together** is the pragmatic combination: SVM confirms the signal is real; GAT explains it in biology-friendly terms.

## 4. How to interpret small differences

- 0.003 AUC between SVM and GAT is within the fold-to-fold noise for 5-fold CV with 50 normals. Do **not** make a predictive-performance claim either way.
- The 4–5 point gap in balanced accuracy between GAT and SVM (0.945 vs 0.983) is more meaningful and reflects the GNN's weaker handling of the tiny normal class. Could likely be closed with more epochs, SAGPooling, or graph augmentation — deferred as "future work".

## 5. Summary

The copper proteome carries a strong, learnable tumor-vs-normal signal in LIHC. Every reasonable classifier exceeds 0.98 AUC. For a 54-gene problem, the GNN's value is interpretability — specifically GAT attention, which reproduces curated Cu-handling edges and nominates ATOX1 / ATP7B / CP / SOD1 / MT-CO1/CO2 as the axis-level drivers. Tabular models are a fair **scientific baseline** and should be kept as a sanity check for any future graph-based extension (e.g. multi-omics, subtype prediction, survival).

See `real_gnn_summary.md` and `real_baseline_summary.md` for per-block details, and `lihc_biological_interpretation.md` for the biology-facing read-out.
