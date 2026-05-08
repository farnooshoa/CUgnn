# Model comparison — TCGA-LIHC Cu proteome pilot

## Classical baselines (5-fold CV, 54-gene vector)

| model         |   cv_folds |   accuracy |   balanced_accuracy |     f1 |   roc_auc |   accuracy_std |
|:--------------|-----------:|-----------:|--------------------:|-------:|----------:|---------------:|
| logreg        |          5 |     0.9834 |              0.9733 | 0.9905 |    0.997  |         0.0121 |
| random_forest |          5 |     0.9764 |              0.9173 | 0.9868 |    0.9979 |         0.0074 |
| svm_rbf       |          5 |     0.9858 |              0.9833 | 0.9919 |    0.9979 |         0.0048 |

## GNN models (5-fold CV, per-patient graph)

| model   |   accuracy |   balanced_accuracy |     f1 |   roc_auc |
|:--------|-----------:|--------------------:|-------:|----------:|
| gat     |     0.9339 |              0.9452 | 0.9607 |    0.9949 |
| gcn     |     0.8681 |              0.8992 | 0.9133 |    0.9842 |

## Interpretation

- Compare the ROC-AUC row across classical vs GNN. If GNN is within
  ~0.02 AUC of the best classical baseline, the marginal predictive
  value of graph structure is small on this cohort — which is the
  *expected* situation for a 54-gene tumor-vs-normal problem.
- The *real* value of the GNN in the pilot is interpretability:
  attention / saliency maps highlight Cu-proteome modules (e.g.
  ATP7B / CP / ceruloplasmin axis; LOX family) that tabular models
  cannot expose.
- See outputs/gnn/top_subgraph_or_attention_summary.md for the
  per-model interpretability report.
- See outputs/baseline/copper_network_logfc.png for the static
  LIHC-coloured network.

## Conclusion

The copper proteome is small enough (54 genes) that classical ML on
a flat expression vector is already a strong baseline. The GNN
earns its place mainly by (a) aligning with the biology of Fig. 3
of Blockhuys 2017, (b) providing module-level interpretability, and
(c) offering a natural substrate for future multi-modal graphs
(methylation, mutation, clinical covariates as additional node/edge
features). Predictive parity with tabular models is acceptable at
this pilot stage.