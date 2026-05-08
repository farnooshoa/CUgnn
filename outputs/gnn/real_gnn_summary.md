# Real-Data GNN Summary — TCGA-LIHC Cu Proteome

Graph-classification experiments on real TCGA-LIHC (424 samples, 54 Cu nodes, functional graph = 60 edges).

## 1. Cross-validation metrics (real data only, 5-fold stratified)

| model | mean accuracy | balanced accuracy | F1 | ROC-AUC |
|---|---:|---:|---:|---:|
| **GCN** | 0.868 | 0.899 | 0.913 | 0.984 |
| **GAT** | **0.934** | **0.945** | **0.961** | **0.995** |

Values rounded; per-fold numbers in `model_metrics.csv`, JSON summary in `model_metrics_summary.json`.

GCN shows one weak fold (fold 2, accuracy 0.635) that drags its mean down — the remaining four folds are ≥ 0.90 accuracy. GAT is consistent across folds (accuracy range 0.87–0.98).

## 2. Class imbalance

Tumor : Normal ≈ 374 : 50 (7.5 : 1). Handled with:
- **Stratified 5-fold CV** — each fold sees ≈ 10 normals and ≈ 75 tumors.
- **Class-weighted cross-entropy** — weights = `1/count` per class, rebalancing gradient contribution.

With these in place, balanced accuracy is within 1–3 points of plain accuracy — i.e., the normal class is not being ignored. Without class weighting GCN collapsed to the majority class on fold 2 (observed before adding weights; not shown).

## 3. Which model performed better

**GAT > GCN on every metric** (+0.007 AUC, +0.046 balanced accuracy, +0.048 F1). The attention mechanism helps in a sparse graph (average degree ≈ 2.2) where uniform message passing is easily dominated by a few high-degree hubs. With 54 nodes and a relatively simple tumor-vs-normal separation, the gains are not dramatic — but GAT is the right choice for this graph.

Best model saved as the one used for final embeddings and interpretability outputs (see below).

## 4. Which copper genes were most important

From `node_importance.csv` (saliency = mean absolute gradient of the tumor logit w.r.t. each node's features, averaged over 50 graphs):

| rank | gene | importance |
|---|---|---:|
| 1 | **CP** | 1.77 |
| 2 | **SOD1** | 1.14 |
| 3 | **MT-CO1** | 1.09 |
| 4 | **ATOX1** | 1.05 |
| 5 | **MAP2K1** | 0.99 |
| 6 | **ALB** | 0.97 |
| 7 | LOXL4 | 0.89 |
| 8 | MT-CO2 | 0.83 |
| 9 | **LOX** | 0.81 |
| 10 | **ATP7B** | 0.71 |

The top-10 list overlaps substantially with the differential-expression top-10 (CP, SOD1, MT-CO1/2, ALB, LOX, MAP2K1 all appear in both). Two new entries — **ATOX1** and **ATP7B** — rank high in GNN saliency despite having relatively modest log2FC, suggesting the model is using their connectivity (both sit at the centre of the Cu-secretion module) to amplify weaker expression signals.

## 5. Does attention/saliency align with known copper biology?

Top attention-weighted edges from `top_attention_edges.csv` (excluding trivial self-loops on isolated nodes ENOX1, ENOX2, PAM, PARK7, HEPH, HEPHL1, F5, CUTC which show up as artefacts of the attention normalisation at degree-0 nodes):

| source → target | attention | known Cu biology |
|---|---:|---|
| CCS ↔ SLC31A1 | 41.9 | Cu enters via SLC31A1 → chaperoned by CCS to SOD1 |
| ATP7B ↔ AFP | 39.3 | Both in the Cu-secretory/hepatocyte axis; AFP is a classical HCC marker |
| SOD3 ↔ ATOX1 | 34.5 | **Paper-identified interaction**: ATOX1 is a TF for SOD3 |
| PARK7 ↔ SOD1 | 34.0 | Both cytosolic redox proteins |
| S100A12 ↔ S100A13, S100B, S100A5 | ~34 | S100 family (Ca/Cu-binding, EF-hand) |
| MT-CO2 ↔ SCO2 | 31.2 | SCO2 delivers Cu to COX2 in mitochondria |
| SOD1 ↔ CCS | 29.1 | Cu chaperone → SOD1 |
| ATP7B ↔ COMMD1 | 28.6 | COMMD1 regulates ATP7B stability |
| ATP7B ↔ ATP7A | 27.6 | Paralogous Cu-ATPases |

This list is strikingly consistent with the Cu-handling pathway described in Blockhuys 2017 (Fig. 1 legend + Introduction). The model is not just learning a discriminative signal — it is learning **through biologically correct edges**. The SOD3↔ATOX1 edge in particular reflects the very interaction the paper proposes (ATOX1 as transcription factor for SOD3).

## 6. Does the graph add interpretability beyond tabular models?

**Yes — the interpretability is the main win, not the accuracy.**

- Tabular models (logreg / RF / SVM) give a per-feature coefficient / importance, but no structure. Their "top features" list is CP, ALB, DBH, MAP2K1 — just the largest-log2FC genes.
- The GNN's saliency list adds ATOX1 and ATP7B at high ranks despite modest log2FC. Inspecting the attention graph shows these genes sit at the centre of highly-attended Cu-secretion subgraphs — information the tabular model has no way to express.
- The GAT attention map reproduces curated Cu-handling edges (CCS-SLC31A1, SOD3-ATOX1, COMMD1-ATP7B, MT-CO2-SCO2). This is a **sanity check**: the model trusts biology, not just statistics.
- The top-ranked modules (ECM LOX family / mitochondrial COX / hepatocyte secretome) translate directly into three concrete biological hypotheses for LIHC follow-up.

For a 54-gene tumor-vs-normal problem, the GNN does not beat SVM on raw AUC (0.995 vs 0.998). But it provides a **biologically anchored explanation** that SVM cannot.

## Files in this step
- `model_metrics.csv`, `model_metrics_summary.json` — per-fold + aggregate metrics
- `graph_embedding_umap.png` — UMAP of graph embeddings, tumor vs normal
- `node_importance.csv` — per-gene saliency (all 54 genes)
- `top_attention_edges.csv` — top-50 GAT attention edges
- `top_subgraph_or_attention_summary.md` — short human-readable interpretability report
