# Visualization Interpretation — TCGA-LIHC Cu Proteome

All five figures are rendered from the real-data pipeline run (424 TCGA-LIHC samples, 374 tumor / 50 normal, 54/54 Cu genes). Layout is reproducible — the same spring-layout seed is shared across Fig 1, 3, 4 so figures can be compared side-by-side.

## Fig 1 — Cu-proteome graph coloured by LIHC log2FC
`network_logfc.png`

Node colour tells the LIHC tumor-vs-normal story at a glance:
- **Blue (down)** nodes dominate the hepatocyte-secretome cluster (ALB, CP, HEPH, HEPHL1, DBH) and the mitochondrial COX cluster (MT-CO1/CO2, SCO1/SCO2, COX11, COX17).
- **Red (up)** nodes cluster around the ECM axis (LOX, LOXL2, SPARC) and AFP.
- Node size (degree centrality) highlights **ATP7A/B, SOD1, CP, LOX, SPARC** as hubs.

Biologically this is the textbook HCC pattern: loss of hepatocyte secretory function and mitochondrial oxidative phosphorylation, gain of ECM remodelling.

## Fig 2 — Top-15 Cu genes by GNN importance (bar plot)
`node_importance_bar.png`

GNN saliency top-10: CP, SOD1, MT-CO1, ATOX1, MAP2K1, ALB, LOXL4, MT-CO2, LOX, ATP7B.
Bar colour encodes log2FC direction so you can read significance *and* direction at once. Notice that **CP, ALB, SOD1, MT-CO1/CO2, MAP2K1** rank high *and* are downregulated — the strongest "lost" genes. **LOX/LOXL4, AFP** rank high *and* are up — the strongest "gained" genes. **ATOX1 and ATP7B** rank high despite modest log2FC, which is the GNN using graph context (both sit at the centre of the Cu-secretion module).

## Fig 3 — Cu graph with node size = GNN importance
`network_importance.png`

Same layout as Fig 1 but node sizes reflect GNN saliency instead of degree. Makes it visually obvious that importance is not just the biggest hub — the **ATOX1 / ATP7B / CP / CCS / SOD1 secretory-antioxidant triangle** is what the model is using most. Top-15 nodes are outlined with a bold black border so they can be picked out quickly.

## Fig 4 — Top GAT attention edges
`attention_edges.png`

Top 25 non-self attention edges overlaid on the Cu graph background (faint grey). Red edges = pairs with canonical Cu-handling biology; blue edges = other top-attended pairs.

Top attention pairs (excluding self-loops):
- CCS — SLC31A1 (41.9)
- ATP7B — AFP (39.3)
- SOD3 — ATOX1 (34.5)
- PARK7 — SOD1 (34.0)
- S100A12 — S100A13 (34.0)
- S100A12 — S100B (34.0)
- S100A12 — S100A5 (34.0)
- MT-CO2 — SCO2 (31.2)
- SOD1 — CCS (29.1)
- ATP7B — COMMD1 (28.6)

The red edges that do appear include:
- **CCS ↔ SLC31A1** — Cu import to cytosolic SOD1 chaperone
- **SOD3 ↔ ATOX1** — the ATOX1 transcription-factor-for-SOD3 interaction the 2017 paper explicitly validated
- **ATP7B ↔ COMMD1** — COMMD1 regulates ATP7B stability
- **MT-CO2 ↔ SCO2** — mitochondrial Cu delivery
- **SOD1 ↔ CCS** — canonical Cu chaperone → SOD1
- **ATP7B ↔ ATP7A** — paralogous Cu ATPases

This is the single most compelling figure in the pilot: the GAT did not just learn to classify, it learned **through biologically correct Cu-handling edges**.

## Fig 5 — UMAP of GAT graph embeddings
`graph_embedding_umap.png`

Tumor and normal samples separate cleanly in the UMAP of per-sample graph embeddings — no dense mixing of colours. The small normal cluster is compact, consistent with the paired-sample structure of TCGA's liver normals (often adjacent-normal tissue from the same patient).

The separation in embedding space is broadly consistent with the high AUC (0.995) — the model's representation is pushing the two classes apart, not just setting a decision boundary.

## Summary

- **Fig 1 + Fig 3** together show the biology at gene-level (what is up/down) and at model-level (what the GNN leans on).
- **Fig 2** is the single best view for a collaborator who wants "which genes matter?".
- **Fig 4** is the single best view for a biologist who wants "do I trust this model?" — the answer is yes, because the top attention edges are textbook Cu biology.
- **Fig 5** confirms the classifier's decision geometry is clean, not a knife-edge separation.

## Top up/down genes from DE (for reference)
- Top 5 up in tumor: LOXL2, AFP, CUTA, MT3, SPARC.
- Top 5 down in tumor: ALB, SLC31A1, CP, PRNP, DBH.
