# Figure 3 — Notes for Graph Construction

Fig. 3 of Blockhuys et al. (2017) is the most graph-like figure in the paper and is the natural inspiration for how we build a graph over the Cu proteome for LIHC. This note records what the figure actually is, what it is *not*, and how we translate it into our own graph.

---

## How Fig. 3 was generated — conceptually

1. Start from the 8 gene clusters in the Fig. 2 heat-map (Ward's-D2 on Pearson-correlation distance of tumour/normal log2FC across 18 cancers).
2. For each cluster, dump the cluster's gene list into **Cytoscape + GeneMANIA App**.
3. GeneMANIA queries its prior knowledge base and pulls in **first-level partner genes** (any gene known to interact with at least one seed gene). These appear as grey nodes.
4. Three evidence channels are activated at once:
   - **Physical (protein–protein) interactions** — red edges
   - **Co-expression** — purple edges
   - **Genetic interactions** — green edges
5. A weighting scheme is applied: node size is inversely proportional to gene rank in GeneMANIA's result list, so "central" partners in the evidence base draw as larger circles.
6. Seed nodes are then re-coloured by their log2FC direction **in breast cancer only**:
   - red = up in breast tumour
   - blue = down in breast tumour
   - black = no change in breast tumour
   - grey = not a seed (GeneMANIA first-level partner)

In other words, the network topology is **prior biological knowledge** (GeneMANIA) and the node annotation is **breast-cancer expression direction**.

---

## Interaction types used (and what they mean)

| Edge type | GeneMANIA source | What it tells you |
|---|---|---|
| Physical | literature-curated PPI (BioGRID, IntAct, MINT, etc.) | Direct protein–protein interaction |
| Co-expression | pairwise gene correlations in public expression datasets | Genes regulated in coordinated fashion across conditions |
| Genetic | synthetic-lethal / epistasis / genetic-interaction screens | Functional dependence without requiring physical contact |

Fig. 3 does **not** include pathway-membership edges or literature co-mention edges.

---

## Why Fig. 3 is useful for our LIHC graph

- It is proof that domain experts already think of the Cu proteome as a **graph of genes with multi-evidence edges**. Our GNN input format (adjacency matrix over gene nodes with evidence-weighted edges) is therefore biologically honest, not arbitrary.
- It tells us which **edge channels** are worth building: physical, co-expression, genetic. Any modern curated source (STRING, BioGRID, GeneMANIA, HumanNet) exposes the same three broadly.
- It reveals empirically that **clusters 1, 3, 8 are dense; clusters 2, 4 are sparse**, and that seed genes are highly re-used as first-level partners across clusters. This means the right structure for a GNN is **one merged 54-node graph**, not eight disconnected subgraphs.
- It defines the **interpretability bar** we should hit: our trained GNN's important-node/edge output should be comparable with (not identical to) the Fig. 3 clusters.

---

## Why Fig. 3 should NOT be copied as-is into the GNN

1. **Node set is inflated**: Fig. 3 includes first-level partner genes (grey) that are not in the 54-gene Cu proteome. Our graph sticks to the 54 curated genes so that node identity is fixed across samples and interpretation stays tied to Cu biology. Adding arbitrary first-level partners balloons the node count without controllable selection criteria.
2. **Node annotation is breast-specific**: the red/blue/black coloring in Fig. 3 encodes a breast-cancer expression direction. For LIHC we must recompute log2FC and p-values from TCGA-LIHC — we cannot simply reuse the breast-derived colors.
3. **Edges are a snapshot**: GeneMANIA v2017 evidence is different from today's STRING v12. We should rebuild edges from a current source and record a version string.
4. **Edge weighting is qualitative**: Fig. 3 does not expose numerical edge weights. A GCN/GAT can use binary adjacency, but if we want weighted graphs we need to go to the underlying database, not to Fig. 3.
5. **Per-cluster subgraphs are a visualisation choice, not a model choice**: the eight subgraphs are plotted separately for readability. The underlying data is one connected interaction graph.
6. **Cancer context is not encoded in the graph itself**: in Fig. 3 cancer state is only a node *color*, not an edge or a separate graph. In our setup cancer state is the **graph-level label** (tumour vs normal), which is a different learning problem from simply staring at the colors.

---

## Translation into our LIHC graph (summary)

- **Nodes**: the 54 genes from `copper_gene_list.csv`. No first-level partners.
- **Edges**: STRING v12 "combined_score" restricted to the 54-gene subnetwork as `ppi_graph`; optional `functional_graph` from BioGRID or GeneMANIA; `coexpression_graph` computed from TCGA-LIHC expression as a secondary, exploratory variant only.
- **Edge types**: stored as a discrete `edge_type` attribute (`physical`, `coexpression`, `genetic`) so we can later run a relational GNN if we want.
- **Node features (per patient graph)**: TCGA-LIHC normalised expression + optional z-score + functional-category one-hot. Graph-level label is tumour vs normal.
- **Interpretability target**: a module / attention / GNNExplainer report that we can lay next to Fig. 3 and ask "do the same Cu clusters come up in LIHC as in the paper's cross-cancer analysis?"
