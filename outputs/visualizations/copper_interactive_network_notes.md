# Interactive Copper Network — Legend & Notes

Open [`copper_interactive_network.html`](copper_interactive_network.html) in any modern browser. No server needed — it is a single self-contained HTML file that pulls Cytoscape.js from CDN.

## Layout

The initial layout is computed with **fCoSE** — a force-directed algorithm tuned for biological networks. After load you can drag any node to rearrange it, scroll to zoom, and click-drag the background to pan.

## Nodes — 54 Cu proteome genes

| visual | meaning |
|---|---|
| **size** | GNN saliency (`outputs/gnn/node_importance.csv`). Bigger = more influential on the tumor/normal decision. Top-10 genes have a thick black border. |
| **colour** | TCGA-LIHC log2FC from `copper_de_results.csv`. Red = up in tumor, blue = down, near-white = no change. Colour scale is clipped at ±3. |
| **shape** | functional category: circle = transporter, rounded square = enzyme, triangle = other/unknown. |
| **tooltip (hover)** | gene name, log2FC, adj p (BH), GNN importance, category, module id, subcellular localization, and the paper-derived notes field. |

## Edges

| visual | meaning |
|---|---|
| **colour** | interaction type — red = physical / PPI, purple = co-expression, green = genetic, grey = shared subcellular compartment, blue dashed = attention-only (GAT top edge that is not in the fixed topology) |
| **width** | GAT attention sum (`outputs/gnn/top_attention_edges.csv`). Thicker = higher attention. Edges with no attention data get default thin width. |
| **bold red** | canonical Cu-handling pair (CCS–SLC31A1, SOD3–ATOX1, ATP7B–COMMD1, MT-CO2–SCO2, …). Toggleable via the "Canonical Cu-biology edges" checkbox. |
| **tooltip (hover)** | source–target pair, edge type, GAT attention, prior weight, top-attention flag. |

## Controls

- **Search box** — type a gene symbol or prefix to fade everything else.
- **Fit view / Reset highlight** — self-explanatory.
- **Highlight toggles** — turn on/off the top-10 important nodes, top-15 attention edges, canonical Cu-biology edges.
- **Edge type toggles** — show/hide edges by interaction type (useful to see the attention-only blue edges on their own).
- **Click a node** — pins its tooltip and fades everything outside its 1-hop neighbourhood.
- **Click an edge** — pins its tooltip and fades everything except its two endpoints.
- **Click the background** — unpins.

## What to look for

1. The **blue cluster** in the lower-right (ALB, CP, DBH, MT-CO1/CO2, SOD1) is the hepatocyte-secretome + mitochondrial-COX axis, strongly downregulated.
2. The **red cluster** (LOX, LOXL2, SPARC, AFP) is the ECM-remodelling + HCC-marker axis, upregulated.
3. Switch OFF all edge types except "attention-only" to see **purely what the GAT learned** — the canonical Cu-handling pairs (CCS↔SLC31A1, SOD3↔ATOX1, ATP7B↔COMMD1, MT-CO2↔SCO2) should dominate.
4. Click **ATOX1** or **ATP7B** — both have modest log2FC but high GNN importance because their 1-hop neighbourhood is central to Cu homeostasis.
5. Toggle off the grey "shared compartment" edges to declutter the graph and focus on biology-derived interactions.

## Stats
- Nodes: 54 (of which 44 significant at BH<0.05)
- Edges: 60 (including 0 attention-only)
- Top-10 important genes: highlighted with thick black border
- Canonical Cu-biology edges: red glow when the "Canonical Cu-biology edges" toggle is on
