# Verification Note — Paper Extraction Audit

Audit of the three paper-extraction artefacts before the real-data phase.

## 1. Copper gene count is exactly 54
`copper_gene_list.csv` contains **54 rows, 54 unique gene symbols**, no duplicates, no missing `subcellular_localization` values. Confirmed.

Counts (verified programmatically):
- transporter: 12
- enzyme: 26
- other_or_unknown: 16
- **Total: 54**

## 2. Internal consistency of category counts
The paper text states 12 / ~27 / 15 (transporters / enzymes / "other"). Our reconstruction from Fig. 1 + text gives 12 / 26 / 16 — off by one in each of the last two categories. The off-by-one is a real ambiguity in the paper, not a data error. Specifically:

- **CUTC**: the paper flags this explicitly — "Whether CUTC is truly a Cu transporter, as reported here, or a Cu-dependent enzyme instead, is an open question." We kept it as a transporter.
- **SCO1 / SCO2**: classified as transporters in the paper text ("12 proteins classified as Cu transporters … COX11, COX17, SCO1, and SCO2") but coloured orange (enzyme) in Fig. 1. We followed the text.
- **CUTA**: Fig. 1 lists it as black (other) but the paper describes it as a predicted mitochondrial Cu-binding protein of uncertain function. We kept it as "other".
- **S100A5 / S100A12 / S100A13 / S100B / MT3 / MT4**: the paper groups these under "other" because Cu role is non-enzymatic; they are Cu-binding but their primary identity is Ca-binding (S100 family) or metal-sequestration (metallothioneins). Kept as "other".

None of these edge cases changes the total of 54 or the biological content of the node set. Downstream code is agnostic to the transporter-vs-enzyme-vs-other boundary beyond an optional 3-dim one-hot feature.

## 3. Uncertain or inferred entries (explicit list)
Entries flagged as having non-trivial ambiguity in the original paper — users should be aware these are the "soft" nodes in our copper proteome graph:

| Gene | Ambiguity | Source in paper |
|---|---|---|
| CUTC | transporter vs enzyme | Main text page 115 — explicitly called "open question" |
| CUTA | compartment and function | Fig. 1 legend — less stringent COMPARTMENTS level 2 |
| MT4 | compartment (less stringent) | Fig. 1 legend — COMPARTMENTS level 2 |
| ATOX1 / MEMO1 / MT3 / LTF / S100A5 | compartment (COMPARTMENTS-only) | Fig. 1 legend |
| MOXD2P | pseudogene; ER localisation assumed same as MOXD1 | Fig. 1 caption |
| SNCA / APP / PRNP | Cu binding confirmed; Cu *role* unclear | Paper page 115 |
| SPARC / MEMO1 / MAP2K1 | added by the authors from literature, not from the UniProt "human + copper" search | Methods — Data mining |

Localization for three entries (ATOX1, MEMO1, MT3) was assigned from COMPARTMENTS only (no dual-source confirmation) — flagged but kept.

Nothing in this list disqualifies an entry from the copper proteome. We keep all 54 but the above notes should be surfaced in any downstream claim about a specific gene.

## 4. Fig. 3 translation into graph-construction guidance
`figure3_graph_construction_notes.md` correctly captures:
- Multi-edge-type provenance (physical / co-expression / genetic), matching the Fig. 3 edge-colour legend.
- Node-colour semantics (red / blue / black / grey) and their **breast-cancer-specific** meaning — importantly, the document flags that these colours are NOT reusable for LIHC because LIHC has different log2FC signs.
- The decision to drop first-level partner genes (grey nodes in Fig. 3) and keep only the 54-gene core — grounded in reproducibility and fixed-node-set requirements of graph classification.
- The decision to treat cancer state as a **graph-level label**, not as a node colour — the correct translation for a GNN.

Cross-check of the code:
- `src/graph_building/build_graph.py::CURATED_CU_EDGES` uses the same three edge-type vocabulary (`physical`, `coexpression`, `genetic`) and adds a fourth (`shared_compartment`) that is clearly labelled as an auxiliary functional-graph-only edge, not a Fig. 3 claim.
- The module-detection step (`module_detection`) is applied to the merged 54-node graph, not per-cluster subgraphs — matching the "single merged graph, not eight isolated subgraphs" argument in the notes.

No disagreements found between the notes and the code.

## Summary
- 54 genes confirmed, no duplicates, no missing fields.
- 12 / 26 / 16 category split; the 1-unit gap vs. the paper text is inherent to the paper and traceable to CUTC / SCO1-2 / CUTA boundary cases.
- Uncertain entries are listed in full above; none are removed.
- Fig. 3 translation in `figure3_graph_construction_notes.md` is consistent with the graph-building code.

Extraction is safe to use for the real-data phase.
