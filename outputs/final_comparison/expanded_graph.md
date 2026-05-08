# Expanded Graph — 54 Cu Genes + 1-Hop STRING Partners

## Rationale
The 54-node Cu graph may be too narrow for the GNN to show its full
advantage. Expanding with first-order STRING partners adds functional
context (co-pathway genes that are not themselves Cu-binders) and gives
message passing more structure to exploit.

## Construction
- **Seeds**: 54 Cu proteome genes
- **Query**: STRING v12 `network` endpoint with `add_nodes=100`, `min_score=700`
- **Final graph**: 153 nodes, 820 edges
- **Partner genes** (non-Cu): 99
- **Expression re-aggregated** from the raw STAR-Counts TSVs in `data/gdc_raw/`
  (same 424 TCGA-LIHC samples, same log2(FPKM-UQ+1) scale).

## Results (GAT, 5-fold StratifiedGroupKFold)

| task | 54-node graph | **expanded graph** (153 nodes) | Δ |
|---|---:|---:|---:|
| Tumor vs Normal        | 0.993 | **0.989** | -0.004 |
| Stage I/II vs III/IV   | 0.668 | **0.684** | +0.016 |
| 3-year Overall Survival | 0.692 | **0.695** | +0.003 |

## Interpretation

- **Tumor vs Normal**: already saturated at 54 nodes; expansion should
  not help and the Δ reflects noise only.
- **Stage and Survival**: these are the interesting tasks. A positive
  Δ ≥ 0.02 would confirm that graph scale matters for harder clinical
  endpoints. A flat or negative Δ would suggest that Cu biology is a
  complete enough feature space for these tasks and adding
  non-Cu partners adds noise.
- Either outcome is informative — it tells us whether the right
  inductive bias for Cu-biology problems is "narrow and curated"
  (54 nodes) or "wide and connected" (153 nodes).

## Caveats
- Expression for partner genes is on the same platform and distribution
  as Cu genes (same TCGA-LIHC cohort); no additional batch effects.
- STRING's `add_nodes=100` is a convenience parameter — for a final
  analysis it would be cleaner to compute the 1-hop closure explicitly
  and enforce a fixed node budget.
- The graph includes **attention-level cross-talk** to non-Cu partners,
  but the per-node features (expression only) are identical in kind to
  the Cu genes — no partner-specific prior is injected.

## Files produced
- `data/string_v12_expanded_edges.tsv` — expanded STRING edges
- `data/lihc_expression_expanded.tsv` — expanded expression matrix
- `outputs/final_comparison/expanded_graph.md` — this document
- `outputs/final_comparison/expanded_graph_metrics.csv` — per-task comparison
