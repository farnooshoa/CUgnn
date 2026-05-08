# STRING v12 Edges — Replacing the Curated Cu-Interaction Fallback

## Source
- STRING v12 via REST API (`https://string-db.org/api/tsv/network`)
- Species: Homo sapiens (9606)
- Input: the 54 HGNC symbols from `copper_gene_list.csv`
- `required_score = 700` (**high-confidence**; STRING default "high" threshold)
- Downloaded edge list saved as `data/string_v12_copper_edges.tsv`

## Graph statistics

| metric | curated functional (fallback) | STRING v12 HC |
|---|---:|---:|
| nodes | 54 | 54 |
| edges | 60 | **77** |
| isolated nodes | 3 | 18 |
| edge weights | binary / 0.5 | continuous 0.70–1.00 (STRING combined_score / 1000) |
| provenance | manual + paper-derived | STRING v12 multi-evidence aggregation |

## Head-to-head performance (GAT, 5-fold StratifiedGroupKFold)

### Task 1 — Tumor vs Normal (n=424)

| graph | ROC-AUC | balanced acc |
|---|---:|---:|
| curated functional | 0.993 | 0.913 |
| **STRING v12 HC** | **0.984** | **0.909** |
| Δ (STRING − curated) | -0.009 | -0.004 |

### Task 2 — Stage I/II vs III/IV (n=349 tumors)

| graph | ROC-AUC | balanced acc |
|---|---:|---:|
| curated functional | 0.668 | 0.525 |
| **STRING v12 HC** | **0.690** | **0.579** |
| Δ (STRING − curated) | +0.022 | +0.054 |

## Interpretation

- On tumor-vs-normal the task is saturated — any sparse Cu graph gives AUC
  ≈ 0.99 (see `graph_ablation.md`). Moving to STRING changes neither the
  verdict nor the interpretability-only framing.
- On the harder stage task, STRING's denser + evidence-weighted edges may
  or may not help. Compare Δ carefully against the ~0.05 fold-to-fold std
  in `stage_task_metrics.csv`.
- Regardless of AUC, STRING is the more **defensible** choice for a
  real manuscript: pinned database version, documented scoring, reviewer
  immediately understands the provenance. The curated fallback was always
  a pilot-phase shortcut.

## Practical notes for downstream use

To switch the rest of the pipeline to the STRING graph:
```python
import pandas as pd
from src.graph_building import build_ppi_graph
from src.preprocessing import load_lihc_dataset

ds = load_lihc_dataset(require_real=True)
genes = ds.expression.index.tolist()
edges = pd.read_csv("data/string_v12_copper_edges.tsv", sep="\t")
G = build_ppi_graph(genes, external_edges=edges)   # uses score column
```

`build_ppi_graph()` already accepts a pre-filtered `external_edges` argument
(columns: `source`, `target`, `score`) — no code changes are needed beyond
passing the STRING TSV.

## Files produced
- `data/string_v12_copper_edges.tsv` — raw STRING download
- `outputs/final_comparison/string_v12.md` — this document
- `outputs/final_comparison/string_v12_metrics.csv` — per-task summary
