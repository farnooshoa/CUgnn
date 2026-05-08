"""Replace curated Cu edges with STRING v12 Cu-subnetwork.

Pulls the STRING network for the 54 Cu genes (homo sapiens, species 9606),
filters to combined_score >= 700 (high-confidence), then reruns the
tumor-vs-normal and stage-classification experiments on the STRING edge set
for comparison with the curated functional graph.

Outputs:
  data/string_v12_copper_edges.tsv
  outputs/final_comparison/string_v12.md
  outputs/final_comparison/string_v12_metrics.csv
"""
from __future__ import annotations
import sys
import io
import time
from pathlib import Path
import numpy as np
import pandas as pd
import requests
import networkx as nx

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.preprocessing import load_lihc_dataset
from src.graph_building import build_ppi_graph
from src.gnn_models.evaluate import run_gnn_grouped_cv
from src.gnn_models.train import TrainConfig
from src.utils import RANDOM_SEED

OUT = ROOT / "outputs" / "final_comparison"


def fetch_string_network(genes: list[str], min_score: int = 700) -> pd.DataFrame:
    """Return a DataFrame with columns source, target, score (0-1 from combined_score/1000)."""
    url = "https://string-db.org/api/tsv/network"
    params = {
        "identifiers": "\r".join(genes),
        "species": 9606,
        "caller_identity": "CUgnn-pilot",
        "required_score": min_score,  # 700 = high-confidence
    }
    r = requests.post(url, data=params, timeout=120)
    r.raise_for_status()
    df = pd.read_csv(io.StringIO(r.text), sep="\t")
    print(f"[string] API returned {len(df)} interactions at min_score={min_score}")
    if df.empty:
        return df
    df["source"] = df["preferredName_A"].str.upper()
    df["target"] = df["preferredName_B"].str.upper()
    df["score"] = df["score"].astype(float)
    return df[["source", "target", "score"]]


def build_string_graph(genes: list[str], edges: pd.DataFrame) -> nx.Graph:
    gene_set = set(genes)
    G = nx.Graph()
    G.add_nodes_from(genes)
    for row in edges.itertuples(index=False):
        s, t = row.source, row.target
        if s in gene_set and t in gene_set and s != t:
            if G.has_edge(s, t):
                G[s][t]["weight"] = max(G[s][t]["weight"], float(row.score))
            else:
                G.add_edge(s, t, weight=float(row.score), edge_type="physical")
    return G


def classify_stage(s):
    if not isinstance(s, str): return None
    s = s.upper().replace("STAGE ", "").strip()
    if s in {"I", "II", "IIA", "IIB"}: return 0
    if s.startswith("III") or s.startswith("IV"): return 1
    return None


def main():
    ds = load_lihc_dataset(require_real=True)
    genes = ds.expression.index.tolist()

    print(f"[string] fetching STRING v12 network for {len(genes)} Cu genes")
    edges = fetch_string_network(genes, min_score=700)
    if edges.empty:
        raise RuntimeError("STRING returned no edges at min_score=700")
    edges.to_csv(ROOT / "data" / "string_v12_copper_edges.tsv", sep="\t", index=False)

    G_string = build_string_graph(genes, edges)
    print(f"[string] graph: {G_string.number_of_nodes()} nodes, "
          f"{G_string.number_of_edges()} edges (vs 60 in curated functional)")
    isolated = [n for n, d in G_string.degree() if d == 0]
    print(f"[string] isolated nodes: {len(isolated)}")

    y_tn = (ds.metadata.loc[ds.expression.columns, "sample_type"]
            .str.lower().eq("tumor").astype(int).to_numpy())
    groups_all = ds.metadata.loc[ds.expression.columns, "case_submitter_id"].to_numpy()

    rows = []

    # === Tumor vs Normal on STRING graph ===
    print("\n[string] tumor-vs-normal with STRING graph (GAT)")
    cfg = TrainConfig(model="gat", epochs=60)
    sum_gat, _ = run_gnn_grouped_cv(ds, G_string, y_tn, groups_all, cfg, n_splits=5)
    r = sum_gat.iloc[0].to_dict()
    rows.append({
        "task": "tumor_vs_normal", "graph": "string_v12_hc", "model": "gat",
        "n_edges": G_string.number_of_edges(),
        "roc_auc": r["roc_auc_mean"], "roc_auc_std": r["roc_auc_std"],
        "balanced_accuracy": r["balanced_accuracy_mean"], "f1": r["f1_mean"],
    })

    # === Stage classification on STRING graph ===
    md = ds.metadata.loc[ds.expression.columns]
    stage_bin = md["stage"].map(classify_stage)
    keep = (md["sample_type"] == "Tumor") & stage_bin.notna()
    idx = np.where(keep.to_numpy())[0]
    y_stage = stage_bin[keep].astype(int).to_numpy()
    groups_stage = md.loc[keep, "case_submitter_id"].to_numpy()

    print(f"\n[string] stage classification with STRING graph (GAT), n={len(y_stage)}")
    cfg2 = TrainConfig(model="gat", epochs=80)
    sum_gat_s, _ = run_gnn_grouped_cv(ds, G_string, y_stage, groups_stage, cfg2,
                                       n_splits=5, sample_indices=idx)
    r = sum_gat_s.iloc[0].to_dict()
    rows.append({
        "task": "stage_early_vs_late", "graph": "string_v12_hc", "model": "gat",
        "n_edges": G_string.number_of_edges(),
        "roc_auc": r["roc_auc_mean"], "roc_auc_std": r["roc_auc_std"],
        "balanced_accuracy": r["balanced_accuracy_mean"], "f1": r["f1_mean"],
    })

    df = pd.DataFrame(rows).round(4)
    df.to_csv(OUT / "string_v12_metrics.csv", index=False)
    print("\n" + df.to_string(index=False))

    # Comparison: curated-functional numbers from previous runs
    tn_curated_auc = 0.993           # from leakage audit gnn_groupkfold
    tn_curated_balacc = 0.913
    stage_curated_auc = 0.668        # from stage_task
    stage_curated_balacc = 0.525

    string_tn = df[df["task"] == "tumor_vs_normal"].iloc[0]
    string_stage = df[df["task"] == "stage_early_vs_late"].iloc[0]

    (OUT / "string_v12.md").write_text(f"""# STRING v12 Edges — Replacing the Curated Cu-Interaction Fallback

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
| edges | 60 | **{G_string.number_of_edges()}** |
| isolated nodes | 3 | {len(isolated)} |
| edge weights | binary / 0.5 | continuous 0.70–1.00 (STRING combined_score / 1000) |
| provenance | manual + paper-derived | STRING v12 multi-evidence aggregation |

## Head-to-head performance (GAT, 5-fold StratifiedGroupKFold)

### Task 1 — Tumor vs Normal (n=424)

| graph | ROC-AUC | balanced acc |
|---|---:|---:|
| curated functional | {tn_curated_auc:.3f} | {tn_curated_balacc:.3f} |
| **STRING v12 HC** | **{string_tn['roc_auc']:.3f}** | **{string_tn['balanced_accuracy']:.3f}** |
| Δ (STRING − curated) | {string_tn['roc_auc'] - tn_curated_auc:+.3f} | {string_tn['balanced_accuracy'] - tn_curated_balacc:+.3f} |

### Task 2 — Stage I/II vs III/IV (n=349 tumors)

| graph | ROC-AUC | balanced acc |
|---|---:|---:|
| curated functional | {stage_curated_auc:.3f} | {stage_curated_balacc:.3f} |
| **STRING v12 HC** | **{string_stage['roc_auc']:.3f}** | **{string_stage['balanced_accuracy']:.3f}** |
| Δ (STRING − curated) | {string_stage['roc_auc'] - stage_curated_auc:+.3f} | {string_stage['balanced_accuracy'] - stage_curated_balacc:+.3f} |

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
edges = pd.read_csv("data/string_v12_copper_edges.tsv", sep="\\t")
G = build_ppi_graph(genes, external_edges=edges)   # uses score column
```

`build_ppi_graph()` already accepts a pre-filtered `external_edges` argument
(columns: `source`, `target`, `score`) — no code changes are needed beyond
passing the STRING TSV.

## Files produced
- `data/string_v12_copper_edges.tsv` — raw STRING download
- `outputs/final_comparison/string_v12.md` — this document
- `outputs/final_comparison/string_v12_metrics.csv` — per-task summary
""")
    print(f"\n[string] wrote {OUT/'string_v12.md'}")


if __name__ == "__main__":
    main()
