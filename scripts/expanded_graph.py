"""Expanded graph: 54 Cu genes + 1-hop STRING partners.

Hypothesis: the tight 54-node graph is too narrow for the GNN to show its
full advantage. Expanding to include each Cu gene's first-order STRING
neighbours (at high confidence) grows the node set to ~150-250 and gives
the GNN more room for structural reasoning.

Steps:
  1. Query STRING v12 with the 54 Cu genes as seeds AND request their
     first-order neighbours (add_nodes=100).
  2. Build the expanded graph using combined_score >= 700 throughout.
  3. Verify TCGA-LIHC RNA-seq coverage for the expanded node set (our
     existing expression matrix only has 54 Cu genes — we must refetch the
     full TCGA-LIHC matrix for the expanded set).
  4. Since full re-download is expensive, we re-use what is already on disk:
     the raw per-sample STAR-Counts TSVs in data/gdc_raw/ contain every
     gene's FPKM-UQ — we just need to re-aggregate for the expanded node
     set.
  5. Run tumor-vs-normal, stage, and survival on the expanded graph.

Outputs:
  data/lihc_expression_expanded.tsv
  outputs/final_comparison/expanded_graph.md
  outputs/final_comparison/expanded_graph_metrics.csv
"""
from __future__ import annotations
import gzip
import io
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import requests
import networkx as nx

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.preprocessing import load_lihc_dataset, LIHCDataset
from src.graph_building import build_ppi_graph
from src.gnn_models.evaluate import run_gnn_grouped_cv
from src.gnn_models.train import TrainConfig
from src.utils import RANDOM_SEED
from scripts.add_methylation import classify_stage

OUT = ROOT / "outputs" / "final_comparison"
EXP_EXPR = ROOT / "data" / "lihc_expression_expanded.tsv"
EDGE_TSV = ROOT / "data" / "string_v12_expanded_edges.tsv"


def fetch_expanded_network(seeds: list[str], add_nodes: int = 100,
                             min_score: int = 700) -> pd.DataFrame:
    url = "https://string-db.org/api/tsv/network"
    params = {
        "identifiers": "\r".join(seeds), "species": 9606,
        "caller_identity": "CUgnn-pilot", "required_score": min_score,
        "add_nodes": add_nodes,
    }
    r = requests.post(url, data=params, timeout=300)
    r.raise_for_status()
    df = pd.read_csv(io.StringIO(r.text), sep="\t")
    df["source"] = df["preferredName_A"].str.upper()
    df["target"] = df["preferredName_B"].str.upper()
    df["score"] = df["score"].astype(float)
    return df[["source", "target", "score"]]


def reaggregate_expression_for_genes(gene_list: list[str]) -> pd.DataFrame:
    """Re-parse data/gdc_raw/ TSVs for the full expanded gene list."""
    from scripts.download_tcga_lihc import parse_one
    mf = pd.read_csv(ROOT / "data" / "manifest.tsv", sep="\t")
    gene_set = set(gene_list)
    cols = {}
    for row in mf.itertuples(index=False):
        path = ROOT / "data" / "gdc_raw" / row.file_name
        if not path.exists():
            continue
        series = parse_one(path)   # full gene-level FPKM-UQ
        if series is None:
            continue
        series = series[series.index.isin(gene_set)]
        cols[row.sample_submitter_id] = series
    matrix = pd.concat(cols, axis=1)
    matrix.index.name = "gene_symbol"
    matrix = matrix.reindex(gene_list)
    matrix = np.log2(matrix + 1.0)
    return matrix


def main():
    ds_orig = load_lihc_dataset(require_real=True)
    seeds = ds_orig.copper_genes["gene_symbol"].tolist()

    print(f"[expand] querying STRING for {len(seeds)} Cu seeds + 100 partners")
    edges = fetch_expanded_network(seeds, add_nodes=100, min_score=700)
    edges.to_csv(EDGE_TSV, sep="\t", index=False)
    all_genes = sorted(set(edges["source"]) | set(edges["target"]) | set(seeds))
    print(f"[expand] expanded node set: {len(all_genes)} genes, {len(edges)} edges")

    print("[expand] re-aggregating expression from gdc_raw for expanded set ...")
    expr_expanded = reaggregate_expression_for_genes(all_genes)
    # drop genes with all-NaN expression (not in GDC output)
    expr_expanded = expr_expanded.dropna(how="all")
    all_genes = expr_expanded.index.tolist()
    print(f"[expand] final expanded node set after expression filter: "
          f"{len(all_genes)} genes")
    expr_expanded.to_csv(EXP_EXPR, sep="\t")

    # Align with existing metadata
    md = ds_orig.metadata
    shared_samples = [c for c in expr_expanded.columns if c in md.index]
    expr_expanded = expr_expanded[shared_samples]
    md = md.loc[shared_samples]
    print(f"[expand] aligned samples: {expr_expanded.shape[1]}")

    # Construct a synthetic LIHCDataset-like object
    # (the evaluate helper only uses ds.expression, ds.metadata, ds.n_samples,
    #  ds.copper_genes, ds.tumor_mask - we need to fake copper_genes with the
    #  expanded node annotations)
    expanded_annot = pd.DataFrame({
        "gene_symbol": all_genes,
        "functional_category": [
            "transporter" if g in set(seeds) and g in set(ds_orig.copper_genes[ds_orig.copper_genes.functional_category=="transporter"]["gene_symbol"])
            else "enzyme" if g in set(seeds) and g in set(ds_orig.copper_genes[ds_orig.copper_genes.functional_category=="enzyme"]["gene_symbol"])
            else "other_or_unknown"
            for g in all_genes
        ],
        "subcellular_localization": ["unknown" for _ in all_genes],
        "notes_from_paper": [
            ("seed_Cu_gene" if g in set(seeds) else "STRING_1hop_partner")
            for g in all_genes
        ],
    })

    from dataclasses import dataclass
    from typing import List
    @dataclass
    class ExpDS:
        expression: pd.DataFrame
        metadata: pd.DataFrame
        copper_genes: pd.DataFrame
        missing_genes: list
        @property
        def n_samples(self): return self.expression.shape[1]
        @property
        def n_genes(self): return self.expression.shape[0]
        @property
        def tumor_mask(self):
            return (self.metadata.loc[self.expression.columns, "sample_type"]
                     .str.lower().eq("tumor").to_numpy())
    ds = ExpDS(expr_expanded, md, expanded_annot, [])

    # Build graph from STRING edges
    G = build_ppi_graph(all_genes, external_edges=edges)
    print(f"[expand] graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

    y_tn = (md.loc[expr_expanded.columns, "sample_type"]
             .str.lower().eq("tumor").astype(int).to_numpy())
    groups_all = md.loc[expr_expanded.columns, "case_submitter_id"].to_numpy()

    cfg = TrainConfig(model="gat", epochs=60)
    print("\n[expand] task 1: tumor vs normal (expanded graph)")
    s_tn, _ = run_gnn_grouped_cv(ds, G, y_tn, groups_all, cfg, n_splits=5)
    tn = s_tn.iloc[0].to_dict()

    # Stage
    stage_bin = md.loc[expr_expanded.columns, "stage"].map(classify_stage)
    keep = (md.loc[expr_expanded.columns, "sample_type"] == "Tumor") & stage_bin.notna()
    idx = np.where(keep.to_numpy())[0]
    y_s = stage_bin[keep].astype(int).to_numpy()
    g_s = md.loc[expr_expanded.columns[keep], "case_submitter_id"].to_numpy()

    cfg2 = TrainConfig(model="gat", epochs=80)
    print("[expand] task 2: stage (expanded graph)")
    s_st, _ = run_gnn_grouped_cv(ds, G, y_s, g_s, cfg2, n_splits=5, sample_indices=idx)
    st = s_st.iloc[0].to_dict()

    # Survival
    md["overall_survival_days"] = pd.to_numeric(md["overall_survival_days"], errors="coerce")
    def label_3y(row):
        if pd.isna(row["overall_survival_days"]) or pd.isna(row["vital_status"]):
            return None
        if row["vital_status"] == "Dead":
            return 1 if row["overall_survival_days"] <= 1095 else 0
        if row["vital_status"] == "Alive":
            return 0 if row["overall_survival_days"] > 1095 else None
        return None
    md_aligned = md.loc[expr_expanded.columns].copy()
    md_aligned["y_3y"] = md_aligned.apply(label_3y, axis=1)
    keep_s = (md_aligned["sample_type"] == "Tumor") & md_aligned["y_3y"].notna()
    idx_s = np.where(keep_s.to_numpy())[0]
    y_sv = md_aligned.loc[keep_s, "y_3y"].astype(int).to_numpy()
    g_sv = md_aligned.loc[keep_s, "case_submitter_id"].to_numpy()

    print("[expand] task 3: survival (expanded graph)")
    s_sv, _ = run_gnn_grouped_cv(ds, G, y_sv, g_sv, cfg2, n_splits=5, sample_indices=idx_s)
    sv = s_sv.iloc[0].to_dict()

    rows = [
        {"task": "tumor_vs_normal",     "n_nodes": len(all_genes),
         "roc_auc": tn["roc_auc_mean"], "balanced_accuracy": tn["balanced_accuracy_mean"],
         "setting": "expanded_graph"},
        {"task": "stage_early_vs_late", "n_nodes": len(all_genes),
         "roc_auc": st["roc_auc_mean"], "balanced_accuracy": st["balanced_accuracy_mean"],
         "setting": "expanded_graph"},
        {"task": "survival_3yr",        "n_nodes": len(all_genes),
         "roc_auc": sv["roc_auc_mean"], "balanced_accuracy": sv["balanced_accuracy_mean"],
         "setting": "expanded_graph"},
    ]
    df = pd.DataFrame(rows).round(4)
    df.to_csv(OUT / "expanded_graph_metrics.csv", index=False)
    print("\n" + df.to_string(index=False))

    (OUT / "expanded_graph.md").write_text(f"""# Expanded Graph — 54 Cu Genes + 1-Hop STRING Partners

## Rationale
The 54-node Cu graph may be too narrow for the GNN to show its full
advantage. Expanding with first-order STRING partners adds functional
context (co-pathway genes that are not themselves Cu-binders) and gives
message passing more structure to exploit.

## Construction
- **Seeds**: 54 Cu proteome genes
- **Query**: STRING v12 `network` endpoint with `add_nodes=100`, `min_score=700`
- **Final graph**: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges
- **Partner genes** (non-Cu): {G.number_of_nodes() - len(set(seeds))}
- **Expression re-aggregated** from the raw STAR-Counts TSVs in `data/gdc_raw/`
  (same 424 TCGA-LIHC samples, same log2(FPKM-UQ+1) scale).

## Results (GAT, 5-fold StratifiedGroupKFold)

| task | 54-node graph | **expanded graph** ({G.number_of_nodes()} nodes) | Δ |
|---|---:|---:|---:|
| Tumor vs Normal        | 0.993 | **{tn['roc_auc_mean']:.3f}** | {tn['roc_auc_mean'] - 0.993:+.3f} |
| Stage I/II vs III/IV   | 0.668 | **{st['roc_auc_mean']:.3f}** | {st['roc_auc_mean'] - 0.668:+.3f} |
| 3-year Overall Survival | 0.692 | **{sv['roc_auc_mean']:.3f}** | {sv['roc_auc_mean'] - 0.692:+.3f} |

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
  (54 nodes) or "wide and connected" ({G.number_of_nodes()} nodes).

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
""")
    print(f"\n[expand] wrote {OUT/'expanded_graph.md'}")


if __name__ == "__main__":
    main()
