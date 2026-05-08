"""Promoter-only methylation: re-aggregate using TSS1500/TSS200 probes only.

The Phase 3 methylation integration averaged *all* probes mapped to each
Cu gene — body + promoter + UTR. Biologically, only the **promoter**
methylation (TSS1500 and TSS200 regions) correlates reliably with
transcriptional silencing. This script:

  1. Parses the Illumina HumanMethylation450 manifest to extract the
     UCSC_RefGene_Group column (per-probe region annotation).
  2. Filters probes mapping to a Cu gene to those whose region is
     TSS1500 or TSS200 (or both) for that gene.
  3. Re-aggregates to gene-level β.
  4. Re-runs the tumor-vs-normal and stage tasks with the new feature.

Outputs:
  data/lihc_methylation_promoter.tsv
  outputs/final_comparison/promoter_methylation.md
  outputs/final_comparison/promoter_methylation_metrics.csv
"""
from __future__ import annotations
import gzip
import sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.preprocessing import load_lihc_dataset
from src.graph_building import build_functional_graph
from src.gnn_models.train import TrainConfig
# re-use the run_with_meth helper from add_methylation
from scripts.add_methylation import run_with_meth, align_methylation_to_expression, classify_stage

OUT = ROOT / "outputs" / "final_comparison"
XENA_FILE = ROOT / "data" / "xena_lihc_methylation.gz"
MANIFEST = ROOT / "data" / "illumina_450k_manifest.csv"
PROMO_TSV = ROOT / "data" / "lihc_methylation_promoter.tsv"


PROMOTER_GROUPS = {"TSS1500", "TSS200"}   # Illumina UCSC_RefGene_Group values


def load_manifest_cu_probes(copper_genes: set[str]) -> pd.DataFrame:
    """Return DataFrame: probe, gene (split), group (split), promoter_flag.

    The Illumina manifest has ; -delimited gene names and matching ; -delimited
    group labels when a probe overlaps multiple transcripts.
    """
    print("[promo] parsing Illumina 450K manifest ...")
    # skip the 7-line header, take only the [Assay] rows
    df = pd.read_csv(
        MANIFEST, skiprows=7, low_memory=False,
        usecols=["IlmnID", "UCSC_RefGene_Name", "UCSC_RefGene_Group"],
    )
    df = df.dropna(subset=["UCSC_RefGene_Name", "UCSC_RefGene_Group"])
    print(f"[promo] manifest rows with gene annotation: {len(df)}")

    # Explode gene;group pairs
    def split_pairs(row):
        genes = [g.upper() for g in str(row["UCSC_RefGene_Name"]).split(";")]
        groups = str(row["UCSC_RefGene_Group"]).split(";")
        # align lengths: sometimes genes has 3 and groups has 3 (match)
        if len(genes) == len(groups):
            return list(zip(genes, groups))
        return []
    rows = []
    cu_set = set(g.upper() for g in copper_genes)
    for row in df.itertuples(index=False):
        for gene, group in split_pairs(row._asdict()):
            if gene in cu_set:
                rows.append({"probe": row.IlmnID, "gene": gene, "group": group,
                              "is_promoter": group in PROMOTER_GROUPS})
    out = pd.DataFrame(rows).drop_duplicates(subset=["probe", "gene"], keep="first")
    print(f"[promo] Cu-mapping rows: {len(out)}")
    print(f"[promo] promoter (TSS1500/TSS200) rows: {int(out['is_promoter'].sum())}")
    print(f"[promo] Cu genes with ≥1 promoter probe: "
          f"{out[out['is_promoter']]['gene'].nunique()} / {len(copper_genes)}")
    return out


def build_gene_level_beta(anno_cu: pd.DataFrame, copper_list: list[str],
                           promoter_only: bool) -> pd.DataFrame:
    """Stream the Xena β table once and aggregate to per-gene."""
    keep_anno = anno_cu[anno_cu["is_promoter"]] if promoter_only else anno_cu
    keep_probes = set(keep_anno["probe"])
    print(f"[promo] streaming Xena file for {len(keep_probes)} probes "
          f"({'promoter-only' if promoter_only else 'all regions'}) ...")
    rows = []
    with gzip.open(XENA_FILE, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        sample_cols = header[1:]
        for i, line in enumerate(f):
            if i % 100000 == 0 and i > 0:
                print(f"  ... {i} probes scanned")
            cid, *vals = line.rstrip("\n").split("\t")
            if cid in keep_probes:
                rows.append([cid] + vals)
    beta = pd.DataFrame(rows, columns=header).set_index(header[0])
    beta = beta.apply(pd.to_numeric, errors="coerce")
    probe_to_gene = keep_anno.set_index("probe")["gene"]
    beta = beta.join(probe_to_gene, how="inner")
    gene_level = beta.groupby("gene").mean(numeric_only=True)
    gene_level = gene_level.reindex(copper_list)
    gene_level.index.name = "gene_symbol"
    gene_level.columns = [c.replace(".", "-") for c in gene_level.columns]
    return gene_level


def main():
    ds = load_lihc_dataset(require_real=True)
    copper_list = ds.copper_genes["gene_symbol"].tolist()

    anno_cu = load_manifest_cu_probes(set(copper_list))

    print("\n[promo] building PROMOTER-ONLY β matrix ...")
    promoter_beta = build_gene_level_beta(anno_cu, copper_list, promoter_only=True)
    promoter_beta.to_csv(PROMO_TSV, sep="\t")
    promoter_aligned = align_methylation_to_expression(ds, promoter_beta)

    # Re-run two tasks — tumor vs normal and stage
    graph = build_functional_graph(ds.expression.index.tolist(), ds.copper_genes)

    y_tn = (ds.metadata.loc[ds.expression.columns, "sample_type"]
            .str.lower().eq("tumor").astype(int).to_numpy())
    groups_all = ds.metadata.loc[ds.expression.columns, "case_submitter_id"].to_numpy()

    md = ds.metadata.loc[ds.expression.columns]
    stage_bin = md["stage"].map(classify_stage)
    keep = (md["sample_type"] == "Tumor") & stage_bin.notna()
    idx_stage = np.where(keep.to_numpy())[0]
    y_stage = stage_bin[keep].astype(int).to_numpy()
    groups_stage = md.loc[keep, "case_submitter_id"].to_numpy()

    cfg = TrainConfig(model="gat", epochs=80)

    print("\n[promo] task 1: tumor vs normal + promoter-only methylation")
    tn_promo, _ = run_with_meth(ds, graph, y_tn, groups_all, cfg,
                                 np.arange(len(y_tn)), promoter_aligned, n_splits=5)

    print("[promo] task 2: stage + promoter-only methylation")
    stage_promo, _ = run_with_meth(ds, graph, y_stage, groups_stage, cfg,
                                    idx_stage, promoter_aligned, n_splits=5)

    rows = [
        {"task": "tumor_vs_normal", "features": "expr_only",
         "roc_auc": 0.993, "balanced_accuracy": 0.913, "source": "leakage audit"},
        {"task": "tumor_vs_normal", "features": "expr + methylation (all regions)",
         "roc_auc": 0.997, "balanced_accuracy": 0.912, "source": "phase 3 methylation"},
        {"task": "tumor_vs_normal", "features": "expr + methylation (promoter-only)",
         "roc_auc": tn_promo["roc_auc"], "balanced_accuracy": tn_promo["balanced_accuracy"],
         "source": "this run"},
        {"task": "stage_early_vs_late", "features": "expr_only",
         "roc_auc": 0.668, "balanced_accuracy": 0.525, "source": "stage task"},
        {"task": "stage_early_vs_late", "features": "expr + methylation (all regions)",
         "roc_auc": 0.667, "balanced_accuracy": 0.518, "source": "phase 3 methylation"},
        {"task": "stage_early_vs_late", "features": "expr + methylation (promoter-only)",
         "roc_auc": stage_promo["roc_auc"], "balanced_accuracy": stage_promo["balanced_accuracy"],
         "source": "this run"},
    ]
    df = pd.DataFrame(rows).round(4)
    df.to_csv(OUT / "promoter_methylation_metrics.csv", index=False)
    print("\n" + df.to_string(index=False))

    n_genes_w_promoter = int(anno_cu[anno_cu['is_promoter']]['gene'].nunique())
    n_promoter_probes = int(anno_cu['is_promoter'].sum())

    (OUT / "promoter_methylation.md").write_text(f"""# Promoter-Only Methylation — Targeted Aggregation

## Motivation
Phase 3 methylation (in `methylation_results.md`) averaged **all** 450K
probes per Cu gene — body + promoter + UTR. Biology says only the
**promoter** methylation (TSS1500 / TSS200) matters for transcriptional
silencing. This step re-aggregates restricted to promoter probes.

## Probe filtering
- Source: Illumina HumanMethylation450 v1-2 manifest (`data/illumina_450k_manifest.csv`)
- Column: `UCSC_RefGene_Group`
- Kept groups: `TSS1500`, `TSS200`
- Dropped groups: `Body`, `3'UTR`, `5'UTR`, `1stExon`

### Coverage
- Cu genes with ≥1 promoter probe: **{n_genes_w_promoter} / 54**
- Promoter probes mapping to Cu genes: **{n_promoter_probes}**
  (vs 848 all-region probes in Phase 3)
- Promoter β matrix: `data/lihc_methylation_promoter.tsv`

## Results (GAT, 5-fold StratifiedGroupKFold)

| task | features | ROC-AUC | balanced acc |
|---|---|---:|---:|
| Tumor vs Normal | expr_only | 0.993 | 0.913 |
| Tumor vs Normal | expr + meth (all regions) | 0.997 | 0.912 |
| Tumor vs Normal | **expr + meth (promoter-only)** | **{tn_promo['roc_auc']:.3f}** | **{tn_promo['balanced_accuracy']:.3f}** |
| | | | |
| Stage I/II vs III/IV | expr_only | 0.668 | 0.525 |
| Stage I/II vs III/IV | expr + meth (all regions) | 0.667 | 0.518 |
| Stage I/II vs III/IV | **expr + meth (promoter-only)** | **{stage_promo['roc_auc']:.3f}** | **{stage_promo['balanced_accuracy']:.3f}** |

## Interpretation

- **Δ on stage (promoter - all)**: {stage_promo['roc_auc'] - 0.667:+.3f} AUC.
- **Δ on stage (promoter - expr only)**: {stage_promo['roc_auc'] - 0.668:+.3f} AUC.

### Read
- A positive Δ ≥ 0.02 AUC on the stage task would validate the hypothesis
  that promoter-only aggregation is the right choice. If so, Phase 3's null
  result was a methods artefact, not a biological truth.
- A still-flat Δ would suggest that at this scale (54 genes), the Cu
  proteome simply does not carry much stage-related methylation signal,
  regardless of aggregation.

## Caveats
- Promoter-only aggregation loses a lot of probes for some genes. Genes
  like MT-CO1/CO2 have no promoter probes on the 450K array (mitochondrial
  genome) — same gap as in Phase 3.
- Some Cu-gene probes appear in both TSS1500 and TSS200 slots; we keep
  the first row per (probe, gene) pair, which is reasonable but slightly
  arbitrary.
- Gene-level β is still a mean of promoter probes; an alternative is to
  take the *minimum* β (the most-unmethylated probe) as a proxy for
  transcriptional openness.

## Files produced
- `data/lihc_methylation_promoter.tsv` — promoter-only gene-level β
- `outputs/final_comparison/promoter_methylation.md` — this document
- `outputs/final_comparison/promoter_methylation_metrics.csv` — comparison table
""")
    print(f"\n[promo] wrote {OUT/'promoter_methylation.md'}")


if __name__ == "__main__":
    main()
