# Promoter-Only Methylation — Targeted Aggregation

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
- Cu genes with ≥1 promoter probe: **50 / 54**
- Promoter probes mapping to Cu genes: **321**
  (vs 848 all-region probes in Phase 3)
- Promoter β matrix: `data/lihc_methylation_promoter.tsv`

## Results (GAT, 5-fold StratifiedGroupKFold)

| task | features | ROC-AUC | balanced acc |
|---|---|---:|---:|
| Tumor vs Normal | expr_only | 0.993 | 0.913 |
| Tumor vs Normal | expr + meth (all regions) | 0.997 | 0.912 |
| Tumor vs Normal | **expr + meth (promoter-only)** | **0.997** | **0.908** |
| | | | |
| Stage I/II vs III/IV | expr_only | 0.668 | 0.525 |
| Stage I/II vs III/IV | expr + meth (all regions) | 0.667 | 0.518 |
| Stage I/II vs III/IV | **expr + meth (promoter-only)** | **0.673** | **0.599** |

## Interpretation

- **Δ on stage (promoter - all)**: +0.006 AUC.
- **Δ on stage (promoter - expr only)**: +0.005 AUC.

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
