# `data/` — what is and is not in the handoff

## Included (lightweight, ready to use)

| File | Size | Purpose |
|---|---:|---|
| `manifest.tsv` | 104 KB | GDC manifest of all 425 TCGA-LIHC STAR-Counts files. Used by `scripts/download_tcga_lihc.py` to re-pull the raw data. |
| `lihc_expression.tsv` | 427 KB | Preprocessed 58-gene × 424-sample log2(FPKM-UQ+1) matrix. Used directly by the GNN scripts; lets the pipeline run without redownloading raw data. |
| `lihc_expression_expanded.tsv` | 1.1 MB | Expanded gene set (broader than the 58-node Cu proteome) for sensitivity analysis. |
| `lihc_metadata.tsv` | 33 KB | sample_type, stage, vital_status, days_to_death, etc. |
| `gse14520_expression.tsv` | 122 KB | External validation cohort (Roessler 2010, GSE14520), already preprocessed. |
| `gse14520_metadata.tsv` | 25 KB | Matching metadata for the validation cohort. |
| `lihc_methylation.tsv` | 361 KB | DNA methylation per probe, restricted to Cu-proteome promoter regions. Not yet used in the main pipeline. |
| `lihc_methylation_promoter.tsv` | 280 KB | Aggregated promoter-level methylation. |
| `natcomm_agg_per_gene.tsv` | 10 MB | Per-gene aggregated logFC table from the collaborator's 2020 Nat Comm Agilent microarray (8 samples). |
| `string_v12_copper_edges.tsv` | 2 KB | STRING v12 PPI subset for the 58 Cu genes. |
| `string_v12_expanded_edges.tsv` | 14 KB | STRING v12 PPI for the expanded gene set. |
| `EXPECTED_INPUT_FORMAT.md` | 6 KB | Original specification of the expected input file formats. |

## Excluded from the handoff (re-downloadable)

| Item | Original size | How to regenerate |
|---|---:|---|
| `gdc_raw/` (425 STAR-Counts TSVs) | **2.4 GB** | `python ../scripts/download_tcga_lihc.py` (uses `manifest.tsv`). Takes ~30 min on a normal connection. |
| `xena_lihc_methylation.gz` | 377 MB | UCSC Xena: TCGA-LIHC HumanMethylation450 beta-values. Re-download from xenabrowser.net. |
| `illumina_450k_manifest.csv` | 192 MB | Illumina HumanMethylation450 v1.2 manifest. Free download from Illumina support site. |
| `xena_450k_probemap.tsv` | 17 MB | UCSC Xena probe-to-gene mapping. Re-download from xenabrowser.net. |
| `geo_cache/GSE14520_family.soft.gz` | 64 MB | NCBI GEO direct download via `GEOparse` (already-processed `gse14520_expression.tsv` is included, so this raw soft.gz is rarely needed). |

If you only care about the GNN pipeline, you do **not** need to re-download anything — `lihc_expression.tsv` + `lihc_metadata.tsv` are enough.

If you want to rebuild the methylation extension, you will need to re-fetch the Xena methylation file and the Illumina manifest.
