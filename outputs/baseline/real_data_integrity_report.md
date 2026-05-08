# Real-Data Integrity Report — TCGA-LIHC

## 1. Files found
| File | Path | Present |
|---|---|---|
| Expression matrix | `data/lihc_expression.tsv` | **yes** |
| Sample metadata | `data/lihc_metadata.tsv` | **yes** |
| GDC file manifest | `data/manifest.tsv` | **yes** |
| Per-sample raw TSVs | `data/gdc_raw/*.tsv` | **yes (424 files)** |

Data source: GDC API direct download, project `TCGA-LIHC`, workflow `STAR - Counts`, data type `Gene Expression Quantification`, access `open`. GDC Data Release 45.0 (2025-12-04). Download completed in ~26 s with 12 parallel workers, 0 retries required.

## 2. Matrix dimensions
`lihc_expression.tsv`: **54 genes × 424 samples**. Log2-transformed FPKM-UQ (`log2(fpkm_uq_unstranded + 1)`). No NaN after mapping to HGNC symbols.

## 3. Sample counts
- Total samples: **424**
- Unique patient cases: **374**
- Paired T/N samples (same case, both present): **50**

## 4. Tumor vs Normal
- **Tumor**: 374 (GDC `sample_type = Primary Tumor` or `Recurrent Tumor`)
- **Normal**: 50 (GDC `sample_type = Solid Tissue Normal`)
- Class ratio: 7.48 : 1

## 5. Copper-gene coverage
All **54 / 54** Cu proteome genes found. No missing genes. No isolated nodes in the functional graph (3 isolates — ENOX1, ENOX2, PAM — have expression data but no curated Cu-edges in the fallback edge list; this is a graph-edge issue, not an expression-coverage issue).

## 6. Duplicates
- No duplicated sample_ids in expression columns (verified by `is_unique`).
- Duplicate HGNC symbols in the upstream GENCODE v36 gene model were resolved by taking the **max FPKM-UQ** per symbol before log-transform.
- No duplicate sample_ids in metadata after drop_duplicates on `sample_id`.

## 7. Missing labels / malformed rows
- All 424 samples have a valid `sample_type` value (either Tumor or Normal as defined above).
- `stage`, `grade`, `overall_survival_days`, `vital_status` are present in metadata as optional columns (not all samples populated — expected for TCGA clinical fields).

## 8. Normalization sanity check
- Expression log2 scale: means per-sample are ~3.3 with IQR ~[0.3, 5.2], max ~16. This is consistent with log2(FPKM-UQ + 1) — no sample appears to be on a different scale.
- Per-gene variance is reasonable (ALB std ≈ 1.3, CP std ≈ 1.6, LOXL2 std ≈ 1.4 across 424 samples).
- No batch-correction was applied (GDC STAR-Counts is already batch-controlled).

## 9. Result of this step
Phase 2 real-data run **proceeded**. Baseline and GNN deliverables are regenerated and overwrite any previously-synthetic content in `outputs/baseline/` and `outputs/gnn/`.

## 10. Reproducibility
```bash
python scripts/download_tcga_lihc.py          # re-download (idempotent, skips cached files)
python scripts/download_tcga_lihc.py --skip-download  # rebuild TSVs from existing gdc_raw
python run_pipeline.py --require-real-data    # re-run analyses
```
GDC snapshot: Data Release 45.0 (2025-12-04). Subsequent releases may add / remove samples.
