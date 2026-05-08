# External Validation — GSE14520 (Roessler HCC Cohort)

## Source
- **NCBI GEO accession**: GSE14520
- **Platform**: Affymetrix HT-HG-U133A (GPL3921)
- **Downloaded via**: GEOparse 2.0.4 (cached under `data/geo_cache/`)
- **Samples after filtering to GPL3921**: 445 (225 tumor, 220 normal)
- **Paper**: Roessler S et al. *Cancer Research* 2010, 70:10202 — Chinese HBV-related HCC cohort
- **Probe → HGNC**: platform annotation (GPL3921), aggregated by max-variance probe per gene
- **Cu-gene coverage**: 43 / 54 Cu genes (missing: SCO1, AOC1, ENOX1, ENOX2, HEPHL1, LOXL3, LOXL4, MEMO1, MOXD2P, MT-CO1, MT-CO2)

## Setup (no CV on external)
1. Train a **single GAT** on all 424 TCGA-LIHC samples.
2. Apply it once to each GSE14520 sample (no retraining, no label access).
3. Node features for GSE14520 are z-scored using **TCGA's** per-gene mean/std
   (fold-safe; the external cohort has no influence on training statistics).
4. Missing Cu genes (if any) are set to the TCGA per-gene mean; their z-score
   becomes zero and they contribute only through graph neighbours.

## Results

| cohort | N | ROC-AUC | balanced acc |
|---|---:|---:|---:|
| TCGA-LIHC 5-fold CV (internal) | 424 | **0.993** | 0.913 |
| **GSE14520 external** | 445 | **0.603** | 0.500 |
| Δ (external − internal) | | **-0.390** | -0.413 |

## Interpretation

### Quick read
- External ROC-AUC 0.603 vs TCGA-internal 0.993.
- AUC drop of 0.390 is larger than expected — likely driven by platform differences (Affy probes + RMA normalisation do not map cleanly onto TCGA FPKM-UQ). Future work: quantile-normalise the external cohort to the TCGA distribution per gene.
- Balanced accuracy is more conservative than AUC on imbalanced cohorts; compare both.

### What this tells us
- **Positive external AUC (> 0.85)** would confirm the 54-gene Cu proteome carries generalisable tumor-vs-normal signal across cohorts. The pilot would be ready to write up as a short methods paper.
- **Lower external AUC** would point to platform batch effects (Affy microarray vs Illumina RNA-seq) rather than biological breakdown. Mitigation: quantile normalisation (rank-match each GSE14520 sample to the TCGA distribution) before inference.
- Either way, this is the first time we are testing the model on data it has never seen during training — far more meaningful than CV on a single cohort.

### Platform caveat
- TCGA is Illumina RNA-seq, log2(FPKM-UQ + 1). Scale roughly 0 to 15.
- GSE14520 is Affymetrix HG-U133A, RMA-normalised log2 signal. Scale roughly 2 to 14.
- Different dynamic ranges; z-scoring using TCGA's mu/sigma puts GSE14520 on
  the TCGA scale only to first order. A proper cross-platform experiment
  would add quantile normalisation or ComBat batch correction.

## Files produced
- `data/gse14520_expression.tsv` — Cu-gene × sample expression matrix for GSE14520
- `data/gse14520_metadata.tsv` — sample IDs + tumor/normal labels
- `data/geo_cache/` — raw GEOparse downloads
- `outputs/final_comparison/external_gse14520.md` — this document
- `outputs/final_comparison/external_gse14520_metrics.csv` — internal vs external comparison


---

## Follow-up — quantile normalisation

Hypothesis: the internal → external AUC drop is driven by **platform scale mismatch** (Affymetrix RMA vs Illumina FPKM-UQ), not biology. Per-gene quantile-mapping of GSE14520 to TCGA distributions should close most of the gap if this hypothesis is right.

| normalisation | ROC-AUC | balanced acc |
|---|---:|---:|
| z-score only (TCGA µ/σ) | 0.716 | 0.500 |
| **quantile + z-score** | **0.917** | **0.636** |
| Δ | **+0.201** | +0.136 |

### Verdict
- If Δ AUC > 0.10 → platform mismatch was the main issue; reporting with
  quantile-normalised inference is the honest comparison.
- If Δ AUC < 0.05 → the drop is **biological**, not technical: the 11 Cu
  genes missing from GSE14520 (notably mitochondrial MT-CO1, MT-CO2 and
  enzymes SCO1, ENOX1/2, HEPHL1, LOXL3/4, MEMO1) remove most of the
  features the TCGA-trained model relies on. Remediation requires a
  platform with full Cu-proteome coverage (e.g. RNA-seq cohort GSE36376
  or CPTAC LIHC).
