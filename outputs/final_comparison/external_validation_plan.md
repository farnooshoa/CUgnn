# External Validation — Feasibility & Plan (Step 5)

## Status
- **Not executed in this run** — the ICGC DCC REST API was deprecated; the
  current `dcc.icgc.org` and `platform.icgc-argo.org` are JS single-page
  apps and do not expose an automated JSON endpoint.
- What this means for the pilot: external validation is deferrable but
  meaningful; it is the single step that would upgrade the pilot from
  "TCGA-internal" to "generalises" — the standard bar for publication.

## Recommended cohorts

| Source | Cohort | N (HCC) | Data format | How to obtain |
|---|---|---:|---|---|
| ICGC DCC | LIRI-JP (Japan) | 258 | tsv.gz, RNA-seq gene counts + clinical | Manual download from https://dcc.icgc.org/releases/current/Projects/LIRI-JP after GUI selection |
| NCBI GEO | GSE14520 | 488 (China, HBV-related) | Affymetrix microarray (Illumina bead chip) | `GEOparse` or `GEOquery` — fully scriptable |
| NCBI GEO | GSE76427 | 167 | Affymetrix | scriptable |
| NCBI GEO | GSE36376 | 433 (Korea) | Illumina HiSeq | scriptable |
| CPTAC | LIHC proteogenomic | ~165 | proteogenomics portal | PDC API |

GSE14520 is the most commonly-cited external LIHC cohort — pair it with TCGA
for the headline external-generalisation table.

## Minimum pipeline for Step 5 (when ready)

```python
# 1. Download & normalise external cohort
#    Normalize to log2(x+1), subset to 54 Cu genes, remap probes if microarray.
# 2. Align feature space
external = external[ds.expression.index]    # drop non-Cu genes
# 3. Train on TCGA only, predict on external
cfg = TrainConfig(model="gat", epochs=80)
bundle_tcga = build_graph_dataset(ds, graph, zscore_train_mask=np.ones(ds.n_samples, dtype=bool))
# train model fully ...
# 4. Build external bundle with *TCGA-fit* z-score statistics (no leakage)
#    Apply same mean/std learnt on TCGA to external samples.
# 5. Report: external-cohort ROC-AUC, balanced accuracy, overlap of top-10
#    important genes between TCGA and external.
```

## Expected outcome

- **Tumor vs Normal on GSE14520** should stay AUC ≥ 0.95 — same biology,
  same task. If it collapses, it flags a TCGA-specific overfit.
- **Stage classification** is much less transferable (AJCC staging
  conventions differ across cohorts and eras). Expect AUC ≈ 0.60 at best.

## Why we are stopping here for now

External validation requires either:
1. A one-off manual portal download (~30 minutes GUI time) followed by
   scripted processing — practical and cheap.
2. Or a GEO-based automation via `GEOparse` — fully scriptable in ~30 min.

Both are feasible; neither is a code problem but a data-logistics one. This
document captures the plan so Phase 3 can pick it up cleanly.
