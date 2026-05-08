# Missing Real Data — Action Items for Phase 2 Rerun

Phase 2 cannot start until two files are added to `data/`. This document lists the minimal steps to unblock the pipeline.

---

## 1. Files required

```
data/lihc_expression.tsv   # MISSING
data/lihc_metadata.tsv     # MISSING
```

Full specifications (column names, value conventions, example rows) are in `data/EXPECTED_INPUT_FORMAT.md`. Please read that file once before producing the inputs.

## 2. Suggested data sources

Any of these will work; pick whichever you already have authorisation to use:

| Source | Endpoint | Notes |
|---|---|---|
| GDC Data Portal | https://portal.gdc.cancer.gov (Project: TCGA-LIHC, Workflow: STAR - Counts or HTSeq - FPKM-UQ) | Matches the upper-quartile FPKM used in Blockhuys 2017 most closely |
| Xena Toolbox | https://xenabrowser.net/ -> TCGA-LIHC cohort | Prenormalised log2(FPKM-UQ+1); fastest for a pilot |
| `TCGAbiolinks` (R) | `GDCquery(project="TCGA-LIHC", data.category="Transcriptome Profiling", data.type="Gene Expression Quantification", workflow.type="STAR - Counts")` | Reproducible, version-pinned |
| `pyTCGA` / `TCGAutils` | python | Less common, but viable |

For the pilot, Xena's prepackaged TSV is the cheapest path.

## 3. Minimal preprocessing before the rerun

A. **Gene identifier mapping**
   - If the download uses Ensembl IDs, map to HGNC symbols (biomaRt / `mygene` / pyensembl) before saving the TSV.
   - Upper-case gene symbols. Drop any row whose symbol is empty or contains a period (deprecated Ensembl noise).
   - Deduplicate by symbol (keep max-variance row if the source has multiple per symbol).

B. **Log transformation**
   - If values are raw counts or linear FPKM/TPM, apply `log2(x + 1)` before writing.

C. **Sample-type labels**
   - TCGA barcodes encode sample type in positions 14–15:
     - `01`–`09` → Tumor
     - `10`–`19` → Normal
   - Produce `sample_type` in the metadata accordingly.

D. **Column layout**
   - Expression: first column `gene_symbol`, remaining columns = sample IDs matching metadata.
   - Metadata: `sample_id` + `sample_type` are required; `stage`, `grade`, `overall_survival_days`, `vital_status` are optional but recommended for the optional Phase 2 extension tasks.

A reference one-liner to produce the two files from a Xena-style download:
```bash
python scripts/prepare_tcga_lihc.py   # user-supplied; not committed to this repo
# writes data/lihc_expression.tsv and data/lihc_metadata.tsv
```
(The repo does not include a downloader because the pipeline should stay agnostic to the source — produce the two TSVs in whatever way is convenient.)

## 4. Sanity checklist before rerun

```bash
python -c "
from src.preprocessing import load_lihc_dataset
ds = load_lihc_dataset(require_real=True)
print('samples:', ds.n_samples, 'genes:', ds.n_genes)
print('tumor:', int(ds.tumor_mask.sum()),
      'normal:', int((~ds.tumor_mask).sum()))
print('missing Cu genes:', ds.missing_genes)
"
```

Expected for a complete TCGA-LIHC dump:
- samples ≈ 370–424 (cohort varies slightly by snapshot)
- genes: 54
- tumor ≈ 370, normal ≈ 50
- missing Cu genes: should be `[]` or a short list (≤ 3)

If the script raises `FileNotFoundError`, the filenames are wrong. If it prints 0 samples matched, your `sample_id`s do not match the expression column headers.

## 5. Exact command to rerun the full Phase 2 pipeline

After the two files are in place:

```bash
# all three phases — baselines + GNN + comparison summary
python run_pipeline.py --require-real-data

# or split, if you want to inspect intermediate state
python run_pipeline.py baselines --require-real-data
python run_pipeline.py gnn       --require-real-data
python run_pipeline.py compare   --require-real-data
```

The `--require-real-data` flag prevents the synthetic fallback. The pipeline will:
1. Regenerate `outputs/baseline/gene_coverage_report.md`, `copper_de_results.csv`, `copper_heatmap.png`, `pca_scatter.png`, `umap_scatter.png`, `copper_network_logfc.png`, `copper_modules.csv`, `module_summary.md`.
2. Train GCN and GAT; regenerate `outputs/gnn/model_metrics.csv`, `graph_embedding_umap.png`, `node_importance.csv`, `top_subgraph_or_attention_summary.md`.
3. Rewrite `outputs/final_comparison/model_comparison_summary.md`.

The remaining Phase 2 deliverables that are specifically hand-authored summaries (Step-5 `real_gnn_summary.md`, Step-4 `real_baseline_summary.md`, Step-6 `real_model_comparison_summary.md`, Step-7 `lihc_biological_interpretation.md`, Step-8 `collaboration_ready_summary.md`) will need to be produced from the regenerated artefacts — either by rerunning this assistant on the populated outputs, or manually. They are **not** auto-generated because they require biological judgment that should not be inferred from synthetic data.

## 6. If you would like the assistant to produce the narrative summaries

Once the two TSVs are present and `python run_pipeline.py --require-real-data` has run successfully, re-prompt with the Phase 2 instructions and the assistant will read the regenerated artefacts and write:

- `outputs/baseline/real_baseline_summary.md`
- `outputs/gnn/real_gnn_summary.md`
- `outputs/final_comparison/real_model_comparison_summary.md`
- `outputs/final_comparison/lihc_biological_interpretation.md`
- `outputs/final_comparison/collaboration_ready_summary.md`

These files are deliberately not produced in the current run, because deriving them from the synthetic fallback would be misleading.
