# Expected Input Format — TCGA-LIHC Phase 2

This pipeline consumes two files in this directory:

```
data/lihc_expression.tsv
data/lihc_metadata.tsv
```

Both are tab-separated. UTF-8. No byte-order mark. No trailing blank columns.

---

## 1. `data/lihc_expression.tsv`

**Orientation**: rows = genes, columns = samples.
The first column is the gene identifier; all remaining columns are sample IDs.

**First column** (gene id):
- header name: `gene_symbol`
- values: HGNC gene symbols, upper-case, unquoted (e.g. `ATP7A`, `MT-CO1`, `SLC31A1`). **Do not use Ensembl IDs as the primary identifier** — if your source uses Ensembl, please map to HGNC first; otherwise the Cu-proteome subset step will drop everything.
- duplicate symbols: the loader keeps the first occurrence and drops later ones, but prefer a deduplicated matrix.

**All other columns** (sample ids):
- one column per TCGA sample
- header = TCGA barcode or any stable sample id (e.g. `TCGA-DD-A4NS-01A-11R-A27V-07`)
- sample IDs must also appear in `lihc_metadata.tsv` — any sample in the expression matrix without a metadata row will be **dropped silently** during loading (and vice versa).

**Values**:
- Normalised expression on a log-transformed scale. Any of:
  - `log2(FPKM + 1)`
  - `log2(TPM + 1)`
  - `log2(upper-quartile FPKM + 1)` (closest to the 2017 paper)
  - DESeq2 `rlog` or `vst`
- Do **not** provide raw counts without transformation — the DE, PCA, heat-map, and GNN all assume near-Gaussian, log-scale values.
- Missing values: use empty cell or `NA`. The loader does not impute.

### Example snippet
```tsv
gene_symbol	TCGA-DD-A4NS-01A	TCGA-DD-A4NS-11A	TCGA-BC-A10X-01A	TCGA-BC-A10X-11A
ATP7A	6.12	5.84	6.31	5.47
ATP7B	8.93	9.14	8.51	9.02
SLC31A1	7.21	6.98	7.44	7.02
ATOX1	9.05	8.41	9.21	8.33
CP	13.42	12.80	13.11	12.45
LOX	8.77	6.90	9.02	7.12
```

(The example shows paired tumor `-01A` and normal `-11A` samples in the TCGA barcode convention.)

---

## 2. `data/lihc_metadata.tsv`

**Orientation**: rows = samples, columns = annotations. One row per sample.

**Required columns**:
- `sample_id` — must exactly match the column headers in `lihc_expression.tsv` (case and whitespace sensitive after a `.strip()`).
- `sample_type` — values must be **`Tumor`** or **`Normal`** (case-insensitive; the loader title-cases them on read). This is the column used as the tumor-vs-normal label.

**Optional columns** (used only if present, never required):
- `stage` — AJCC stage string, e.g. `I`, `II`, `III`, `IV`, or `IIIA`. Used later for the optional stage-classification task.
- `grade` — histological grade (`G1`, `G2`, …).
- `overall_survival_days` — integer days from diagnosis.
- `vital_status` — `Alive` / `Dead`.
- `gender`, `age_at_diagnosis`, `race` — demographic covariates.

Any other columns are passed through but ignored.

### Example snippet
```tsv
sample_id	sample_type	stage	grade	overall_survival_days	vital_status
TCGA-DD-A4NS-01A	Tumor	II	G2	1245	Alive
TCGA-DD-A4NS-11A	Normal				
TCGA-BC-A10X-01A	Tumor	IIIA	G3	412	Dead
TCGA-BC-A10X-11A	Normal				
```

Empty cells for `Normal` samples are acceptable.

---

## 3. Matching between the two files

- The loader computes `shared = expression.columns ∩ metadata.sample_id` and restricts both to this intersection.
- After intersection, samples are aligned positionally by this intersection, then the 54 Cu proteome genes are subset out of the expression matrix.
- Missing Cu genes are recorded in `outputs/baseline/gene_coverage_report.md` and listed by name.
- The pipeline does not fail on missing Cu genes — it uses whatever subset is present (most TCGA LIHC expression matrices include all 54).

---

## 4. Preprocessing / normalization assumptions

- Expression values are assumed to be log-scale already. If you pass linear-scale FPKM, the tumor-vs-normal log2FC will be mis-scaled and the heat-map will saturate.
- No per-sample normalization is applied by the pipeline — assume the upstream tool (GDC upper-quartile normalisation, DESeq2 vst, etc.) has already handled it.
- No batch correction is applied (TCGA-LIHC is already relatively batch-controlled; adding ComBat here would risk over-correcting the tumor-vs-normal signal).
- Gene symbols are upper-cased and whitespace-trimmed on read.

---

## 5. Missing Cu genes

If a Cu proteome gene is not present in the expression file:
1. It is recorded in `outputs/baseline/gene_coverage_report.md`.
2. It becomes an **isolated node** in any graph that references it, so downstream message-passing ignores it.
3. It is still a node of the graph — the node set is always all 54 genes so that graph topology is identical across patients.
4. For any patient feature vector the missing gene is given `NaN` expression; the loader currently **requires** all subset genes to be present (raises `KeyError` otherwise). If you expect missing Cu genes, pre-fill them with `NA` rows in the expression file and the loader will drop them.

A sensible workflow is:
1. Run a one-liner: `python -c "from src.preprocessing import load_lihc_dataset; ds = load_lihc_dataset(); print(len(ds.missing_genes), ds.missing_genes)"` to see missing genes.
2. Decide whether to (a) accept isolated nodes or (b) add rows of `NA` for the missing genes.

---

## 6. Quick sanity test

Before running the full pipeline, run:

```bash
python -c "
from src.preprocessing import load_lihc_dataset
ds = load_lihc_dataset()
print('samples:', ds.n_samples, 'genes:', ds.n_genes)
print('tumor:', int(ds.tumor_mask.sum()),
      'normal:', int((~ds.tumor_mask).sum()))
print('missing Cu genes:', ds.missing_genes)
"
```

Expected output for TCGA-LIHC (~371 samples, ~50 normal) with all 54 Cu genes present:

```
samples: 421
genes: 54
tumor: 371 normal: 50
missing Cu genes: []
```

If the script instead prints the synthetic-fallback banner, your files are not being detected — check that the filenames are **exactly** `lihc_expression.tsv` and `lihc_metadata.tsv` and live in `./data/`.
