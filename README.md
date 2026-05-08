# TCGA-LIHC Copper Proteome Graph Learning — Pilot

Preliminary pipeline that maps TCGA-LIHC RNA-seq onto the 54-gene human
copper proteome defined by Blockhuys *et al.* (*Metallomics*, 2017,
9, 112–123) and explores whether tumor vs normal liver samples can be
distinguished with biologically informed graph learning.

This is a **pilot**: the goals are feasibility, interpretability, and
artefact generation — not a final production model.

---

## What's in this repo

| Path | Content |
|---|---|
| `2017-*.pdf` | The source paper |
| `outputs/paper_2017_extraction/` | Structured extraction of the paper: summary, figures, 54-gene CSV, key-findings JSON, Fig. 3 notes |
| `data/` | Drop `lihc_expression.tsv` and `lihc_metadata.tsv` here |
| `src/preprocessing/` | Expression + metadata loading, Cu-gene subset, synthetic fallback |
| `src/graph_building/` | PPI / functional / co-expression graph builders over the 54 nodes |
| `src/baseline_models/` | DE, heatmap, PCA, UMAP/t-SNE, static network, module detection, classical ML |
| `src/gnn_models/` | GCN + GAT graph classifiers, training, embedding viz, saliency, attention |
| `outputs/baseline/` | Baseline artefacts |
| `outputs/gnn/` | GNN metrics, embeddings, node importance |
| `outputs/final_comparison/` | Classical vs GNN comparison |
| `run_pipeline.py` | Top-level driver |

## Expected input

```
data/lihc_expression.tsv   # rows=gene_symbol, cols=sample_id, normalised FPKM / TPM / log2
data/lihc_metadata.tsv     # columns: sample_id, sample_type ('Tumor' / 'Normal'),
                           # optional: stage, grade, overall_survival_days, vital_status
```

If these files are absent, the pipeline falls back to a **synthetic demo**
dataset (120 tumor + 30 normal samples with a planted tumor signal). This
lets you verify the full pipeline runs end-to-end before plugging in real
TCGA-LIHC data. **Do not interpret synthetic outputs biologically.**

## Quick start

```bash
pip install -r requirements.txt

# Full pipeline (baselines + GNN + comparison summary)
python run_pipeline.py

# Just the baseline block
python run_pipeline.py baselines

# Just GNN training
python run_pipeline.py gnn
```

Outputs land in `outputs/`. Re-running is idempotent.

## Design choices

- **Nodes**: the 54 Cu-proteome genes from Blockhuys 2017 Fig. 1. Fixed
  across all patient graphs.
- **Edges**: curated Cu-homeostasis edges (ATOX1↔ATP7A/B, CCS↔SOD1,
  COX17↔SCO1/2↔MT-CO1/2, LOX family, S100 family, ...) plus optional
  shared-compartment edges. A co-expression variant is also generated
  for exploration. Plug in a STRING / BioGRID export if available.
- **Node features (per patient graph)**: expression, optional z-score,
  and a 3-dimensional functional-category one-hot (transporter / enzyme
  / other). One graph per patient -> **graph classification**.
- **Task**: tumor vs normal. Stratified 5-fold CV, class weights to
  handle the normal-sample scarcity typical of TCGA.
- **Baselines**: logistic regression, random forest, SVM on the same
  54-gene feature vector. Include them *before* claiming GNN wins.
- **Interpretability**: node saliency (all models), attention weights
  (GAT). Output is a small top-genes table that can be laid next to
  Fig. 3 of the paper for sanity-checking.

## Caveats

- RNA != protein. See the paper's Fig. 4 for gene-level mRNA/protein
  correlations; genes with low r (e.g. ENOX2, COMMD1, AOC2, ALB, F5)
  need protein-level validation before any strong claim.
- The curated edge list is a reasonable starting point but is *not* a
  substitute for a proper STRING / BioGRID / GeneMANIA dump at a
  pinned version.
- With 54 nodes and a ~150-sample LIHC cohort, the GNN and a classical
  model will often perform similarly on tumor-vs-normal. Parity is the
  expected outcome; the GNN earns its place on interpretability and
  extensibility, not on marginal accuracy gains.
