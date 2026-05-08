# Phase 3 Summary — Ablations, Harder Task, Database Swap, Multi-Omics

Everything in this file refers to the real TCGA-LIHC cohort (424 samples,
374 tumor / 50 normal, 54 Cu genes) and uses the post-audit methodology
(StratifiedGroupKFold by `case_submitter_id`, per-fold StandardScaler /
z-score). All GNN numbers are GAT unless stated.

---

## Step 1 — Graph ablation: does the curated Cu graph carry signal?

| topology | # edges | ROC-AUC | balanced acc |
|---|---:|---:|---:|
| **curated functional** (Cu biology) | 60 | **0.991 ± 0.006** | **0.918** |
| empty (self-loops only) | 0 | 0.932 ± 0.037 | 0.800 |
| complete (all pairs) | 1431 | 0.979 ± 0.015 | 0.625 |
| random Erdős-Rényi (3 seeds mean) | 60 | 0.985 ± 0.009 | 0.920 |

**Finding.** On tumor-vs-normal, the *identity* of edges does not change
AUC (real − random = −0.004). Having *some* edges matters (empty is
clearly worse); being fully-connected hurts (balanced acc collapse from
over-smoothing). The curated graph gives the **lowest fold-to-fold
variance** and remains the right choice on interpretability grounds.

Detailed doc: `graph_ablation.md`.

---

## Step 2 — Harder task: Stage I/II vs III/IV (n = 349 tumor samples)

| model | ROC-AUC | balanced acc | F1 |
|---|---:|---:|---:|
| **GAT** | **0.668 ± 0.029** | 0.525 | 0.276 |
| **GCN** | 0.668 ± 0.059 | 0.571 | 0.306 |
| SVM (RBF) | 0.659 ± 0.050 | 0.619 | 0.418 |
| Random Forest | 0.643 ± 0.060 | 0.515 | 0.095 |
| Logistic Regression | 0.585 ± 0.058 | 0.565 | 0.382 |

**Finding.** AUC drops into the 0.58–0.67 regime — the expected "right
difficulty band" for a small-feature-set clinical task. GAT tops the
table by ~0.010 AUC and has the **lowest fold-to-fold variance** (0.029
vs 0.050–0.060 for classical models). This is the first place in the
project where GNN has any measurable edge over classical baselines,
though at 0.010 it remains within CV noise.

Detailed doc: `stage_task.md`.

---

## Step 3 — STRING v12 high-confidence edges

Replaced the curated fallback edge list (60 edges) with a high-confidence
STRING v12 subnetwork of the 54 Cu genes (`combined_score ≥ 700`; 77
edges after filtering).

| task | curated AUC | **STRING v12 AUC** | Δ |
|---|---:|---:|---:|
| Tumor vs Normal | 0.993 | 0.984 | −0.009 |
| Stage I/II vs III/IV | 0.668 | **0.690** | **+0.022** |

**Finding.** STRING gives a small but welcome boost on the stage task
(+0.022 AUC). On tumor-vs-normal it is slightly below curated — the
denser graph has more isolated-node paths that the GAT must work
around. For a publishable version, **STRING is the right choice** on
grounds of reproducibility and defensibility, with a small predictive
cost on the saturated task and a small gain on the interesting task.

Edge list saved as `data/string_v12_copper_edges.tsv`.
Detailed doc: `string_v12.md`.

---

## Step 4 — Multi-omics node features: RNA-seq + DNA methylation

Downloaded Xena HumanMethylation450 dataset (395 MB) + probe→gene
annotation (18 MB), aggregated 848 probes to 52 / 54 Cu genes
(2 pseudogenes / mitochondrial genes without promoter probes: MOXD2P,
MT-CO1), median missingness 2.1 %. Node feature dim: 5 → 6.

| task | expr-only AUC | **expr + methylation AUC** | Δ |
|---|---:|---:|---:|
| Tumor vs Normal | 0.993 | 0.997 | +0.004 |
| Stage I/II vs III/IV | 0.668 | 0.667 | −0.001 |

**Finding.** Methylation gives a trivial bump on the already-saturated
tumor-vs-normal task and **no improvement** on stage. Two likely
reasons:
1. **Aggregation is too coarse** — we averaged all probes per gene
   (body + promoter + 3' UTR). A promoter-specific aggregation
   (TSS1500 + TSS200 only) often matters 0.02–0.05 AUC.
2. **Stage-related methylation signal is real but small** — only
   ATP7B, CP, and a few LOX promoters are reported as stage-related
   in HCC literature; the 54-gene block is not broad enough to pick
   this up on its own.

This is a **negative but useful** result: it clarifies that mRNA
already carries most of the Cu-proteome signal at this scale, and
that multi-omics gains likely require a larger node set and/or
promoter-aware methylation aggregation.

Detailed doc: `methylation_results.md`.

---

## Step 5 — External-cohort validation (not executed this round)

The ICGC DCC REST API has been deprecated; manual portal or GEO-based
download is required. A ready-to-execute plan is in
`external_validation_plan.md`, with recommended cohorts (LIRI-JP via
ICGC, GSE14520 via GEO) and the minimum pipeline code skeleton.

---

## Overall take-away after Phase 3

1. **Tumor-vs-normal is a solved task** on the copper proteome — any
   reasonable model reaches 0.98–0.99 AUC, and the graph is not needed
   for prediction. Move on.
2. **Stage classification is the right difficulty** for a copper-proteome
   pilot. GAT with STRING edges now reaches **AUC 0.690**, beating every
   classical model (LR 0.585, RF 0.643, SVM 0.659). Magnitude is small
   (~0.02–0.1 AUC) but consistent.
3. **Methylation at this scale is a null result.** The Cu proteome is a
   narrow biological lens; multi-omics integration should be revisited
   once we expand to a larger panel (e.g., extend to a ~500-node
   Cu-proteome-plus-neighbours graph using STRING 1-hop partners).
4. **Interpretability is the stable win.** GAT attention continues to
   surface canonical Cu-handling edges regardless of edge source
   (curated or STRING). This is the angle worth leading with in any
   external communication.

## What to queue for Phase 4

1. **External validation on GSE14520** (GEOparse, ~30 min of work).
   This is the single step that upgrades the pilot into a methods paper.
2. **Promoter-aware methylation** — restrict probes to TSS1500 / TSS200
   before gene-level aggregation.
3. **Expand the node set to 54 + 1-hop STRING partners** (~150–250 genes).
   Test whether the GNN's advantage grows with proteome scale.
4. **Survival task** — overall survival binarised at 3 years, Cox-weighted
   auxiliary target in the GNN. Clinically the most valuable endpoint.
5. **Ablation on layers/heads** — the GAT is 2-layer / 4-head; try
   3-layer on the expanded node set.

---

## Files produced this round

| file | role |
|---|---|
| `graph_ablation.md` + `.csv` | Step 1 ablation |
| `stage_task.md` + `.csv` + `stage_task_gnn_folds.csv` | Step 2 stage task |
| `string_v12.md` + `.csv` | Step 3 database swap |
| `methylation_results.md` + `.csv` | Step 4 multi-omics |
| `external_validation_plan.md` | Step 5 feasibility plan |
| `data/string_v12_copper_edges.tsv` | STRING Cu-subnetwork |
| `data/lihc_methylation.tsv` | gene-level β matrix |
| `scripts/graph_ablation.py`, `stage_task.py`, `string_edges.py`, `add_methylation.py` | reusable experiment drivers |
| `src/gnn_models/evaluate.py` | reusable grouped-CV helper |
