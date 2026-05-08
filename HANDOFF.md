# Handoff — TCGA-LIHC Copper Proteome Graph-Learning Pilot

**Snapshot date**: 2026-04-22
**Project root (zipped)**: `handoff_CUgnn_2026-04-22.zip`
**Original author**: Ruiheng
**Collaborator**: Dr. Schwartz-Duval

This document is the handoff briefing. Read it once, then dive into the
referenced files. Everything here is grounded in the actual code and outputs
in this folder — no external state required.

---

## 1. What this project does — one paragraph

A pilot graph-learning pipeline that maps TCGA-LIHC bulk RNA-seq onto the
54-protein human copper proteome (Blockhuys 2017 *Metallomics*) plus 4
representative histones (Attar 2020 *Science*, H3-H4 Cu reductase). Each
patient is represented as the same 58-node graph with a fixed curated
edge set; node features are log2(FPKM-UQ+1) expression + z-score +
functional-category one-hot. A GAT graph classifier predicts (a) tumor
vs normal, (b) stage I/II vs III/IV, (c) 3-year overall survival.
Interpretability comes from gradient saliency on nodes and learned
attention on edges. The pilot has been extended to a cartoon-cell figure
that overlays an extracellular **Au compounds** node with literature-
backed (solid) and provisional (dashed) interactions.

## 2. Where things live

| Path | What it is |
|---|---|
| `README.md` | Original repo readme. Quick-start commands, expected inputs. |
| `requirements.txt` | Python deps. PyG + torch + sklearn + pandas. |
| `run_pipeline.py` | Top-level driver: baselines → GNN → comparison. |
| `src/preprocessing/` | Loads `data/lihc_expression.tsv` + metadata; subsets to 58 Cu genes; falls back to synthetic demo if data is missing. |
| `src/graph_building/build_graph.py` | `CURATED_CU_EDGES` (≈75 hand-curated Cu-biology edges). Three graph variants: PPI / functional / co-expression. **Functional graph is the production graph.** Histone edges (Attar 2020) are at the bottom of the constant. |
| `src/baseline_models/` | DE, heatmap, PCA, UMAP/t-SNE, static network, module detection, LR/RF/SVM tabular baselines. |
| `src/gnn_models/dataset.py` | Builds per-patient PyG `Data` objects (5-dim node features, shared edge_index). |
| `src/gnn_models/models.py` | GCN + GAT classifier definitions. |
| `src/gnn_models/train.py` | Training loop, CV, **gradient saliency** (`node_importance_from_gradient`, line 205), **attention readout** (`attention_summary`, line 228). |
| `scripts/histone_rerun.py` | The script that produced the histone-augmented results (saliency + attention CSVs in `outputs/final_comparison/`). |
| `scripts/make_cell_cartoon_network.py` | Builds the cartoon-cell HTML. Has Save/Load layout, Au-compounds node, citation tooltips. |
| `scripts/make_visualizations.py` | Static figures: top-15 saliency bar (`fig2_node_importance_bar`), network with importance overlay, etc. |
| `scripts/download_tcga_lihc.py` | Re-downloads STAR-Counts from GDC using `data/manifest.tsv`. |
| `outputs/paper_2017_extraction/copper_gene_list.csv` | The 58-gene node list. |
| `outputs/baseline/` | DE results, heatmap, modules, PCA. |
| `outputs/gnn/` | 54-node baseline GNN outputs (node importance, attention edges). |
| `outputs/final_comparison/histone_*` | **58-node (with histones) outputs.** This is the current production set. |
| `outputs/visualizations/` | Static figures + the interactive HTMLs. |
| `outputs/email_package/` | Files sent to collaborator on round 1 (cartoon HTML + histone_results.md + natcomm_qc.md + 2 figures). |
| `outputs/reply_4_22/` | **Files prepared for the round-2 reply** (updated cartoon HTML + glossary.md + email_body.md). Not yet sent at handoff time. |

## 3. The 58 nodes and the graph

- 54 Cu proteome genes from Blockhuys 2017 + 4 histones (H3-3A, H3-3B, H3C1, H4C1) from Attar 2020.
- Functional categories: 30 enzyme, 12 transporter, 16 other_or_unknown.
- ~75 curated edges (`CURATED_CU_EDGES`): physical / coexpression / genetic types from the source papers' figures and known Cu biology.
- ~10–20 auto-added shared-compartment edges (from `subcellular_localization` field, only when 2 ≤ compartment members ≤ 8).
- All patients share the same topology; only node features (expression) vary.

## 4. Headline results (58-node, 5-fold StratifiedGroupKFold)

| Task | GAT ROC-AUC | Comparison |
|---|---:|---|
| Tumor vs Normal | 0.994 | SVM tabular baseline 0.998 (task is easy) |
| Stage I/II vs III/IV | 0.679 | +0.011 over 54-node |
| 3-year Overall Survival | 0.732 | **+0.040 over 54-node — main interpretability gain** |

Top saliency hits: HSPA6, ATF3 / heat-shock, ZFAND2A (zinc finger), H3-3B (rank 7), HMOX1, etc. ATP7B and ATOX1 sit in the top-10 attention-edge band; SOD3 ↔ ATOX1 (the breast TF interaction the collaborator validated previously) and CCS ↔ SLC31A1 are recovered.

## 5. Datasets

**TCGA-LIHC (training cohort, 374 tumor + 50 normal)**
- Raw STAR-Counts: **NOT included in this handoff** (2.4 GB).
- `data/manifest.tsv` (104 KB) lists all 425 GDC files. Re-download with `python scripts/download_tcga_lihc.py`.
- `data/lihc_expression.tsv` (preprocessed 58-gene × 424-sample log2 matrix) is included so the GNN scripts run without redownloading.
- `data/lihc_metadata.tsv` carries sample_type, stage, vital_status, days_to_death.

**GSE14520 (external validation cohort)**
- `data/gse14520_expression.tsv` + `data/gse14520_metadata.tsv` already preprocessed.
- Originals from GEO via `geo_cache/` (excluded from handoff, re-fetchable).

**2020 Nat Comm transcriptomics (collaborator's data)**
- `data/natcomm_agg_per_gene.tsv` — per-gene aggregated logFC from his Agilent 2-colour microarray (8 samples: 4 condition A vs 4 condition C).
- Raw `Raw Nat Comm Transcriptomics.xlsx` at the project root.
- QC report: `outputs/email_package/natcomm_qc.md` (3 Cu hits vs 82 Zn hits; 27.3× ratio; supports collaborator's claim of Au-induced Zn-finger displacement).

## 6. Collaborator interaction history

**Round 1** — sent on (date when first email_package was zipped). Files:
- `outputs/email_package/copper_cell_cartoon_network.html`
- `outputs/email_package/histone_results.md`
- `outputs/email_package/natcomm_qc.md`
- `outputs/email_package/natcomm_cu_vs_zn.png`
- `outputs/email_package/natcomm_cu_proteome_heatmap.png`

Email draft: `outputs/final_comparison/email_draft_followup.md`.

**Round 2 — in flight at handoff**.
Collaborator's reply (paraphrased) requested:
1. Save layout button on the cartoon (so he can drag and keep arrangement) ✓ added
2. Add ENOX2 + PRNP as nodes ⚠️ they were already in the 58-node graph; ranks reported instead
3. Add an "Au compounds" symbol with edges to ATP7A/B, ATOX1, LOX, histones (H1/H3/H4/H2A/H2B), PRNP, ENOX2 ✓ added (with literature/unpublished tier)
4. Send 4 PDFs as Au-interaction references ✓ extracted to `outputs/reply_4_22/email_body.md` table
5. Video call request ✓ accepted in draft email
6. Mentioned several field-specific terms unfamiliar ✓ wrote `outputs/reply_4_22/glossary.md`

Files prepared for the reply, **not yet sent**:
- `outputs/reply_4_22/email_body.md` — draft email
- `outputs/reply_4_22/copper_cell_cartoon_network.html` — updated cartoon (Au node + Save/Load)
- `outputs/reply_4_22/glossary.md` — 1-page ML/stats glossary

Two open questions for collaborator (flagged in the email):
1. **LOX naming ambiguity** — our LOX = lysyl oxidase; the 2011 paper's LOX = lipoxygenase. Confirm which one he meant.
2. **SLC31A1 (CTR1)** — explicitly ruled out as Au uptake route in Spreckelmeyer 2018; the cartoon does not draw it. Confirm scope.

Awaiting his unpublished data on H1/H2A/H2B + PRNP + ENOX2 Au interaction to upgrade dashed edges to solid. He also mentioned wanting H1 / H2A / H2B added — they are **not** in the current 58-node graph; would need new node additions and likely retraining.

## 7. The Au compounds node — design rationale

`scripts/make_cell_cartoon_network.py:61` defines `AU_NODE` and `AU_EDGES`. Important:

- The Au node is a **visualization-only** node. It is **not** part of the trained graph (`build_functional_graph` does not add it). It is appended to the cartoon's nodes list at the end of `load_payload()`.
- Edges have a `confidence` field: `"lit"` (solid line) or `"unpublished"` (dashed line).
- Tooltips show citation + 1-sentence mechanism.

If you decide to bring Au into the trained model (per a future grant aim), you'll need to:
1. Decide what "node features" the Au node should carry (it has no expression — placeholder zeros today).
2. Add the Au edges to `CURATED_CU_EDGES` with a new edge type.
3. Augment `_categorical_features` to include an `au_compound` category.
4. Retrain. Document everything as a new pipeline branch.

## 8. Saliency + attention — how the numbers are produced

- **Saliency** (`train.py:205`): `|∂ logit_tumor / ∂ x|` aggregated over 5 features and averaged over the first 50 samples. Single seed (RANDOM_SEED=42), single model. **Point estimate, no uncertainty quantification.**
- **Attention** (`train.py:228`): GAT last-layer α, averaged across 8 heads, summed across 50 samples. Reported as `attention_sum` per directed edge.
- Reports use **rank**, not raw value, when communicating to the collaborator (raw values have mixed units and are hard to interpret in absolute terms).

## 9. Pending / next-step ideas (not yet started)

1. **Signed saliency**. Drop the `np.abs()` on the gradient and report the signed value of the expression-feature gradient. Lets you compare model direction against log2FC direction (textbook biology check). 5-line code change in `train.py`.
2. **Multi-seed bootstrap** for saliency CIs. Run `histone_rerun.py` with 5–10 seeds, report mean ± std and rank IQR. ~5–10× compute, still cheap.
3. **GSE14520 external validation**. Pull the 58-gene subset from `data/gse14520_expression.tsv`, score with the trained TCGA model (no retraining), report transferred AUCs. This was on the round-1 to-do list and never finished.
4. **Zn proteome extension** for the grant aim. Curate ~60–80 Zn-binding genes (UniProt zinc-binding ∪ known ZnF-TF hubs), build a parallel pipeline, validate on the 2020 Nat Comm cohort (which is Zn-loud per QC). Ideal next push for the metallome-GNN grant story.
5. **H1 / H2A / H2B addition** if the collaborator's unpublished data justifies it. Pick representative genes, add to `copper_gene_list.csv`, define edges in `CURATED_CU_EDGES`, retrain.
6. **Pin STRING v12 subnetwork** as the PPI source. Currently the fallback uses `CURATED_CU_EDGES`; `data/string_v12_copper_edges.tsv` is already downloaded but not wired in.
7. **Patient-level GroupKFold** for tumor-vs-normal CV. The `StratifiedGroupKFold` is already used for stage and survival; double-check tumor/normal split is consistent.

## 10. How to run

```bash
# install
pip install -r requirements.txt

# preprocessing → baselines → GNN → comparison
python run_pipeline.py

# rerun the 58-node histone variant (regenerates outputs/final_comparison/histone_*)
python scripts/histone_rerun.py

# regenerate the cartoon HTML (after editing AU_NODE / MANUAL_POS)
python scripts/make_cell_cartoon_network.py

# regenerate static figures
python scripts/make_visualizations.py
```

If `data/lihc_expression.tsv` is missing the pipeline auto-falls-back to a synthetic demo. **Do not interpret synthetic outputs.** To re-fetch real TCGA:

```bash
python scripts/download_tcga_lihc.py    # uses data/manifest.tsv
```

## 11. Known gotchas

- `RANDOM_SEED` is hard-coded in `src/utils/`. Saliency / attention are deterministic given a fixed seed.
- z-score features in CV use train-fold stats only when `zscore_train_mask` is passed — make sure to pass it; the legacy default uses all-sample stats and leaks slightly.
- The cartoon `MANUAL_POS` dictionary in `make_cell_cartoon_network.py:80` defines node positions. To update from a saved layout JSON: copy the `positions` field from the JSON, replace the dict values. (No auto-bake script yet; could be added in 10 lines.)
- `outputs/cu_followup_package.zip` is the round-1 zip the user previously generated. Round-2 zip would be from `outputs/reply_4_22/`.

## 12. Contact / context for successor

If picking this up cold:

1. Read `README.md` for the original framing.
2. Read this `HANDOFF.md` for current state.
3. Open `outputs/final_comparison/histone_results.md` and `outputs/email_package/natcomm_qc.md` for the two main reports already shared with the collaborator.
4. Open `outputs/visualizations/copper_cell_cartoon_network.html` in a browser to see the cartoon.
5. Read `outputs/reply_4_22/email_body.md` to see what's next to send.

The collaborator's grant narrative direction is **multi-metallome GNN**: Cu pipeline shown, Zn extension as the planned proof of concept, Au as a perturbation overlay. Keep this in mind when prioritising next steps.
