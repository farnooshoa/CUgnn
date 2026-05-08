# Real-Data Baseline Summary — TCGA-LIHC Cu Proteome

Source: GDC TCGA-LIHC STAR-Counts (Data Release 45.0). 424 samples (374 tumor / 50 normal), 54/54 Cu proteome genes covered. Values = `log2(FPKM-UQ + 1)`.

## 1. Strongest tumor-vs-normal copper changes

Top 15 most significant (by Benjamini–Hochberg adjusted p) from `copper_de_results.csv`:

| gene | log2FC | adj p (BH) | direction |
|---|---:|---:|---|
| ALB | **-2.38** | 1.2e-45 | down |
| SLC31A1 | **-1.08** | 2.6e-34 | down |
| CP | **-2.06** | 6.8e-32 | down |
| PRNP | -1.19 | 3.6e-31 | down |
| DBH | **-2.84** | 2.1e-28 | down |
| LOXL2 | **+1.01** | 3.6e-27 | up |
| SOD1 | -0.82 | 1.1e-18 | down |
| AFP | **+1.73** | 1.7e-17 | up |
| MAP2K1 | -1.02 | 1.9e-17 | down |
| AOC3 | -0.71 | 1.6e-16 | down |
| MT-CO1 | -1.14 | 2.6e-16 | down |
| SCO1 | -0.35 | 3.4e-16 | down |
| APP | -0.59 | 3.1e-15 | down |
| CUTA | +0.59 | 7.5e-15 | up |
| MT-CO2 | -0.92 | 7.2e-14 | down |

Every one of the 54 Cu genes is significant at BH < 0.05 — the copper proteome as a whole is systematically dysregulated in LIHC.

**Biology** (see `lihc_biological_interpretation.md` for deeper discussion):
- **AFP↑** and **ALB↓** are textbook HCC markers — the hepatocyte fetal/adult switch.
- **CP↓** is consistent with loss of hepatocyte secretory function.
- **Mitochondrial COX (MT-CO1/CO2) + SCO1/SCO2 ↓** is consistent with the Warburg shift.
- **LOX / LOXL2 / SPARC ↑** mirrors the ECM-remodelling signature the 2017 paper flags as broadly cancer-associated.
- **SLC31A1↓, ATP7B unchanged, ATOX1 ≈ baseline** paints a Cu-import-down / Cu-export-preserved picture — opposite to the breast-cancer pattern in the paper (import↑, ATP7B↑).

## 2. Does the copper gene space separate tumor vs normal?

Yes, **cleanly**. PCA of the 54-gene z-score matrix splits tumor and normal along PC1 (see `pca_scatter.png`). t-SNE (`umap_scatter.png`) recovers two well-defined clusters with almost no overlap. All three classical ML models achieve balanced accuracy ≥ 0.92 and AUC ≥ 0.997 on 5-fold CV (`classical_model_metrics.csv`).

## 3. Which nodes look central in the Cu graph?

Centrality on the **functional graph** (curated Cu edges + shared-compartment edges, 54 nodes / 60 edges):

Top-degree nodes (hand-counted from the functional graph edge list):
- **ATP7B, ATP7A** — Golgi secretory-path ATPases
- **SOD1** — bridges CCS and SOD3 (cytosolic antioxidant axis)
- **CP** — extracellular Cu hub, connected to HEPH, HEPHL1, ATP7B
- **LOX, SPARC** — ECM-remodelling hub
- **COX17** — mitochondrial Cu chaperone; connects to SCO1/SCO2/COX11

See `copper_network_logfc.png` for the visualisation with node size = centrality and colour = LIHC log2FC.

## 4. Which modules appear most relevant?

14 modules detected by greedy modularity; the 6 most Cu-relevant (by gene set) are:

| module | genes | interpretation |
|---|---|---|
| 0 | ATOX1, ATP7A, ATP7B, SLC31A2, AFP, APP, PRNP, SNCA | Cu secretory-path axis + amyloidogenic passengers + AFP (HCC marker) |
| 1 | LOX, LOXL1, LOXL2, LOXL3, LOXL4, SPARC, GPC1 | **ECM remodelling / fibrosis / metastasis** |
| 2 | COX11, COX17, MT-CO1, MT-CO2, SCO1, SCO2 | **Mitochondrial Cu / cytochrome c oxidase assembly** |
| 3 | ALB, CP, HEPH, HEPHL1, LTF, F5 | **Hepatocyte secretome + plasma Cu/Fe delivery** |
| 4 | MAP2K1, MEMO1, S100A5, S100A12, S100A13, S100B | MAPK / EMT / S100 family |
| 5 | CCS, PARK7, SLC31A1, SOD1, SOD3 | Cu import + cytosolic antioxidant axis |

Three of these (1, 2, 3) are directly interpretable for LIHC:
- module 3 is the **hepatocyte-specific** Cu/plasma-protein secretion axis, strongly downregulated in tumor (ALB, CP both down).
- module 2 is **mitochondrial Cu assembly**, uniformly downregulated — consistent with the metabolic reprogramming expected in HCC.
- module 1 is the **LOX/SPARC ECM axis**, upregulated — matches the paper's pan-cancer observation.

## 5. Ranked list of top changed Cu genes (|log2FC|, sig only)

| rank | gene | log2FC | direction |
|---|---|---:|---|
| 1 | DBH | -2.84 | down |
| 2 | ALB | -2.38 | down |
| 3 | CP | -2.06 | down |
| 4 | AFP | +1.73 | up |
| 5 | S100A12 | -1.37 | down |
| 6 | PRNP | -1.19 | down |
| 7 | MT-CO1 | -1.14 | down |
| 8 | SLC31A1 | -1.08 | down |
| 9 | MAP2K1 | -1.02 | down |
| 10 | LOXL2 | +1.01 | up |
| 11 | MT-CO2 | -0.92 | down |
| 12 | SOD1 | -0.82 | down |
| 13 | SPARC | +0.75 | up |
| 14 | AOC3 | -0.71 | down |
| 15 | LOX | +0.61 | up |

Full list with adjusted p-values: `copper_de_results.csv`.

## Files in this step
- `gene_coverage_report.md` — 54/54 coverage, no missing genes
- `copper_de_results.csv` — full DE table
- `copper_heatmap.png` — z-score clustermap with Tumor/Normal annotation bar
- `pca_scatter.png`, `umap_scatter.png` — dimensionality-reduction views
- `copper_network_logfc.png` — functional network coloured by LIHC log2FC
- `copper_modules.csv`, `module_summary.md` — module detection
- `graphs/*` — edge lists + adjacency matrices for all three graph variants
