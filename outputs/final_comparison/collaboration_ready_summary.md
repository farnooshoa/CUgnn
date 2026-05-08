# Collaboration-Ready Summary — TCGA-LIHC Copper Proteome Pilot

## Data used
- **TCGA-LIHC RNA-seq** (GDC STAR-Counts workflow, Data Release 45.0, Dec 2025).
- **424 samples** (374 tumor, 50 normal), expression values as `log2(FPKM-UQ + 1)`.
- **54-gene node set** from Blockhuys *et al.* 2017, *Metallomics* — the complete human copper proteome. All 54 genes present in the expression matrix.

## How the graph was defined
- Nodes: the 54 Cu-binding proteins (12 transporters / 26 enzymes / 16 other).
- Edges: curated Cu-homeostasis interactions (ATOX1↔ATP7A/B, CCS↔SOD1, COX17↔SCO1/2↔MT-CO1/2, LOX family, S100 family, ATP7B↔COMMD1, …) + shared-subcellular-compartment edges. Average degree ≈ 2.2.
- One graph per patient; same topology, patient-specific node features (expression + z-score + category one-hot).
- Task: graph-level classification, tumor vs normal.

## Is the copper-proteome concept useful in LIHC?

**Yes.** Every classifier we tried reaches AUC ≥ 0.98 on the 54-gene feature space, indicating the copper proteome carries a clear tumor-vs-normal signal in liver cancer. More importantly, the specific genes the models highlight map onto three well-known LIHC axes: the hepatocyte secretory program (ALB↓, CP↓, AFP↑), mitochondrial Cu-dependent oxidative phosphorylation (MT-CO1/CO2↓), and the ECM-remodelling LOX family (LOX↑, LOXL2↑, SPARC↑). GAT attention reproduces textbook Cu-handling edges (CCS–SLC31A1, SOD3–ATOX1, ATP7B–COMMD1, MT-CO2–SCO2) — evidence that the model is learning through biology, not just statistics.

## Top genes / modules worth discussing

1. **CP (ceruloplasmin)** — log2FC −2.06, adj p < 10⁻³¹. The largest hepatocyte-specific Cu signal in the cohort.
2. **AFP** — log2FC +1.73. Positive control for HCC biology; validates the approach.
3. **ATP7B** — high GNN centrality + saliency despite modest mean shift. Wilson's-disease gene; easy hook for a clinical collaborator.
4. **ATOX1** — Cu chaperone; high GNN saliency; the 2017 paper's flagship protein in breast. Worth an IHC follow-up in LIHC tissue.
5. **SLC31A1** — log2FC −1.08. **Down in LIHC, opposite to the breast pattern described in the 2017 paper**. Candidate organ-specific Cu-flux switch.
6. **Mitochondrial COX module** (MT-CO1, MT-CO2, SCO1, SCO2, COX11, COX17) — uniformly downregulated. Links Cu biology to the Warburg effect in HCC.
7. **LOX / LOXL2 / SPARC ECM module** — uniformly upregulated. Consistent with the paper's cross-cancer ECM-remodelling finding.

## Is this promising as preliminary data?

**Yes, for three reasons:**
- The biology is legible (AFP↑, ALB↓, LOX↑, COX↓) — a biologist can immediately see what the model is doing.
- The graph + attention + saliency outputs provide a natural interpretability story that a tabular model cannot.
- The 54-gene node set is small enough to extend gracefully to multi-modal node features (methylation, mutation, protein if available) and to clinical-covariate-based edge types.

**Caveats:**
- 50-normal sample size is small; all tumor-vs-normal metrics should be read as "signal is easy to detect", not "model is production-ready".
- mRNA ≠ protein; the paper itself flags ~22% of Cu-binders with poor mRNA/protein correlation.
- The curated edge list is a reasonable starting point but is NOT a substitute for a pinned STRING / BioGRID release. Switching to a curated database is a mechanical next step.

## Practical next steps

1. **Replace the fallback edge list** with a pinned STRING v12 Cu-proteome subnetwork; rerun and confirm top-genes/edges are stable.
2. **Add node features beyond expression**: the obvious wins are DNA methylation (β-values per gene) and mutation status (binary). LIHC has TCGA-level methylation and maf data; layering these as additional node features should improve the GNN's advantage over SVM.
3. **Move from tumor-vs-normal to clinically useful tasks**: early vs late stage, BCLC stage, and 3-year survival. These are harder and are where graph structure is more likely to add predictive value.
4. **Validate 2–3 candidate genes at the protein level**: ATOX1, ATP7B, CP, LOX are the obvious first targets for TMA-based IHC in a small retrospective cohort — a direct parallel to Blockhuys Fig. 5–6.
5. **Cross-cancer comparison**: rerun the same 54-node graph on TCGA-BRCA and TCGA-BLCA (breast and bladder — one up-cluster, one down-cluster in Blockhuys 2017 Fig. 2) and compare module saliency. If the LIHC SLC31A1-down vs BRCA SLC31A1-up finding replicates, we have an interesting organ-specific copper-flux story worth a focused methods paper.
6. **Consider a relational GNN**: with multiple edge types (physical / coexpression / genetic / shared-compartment), an R-GCN or CompGCN would exploit structure the current GCN/GAT does not. Low priority until Step 2 improves node features.

## One-line status
Preliminary data support the hypothesis that the Blockhuys 2017 copper proteome is a useful, compact, biologically interpretable feature set for LIHC — with clear next-step experiments that would turn this pilot into a short methods paper or grant aim.
