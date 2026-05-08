# LIHC Biological Interpretation — Cu Proteome (Exploratory)

**Tone**: hypothesis-generating. Everything below is consistent with the data in `outputs/baseline/` and `outputs/gnn/` but needs orthogonal validation before any strong claim.

## 1. Which Cu-associated genes appear most dysregulated in LIHC?

**Down in tumor (absolute log2FC):**
- **DBH (-2.84)**, **ALB (-2.38)**, **CP (-2.06)**, **MT-CO1 (-1.14)**, **SLC31A1 (-1.08)**, **MAP2K1 (-1.02)**, **MT-CO2 (-0.92)**, **SOD1 (-0.82)**, **AOC3 (-0.71)**.

**Up in tumor:**
- **AFP (+1.73)**, **LOXL2 (+1.01)**, **SPARC (+0.75)**, **LOX (+0.61)**, **CUTA (+0.59)**, **MT3 (+0.36)**.

The pattern is dominated by loss of hepatocyte-specific secretory proteins (ALB, CP, DBH), collapse of mitochondrial Cu-dependent oxidative phosphorylation (MT-CO1, MT-CO2, SCO1), and gain of ECM-remodelling / metastatic signal (LOX, LOXL2, SPARC). AFP follows the textbook HCC fetal/adult switch.

Both GNN saliency (top 10: CP, SOD1, MT-CO1, ATOX1, MAP2K1, ALB, LOXL4, MT-CO2, LOX, ATP7B) and classical logistic-regression coefficients rank the same core set of genes. ATOX1 and ATP7B are pulled up by the GNN despite modest log2FC because they sit at the centre of highly-attended Cu-secretion subgraphs.

## 2. Which Cu-associated subgraphs/modules appear most important?

From greedy modularity on the functional graph (see `outputs/baseline/module_summary.md`), three modules stand out:

**(a) Hepatocyte Cu/Fe secretome — module 3** (ALB, CP, HEPH, HEPHL1, LTF, F5).
Strongly downregulated. All six genes are significant at BH < 0.001. This is the most biologically intuitive finding — HCC hepatocytes lose the adult secretory program.

**(b) Mitochondrial Cu / COX assembly — module 2** (COX11, COX17, MT-CO1, MT-CO2, SCO1, SCO2).
Uniformly downregulated. Consistent with the Warburg-style metabolic shift and with reduced Cu delivery to cytochrome c oxidase. Mechanistically interesting because the Cu import is also blunted (SLC31A1↓), i.e. upstream and downstream of the mitochondrial axis move together.

**(c) ECM remodelling / LOX axis — module 1** (LOX, LOXL1-4, SPARC, GPC1).
LOX and LOXL2 upregulated with high statistical significance; SPARC up; GPC1 down. This mirrors Blockhuys 2017's pan-cancer finding (LOX, LOXL1-2, SPARC frequently up) — the liver is not an exception. The LOX axis is known to drive premetastatic niche formation in other cancers.

Top GAT attention edges point to a fourth, cross-module axis: **CCS–SLC31A1–SOD1–SOD3 cytosolic antioxidant handover** and the **ATP7B–COMMD1–ATOX1 secretory triad**. These did not emerge as unified modules in the greedy-modularity output because their edges span the secretory and cytoplasmic partitions; attention captures the cross-partition signal better than community detection.

## 3. Does the LIHC pattern look organ-specific?

Partly yes. Comparing against the Blockhuys 2017 pan-cancer Fig. 2 and the breast-specific Section IV:

| Observation | LIHC | Blockhuys 2017 breast |
|---|---|---|
| SLC31A1 (Cu import) | **down** | up |
| ATP7B (Wilson) | unchanged | up |
| ATOX1 | near-baseline (mildly up in IHC across multiple cancers per Fig S1) | up (mRNA + IHC) |
| COX17 / SCO2 | **down** | up |
| LOX / LOXL1-2 | **up** | up |
| SPARC | **up** | up |
| ALB | **down** | context-specific |
| CP | **down** | context-specific |
| AFP | **up** | context-specific (and notable as an HCC marker) |

The **Cu-import / mitochondrial-delivery axis moves in opposite directions in LIHC vs breast** (down in LIHC, up in breast). This is a potentially organ-specific finding: the 2017 paper's explicit interpretation in breast was "increased Cu flow via SLC31A1 → ATP7B → mitochondria", whereas LIHC appears to be doing the opposite. If confirmed, it suggests LIHC does not depend on the same Cu-trafficking ramp-up that breast cancer does.

The **ECM / LOX axis moves up in both LIHC and breast** — likely a more universal cancer mechanism.

## 4. Findings consistent with the 2017 paper's cross-cancer logic

- Liver sits in the paper's **group 2** cancer cluster (soft tissue, kidney, pancreas, bile duct, liver) — broadly a Cu-downregulation cluster in Fig. 2. Our LIHC DE pattern (net downregulation across many nodes, with specific exceptions for LOX / LOXL2 / SPARC / AFP) is consistent with that position.
- The paper states SPARC, LOXL1-2 and ENOX2 are broadly upregulated across cancers. We see LOXL1 only marginally but LOXL2 and SPARC clearly upregulated, consistent.
- The paper states S100B, PRNP, ENOX1, SOD3 are broadly downregulated. We see PRNP strongly down (adj p 3.6e-31), SOD3 and ENOX1 modestly changed — partially consistent.
- Module analysis in our LIHC graph recovers the **Cu secretory axis** (ATOX1-ATP7A-ATP7B) that Fig. 1 of the paper lays out as the flagship copper-handling pathway.

## 5. Findings that are exploratory and need validation

- **ATOX1 as a liver-specific driver**: GNN saliency ranks ATOX1 fourth, but the log2FC in our cohort is modest. The paper's breast-cancer protein-level validation (Fig. 5–6) should be repeated in liver tissue before any claim.
- **ATP7B as a high-value node**: centrality + saliency both rank it high, but its mean expression does not shift dramatically. Needs a Wilson-disease-orthogonal line of evidence (functional Cu-efflux assay) before drawing conclusions.
- **MAP2K1 (MEK1)**: strongly downregulated in LIHC (log2FC -1.02). The paper lists MAP2K1 as Cu-dependent and tumor-growth-associated. Our direction is opposite to the default expectation of a growth driver; this may reflect LIHC's lower proliferation-dependent mitogen demand or a Cu-independent MEK downregulation. Exploratory.
- **AFP ↔ ATP7B attention edge**: the model gives this pair high attention, but the biological rationale is indirect (both are liver-secretory-path proteins). Worth noting, not worth claiming as a mechanistic coupling.
- Every result uses mRNA as a proxy for protein. The paper's Fig. 4 shows only 78% of Cu-binders correlate well between mRNA and protein — LIHC-specific protein measurements are still needed.

## 6. Which genes are most reasonable to highlight in a preliminary collaboration discussion?

A defensible short-list to open the conversation:

1. **CP (ceruloplasmin)** — largest-effect hepatocyte-specific Cu gene, down 2 log2 in tumor. Well-known plasma Cu carrier, easy to measure and interpret.
2. **AFP** — canonical HCC biomarker and a Cu-binding member of the paper's list. Sanity check that the approach is picking up real HCC biology.
3. **ATP7B** — Wilson's disease gene, central to liver Cu export; ranks high in GNN centrality despite modest log2FC.
4. **ATOX1** — Cu chaperone that the paper spent Fig. 5–6 validating in breast; GNN saliency suggests it also matters in LIHC.
5. **LOX + LOXL2** — ECM/fibrosis axis; upregulated; strong prior literature in HCC and metastasis.
6. **Mitochondrial COX module (MT-CO1, MT-CO2, SCO1/SCO2)** — treat as a single module, strongly downregulated; links Cu biology to the Warburg shift.
7. **SLC31A1** — the apical Cu importer; downregulated here, **opposite to breast**. Potential organ-specific switch.

Frame the discussion as "the copper proteome behaves differently in LIHC than in breast; here are six concrete axes to probe next."
