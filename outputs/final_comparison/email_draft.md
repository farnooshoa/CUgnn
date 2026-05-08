# Draft email — for review before sending

**To:** Dr. Schwartz-Duval *(please confirm spelling / salutation)*
**Subject:** Pilot results — TCGA-LIHC copper proteome graph model
**Attachments to consider:** `project_report.md`, `copper_interactive_network.html`, `collaboration_ready_summary.md`

---

Dear Dr. Schwartz-Duval,

I wanted to share a short update on the TCGA-LIHC copper-proteome pilot that builds on Blockhuys et al. 2017 (*Metallomics*, 9, 112–123).

**What was done.** I downloaded the full TCGA-LIHC STAR-Counts cohort from GDC (424 samples: 374 tumor, 50 matched normal) and subset expression to the 54-protein human copper proteome defined in the paper (100 % coverage). Each patient is represented as a graph: the 54 proteins are nodes with shared fixed topology (curated Cu-handling edges — ATOX1↔ATP7A/B, CCS↔SOD1, COX17↔SCO1/2↔MT-CO1/2, LOX family, etc. — plus shared-compartment edges from the paper's Fig. 1). Node features are log2(FPKM-UQ+1) expression, z-score, and functional category. I trained GCN and GAT graph classifiers against logistic regression, random forest, and SVM as tabular baselines.

**Key findings.**
1. **Signal is strong and LIHC-consistent.** AFP↑, ALB↓, CP↓, mitochondrial COX↓ (Warburg-like), LOX/LOXL2/SPARC↑ — all textbook HCC biology. 5-fold CV ROC-AUC is 0.995 (GAT), 0.998 (SVM).
2. **The GNN earns its place on interpretability, not accuracy.** GAT attention surfaces canonical Cu-handling edges — CCS↔SLC31A1, SOD3↔ATOX1, ATP7B↔COMMD1, MT-CO2↔SCO2 — i.e. the model passes messages through real biology, not statistical shortcuts. SOD3↔ATOX1 in particular is the transcription-factor interaction you and colleagues validated in breast.
3. **A potentially organ-specific finding.** SLC31A1 is downregulated in LIHC (log2FC −1.08) — opposite to the upregulation the 2017 paper reports in breast. ATP7B and ATOX1 rank high in GNN importance despite modest log2FC because they sit at the centre of the Cu-secretory module.

**Caveats I want to flag honestly.** The tumor-vs-normal task is easy (single-gene DBH already reaches 0.96 AUC); the current CV does not group folds by patient, which likely gives a mildly optimistic estimate; all features are mRNA-level and 22 % of Cu-binders had poor mRNA/protein correlation in the paper's Fig. 4. The pilot is therefore best read as **feasibility + hypothesis generation**, not a final methods paper.

**Proposed next steps.** (1) Tighten the evaluation with patient-level GroupKFold and a permutation-label sanity check; (2) swap the fallback edge list for a pinned STRING v12 subnetwork; (3) move to harder tasks (early-vs-late stage, overall survival); (4) identify two or three candidates for IHC validation — ATP7B, ATOX1, and CP are the most defensible first targets given both the DE and GNN rankings.

A self-contained interactive network visualisation is attached (`copper_interactive_network.html`, opens in any browser — drag/hover/zoom). The full methods and results are in `project_report.md`; the one-page collaborator summary is `collaboration_ready_summary.md`.

Happy to walk through any of this when you have time — either a short call or asynchronous comments on the report would both be useful.

Best regards,
*[Your name]*
