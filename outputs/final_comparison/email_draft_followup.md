# Draft follow-up email

**To:** Dr. Schwartz-Duval
**Subject:** Follow-up on histones, Nat Comm data, cartoon cell
**Attachments:** `copper_cell_cartoon_network.html`, `histone_results.md`, `natcomm_qc.md`

---

Dear Dr. Schwartz-Duval,

Quick update on the five items from your last email.

**Histones.** I added H3-3A, H3-3B, H3C1, H4C1 and the interactions from Attar et al. 2020 (tetramer assembly, Cu hand-off to ATOX1, the H3H113N phenocopy of ctr1Δ, and the mitochondrial respiration defect). On the 58-node graph, 3-year survival ROC-AUC went from 0.692 to 0.732. H3-3B ranks 7 of 58 in GAT saliency, and all 11 paper-proposed histone edges appear in the attention readout. Details in `histone_results.md`.

**Nat Comm transcriptomics.** Your observation holds. Within the Cu proteome only 3 genes pass adj.P < 0.05 (PRNP, ATP7A, ENOX2). The zinc-annotated set has 82 significant hits out of 574 detected, a 27× ratio. The top zinc hits are almost all ZnF transcription factors, and the heat-shock markers (HSPA6, ATF3, HMOX1, DDIT3) are very strong. Consistent with Au-induced zinc displacement from zinc-finger proteins. The cohort is too small to train on (n = 4 vs 4, 29 of 58 Cu genes on the array), but it works well as a perturbation validation cohort.

**Cartoon cell.** Attached as `copper_cell_cartoon_network.html`. All 58 genes placed by subcellular localization over a simple cell drawing, in the same visual style as your slides. Node color is LIHC tumor-vs-normal log2FC, node size is GAT saliency, edge width is GAT attention. The cartoon scales with zoom, nodes are draggable, and Reset restores the layout.

**MetalPDB.** Kept on hold per your note.

**Proposed next step.** Given the Cu-quiet / Zn-loud contrast in your data, the cleanest next move is to curate a compact zinc proteome (60 to 80 genes) and rerun the same three tasks on TCGA-LIHC, then read GAT attention against your A vs C differentials as an orthogonal check. That gives us one story: metallome GNN framework, shown on Cu with TCGA-LIHC, validated on Zn with your Nat Comm cohort.

Alternatively I can first close the loop on the 58-node GSE14520 external validation before moving to zinc. Either order works. Let me know which you prefer.

Best regards,
*[Your name]*
