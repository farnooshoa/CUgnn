# Important Figures — Blockhuys et al. 2017

Colors are described verbatim from the figure legends and the main text. Where the paper makes a claim about what a color "means", that claim is preserved here. Anything inferred is flagged as *interpretation*.

---

## Fig. 1 — Cellular localisation of the 54 Cu-binding proteins

### What it shows
A human cell cartoon with 10 compartments — cell membrane, intracellular vesicles, cytoskeleton, mitochondrion, endoplasmic reticulum, cytoplasm, nucleus, Golgi apparatus, extracellular space — each labelled with the Cu-binding proteins the authors assigned to that compartment. Assignments come from Genecards "Localization" with UniProtKB + COMPARTMENTS confidence 4–5; ATOX1, MEMO1, MT3, LTF, S100A5 use COMPARTMENTS-only assignment; CUTA and MT4 use COMPARTMENTS confidence 2 (less stringent). MOXD2P (a pseudogene of MOXD1) is assumed to co-localise with MOXD1 in the ER.

Several genes appear in multiple compartments (e.g. ATP7A in cell membrane + Golgi; PRNP in mitochondrion, cytoplasm, nucleus, Golgi, extracellular; SNCA in five compartments).

### Color meaning (from legend)
- **Blue** = transporter
- **Orange** = enzyme
- **Black** = protein with other or unknown function

### Why it matters
- Establishes the **54-gene "human copper proteome"** that the rest of the paper — and our LIHC project — uses as the node set.
- Provides two clean, machine-usable per-node features: **functional category** (blue/orange/black) and **subcellular compartment**.
- Makes it visually obvious that Cu flux is distributed: no single compartment owns Cu, so a graph with cross-compartment edges is biologically meaningful.

---

## Fig. 2 — Expression heat map of the Cu proteome across 18 cancer types

### What it shows
A matrix of log2(tumor / normal) for each of the 54 Cu-binding genes (rows, y-axis) across 18 TCGA tumour types (columns, x-axis): thymus, head & neck, esophagus, adrenal gland, bladder, stomach, soft tissue, kidney, pancreas, bile duct, **liver**, prostate, cervix, breast, uterus, thyroid, colorectal, lung.

Row and column dendrograms come from unsupervised hierarchical clustering with Ward's-D2 linkage on a Pearson-correlation dissimilarity matrix. The dendrograms split:
- Genes into **8 clusters** (cluster 1: AOC1→DBH; 2: MT4→GPC1; 3: CUTA→AFP; 4: CP→ALB; 5: S100A5→S100A13; 6: ENOX2→SNCA; 7: LOX→TYR; 8: LOXL1→AOC3).
- Cancer types into **3 groups**: group 1 (thymus, head & neck, esophagus, adrenal, bladder, stomach); group 2 (soft tissue, kidney, pancreas, bile duct, **liver**); group 3 (prostate, cervix, breast, uterus, thyroid, colorectal, lung).

### Color meaning (from legend)
- **Red** = upregulation in tumor vs normal (log2FC > 0)
- **Blue** = downregulation (log2FC < 0)
- **White** = no change
- **Asterisk (*)** = |log2FC| ≥ 0.4 AND t-test p < 0.05 (considered significant)
- Color scale runs approximately from −10 to +10 log2FC.

### Why it matters
- Makes it visible that Cu-proteome dysregulation is **not uniform** — LOX/LOXL1-2, SPARC, ENOX2 are broadly up; S100B, PRNP, ENOX1, SNCA, SOD3, HEPHL1, AOC3 are broadly down.
- The 8 gene clusters in Fig. 2 become the nodes that Fig. 3 expands into interaction networks — the clustering drives graph construction.
- Liver's membership in group 2 — closest to soft tissue, kidney, pancreas, bile duct — is a direct LIHC-relevant finding we can cross-reference later.

---

## Fig. 3 — Network analysis of the 8 gene clusters

### What it shows
Eight panels (one per cluster from Fig. 2). Each panel is a subnetwork generated in Cytoscape with the GeneMANIA App, seeded with the genes of that cluster; first-level partners are allowed in (they appear as grey circles). Networks are layered on physical + co-expression + genetic interactions.

### Color meaning (from legend)
**Node colors (applied to the seed Cu-binding genes, grounded in breast cancer):**
- **Red** = upregulated in breast cancer
- **Blue** = downregulated in breast cancer
- **Black** = no change in breast cancer
- **Grey** = first-level interaction partner (not a Cu-binding seed)

**Edge colors:**
- **Red** = physical (protein–protein) interactions
- **Purple** = co-expression
- **Green** = genetic interactions

### What "close proximity" means
GeneMANIA's layout weights inverse to gene rank in the result list; genes that are strongly connected end up visually close. The paper uses this proximity qualitatively: tight clusters (1, 3, 8) = "closely related in terms of known functional parameters"; spread clusters (2, 4) = less connected.

### Why it matters for graph construction
- Fig. 3 is a **multi-edge-type interaction graph over Cu-binding proteins + their first-degree neighbours**. That is essentially the input format modern GNNs want.
- It makes explicit which **edge types** domain experts consider relevant: physical interaction, co-expression, genetic interaction. This motivates using STRING/BioGRID/GeneMANIA-style evidence channels when we build edges for LIHC.
- It shows that Cu-proteome clusters are **highly inter-connected** — a gene assigned by expression to one cluster often has first-level partners in several other clusters. This suggests a single merged graph over all 54 genes (not 8 isolated subgraphs) is the right scale.

See `figure3_graph_construction_notes.md` for how we translate this into our LIHC graph.

---

## Fig. 4 — Proteogenomic correlations in breast-cancer PAM50 subtypes

### What it shows
Per-gene Pearson correlation between mRNA (RNA-seq) and protein (MS/MS) levels for 36 of the 54 Cu-binding proteins in 29 luminal A, 33 luminal B, 18 HER2-enriched and 25 basal-like breast-cancer samples from Mertins et al. 2016. Rows ordered by correlation coefficient (high at top).

### Color meaning
- **Left column (Pearson's r)** colour blocks: yellow = high correlation, yellow/green bands = moderate, red = low/no correlation (legend: high / moderate / low / no correlation).
- **Main heat map** uses a single per-row scaled blue → red colour (row min → row max) for RNA-seq / Protein abundance — not a cross-gene absolute scale.
- Grey/white square = missing value.

### Why it matters
- Calibrates RNA-based inference: 28/36 (78%) proteins show high/moderate correlation between mRNA and protein, which is why the authors (and we) treat RNA-seq as a reasonable first proxy.
- Flags genes (e.g. ENOX2, COMMD1, AOC2, F5, ALB — low r) where mRNA-level predictions may not translate to protein — useful when interpreting our GNN's top features.

---

## Fig. 5 & Fig. 6 — ATOX1 validation in breast tissue

### What Fig. 5 shows
Two IHC micrographs of ductal breast tissue — panel A = benign, panel B = cancerous — stained with an anti-ATOX1 antibody. Scale bar 100 μm. Illustrative, no quantification.

### What Fig. 6 shows
- Panel A: reference TMA cores scored as negative / weak / moderate ATOX1 staining.
- Panel B: bar chart of 67 cancerous + 4 normal breast-tissue TMA cores, binned by ATOX1 staining intensity (negative / weak / moderate) and PAM50 subtype (TNBC, HER2, luminal A, luminal B). Luminal A has the largest "negative" bar (17); luminal B has 12 negative, 9 weak, 1 moderate; HER2 has 3 negative, 4 weak; TNBC has 8, 4, 1. Fisher's exact test between subtypes and between tumour vs. normal was **not significant** given the small and unmatched sample size.

### Why they matter
- Converts a genomic observation (Fig. 2: ATOX1 up in breast) into a validated protein-level observation in real patient tissue — proof-of-concept for the paper's workflow.
- Confirms cytoplasmic + **nuclear** ATOX1 signal, consistent with its reported transcription-factor-like function.
- For our LIHC project, these figures are a **template**: once our pipeline nominates a candidate Cu gene driving LIHC classification, IHC (or an equivalent validation) is the appropriate next step.
