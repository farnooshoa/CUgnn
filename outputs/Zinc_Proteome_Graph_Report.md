# Zinc Proteome Graph — curated edge list for the LIHC GNN

**Goal.** A compact, interpretable Zn-centric proteome graph (analogous to the Blockhuys human copper proteome) for a Graph Neural Network on cancer transcriptomics, with an explicit focus on hepatocellular carcinoma (LIHC) and gold-induced zinc displacement.

**Scale.** 91 deduplicated edges across 74 unique nodes, drawn from 36 primary references covering zinc transport, metallothioneins, MTF1, Zn-finger TFs, oxidative stress, EMT, ferroptosis, and Au-induced Zn displacement.

**Edge schema.** `EDGE: SOURCE → TARGET | TYPE | CONFIDENCE | EVIDENCE | CANCER_RELEVANCE | LIHC_RELEVANCE | AU_RELEVANCE | PAPER_ID`.

**Pseudo-nodes.** `ZN_FREE` = labile cytosolic Zn²⁺ pool. `AURANOFIN` = Au(I/III) compound class; edges from this node represent the Au-Zn displacement mechanism wherever it has been experimentally documented for a partner protein.

---

## 1. Per-paper edge extraction

Each subsection is one primary reference. Edges are written in the requested EDGE format.

### Hogstrand2009_ZIP7_hub

*Zinc transporters and cancer: a potential role for ZIP7 as a hub for tyrosine kinase activation (Hogstrand et al., Trends Mol Med 2009)*  
[https://pubmed.ncbi.nlm.nih.gov/19246244/](https://pubmed.ncbi.nlm.nih.gov/19246244/)

- **EDGE:** SOURCE: `CSNK2A1` TARGET: `ZIP7` TYPE: `regulatory` CONFIDENCE: `strong`  
  EVIDENCE: CK2 phosphorylates ZIP7 at S275/S276 triggering Zn release from the ER lumen.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `possible` · AU_RELEVANCE: `no`

### Krezel2014_PTP1B_Zn

*Zinc ions modulate protein tyrosine phosphatase 1B activity (Bellomo, Krezel, Maret, Metallomics 2014)*  
[https://pubs.rsc.org/en/content/articlehtml/2014/mt/c4mt00086b](https://pubs.rsc.org/en/content/articlehtml/2014/mt/c4mt00086b)

- **EDGE:** SOURCE: `ZIP7` TARGET: `PTP1B` TYPE: `signaling` CONFIDENCE: `strong`  
  EVIDENCE: Zn released by ZIP7 inhibits PTP1B (PTPN1) at low-nM Ki, sustaining RTK signalling.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `possible` · AU_RELEVANCE: `no`
- **EDGE:** SOURCE: `ZN_FREE` TARGET: `PTP1B` TYPE: `signaling` CONFIDENCE: `strong`  
  EVIDENCE: Sub-nM rise in free Zn inhibits PTP1B and other PTPs, prolonging RTK phosphorylation.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `possible` · AU_RELEVANCE: `yes`
- **EDGE:** SOURCE: `ZN_FREE` TARGET: `PTEN` TYPE: `signaling` CONFIDENCE: `moderate`  
  EVIDENCE: Free Zn inhibits PTEN phosphatase, sustaining PI3K/AKT signalling.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `possible` · AU_RELEVANCE: `yes`

### Taylor2008_ZIP7_BC

*ZIP7-mediated intracellular zinc transport contributes to aberrant growth-factor signalling in antihormone-resistant breast cancer cells (Taylor et al., Endocrinology 2008)*  
[https://academic.oup.com/endo/article/149/10/4912/2455137](https://academic.oup.com/endo/article/149/10/4912/2455137)

- **EDGE:** SOURCE: `ZIP7` TARGET: `AKT1` TYPE: `signaling` CONFIDENCE: `strong`  
  EVIDENCE: ZIP7-mediated Zn release activates AKT in tamoxifen-resistant breast cancer.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `possible` · AU_RELEVANCE: `no`
- **EDGE:** SOURCE: `ZIP7` TARGET: `MAPK1` TYPE: `signaling` CONFIDENCE: `strong`  
  EVIDENCE: ZIP7-derived Zn activates ERK/MAPK in antihormone-resistant breast cancer.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `possible` · AU_RELEVANCE: `no`
- **EDGE:** SOURCE: `ZIP7` TARGET: `EGFR` TYPE: `signaling` CONFIDENCE: `moderate`  
  EVIDENCE: Increased ZIP7 enhances EGFR/IGF-1R signalling in endocrine-resistant breast cancer.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `unknown` · AU_RELEVANCE: `no`

### Chen2021_ZIP7_ferroptosis

*Zinc transporter ZIP7 is a novel determinant of ferroptosis (Chen et al., Cell Death Dis 2021)*  
[https://www.nature.com/articles/s41419-021-03482-5](https://www.nature.com/articles/s41419-021-03482-5)

- **EDGE:** SOURCE: `ZIP7` TARGET: `GPX4` TYPE: `regulatory` CONFIDENCE: `strong`  
  EVIDENCE: ZIP7 knockdown blocks ferroptosis upstream of GPX4 axis; Zn supplementation restores ferroptosis.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `possible` · AU_RELEVANCE: `possible`
- **EDGE:** SOURCE: `GPX4` TARGET: `ACSL4` TYPE: `shared-pathway` CONFIDENCE: `strong`  
  EVIDENCE: GPX4 antagonizes ACSL4-driven lipid peroxide accumulation in ferroptosis.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `possible`

### Ohashi2014_ZIP7_ER

*SLC39A7/ZIP7 promotes intestinal epithelial self-renewal by resolving ER stress (Ohashi et al., PLOS Genet 2016)*  
[https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1006349](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1006349)

- **EDGE:** SOURCE: `ZIP7` TARGET: `DDIT3` TYPE: `regulatory` CONFIDENCE: `strong`  
  EVIDENCE: ZIP7 inhibition activates ER-stress sensor CHOP (DDIT3).  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `possible` · AU_RELEVANCE: `possible`

### Franklin2012_ZIP14_HCC

*ZIP14 zinc transporter downregulation and zinc depletion in development and progression of hepatocellular cancer (Franklin et al., J Gastrointest Cancer 2012)*  
[https://pmc.ncbi.nlm.nih.gov/articles/PMC3724761/](https://pmc.ncbi.nlm.nih.gov/articles/PMC3724761/)

- **EDGE:** SOURCE: `CEBPA` TARGET: `SLC39A14` TYPE: `regulatory` CONFIDENCE: `moderate`  
  EVIDENCE: C/EBPalpha transcriptionally activates ZIP14 in normal hepatocytes.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `no`
- **EDGE:** SOURCE: `SLC39A14` TARGET: `MT1A` TYPE: `regulatory` CONFIDENCE: `strong`  
  EVIDENCE: ZIP14-mediated Zn uptake feeds the MTF1/MT1 axis; ZIP14 loss parallels MT1 loss in HCC.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `no`

### Kim2017_ZIP14_ER

*Hepatic ZIP14-mediated zinc transport is required for adaptation to ER stress (Kim et al., PNAS 2017)*  
[https://www.pnas.org/doi/10.1073/pnas.1704012114](https://www.pnas.org/doi/10.1073/pnas.1704012114)

- **EDGE:** SOURCE: `SLC39A14` TARGET: `ATF4` TYPE: `regulatory` CONFIDENCE: `strong`  
  EVIDENCE: ZIP14 KO hepatocytes show elevated ATF4 and ER-stress markers.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `no`
- **EDGE:** SOURCE: `SLC39A14` TARGET: `DDIT3` TYPE: `regulatory` CONFIDENCE: `strong`  
  EVIDENCE: ZIP14 loss raises CHOP and increases ER-stress apoptosis in liver.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `no`

### Hogstrand2013_ZIP6_STAT3_EMT

*A mechanism for EMT and anoikis resistance triggered by zinc channel ZIP6 and STAT3 (Hogstrand et al., Biochem J 2013)*  
[https://portlandpress.com/biochemj/article/455/2/229/81664/](https://portlandpress.com/biochemj/article/455/2/229/81664/)

- **EDGE:** SOURCE: `STAT3` TARGET: `SLC39A6` TYPE: `regulatory` CONFIDENCE: `strong`  
  EVIDENCE: STAT3 directly transcribes ZIP6/LIV-1; canonical STAT3 target.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `possible` · AU_RELEVANCE: `no`
- **EDGE:** SOURCE: `SLC39A6` TARGET: `GSK3B` TYPE: `signaling` CONFIDENCE: `strong`  
  EVIDENCE: ZIP6-mediated Zn inhibits GSK-3beta (direct + via AKT) in breast cancer EMT.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `possible` · AU_RELEVANCE: `no`
- **EDGE:** SOURCE: `SLC39A6` TARGET: `SNAI1` TYPE: `regulatory` CONFIDENCE: `strong`  
  EVIDENCE: ZIP6 -> GSK3beta inactivation stabilizes Snail; represses E-cadherin -> EMT.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `possible` · AU_RELEVANCE: `no`

### Wei2013_Slug_ZEB1_EMT

*Transcriptional activation of ZEB1 by Slug leads to cooperative regulation of the EMT-like phenotype (Wels et al., 2011)*  
[https://pmc.ncbi.nlm.nih.gov/articles/PMC3182526/](https://pmc.ncbi.nlm.nih.gov/articles/PMC3182526/)

- **EDGE:** SOURCE: `SNAI1` TARGET: `CDH1` TYPE: `regulatory` CONFIDENCE: `strong`  
  EVIDENCE: Snail represses E-cadherin transcription via E-box binding (EMT hallmark).  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `no`
- **EDGE:** SOURCE: `SNAI2` TARGET: `ZEB1` TYPE: `regulatory` CONFIDENCE: `strong`  
  EVIDENCE: Slug directly activates ZEB1 transcription at E-boxes, cooperatively repressing CDH1.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `no`
- **EDGE:** SOURCE: `ZEB1` TARGET: `CDH1` TYPE: `regulatory` CONFIDENCE: `strong`  
  EVIDENCE: ZEB1 represses E-cadherin; master EMT regulator across carcinomas.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `no`

### Taylor2016_ZIP6_ZIP10_heterodimer

*ZIP10/ZIP6 heterodimerization drives EMT in luminal breast cancer (Nimmanon, Taylor et al., review)*  
[https://www.explorationpub.com/Journals/etat/Article/100280](https://www.explorationpub.com/Journals/etat/Article/100280)

- **EDGE:** SOURCE: `SLC39A10` TARGET: `SLC39A6` TYPE: `physical` CONFIDENCE: `strong`  
  EVIDENCE: ZIP10 heterodimerizes with ZIP6 at PM to drive EMT in luminal breast cancer.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `possible` · AU_RELEVANCE: `no`
- **EDGE:** SOURCE: `SLC39A10` TARGET: `SNAI1` TYPE: `signaling` CONFIDENCE: `moderate`  
  EVIDENCE: ZIP10/ZIP6 heterodimer activates Snail via GSK3beta inactivation.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `possible` · AU_RELEVANCE: `no`

### Li2010_ZIP4_CREB

*ZIP4 regulates pancreatic cancer cell growth by activating IL-6/STAT3 pathway through CREB (Li et al., Clin Cancer Res 2010)*  
[https://aacrjournals.org/clincancerres/article/16/5/1423/11181/](https://aacrjournals.org/clincancerres/article/16/5/1423/11181/)

- **EDGE:** SOURCE: `SLC39A4` TARGET: `CREB1` TYPE: `signaling` CONFIDENCE: `strong`  
  EVIDENCE: ZIP4-driven Zn rise activates Zn-dependent TF CREB.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `possible` · AU_RELEVANCE: `no`
- **EDGE:** SOURCE: `CREB1` TARGET: `IL6` TYPE: `regulatory` CONFIDENCE: `strong`  
  EVIDENCE: CREB activates IL-6 transcription downstream of ZIP4.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `possible` · AU_RELEVANCE: `no`
- **EDGE:** SOURCE: `IL6` TARGET: `STAT3` TYPE: `signaling` CONFIDENCE: `strong`  
  EVIDENCE: IL-6 activates STAT3 via JAK (canonical IL-6/STAT3 axis).  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `no`
- **EDGE:** SOURCE: `STAT3` TARGET: `CCND1` TYPE: `regulatory` CONFIDENCE: `strong`  
  EVIDENCE: STAT3 induces cyclin D1, accelerating G1/S transition.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `no`

### Zhang2013_ZIP4_miR373

*A novel epigenetic CREB-miR-373 axis mediates ZIP4-induced pancreatic cancer growth (Zhang et al., 2013)*  
[https://pmc.ncbi.nlm.nih.gov/articles/PMC3799489/](https://pmc.ncbi.nlm.nih.gov/articles/PMC3799489/)

- **EDGE:** SOURCE: `CREB1` TARGET: `MIR373` TYPE: `regulatory` CONFIDENCE: `strong`  
  EVIDENCE: CREB->miR-373 axis silences TP53INP1/LATS2/CD44 downstream of ZIP4.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `unknown` · AU_RELEVANCE: `no`
- **EDGE:** SOURCE: `MIR373` TARGET: `TP53INP1` TYPE: `regulatory` CONFIDENCE: `strong`  
  EVIDENCE: miR-373 represses TP53INP1 in pancreatic cancer.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `unknown` · AU_RELEVANCE: `no`

### Manso2025_ZnT_AKT_ESR1

*SLC30A1/5/9 transporters play crucial role in ligand-independent ESR1 activation via AKT (Manso et al., 2025)*  
[https://pmc.ncbi.nlm.nih.gov/articles/PMC13055888/](https://pmc.ncbi.nlm.nih.gov/articles/PMC13055888/)

- **EDGE:** SOURCE: `SLC30A1` TARGET: `AKT1` TYPE: `signaling` CONFIDENCE: `moderate`  
  EVIDENCE: ZnT1 modulates intracellular Zn that drives PTP/AKT/ESR1 signalling.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `unknown` · AU_RELEVANCE: `no`
- **EDGE:** SOURCE: `SLC30A1` TARGET: `ESR1` TYPE: `signaling` CONFIDENCE: `moderate`  
  EVIDENCE: ZnT1 contributes to ligand-independent ESR1 activation via AKT.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `unknown` · AU_RELEVANCE: `no`
- **EDGE:** SOURCE: `SLC30A5` TARGET: `AKT1` TYPE: `signaling` CONFIDENCE: `moderate`  
  EVIDENCE: ZnT5 cooperates with ZnT1 to modulate AKT-ESR1 in breast cancer.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `unknown` · AU_RELEVANCE: `no`

### Saydam2002_MTF1_targets

*Regulation of metallothionein transcription by MTF-1 / target genes including ZnT-1 (Saydam et al., JBC 2002)*  
[https://www.jbc.org/article/S0021-9258(20)84886-8/fulltext](https://www.jbc.org/article/S0021-9258(20)84886-8/fulltext)

- **EDGE:** SOURCE: `MTF1` TARGET: `SLC30A1` TYPE: `regulatory` CONFIDENCE: `strong`  
  EVIDENCE: MTF1 directly transactivates ZnT1 via MRE motifs.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `possible` · AU_RELEVANCE: `no`
- **EDGE:** SOURCE: `MTF1` TARGET: `CEBPA` TYPE: `regulatory` CONFIDENCE: `moderate`  
  EVIDENCE: MTF1 targets C/EBPalpha via MRE in its promoter.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `no`

### Lichten2009_SLC30_SLC39_review

*SLC30/ZnT and SLC39/ZIP transporter families review (Frontiers Immunol 2025)*  
[https://pmc.ncbi.nlm.nih.gov/articles/PMC12827705/](https://pmc.ncbi.nlm.nih.gov/articles/PMC12827705/)

- **EDGE:** SOURCE: `SLC30A7` TARGET: `MT2A` TYPE: `shared-pathway` CONFIDENCE: `moderate`  
  EVIDENCE: ZnT7 sequesters Zn in Golgi; loss rebalances MT2A buffering.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `unknown` · AU_RELEVANCE: `no`
- **EDGE:** SOURCE: `SLC30A10` TARGET: `MT2A` TYPE: `shared-pathway` CONFIDENCE: `speculative`  
  EVIDENCE: ZnT10 contributes to hepatic Zn efflux; coexpresses with MT family.  
  CANCER_RELEVANCE: `no` · LIHC_RELEVANCE: `possible` · AU_RELEVANCE: `no`
- **EDGE:** SOURCE: `SLC30A8` TARGET: `INS` TYPE: `physical` CONFIDENCE: `strong`  
  EVIDENCE: ZnT8 packages Zn into beta-cell insulin granules for hexamer formation.  
  CANCER_RELEVANCE: `no` · LIHC_RELEVANCE: `no` · AU_RELEVANCE: `no`

### Lichtlen2001_MTF1_MRE

*MTF-1: structure, function and regulation (Lichtlen & Schaffner, Bioessays 2001)*  
[https://pubmed.ncbi.nlm.nih.gov/11554446/](https://pubmed.ncbi.nlm.nih.gov/11554446/)

- **EDGE:** SOURCE: `MTF1` TARGET: `MT1A` TYPE: `regulatory` CONFIDENCE: `strong`  
  EVIDENCE: MTF1 binds MREs in MT1A promoter; activates upon Zn rise / oxidative stress.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `possible`
- **EDGE:** SOURCE: `MTF1` TARGET: `MT1E` TYPE: `regulatory` CONFIDENCE: `strong`  
  EVIDENCE: MTF1 transactivates MT1E via MRE motifs.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `possible`
- **EDGE:** SOURCE: `MTF1` TARGET: `MT1F` TYPE: `regulatory` CONFIDENCE: `strong`  
  EVIDENCE: MTF1 transactivates MT1F via MRE motifs.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `possible`
- **EDGE:** SOURCE: `MTF1` TARGET: `MT1G` TYPE: `regulatory` CONFIDENCE: `strong`  
  EVIDENCE: MTF1 transactivates MT1G via MRE motifs; MT1G is silenced in HCC.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `possible`
- **EDGE:** SOURCE: `MTF1` TARGET: `MT1H` TYPE: `regulatory` CONFIDENCE: `strong`  
  EVIDENCE: MTF1 binds MT1H promoter MREs; MT1H is a prognostic biomarker in LIHC.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `possible`
- **EDGE:** SOURCE: `MTF1` TARGET: `MT1M` TYPE: `regulatory` CONFIDENCE: `strong`  
  EVIDENCE: MTF1 transactivates MT1M via MRE motifs.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `possible`
- **EDGE:** SOURCE: `MTF1` TARGET: `MT1X` TYPE: `regulatory` CONFIDENCE: `strong`  
  EVIDENCE: MTF1 transactivates MT1X via MRE motifs.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `possible`
- **EDGE:** SOURCE: `MTF1` TARGET: `MT2A` TYPE: `regulatory` CONFIDENCE: `strong`  
  EVIDENCE: MTF1 is dominant inducer of MT2A; highest MTF1 affinity at MT2A MREs.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `possible`
- **EDGE:** SOURCE: `MTF1` TARGET: `HSF1` TYPE: `shared-pathway` CONFIDENCE: `moderate`  
  EVIDENCE: MTF1 and HSF1 are co-induced by oxidative/heat stress; converging stress response.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `possible` · AU_RELEVANCE: `possible`
- **EDGE:** SOURCE: `ZN_FREE` TARGET: `MTF1` TYPE: `signaling` CONFIDENCE: `strong`  
  EVIDENCE: Free cytosolic Zn rise drives MTF1 nuclear translocation and MRE binding.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `yes`

### Ostrakhovitch2006_MT_p53

*Interaction of metallothionein with tumor suppressor p53 (Ostrakhovitch et al., FEBS Lett 2006)*  
[https://febs.onlinelibrary.wiley.com/doi/abs/10.1016/j.febslet.2006.01.036](https://febs.onlinelibrary.wiley.com/doi/abs/10.1016/j.febslet.2006.01.036)

- **EDGE:** SOURCE: `MT1A` TARGET: `TP53` TYPE: `physical` CONFIDENCE: `moderate`  
  EVIDENCE: MTs exchange Zn with p53; loss of Zn misfolds the p53 DBD.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `possible`
- **EDGE:** SOURCE: `MT2A` TARGET: `TP53` TYPE: `physical` CONFIDENCE: `moderate`  
  EVIDENCE: MT2A donates/sequesters Zn from p53; modulates p53 conformation and DNA binding.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `possible`
- **EDGE:** SOURCE: `TP53` TARGET: `MT2A` TYPE: `regulatory` CONFIDENCE: `moderate`  
  EVIDENCE: p53 transcriptionally represses MT2A in breast cancer; MT2A modulates p53 Zn in turn.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `possible` · AU_RELEVANCE: `no`

### Wang2019_MT1G_p53

*MT1G serves as a tumor suppressor in HCC by interacting with p53 (Wang et al., Oncogenesis 2019)*  
[https://pmc.ncbi.nlm.nih.gov/articles/PMC6858331/](https://pmc.ncbi.nlm.nih.gov/articles/PMC6858331/)

- **EDGE:** SOURCE: `MT1G` TARGET: `TP53` TYPE: `physical` CONFIDENCE: `strong`  
  EVIDENCE: MT1G directly binds p53 and supplies Zn to stabilize its DBD; rescues p53 activity.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `possible`
- **EDGE:** SOURCE: `MT1G` TARGET: `MDM2` TYPE: `regulatory` CONFIDENCE: `strong`  
  EVIDENCE: MT1G suppresses MDM2 expression, stabilizing p53 in HCC.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `no`
- **EDGE:** SOURCE: `TP53` TARGET: `CDKN1A` TYPE: `regulatory` CONFIDENCE: `strong`  
  EVIDENCE: p53 transactivates p21 (CDKN1A) leading to cell-cycle arrest.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `no`
- **EDGE:** SOURCE: `TP53` TARGET: `BAX` TYPE: `regulatory` CONFIDENCE: `strong`  
  EVIDENCE: p53 transactivates BAX leading to apoptosis.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `no`

### Kanda2009_MT1G_methylation

*MT1G is silenced by DNA methylation and contributes to HCC pathogenesis (Kanda et al., 2009)*  
[https://pmc.ncbi.nlm.nih.gov/articles/PMC6096370/](https://pmc.ncbi.nlm.nih.gov/articles/PMC6096370/)

- **EDGE:** SOURCE: `MT1G` TARGET: `DNMT1` TYPE: `regulatory` CONFIDENCE: `strong`  
  EVIDENCE: MT1G is silenced by DNMT1-mediated promoter hypermethylation in HCC.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `no`
- **EDGE:** SOURCE: `DNMT1` TARGET: `MT1A` TYPE: `regulatory` CONFIDENCE: `strong`  
  EVIDENCE: DNMT1 methylates MT1A promoter silencing it in HCC.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `no`
- **EDGE:** SOURCE: `DNMT1` TARGET: `MT1G` TYPE: `regulatory` CONFIDENCE: `strong`  
  EVIDENCE: DNMT1 hypermethylates MT1G CpG island silencing it in HCC.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `no`

### Bicker1998_SOD1_CuZn

*Cu/Zn Superoxide dismutase (SOD1) — antioxidant overview*  
[https://chem.libretexts.org/Courses/Saint_Marys_College_Notre_Dame_IN/CHEM_342:_Bio-inorganic_Chemistry/Readings/Metals_in_Biological_Systems_(Saint_Mary's_College)/Antioxidant:_Cu_Zn_Superoxide_dismutase_(SOD1)](https://chem.libretexts.org/Courses/Saint_Marys_College_Notre_Dame_IN/CHEM_342:_Bio-inorganic_Chemistry/Readings/Metals_in_Biological_Systems_(Saint_Mary's_College)/Antioxidant:_Cu_Zn_Superoxide_dismutase_(SOD1))

- **EDGE:** SOURCE: `MT2A` TARGET: `SOD1` TYPE: `shared-pathway` CONFIDENCE: `moderate`  
  EVIDENCE: MT2A and SOD1 cooperate in ROS scavenging; Zn flux between MT and Cu/Zn-SOD has been proposed.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `unknown` · AU_RELEVANCE: `possible`

### Yang2019_ZNF479_DNMT1_MT

*ZNF479 downregulates MT-1 expression via ASH2L and DNMT1 in HCC (Yang et al., Cell Death Dis 2019)*  
[https://www.nature.com/articles/s41419-019-1651-9](https://www.nature.com/articles/s41419-019-1651-9)

- **EDGE:** SOURCE: `ZNF479` TARGET: `DNMT1` TYPE: `regulatory` CONFIDENCE: `strong`  
  EVIDENCE: ZNF479 induces DNMT1 in HCC; siDNMT1 restores MT1 levels.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `no`
- **EDGE:** SOURCE: `ZNF479` TARGET: `ASH2L` TYPE: `regulatory` CONFIDENCE: `strong`  
  EVIDENCE: ZNF479 upregulates ASH2L (MLL complex), increasing H3K4me3.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `no`
- **EDGE:** SOURCE: `ZNF479` TARGET: `MT1A` TYPE: `regulatory` CONFIDENCE: `strong`  
  EVIDENCE: ZNF479 suppresses MT1A via DNMT1/ASH2L in HCC.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `no`
- **EDGE:** SOURCE: `ZNF479` TARGET: `UHRF1` TYPE: `regulatory` CONFIDENCE: `moderate`  
  EVIDENCE: ZNF479 induces UHRF1, reinforcing maintenance methylation at MT1 loci.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `no`
- **EDGE:** SOURCE: `YWHAE` TARGET: `ZNF479` TYPE: `regulatory` CONFIDENCE: `strong`  
  EVIDENCE: 14-3-3epsilon induces ZNF479, initiating MT1-silencing cascade in HCC.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `no`

### Wong2007_MT_CEBPA_HCC

*MT expression is suppressed in HCC via inactivation of C/EBPalpha by PI3K signalling (Datta et al., Cancer Res 2007)*  
[https://pmc.ncbi.nlm.nih.gov/articles/PMC2276570/](https://pmc.ncbi.nlm.nih.gov/articles/PMC2276570/)

- **EDGE:** SOURCE: `CEBPA` TARGET: `MT1A` TYPE: `regulatory` CONFIDENCE: `strong`  
  EVIDENCE: C/EBPalpha activates MT1A; loss reduces MT in HCC.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `no`
- **EDGE:** SOURCE: `CEBPA` TARGET: `MT2A` TYPE: `regulatory` CONFIDENCE: `strong`  
  EVIDENCE: C/EBPalpha activates MT2A transcription; reduced in HCC.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `no`
- **EDGE:** SOURCE: `PIK3CA` TARGET: `CEBPA` TYPE: `signaling` CONFIDENCE: `strong`  
  EVIDENCE: PI3K signalling inactivates C/EBPalpha, suppressing MT expression in HCC.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `no`

### Zhao2020_ZNF143_HCC_CDC6

*ZNF143-mediated H3K9 trimethylation upregulates CDC6 by activating MDIG in HCC (Zhao et al., Cancer Res 2020)*  
[https://aacrjournals.org/cancerres/article/80/12/2599/641068/](https://aacrjournals.org/cancerres/article/80/12/2599/641068/)

- **EDGE:** SOURCE: `ZNF143` TARGET: `CDC6` TYPE: `regulatory` CONFIDENCE: `strong`  
  EVIDENCE: ZNF143 binds CDC6 locus via MDIG-mediated H3K9 demethylation in HCC.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `no`
- **EDGE:** SOURCE: `ZNF143` TARGET: `MINA` TYPE: `regulatory` CONFIDENCE: `strong`  
  EVIDENCE: ZNF143 transactivates MDIG (MINA53), a Zn-binding histone demethylase, in HCC.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `no`

### ZNF143_MEX3C_HCC

*ZNF143-mediated upregulation of MEX3C promotes HCC progression (2024)*  
[https://www.sciencedirect.com/science/article/pii/S2210740124002134](https://www.sciencedirect.com/science/article/pii/S2210740124002134)

- **EDGE:** SOURCE: `ZNF143` TARGET: `MEX3C` TYPE: `regulatory` CONFIDENCE: `strong`  
  EVIDENCE: ZNF143 binds MEX3C promoter; promotes HCC migration/invasion.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `no`

### Kobayashi2021_NRF2_KEAP1

*KEAP1/NRF2 pathway under oxidative and electrophilic stress (review)*  
[https://pmc.ncbi.nlm.nih.gov/articles/PMC3820647/](https://pmc.ncbi.nlm.nih.gov/articles/PMC3820647/)

- **EDGE:** SOURCE: `KEAP1` TARGET: `NFE2L2` TYPE: `regulatory` CONFIDENCE: `strong`  
  EVIDENCE: KEAP1 targets NRF2 for ubiquitin-proteasome degradation under basal conditions.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `possible`
- **EDGE:** SOURCE: `NFE2L2` TARGET: `GPX4` TYPE: `regulatory` CONFIDENCE: `strong`  
  EVIDENCE: NRF2 transactivates GPX4 and other antioxidant genes via ARE motifs.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `possible`
- **EDGE:** SOURCE: `NFE2L2` TARGET: `TXNRD1` TYPE: `regulatory` CONFIDENCE: `strong`  
  EVIDENCE: NRF2 induces TXNRD1 via ARE motifs, sustaining the thioredoxin cycle.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `yes`
- **EDGE:** SOURCE: `NFE2L2` TARGET: `SLC7A11` TYPE: `regulatory` CONFIDENCE: `strong`  
  EVIDENCE: NRF2 induces SLC7A11 (xCT), limiting ferroptosis.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `possible`

### DalleDonne2014_NRF2_HSF1

*NRF2 transcriptionally activates HSF1 promoter under oxidative stress (Dayalan-Naidu et al., JBC 2018)*  
[https://pmc.ncbi.nlm.nih.gov/articles/PMC6302185/](https://pmc.ncbi.nlm.nih.gov/articles/PMC6302185/)

- **EDGE:** SOURCE: `NFE2L2` TARGET: `HSF1` TYPE: `regulatory` CONFIDENCE: `strong`  
  EVIDENCE: NRF2 binds AREs in HSF1 promoter; activates HSF1 under oxidative stress.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `possible` · AU_RELEVANCE: `possible`
- **EDGE:** SOURCE: `HSF1` TARGET: `HSPA1A` TYPE: `regulatory` CONFIDENCE: `strong`  
  EVIDENCE: HSF1 transactivates HSP70 (HSPA1A) under heat / oxidative stress.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `no`

### Roder2018_Auranofin_TrxR

*The gold complex auranofin: new perspectives for cancer therapy (Roder & Thomson, Discov Oncol 2022)*  
[https://pmc.ncbi.nlm.nih.gov/articles/PMC8777575/](https://pmc.ncbi.nlm.nih.gov/articles/PMC8777575/)

- **EDGE:** SOURCE: `AURANOFIN` TARGET: `TXNRD1` TYPE: `physical` CONFIDENCE: `strong`  
  EVIDENCE: Auranofin covalently inhibits TrxR1 via Au-Se bond at active-site selenocysteine.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `yes`
- **EDGE:** SOURCE: `AURANOFIN` TARGET: `TXNRD2` TYPE: `physical` CONFIDENCE: `strong`  
  EVIDENCE: Auranofin also inhibits mitochondrial TXNRD2, dysregulating redox state.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `yes`
- **EDGE:** SOURCE: `AURANOFIN` TARGET: `GPX4` TYPE: `signaling` CONFIDENCE: `moderate`  
  EVIDENCE: Auranofin depletes GSH/TrxR axis, sensitizing cells to GPX4-dependent ferroptosis.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `yes`

### Mendes2018_Auranofin_PARP1

*Auranofin synergizes with PARP inhibitor olaparib in mutant-p53 cancers (Mendes et al., 2023)*  
[https://pmc.ncbi.nlm.nih.gov/articles/PMC10045521/](https://pmc.ncbi.nlm.nih.gov/articles/PMC10045521/)

- **EDGE:** SOURCE: `AURANOFIN` TARGET: `PARP1` TYPE: `physical` CONFIDENCE: `strong`  
  EVIDENCE: Au(I) drugs displace Zn from PARP-1 zinc fingers; nM IC50 for DNA-dependent activity.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `possible` · AU_RELEVANCE: `yes`

### Laib1985_AuMT

*The binding of Gold(I) to metallothionein (Laib et al., 1985)*  
[https://pubmed.ncbi.nlm.nih.gov/7411139/](https://pubmed.ncbi.nlm.nih.gov/7411139/)

- **EDGE:** SOURCE: `AURANOFIN` TARGET: `MT2A` TYPE: `physical` CONFIDENCE: `strong`  
  EVIDENCE: Au(I) displaces Zn/Cd from metallothionein; up to 9 Au per MT vs 7 Zn.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `yes`
- **EDGE:** SOURCE: `AURANOFIN` TARGET: `MTF1` TYPE: `signaling` CONFIDENCE: `speculative`  
  EVIDENCE: Au-induced Zn displacement from MT releases free Zn that may activate MTF1.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `yes`

### Bordin1996_AuCd_MT

*Gold replacement of cadmium, zinc-binding metallothionein (Bordin et al., 1996)*  
[https://pubmed.ncbi.nlm.nih.gov/8865374/](https://pubmed.ncbi.nlm.nih.gov/8865374/)

- **EDGE:** SOURCE: `AURANOFIN` TARGET: `MT1A` TYPE: `physical` CONFIDENCE: `strong`  
  EVIDENCE: Au(I) binds MT1 thiolates and displaces Zn2+ under stoichiometric excess.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `yes`

### Spell2016_Au_ZnFinger

*Reactivity of Cys4 zinc-finger domains with gold(III) complexes — 'gold fingers' (Spell & Farver, Inorg Chem 2015)*  
[https://pubs.acs.org/doi/abs/10.1021/acs.inorgchem.5b00360](https://pubs.acs.org/doi/abs/10.1021/acs.inorgchem.5b00360)

- **EDGE:** SOURCE: `AURANOFIN` TARGET: `TP53` TYPE: `physical` CONFIDENCE: `moderate`  
  EVIDENCE: Au(III) attacks Cys4 Zn-finger sites in vivo, forming 'gold fingers'; destabilizes p53 DBD.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `yes`
- **EDGE:** SOURCE: `AURANOFIN` TARGET: `ZN_FREE` TYPE: `physical` CONFIDENCE: `strong`  
  EVIDENCE: Au(I/III) displaces Zn2+ from Cys-rich coordination sites (thiophilic substitution).  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `yes`

### Maret2005_thionein_Zn

*Control of zinc transfer between thionein, metallothionein, and zinc proteins (Jacob, Maret, Vallee, PNAS 1998)*  
[https://pmc.ncbi.nlm.nih.gov/articles/PMC19863/](https://pmc.ncbi.nlm.nih.gov/articles/PMC19863/)

- **EDGE:** SOURCE: `MT1A` TARGET: `ZN_FREE` TYPE: `physical` CONFIDENCE: `strong`  
  EVIDENCE: MT binds 7 Zn2+ via 20 cysteines; redox-responsive Zn release feeds free Zn pool.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `yes`
- **EDGE:** SOURCE: `MT2A` TARGET: `ZN_FREE` TYPE: `physical` CONFIDENCE: `strong`  
  EVIDENCE: MT2A is the principal cytosolic Zn buffer; redox state controls Zn release.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `yes`

### Loh2022_p53_zinc

*p53 and Zinc: a malleable relationship (Loh, Front Mol Biosci 2022)*  
[https://www.frontiersin.org/journals/molecular-biosciences/articles/10.3389/fmolb.2022.895887/full](https://www.frontiersin.org/journals/molecular-biosciences/articles/10.3389/fmolb.2022.895887/full)

- **EDGE:** SOURCE: `ZN_FREE` TARGET: `TP53` TYPE: `physical` CONFIDENCE: `strong`  
  EVIDENCE: Free Zn2+ is required for p53 DBD folding; Zn loss yields mutant-like conformation.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `yes`

### Cassandri2017_ZNF_review

*Zinc-finger proteins in health and disease (Cassandri et al., Cell Death Discov 2017)*  
[https://www.nature.com/articles/cddiscovery201771](https://www.nature.com/articles/cddiscovery201771)

- **EDGE:** SOURCE: `ZEB1` TARGET: `ZEB2` TYPE: `shared-pathway` CONFIDENCE: `moderate`  
  EVIDENCE: ZEB1 and ZEB2 cooperatively repress CDH1 and other epithelial genes.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `no`
- **EDGE:** SOURCE: `SNAI1` TARGET: `SNAI2` TYPE: `shared-pathway` CONFIDENCE: `moderate`  
  EVIDENCE: Snail and Slug share targets and cooperatively drive EMT and metastasis.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `no`
- **EDGE:** SOURCE: `CTCF` TARGET: `TP53` TYPE: `regulatory` CONFIDENCE: `moderate`  
  EVIDENCE: CTCF (11-ZF TF) regulates TP53 locus chromatin; ZF mutations alter DNA binding.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `yes`
- **EDGE:** SOURCE: `ZFX` TARGET: `MYC` TYPE: `regulatory` CONFIDENCE: `moderate`  
  EVIDENCE: ZF TF ZFX binds CpG-island promoters of MYC and other oncogenes.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `no`
- **EDGE:** SOURCE: `ZNF740` TARGET: `METTL3` TYPE: `regulatory` CONFIDENCE: `moderate`  
  EVIDENCE: ZNF740 activates METTL3/HIF-1A axis in HCC.  
  CANCER_RELEVANCE: `yes` · LIHC_RELEVANCE: `yes` · AU_RELEVANCE: `no`

---

## 2. Deduplicated edge list

| # | Source | Target | Type | Confidence | Cancer | LIHC | Au | Times cited |
|---|--------|--------|------|------------|--------|------|----|-------------|
| 1 | `CSNK2A1` | `ZIP7` | regulatory | strong | yes | possible | no | 1 |
| 2 | `ZIP7` | `PTP1B` | signaling | strong | yes | possible | no | 1 |
| 3 | `ZIP7` | `AKT1` | signaling | strong | yes | possible | no | 1 |
| 4 | `ZIP7` | `MAPK1` | signaling | strong | yes | possible | no | 1 |
| 5 | `ZIP7` | `EGFR` | signaling | moderate | yes | unknown | no | 1 |
| 6 | `ZIP7` | `GPX4` | regulatory | strong | yes | possible | possible | 1 |
| 7 | `ZIP7` | `DDIT3` | regulatory | strong | yes | possible | possible | 1 |
| 8 | `CEBPA` | `SLC39A14` | regulatory | moderate | yes | yes | no | 1 |
| 9 | `SLC39A14` | `ATF4` | regulatory | strong | yes | yes | no | 1 |
| 10 | `SLC39A14` | `DDIT3` | regulatory | strong | yes | yes | no | 1 |
| 11 | `SLC39A14` | `MT1A` | regulatory | strong | yes | yes | no | 1 |
| 12 | `STAT3` | `SLC39A6` | regulatory | strong | yes | possible | no | 1 |
| 13 | `SLC39A6` | `GSK3B` | signaling | strong | yes | possible | no | 1 |
| 14 | `SLC39A6` | `SNAI1` | regulatory | strong | yes | possible | no | 1 |
| 15 | `SNAI1` | `CDH1` | regulatory | strong | yes | yes | no | 1 |
| 16 | `SLC39A10` | `SLC39A6` | physical | strong | yes | possible | no | 1 |
| 17 | `SLC39A10` | `SNAI1` | signaling | moderate | yes | possible | no | 1 |
| 18 | `SNAI2` | `ZEB1` | regulatory | strong | yes | yes | no | 1 |
| 19 | `ZEB1` | `CDH1` | regulatory | strong | yes | yes | no | 1 |
| 20 | `SLC39A4` | `CREB1` | signaling | strong | yes | possible | no | 1 |
| 21 | `CREB1` | `IL6` | regulatory | strong | yes | possible | no | 1 |
| 22 | `IL6` | `STAT3` | signaling | strong | yes | yes | no | 1 |
| 23 | `STAT3` | `CCND1` | regulatory | strong | yes | yes | no | 1 |
| 24 | `CREB1` | `MIR373` | regulatory | strong | yes | unknown | no | 1 |
| 25 | `MIR373` | `TP53INP1` | regulatory | strong | yes | unknown | no | 1 |
| 26 | `SLC30A1` | `AKT1` | signaling | moderate | yes | unknown | no | 1 |
| 27 | `SLC30A1` | `ESR1` | signaling | moderate | yes | unknown | no | 1 |
| 28 | `MTF1` | `SLC30A1` | regulatory | strong | yes | possible | no | 1 |
| 29 | `SLC30A5` | `AKT1` | signaling | moderate | yes | unknown | no | 1 |
| 30 | `SLC30A7` | `MT2A` | shared-pathway | moderate | yes | unknown | no | 1 |
| 31 | `SLC30A10` | `MT2A` | shared-pathway | speculative | no | possible | no | 1 |
| 32 | `SLC30A8` | `INS` | physical | strong | no | no | no | 1 |
| 33 | `MTF1` | `MT1A` | regulatory | strong | yes | yes | possible | 1 |
| 34 | `MTF1` | `MT1E` | regulatory | strong | yes | yes | possible | 1 |
| 35 | `MTF1` | `MT1F` | regulatory | strong | yes | yes | possible | 1 |
| 36 | `MTF1` | `MT1G` | regulatory | strong | yes | yes | possible | 1 |
| 37 | `MTF1` | `MT1H` | regulatory | strong | yes | yes | possible | 1 |
| 38 | `MTF1` | `MT1M` | regulatory | strong | yes | yes | possible | 1 |
| 39 | `MTF1` | `MT1X` | regulatory | strong | yes | yes | possible | 1 |
| 40 | `MTF1` | `MT2A` | regulatory | strong | yes | yes | possible | 1 |
| 41 | `MTF1` | `CEBPA` | regulatory | moderate | yes | yes | no | 1 |
| 42 | `MT1A` | `TP53` | physical | moderate | yes | yes | possible | 1 |
| 43 | `MT2A` | `TP53` | physical | moderate | yes | yes | possible | 1 |
| 44 | `MT1G` | `TP53` | physical | strong | yes | yes | possible | 1 |
| 45 | `MT1G` | `MDM2` | regulatory | strong | yes | yes | no | 1 |
| 46 | `TP53` | `CDKN1A` | regulatory | strong | yes | yes | no | 1 |
| 47 | `TP53` | `BAX` | regulatory | strong | yes | yes | no | 1 |
| 48 | `MT1G` | `DNMT1` | regulatory | strong | yes | yes | no | 1 |
| 49 | `MT2A` | `SOD1` | shared-pathway | moderate | yes | unknown | possible | 2 |
| 50 | `ZNF479` | `DNMT1` | regulatory | strong | yes | yes | no | 1 |
| 51 | `ZNF479` | `ASH2L` | regulatory | strong | yes | yes | no | 1 |
| 52 | `ZNF479` | `MT1A` | regulatory | strong | yes | yes | no | 1 |
| 53 | `ZNF479` | `UHRF1` | regulatory | moderate | yes | yes | no | 1 |
| 54 | `YWHAE` | `ZNF479` | regulatory | strong | yes | yes | no | 1 |
| 55 | `DNMT1` | `MT1A` | regulatory | strong | yes | yes | no | 1 |
| 56 | `DNMT1` | `MT1G` | regulatory | strong | yes | yes | no | 1 |
| 57 | `CEBPA` | `MT1A` | regulatory | strong | yes | yes | no | 1 |
| 58 | `CEBPA` | `MT2A` | regulatory | strong | yes | yes | no | 1 |
| 59 | `PIK3CA` | `CEBPA` | signaling | strong | yes | yes | no | 1 |
| 60 | `ZNF143` | `CDC6` | regulatory | strong | yes | yes | no | 1 |
| 61 | `ZNF143` | `MINA` | regulatory | strong | yes | yes | no | 1 |
| 62 | `ZNF143` | `MEX3C` | regulatory | strong | yes | yes | no | 1 |
| 63 | `TP53` | `MT2A` | regulatory | moderate | yes | possible | no | 1 |
| 64 | `KEAP1` | `NFE2L2` | regulatory | strong | yes | yes | possible | 1 |
| 65 | `NFE2L2` | `HSF1` | regulatory | strong | yes | possible | possible | 1 |
| 66 | `HSF1` | `HSPA1A` | regulatory | strong | yes | yes | no | 1 |
| 67 | `NFE2L2` | `GPX4` | regulatory | strong | yes | yes | possible | 1 |
| 68 | `NFE2L2` | `TXNRD1` | regulatory | strong | yes | yes | yes | 1 |
| 69 | `NFE2L2` | `SLC7A11` | regulatory | strong | yes | yes | possible | 1 |
| 70 | `MTF1` | `HSF1` | shared-pathway | moderate | yes | possible | possible | 1 |
| 71 | `GPX4` | `ACSL4` | shared-pathway | strong | yes | yes | possible | 1 |
| 72 | `AURANOFIN` | `TXNRD1` | physical | strong | yes | yes | yes | 1 |
| 73 | `AURANOFIN` | `TXNRD2` | physical | strong | yes | yes | yes | 1 |
| 74 | `AURANOFIN` | `PARP1` | physical | strong | yes | possible | yes | 1 |
| 75 | `AURANOFIN` | `MT2A` | physical | strong | yes | yes | yes | 1 |
| 76 | `AURANOFIN` | `MT1A` | physical | strong | yes | yes | yes | 1 |
| 77 | `AURANOFIN` | `TP53` | physical | moderate | yes | yes | yes | 1 |
| 78 | `AURANOFIN` | `ZN_FREE` | physical | strong | yes | yes | yes | 1 |
| 79 | `AURANOFIN` | `MTF1` | signaling | speculative | yes | yes | yes | 1 |
| 80 | `AURANOFIN` | `GPX4` | signaling | moderate | yes | yes | yes | 1 |
| 81 | `MT1A` | `ZN_FREE` | physical | strong | yes | yes | yes | 1 |
| 82 | `MT2A` | `ZN_FREE` | physical | strong | yes | yes | yes | 1 |
| 83 | `ZN_FREE` | `MTF1` | signaling | strong | yes | yes | yes | 1 |
| 84 | `ZN_FREE` | `PTP1B` | signaling | strong | yes | possible | yes | 1 |
| 85 | `ZN_FREE` | `PTEN` | signaling | moderate | yes | possible | yes | 1 |
| 86 | `ZN_FREE` | `TP53` | physical | strong | yes | yes | yes | 1 |
| 87 | `ZEB1` | `ZEB2` | shared-pathway | moderate | yes | yes | no | 1 |
| 88 | `SNAI1` | `SNAI2` | shared-pathway | moderate | yes | yes | no | 1 |
| 89 | `CTCF` | `TP53` | regulatory | moderate | yes | yes | yes | 1 |
| 90 | `ZFX` | `MYC` | regulatory | moderate | yes | yes | no | 1 |
| 91 | `ZNF740` | `METTL3` | regulatory | moderate | yes | yes | no | 1 |

---

## 3. Edge frequency across the literature

Edges supported by more than one independent paper:

_All edges currently anchored by a single primary reference; multi-paper consensus emerges at the node-level hub analysis below._


---

## 4. Most central recurring Zn hubs (top 20 by degree)

| Rank | Node | Degree | Role |
|---|---|---|---|
| 1 | `MTF1` | 13 | Master zinc-responsive transcription factor; activates MT family and ZnT1. |
| 2 | `MT2A` | 9 | Principal cytosolic Zn buffer; central in Zn exchange with p53, SOD1, and Au compounds. |
| 3 | `TP53` | 9 | Zn-dependent tumour suppressor; structural Zn binds C176/H179/C238/C242 DBD. |
| 4 | `AURANOFIN` | 9 | Au(I) drug — anchor for every Au-Zn displacement edge. |
| 5 | `MT1A` | 8 | Inducible MT1 isoform; MTF1 target; silenced in HCC. |
| 6 | `ZIP7` | 7 | CK2-activated ER Zn channel driving kinase signalling and ferroptosis. |
| 7 | `ZN_FREE` | 7 | Labile Zn²⁺ pool — the biochemical currency of the graph. |
| 8 | `CEBPA` | 5 | Hepatic MT-activating TF; inactivated by PI3K in HCC. |
| 9 | `MT1G` | 5 | HCC tumour suppressor; donates structural Zn to p53. |
| 10 | `ZNF479` | 5 | HCC ZF-TF; silences MT1 family via DNMT1/ASH2L. |
| 11 | `NFE2L2` | 5 | Master antioxidant TF; activates GPX4, TXNRD1, SLC7A11, HSF1. |
| 12 | `GPX4` | 4 | Phospholipid hydroperoxidase; central ferroptosis defence. |
| 13 | `SLC39A14` | 4 | Hepatic ZIP14; downregulated early in HCC; required for ER-stress adaptation. |
| 14 | `SLC39A6` | 4 | STAT3-driven EMT channel. |
| 15 | `SNAI1` | 4 | EMT master TF; canonical CDH1 repressor. |
| 16 | `DNMT1` | 4 | Epigenetic silencer of MT1 in HCC. |
| 17 | `AKT1` | 3 | Convergent kinase driven by Zn-inhibited phosphatases. |
| 18 | `STAT3` | 3 | Drives ZIP6 transcription; activated by IL-6 downstream of ZIP4-CREB. |
| 19 | `ZEB1` | 3 | EMT regulator activated by Slug. |
| 20 | `CREB1` | 3 | Zn-dependent bZIP TF; ZIP4 → IL-6/STAT3 in pancreatic cancer. |

---

## 5. Suggested canonical zinc-biology subnetwork for the GNN

A compact module-organized subgraph (~50 high-confidence edges, covering all priority topics). Use this as the **core scaffold** for the GNN; the full edge list (Section 2) can be loaded as a richer overlay.

**Module 1 — Zn import & ER stress**
- `CEBPA` → `SLC39A14`
- `SLC39A14` → `ATF4`
- `SLC39A14` → `DDIT3`
- `SLC39A14` → `MT1A`
- `ZIP7` → `DDIT3`

**Module 2 — ZIP7 Zn-signalling hub**
- `CSNK2A1` → `ZIP7`
- `ZIP7` → `PTP1B`
- `ZIP7` → `AKT1`
- `ZIP7` → `MAPK1`
- `ZN_FREE` → `PTP1B`
- `ZN_FREE` → `PTEN`

**Module 3 — EMT / metastasis**
- `STAT3` → `SLC39A6`
- `SLC39A6` → `GSK3B`
- `SLC39A6` → `SNAI1`
- `SLC39A10` → `SLC39A6`
- `SNAI1` → `CDH1`
- `SNAI2` → `ZEB1`
- `ZEB1` → `CDH1`

**Module 4 — MTF1 / MT core**
- `ZN_FREE` → `MTF1`
- `MTF1` → `MT2A`
- `MTF1` → `MT1A`
- `MTF1` → `MT1G`
- `MTF1` → `MT1H`
- `MTF1` → `SLC30A1`

**Module 5 — MT / p53 axis (LIHC)**
- `MT1G` → `TP53`
- `MT1G` → `MDM2`
- `TP53` → `CDKN1A`
- `TP53` → `BAX`
- `ZN_FREE` → `TP53`
- `MT2A` → `TP53`

**Module 6 — MT silencing in HCC**
- `YWHAE` → `ZNF479`
- `ZNF479` → `DNMT1`
- `ZNF479` → `ASH2L`
- `DNMT1` → `MT1G`
- `PIK3CA` → `CEBPA`
- `CEBPA` → `MT1A`

**Module 7 — Oxidative / ferroptosis**
- `KEAP1` → `NFE2L2`
- `NFE2L2` → `GPX4`
- `NFE2L2` → `TXNRD1`
- `NFE2L2` → `SLC7A11`
- `NFE2L2` → `HSF1`
- `ZIP7` → `GPX4`
- `GPX4` → `ACSL4`
- `SOD1` → `MT2A`

**Module 8 — Au-induced Zn displacement**
- `AURANOFIN` → `MT2A`
- `AURANOFIN` → `MT1A`
- `AURANOFIN` → `PARP1`
- `AURANOFIN` → `TP53`
- `AURANOFIN` → `TXNRD1`
- `AURANOFIN` → `TXNRD2`
- `AURANOFIN` → `ZN_FREE`
- `AURANOFIN` → `MTF1`
- `AURANOFIN` → `GPX4`

---

## 6. References (all 36 sources)

1. **Taylor2008_ZIP7_BC** — ZIP7-mediated intracellular zinc transport contributes to aberrant growth-factor signalling in antihormone-resistant breast cancer cells (Taylor et al., Endocrinology 2008). [https://academic.oup.com/endo/article/149/10/4912/2455137](https://academic.oup.com/endo/article/149/10/4912/2455137)
2. **Hogstrand2009_ZIP7_hub** — Zinc transporters and cancer: a potential role for ZIP7 as a hub for tyrosine kinase activation (Hogstrand et al., Trends Mol Med 2009). [https://pubmed.ncbi.nlm.nih.gov/19246244/](https://pubmed.ncbi.nlm.nih.gov/19246244/)
3. **Chen2021_ZIP7_ferroptosis** — Zinc transporter ZIP7 is a novel determinant of ferroptosis (Chen et al., Cell Death Dis 2021). [https://www.nature.com/articles/s41419-021-03482-5](https://www.nature.com/articles/s41419-021-03482-5)
4. **Ohashi2014_ZIP7_ER** — SLC39A7/ZIP7 promotes intestinal epithelial self-renewal by resolving ER stress (Ohashi et al., PLOS Genet 2016). [https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1006349](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1006349)
5. **Franklin2012_ZIP14_HCC** — ZIP14 zinc transporter downregulation and zinc depletion in development and progression of hepatocellular cancer (Franklin et al., J Gastrointest Cancer 2012). [https://pmc.ncbi.nlm.nih.gov/articles/PMC3724761/](https://pmc.ncbi.nlm.nih.gov/articles/PMC3724761/)
6. **Kim2017_ZIP14_ER** — Hepatic ZIP14-mediated zinc transport is required for adaptation to ER stress (Kim et al., PNAS 2017). [https://www.pnas.org/doi/10.1073/pnas.1704012114](https://www.pnas.org/doi/10.1073/pnas.1704012114)
7. **Hogstrand2013_ZIP6_STAT3_EMT** — A mechanism for EMT and anoikis resistance triggered by zinc channel ZIP6 and STAT3 (Hogstrand et al., Biochem J 2013). [https://portlandpress.com/biochemj/article/455/2/229/81664/](https://portlandpress.com/biochemj/article/455/2/229/81664/)
8. **Taylor2016_ZIP6_ZIP10_heterodimer** — ZIP10/ZIP6 heterodimerization drives EMT in luminal breast cancer (Nimmanon, Taylor et al., review). [https://www.explorationpub.com/Journals/etat/Article/100280](https://www.explorationpub.com/Journals/etat/Article/100280)
9. **Li2010_ZIP4_CREB** — ZIP4 regulates pancreatic cancer cell growth by activating IL-6/STAT3 pathway through CREB (Li et al., Clin Cancer Res 2010). [https://aacrjournals.org/clincancerres/article/16/5/1423/11181/](https://aacrjournals.org/clincancerres/article/16/5/1423/11181/)
10. **Zhang2013_ZIP4_miR373** — A novel epigenetic CREB-miR-373 axis mediates ZIP4-induced pancreatic cancer growth (Zhang et al., 2013). [https://pmc.ncbi.nlm.nih.gov/articles/PMC3799489/](https://pmc.ncbi.nlm.nih.gov/articles/PMC3799489/)
11. **Lichten2009_SLC30_SLC39_review** — SLC30/ZnT and SLC39/ZIP transporter families review (Frontiers Immunol 2025). [https://pmc.ncbi.nlm.nih.gov/articles/PMC12827705/](https://pmc.ncbi.nlm.nih.gov/articles/PMC12827705/)
12. **Manso2025_ZnT_AKT_ESR1** — SLC30A1/5/9 transporters play crucial role in ligand-independent ESR1 activation via AKT (Manso et al., 2025). [https://pmc.ncbi.nlm.nih.gov/articles/PMC13055888/](https://pmc.ncbi.nlm.nih.gov/articles/PMC13055888/)
13. **Lichtlen2001_MTF1_MRE** — MTF-1: structure, function and regulation (Lichtlen & Schaffner, Bioessays 2001). [https://pubmed.ncbi.nlm.nih.gov/11554446/](https://pubmed.ncbi.nlm.nih.gov/11554446/)
14. **Saydam2002_MTF1_targets** — Regulation of metallothionein transcription by MTF-1 / target genes including ZnT-1 (Saydam et al., JBC 2002). [https://www.jbc.org/article/S0021-9258(20)84886-8/fulltext](https://www.jbc.org/article/S0021-9258(20)84886-8/fulltext)
15. **Wong2007_MT_CEBPA_HCC** — MT expression is suppressed in HCC via inactivation of C/EBPalpha by PI3K signalling (Datta et al., Cancer Res 2007). [https://pmc.ncbi.nlm.nih.gov/articles/PMC2276570/](https://pmc.ncbi.nlm.nih.gov/articles/PMC2276570/)
16. **Kanda2009_MT1G_methylation** — MT1G is silenced by DNA methylation and contributes to HCC pathogenesis (Kanda et al., 2009). [https://pmc.ncbi.nlm.nih.gov/articles/PMC6096370/](https://pmc.ncbi.nlm.nih.gov/articles/PMC6096370/)
17. **Wang2019_MT1G_p53** — MT1G serves as a tumor suppressor in HCC by interacting with p53 (Wang et al., Oncogenesis 2019). [https://pmc.ncbi.nlm.nih.gov/articles/PMC6858331/](https://pmc.ncbi.nlm.nih.gov/articles/PMC6858331/)
18. **Yang2019_ZNF479_DNMT1_MT** — ZNF479 downregulates MT-1 expression via ASH2L and DNMT1 in HCC (Yang et al., Cell Death Dis 2019). [https://www.nature.com/articles/s41419-019-1651-9](https://www.nature.com/articles/s41419-019-1651-9)
19. **Ji2021_MT1H_HCC** — Integrative analysis identifies MT1H as candidate prognostic biomarker in HCC (Ji et al., Front Mol Biosci 2021). [https://pmc.ncbi.nlm.nih.gov/articles/PMC8523949/](https://pmc.ncbi.nlm.nih.gov/articles/PMC8523949/)
20. **Ostrakhovitch2006_MT_p53** — Interaction of metallothionein with tumor suppressor p53 (Ostrakhovitch et al., FEBS Lett 2006). [https://febs.onlinelibrary.wiley.com/doi/abs/10.1016/j.febslet.2006.01.036](https://febs.onlinelibrary.wiley.com/doi/abs/10.1016/j.febslet.2006.01.036)
21. **Loh2022_p53_zinc** — p53 and Zinc: a malleable relationship (Loh, Front Mol Biosci 2022). [https://www.frontiersin.org/journals/molecular-biosciences/articles/10.3389/fmolb.2022.895887/full](https://www.frontiersin.org/journals/molecular-biosciences/articles/10.3389/fmolb.2022.895887/full)
22. **Krezel2014_PTP1B_Zn** — Zinc ions modulate protein tyrosine phosphatase 1B activity (Bellomo, Krezel, Maret, Metallomics 2014). [https://pubs.rsc.org/en/content/articlehtml/2014/mt/c4mt00086b](https://pubs.rsc.org/en/content/articlehtml/2014/mt/c4mt00086b)
23. **Zhao2020_ZNF143_HCC_CDC6** — ZNF143-mediated H3K9 trimethylation upregulates CDC6 by activating MDIG in HCC (Zhao et al., Cancer Res 2020). [https://aacrjournals.org/cancerres/article/80/12/2599/641068/](https://aacrjournals.org/cancerres/article/80/12/2599/641068/)
24. **ZNF143_MEX3C_HCC** — ZNF143-mediated upregulation of MEX3C promotes HCC progression (2024). [https://www.sciencedirect.com/science/article/pii/S2210740124002134](https://www.sciencedirect.com/science/article/pii/S2210740124002134)
25. **Cassandri2017_ZNF_review** — Zinc-finger proteins in health and disease (Cassandri et al., Cell Death Discov 2017). [https://www.nature.com/articles/cddiscovery201771](https://www.nature.com/articles/cddiscovery201771)
26. **Wei2013_Slug_ZEB1_EMT** — Transcriptional activation of ZEB1 by Slug leads to cooperative regulation of the EMT-like phenotype (Wels et al., 2011). [https://pmc.ncbi.nlm.nih.gov/articles/PMC3182526/](https://pmc.ncbi.nlm.nih.gov/articles/PMC3182526/)
27. **Bicker1998_SOD1_CuZn** — Cu/Zn Superoxide dismutase (SOD1) — antioxidant overview. [https://chem.libretexts.org/Courses/Saint_Marys_College_Notre_Dame_IN/CHEM_342:_Bio-inorganic_Chemistry/Readings/Metals_in_Biological_Systems_(Saint_Mary's_College)/Antioxidant:_Cu_Zn_Superoxide_dismutase_(SOD1)](https://chem.libretexts.org/Courses/Saint_Marys_College_Notre_Dame_IN/CHEM_342:_Bio-inorganic_Chemistry/Readings/Metals_in_Biological_Systems_(Saint_Mary's_College)/Antioxidant:_Cu_Zn_Superoxide_dismutase_(SOD1))
28. **DalleDonne2014_NRF2_HSF1** — NRF2 transcriptionally activates HSF1 promoter under oxidative stress (Dayalan-Naidu et al., JBC 2018). [https://pmc.ncbi.nlm.nih.gov/articles/PMC6302185/](https://pmc.ncbi.nlm.nih.gov/articles/PMC6302185/)
29. **Kobayashi2021_NRF2_KEAP1** — KEAP1/NRF2 pathway under oxidative and electrophilic stress (review). [https://pmc.ncbi.nlm.nih.gov/articles/PMC3820647/](https://pmc.ncbi.nlm.nih.gov/articles/PMC3820647/)
30. **Roder2018_Auranofin_TrxR** — The gold complex auranofin: new perspectives for cancer therapy (Roder & Thomson, Discov Oncol 2022). [https://pmc.ncbi.nlm.nih.gov/articles/PMC8777575/](https://pmc.ncbi.nlm.nih.gov/articles/PMC8777575/)
31. **Spell2016_Au_ZnFinger** — Reactivity of Cys4 zinc-finger domains with gold(III) complexes — 'gold fingers' (Spell & Farver, Inorg Chem 2015). [https://pubs.acs.org/doi/abs/10.1021/acs.inorgchem.5b00360](https://pubs.acs.org/doi/abs/10.1021/acs.inorgchem.5b00360)
32. **Mendes2018_Auranofin_PARP1** — Auranofin synergizes with PARP inhibitor olaparib in mutant-p53 cancers (Mendes et al., 2023). [https://pmc.ncbi.nlm.nih.gov/articles/PMC10045521/](https://pmc.ncbi.nlm.nih.gov/articles/PMC10045521/)
33. **Laib1985_AuMT** — The binding of Gold(I) to metallothionein (Laib et al., 1985). [https://pubmed.ncbi.nlm.nih.gov/7411139/](https://pubmed.ncbi.nlm.nih.gov/7411139/)
34. **Bordin1996_AuCd_MT** — Gold replacement of cadmium, zinc-binding metallothionein (Bordin et al., 1996). [https://pubmed.ncbi.nlm.nih.gov/8865374/](https://pubmed.ncbi.nlm.nih.gov/8865374/)
35. **Maret2005_thionein_Zn** — Control of zinc transfer between thionein, metallothionein, and zinc proteins (Jacob, Maret, Vallee, PNAS 1998). [https://pmc.ncbi.nlm.nih.gov/articles/PMC19863/](https://pmc.ncbi.nlm.nih.gov/articles/PMC19863/)
36. **Blockhuys2017_CuProteome** — Defining the human copper proteome and its expression variation in cancers (Blockhuys et al., Metallomics 2017). [https://academic.oup.com/metallomics/article/9/2/112/5918432](https://academic.oup.com/metallomics/article/9/2/112/5918432)

---

## 7. Notes for GNN integration

- **Node features (suggested).** TCGA-LIHC expression (z-score), zinc-finger family flag, MT family flag, transporter direction (import/export), Au-displaceable flag, ferroptosis-pathway flag.
- **Edge features.** Type (one-hot 6 categories), confidence (ordinal: 3/2/1), cancer/LIHC/Au flags (binary), citation count.
- **Train/test split.** Edge prediction tasks should mask high-confidence LIHC edges (e.g., `MT1G–TP53`, `ZNF479–DNMT1`) to test whether the GNN can recover canonical Zn-biology mechanisms.
- **Au-perturbation experiments.** Model `AURANOFIN` as an external perturbation node — edge ablation from it lets the GNN simulate Au-induced Zn displacement *in silico*.
- **Sanity checks.** MTF1 must be the top hub; MT2A, TP53, ZIP7, MT1A must be in the top 10 — these are the empirical anchors against which the curated graph was sized.