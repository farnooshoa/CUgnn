# Paper Summary — Defining the Human Copper Proteome (2017)

## Full citation
Blockhuys S., Celauro E., Hildesjö C., Feizi A., Stål O., Fierro-González J. C., Wittung-Stafshede P.
**Defining the human copper proteome and analysis of its expression variation in cancers.**
*Metallomics*, 2017, 9, 112–123. DOI: 10.1039/C6MT00202A.

## Central question
Which human proteins bind or transport copper (Cu), and how are their transcript levels altered across human cancers? The authors aim to (1) curate a comprehensive list of Cu-binding proteins — "the human copper proteome" — and (2) screen its RNA-level behaviour across TCGA cancer types to nominate Cu-dependent proteins that may matter in cancer biology.

## Why copper is relevant to cancer
- Cu is a required cofactor for many enzymes (redox, secretory, mitochondrial respiration).
- Cancer hallmarks that depend on Cu: sustained proliferation (Cu-dependent MAPK/MEK1), angiogenesis (LOX, VEGF via HIF-1), invasion and metastasis (LOX in premetastatic niche; MEMO1, SPARC).
- Serum Cu rises in many cancer patients, while iron/zinc drop.
- Several Cu chaperones (ATOX1, CCS) have newly described non-classical roles, e.g. transcription-factor-like activity linked to proliferation.

## What the authors did
1. Curated all human Cu-binding / Cu-transporting proteins from UniProt, with manual review and literature additions (SPARC, MEMO1, MAP2K1) → **54 genes** ("human copper proteome").
2. Retrieved upper-quartile-normalised FPKM RNA-seq from GDC/TCGA for 25 cancer types; 18 had data for all 54 genes.
3. Computed log2(tumor/normal) and t-test p-values per gene per cancer.
4. Produced a heat map (Fig. 2) and performed Ward's-D2 hierarchical clustering on Pearson-correlation distance, yielding 8 gene clusters and 3 cancer super-groups.
5. Built per-cluster interaction networks in Cytoscape + GeneMANIA using physical, co-expression, and genetic interactions (Fig. 3).
6. Drilled down on breast cancer: integrated with PAM50 subtype proteogenomics (Mertins et al. 2016) and GO-term enrichment via GOrilla/ReviGO.
7. Validated ATOX1 protein upregulation in 67 breast-cancer tissue microarrays by IHC (Figs 5–6).

## Main methods (at a glance)
| Step | Tool |
|---|---|
| Gene curation | UniProt + Genecards + manual |
| Localisation | Genecards "Localization" (UniProt/COMPARTMENTS conf. 4/5) |
| RNA-seq | TCGA/GDC upper-quartile FPKM, log2 fold change, t-test |
| Heatmap/clustering | R 3.3.0, `gplots::heatmap.2`, Ward's-D2 on Pearson distance |
| Networks | Cytoscape + GeneMANIA App (physical, co-expression, genetic) |
| GO enrichment | GOrilla + ReviGO |
| Protein validation | IHC on tissue microarrays (anti-ATOX1 ab) |

## Main findings
- **54 Cu-binding proteins** define the human Cu proteome (<0.5% of the proteome); 12 transporters, ~27 enzymes, ~15 "other / unknown".
- Genes distribute across all major compartments — with enrichment in cytoplasm, nucleus, extracellular space, mitochondrion and Golgi.
- Expression changes are **non-trivial and cancer-type-specific** (Fig. 2): LOX, LOXL1-2, SPARC, ENOX2 frequently up; S100B, PRNP, ENOX1, SNCA, SOD3, HEPHL1, AOC3 frequently down.
- 8 gene clusters and 3 cancer super-groups emerge from clustering. Liver cancer clusters with prostate, cervix, breast, uterus, thyroid, colorectal, lung (group 3).
- Network analysis (Fig. 3) shows clusters 1, 3, and 8 are tightly inter-connected; clusters 2 and 4 are more dispersed. Cu-binding proteins in one heat-map cluster frequently reappear as first-level partners in another network — high inter-connectivity of the copper proteome.
- In breast cancer: five upregulated transporters (SLC31A1, ATOX1, ATP7B, COX17, SCO2) suggest elevated Cu import → ATP7B and Cu delivery to mitochondria. Upregulated enzymes include HEPHL1, TYRP1, LOXL1, LOXL2, MOXD1.
- ATOX1 protein is detectable in cancerous but generally absent in normal breast tissue (67-sample TMA), with strongest staining in luminal subtypes.

## Main conclusions
- A globally usable catalogue of Cu-binding human proteins exists and is tractable for systems-level analysis.
- Cu-binding protein expression varies by tissue of origin and cancer type; there is no single "copper cancer signature".
- Network-level view reveals tightly connected Cu clusters that are likely functional modules worth mechanistic follow-up.
- Breast cancer case study supports integrated Cu-flux remodelling (import ↑, mitochondrial delivery ↑, secretory/ATP7B ↑).
- ATOX1 is a candidate Cu-proteome biomarker in breast cancer with cytoplasmic + nuclear signal consistent with a proposed transcription-factor role.

## Limitations
- RNA ≠ protein; only 28/36 Cu-binders had high/moderate mRNA–protein correlation in breast PAM50 data (Fig. 4). The other ~22% may require direct proteomics.
- TCGA tumour/normal ratios are imbalanced (few normals) and technology-specific.
- Only 18 of 25 TCGA cancers had complete coverage for all 54 genes.
- The list is a lower bound ("list will grow with time") and some entries (SNCA, APP, PRNP, CUTC) have unclear Cu function.
- Breast-cancer IHC cohort (67 + 4 normals) is small; no matched tumour/normal from the same patient; no statistically significant subtype differentiation.

## Why this paper is useful for a TCGA-LIHC graph project
- **Node set**: a ready-made, biologically motivated list of 54 Cu-proteome genes — small enough for whole-cohort graphs, large enough to be interesting.
- **Functional annotations**: transporter / enzyme / other, plus subcellular compartment — both useful as node features and for interpretability.
- **Prior beliefs about edges**: Fig. 3 clusters (Ward's-D2 on tumor/normal log2FC) show that co-regulation, physical interaction, and co-expression all jointly shape Cu-protein modules — motivating multi-evidence edge building (STRING / BioGRID / GeneMANIA).
- **LIHC-specific relevance**: liver sits in Group 3 of the cancer dendrogram; liver cancer is mentioned specifically as having ATOX1 downregulation (Fig. S1) and is tightly linked to Cu metabolism via ceruloplasmin (CP), ATP7B (Wilson's disease), and ferroxidase activity. A copper-proteome view of LIHC is scientifically motivated.
- **Task framing**: the paper frames tumor-vs-normal log2FC on a fixed gene set — directly reusable as a graph-classification formulation.
- **Interpretability template**: Fig. 2 heatmap and Fig. 3 interaction networks are the kinds of artefacts we want our baseline and GNN outputs to echo.
