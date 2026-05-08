# Glossary — terms used in the Cu-proteome pilot materials

A one-page reference for the ML / statistical language in the earlier emails, `histone_results.md`, and `natcomm_qc.md`. Roughly grouped by topic.

## Data and cohorts

- **TCGA-LIHC** — The Cancer Genome Atlas, Liver Hepatocellular Carcinoma cohort: 374 tumor and 50 matched normal RNA-seq samples used for training.
- **GDC STAR-Counts** — the Genomic Data Commons release of gene-level read counts aligned with the STAR aligner. The raw table we subset down to the 58 Cu-proteome genes.
- **FPKM-UQ** — upper-quartile normalized Fragments Per Kilobase per Million reads. A TCGA normalization that makes expression values comparable across samples.
- **Agilent 2-colour microarray** — the platform behind your 2020 Nat Comm transcriptomics. Each spot measures fluorescence in two channels; `logFC = log2(A / C)` is the ratio of condition A vs condition C on that probe.

## Differential expression statistics

- **logFC / log2 fold change** — log2 of (tumor mean / normal mean). +1 means 2× higher in tumor, -1 means 2× lower.
- **p-value** — probability of seeing a difference at least this big by chance if there were no real effect.
- **adj.P.Val / Benjamini-Hochberg / FDR** — a p-value corrected for testing thousands of genes at once. `adj.P < 0.05` means the false discovery rate is controlled at 5%.
- **significant** — in the tables it means `adj.P.Val < 0.05`.

## Graph / GNN vocabulary

- **Node** — one gene (or, in the updated cartoon, the Au compounds panel).
- **Edge** — a connection between two nodes. Types used in this project: physical (protein-protein binding), co-expression, genetic (pathway/epistasis), shared compartment (weak prior), and Au interaction (new).
- **Graph** — the whole set of nodes plus edges. In this pilot every patient is represented as the same 58-node graph, with their own expression values on the nodes.
- **GCN (Graph Convolutional Network)** — a neural network that updates each node by averaging its neighbours' features, repeated a few layers. Simple, strong baseline.
- **GAT (Graph Attention Network)** — same idea, but before averaging, the model learns a weight (attention) for each neighbour. The attention values are the extra interpretability layer we read out.
- **Message passing** — the general term for node-updates-from-neighbours. Both GCN and GAT are message-passing networks.
- **Node features** — the numbers attached to each node. In this pilot: log2 expression, z-scored expression, and a one-hot for functional category.
- **Node importance / saliency** — how much the model's output changes if you perturb that one node's features. Higher saliency means the model relies on that gene more.
- **Attention (for an edge)** — the weight the GAT gave that edge when message-passing. `attention_sum` aggregates across all patients and both directions. Used as a ranking, not an absolute quantity.
- **Canonical edges** — the hand-curated Cu-biology edges from the Blockhuys 2017 paper (ATOX1↔ATP7A/B, CCS↔SOD1, COX17↔SCO1/2, etc.). Drawn darker in the cartoon.

## Evaluation

- **ROC-AUC** — area under the ROC curve. 1.0 is perfect, 0.5 is coin-flip. On tumor-vs-normal the task is easy and almost all models saturate near 1.0; on stage or survival the task is harder and AUC around 0.7 is meaningful.
- **Cross-validation (k-fold CV)** — split the data into k parts, train on k-1 and test on the remaining one, repeat. Used to avoid overfitting.
- **StratifiedGroupKFold** — a safer form of k-fold that (a) keeps class ratios balanced across folds (stratified) and (b) keeps all samples from the same patient or group in the same fold (grouped). Without this, the same patient can leak across train and test and inflate the AUC.
- **Permutation-label sanity check** — shuffle the labels and retrain; a correct model should collapse to AUC ~0.5. Confirms the real AUC is not an artifact.

## Specific to this pilot

- **Cu proteome (58 nodes)** — 54 human Cu-binding / Cu-handling proteins from Blockhuys 2017 *Metallomics* + 4 histone representatives (H3-3A, H3-3B, H3C1, H4C1) added from Attar 2020 *Science*.
- **Shared fixed topology** — every patient's graph has the same 58 nodes and same edges; only the node features (expression values) differ across patients. Lets the GNN learn a single structural prior.
- **Au compounds node** — new extracellular node representing the Au(I) / Au(III) drug panel as a whole. Connecting lines mean "reported or proposed Au-biomolecule interaction" with evidence strength encoded by line style (solid = published, dashed = unpublished).
- **Hub / top-band node** — a node whose removal degrades the model disproportionately; typically sits at the centre of many edges, like ATOX1 or ATP7B in our attention readout.

If any of these need more depth or a worked example during the call, happy to walk through on a whiteboard.
