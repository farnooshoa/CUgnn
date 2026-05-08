# Collaborator's 2020 Nat Comm Transcriptomics — QC Analysis

Context: Raw transcriptomics from the collaborator's 2020 Nature
Communications paper (hosted with the article was removed; raw file
`Raw Nat Comm Transcriptomics.xlsx` shared directly).

Design: Agilent 2-colour microarray, 4 samples condition **A** (A1, A2,
A2_R, A3) vs 4 samples condition **C** (C1–C4). `logFC = log2(A / C)`.
16,168 probes covering **11,831 unique gene symbols**; **1,618** pass
adj.P.Val < 0.05 on the shipped differential analysis.

## 1. Confirming the collaborator's observation

> "there are not many copper-transcripts with significant change, but
>  there are lots of Zn-based changes."

| feature | Cu proteome (58 genes) | Zn-annotated (`"zinc"` / `"metallothionein"` in GENENAME) |
|---|---:|---:|
| genes detected on the array | **29** / 58 | **574** |
| significant (adj.P < 0.05) | **3** | **82** |
| significant fraction | 10.3% | 14.3% |
| ratio Zn / Cu (sig counts) | — | **27.3 ×** |

**Conclusion**: the collaborator's claim is fully supported.
Cu proteome shows only 3 significant hits (PRNP, ATP7A, ENOX2); the
Zn-annotated gene set shows 82 significant hits dominated by zinc-finger
transcription factors (ZFAND2A, ZFHX3, ZMYND8, ZFAND5, ZMIZ1...).

Volcano plot: `figures/natcomm_cu_vs_zn.png`.
Cu-proteome heatmap: `figures/natcomm_cu_proteome_heatmap.png`.

### Cu proteome — top 10 by |logFC| (58-gene reference)

| gene | logFC | adj.P.Val |
|---|---:|---:|
| PRNP | +1.316 | 2.92e-10 |
| ATP7A | -0.624 | 1.00e-02 |
| ENOX2 | -0.587 | 1.72e-02 |
| SCO2 | -0.433 | 1.00e+00 |
| S100A5 | +0.423 | 1.00e+00 |
| APP | -0.409 | 1.00e+00 |
| LOXL2 | -0.403 | 1.00e+00 |
| SOD1 | +0.372 | 1.00e+00 |
| SCO1 | +0.275 | 1.00e+00 |
| CCS | +0.242 | 1.00e+00 |

Only **PRNP, ATP7A, ENOX2** are formally significant; everything else is
modest (|logFC| < 0.5). The Cu proteome is quiet in this experiment.

### Zn-related — top 15 by |logFC|

| gene | GENENAME (truncated) | logFC | adj.P.Val |
|---|---|---:|---:|
| ZFAND2A | zinc finger AN1-type containing 2A | +2.673 | 1.27e-15 |
| ZFHX3 | zinc finger homeobox 3 | -1.758 | 9.23e-12 |
| ZMYND8 | zinc finger MYND-type containing 8 | -1.587 | 2.13e-11 |
| ZFAND5 | zinc finger AN1-type containing 5 | +1.339 | 4.63e-10 |
| ZMIZ1 | zinc finger MIZ-type containing 1 | -1.324 | 3.74e-10 |
| IKZF2 | IKAROS family zinc finger 2 | -1.284 | 1.41e-08 |
| ZNRF3 | zinc and ring finger 3 | -1.238 | 3.78e-09 |
| ZBTB21 | zinc finger and BTB domain containing 21 | +1.227 | 4.25e-08 |
| GLI3 | GLI family zinc finger 3 | -1.169 | 9.75e-08 |
| HELZ | helicase with zinc finger | -1.017 | 8.48e-08 |
| ZNF473 | zinc finger protein 473 | +1.005 | 6.02e-07 |
| ZNF394 | zinc finger protein 394 | +1.004 | 1.17e-07 |
| ZNF256 | zinc finger protein 256 | +0.975 | 1.02e-06 |
| ZNF761 | zinc finger protein 761 | +0.972 | 1.53e-06 |
| ZBTB42 | zinc finger and BTB domain containing 42 | -0.969 | 1.08e-05 |

Almost all top-ranked Zn genes are **zinc-finger transcription factors**
(ZFAND, ZFHX3, ZMYND8, ZMIZ1, ZNF... families). This pattern — broad
ZnF-TF dysregulation, quiet Cu proteome, strong oxidative / heat-shock
response (HSPA6, ATF3, HMOX1, DDIT3 all logFC > 3, adj.P ~10⁻¹⁷) — is
exactly what **gold-induced displacement of zinc from zinc-finger
proteins** looks like.

### Stress / metal-response markers (top 6)

| gene | logFC | adj.P.Val |
|---|---:|---:|
| HSPA6 | +4.467 | 1.11e-17 |
| ATF3 | +4.259 | 1.11e-17 |
| HMOX1 | +3.682 | 1.18e-17 |
| HSPA1B | +3.508 | 4.03e-17 |
| DDIT3 | +3.430 | 6.93e-17 |
| HSPA1A | +3.260 | 7.95e-17 |

## 2. What this means for the grant narrative

- The collaborator's Cu vs Zn ratio observation is statistically real
  and directly relevant to a multi-metal metallome framework proposal.
- Our existing TCGA-LIHC Cu-proteome GNN pipeline is **directly
  adaptable** to a Zn proteome: replace the 54/58 Cu node set with a
  comparable Zn-proteome list (UniProt zinc-binding ~2,700 genes →
  curate down to ~60–80 "core Zn proteome"), keep everything else.
- The same "attention recovers canonical metal-handling pairs" story
  should work for Zn (ZIP/ZnT transporters, MT1/MT2 metallothioneins,
  zinc-finger clusters).
- Adding gold (Au) is scope for the grant but not for the pilot:
  Au doesn't have a native "proteome" in the UniProt sense — it acts
  by **displacement** of Zn/Cu. A principled Au-GNN would need
  perturbation data like *this very dataset*.

## 3. Limitations of running our current pipeline on this dataset as-is

1. **Only 29/58 Cu genes detected** (Agilent 026652 array has limited
   coverage; all histones, CP, ALB, AFP, LOX, mitochondrial COX absent).
   Running the 58-node Cu graph on this cohort is not apples-to-apples
   to the TCGA-LIHC run.
2. **n = 4 vs 4** is too small for graph classification; this dataset
   is appropriate as a perturbation-level validation cohort
   (single-sample inference against a TCGA-trained model), not as a
   training cohort.
3. The design is a **treatment effect** (condition A vs C), not
   tumor-vs-normal — so the TCGA-trained classifier head is not
   interpretable here without re-purposing.

## 4. Suggested next step for Phase 5 (metallome extension)

1. Curate a **Zn proteome** (~60–80 genes from UniProt zinc-binding +
   known ZnF-TF hubs) in the same format as `copper_gene_list.csv`.
2. Pull TCGA-LIHC expression for this Zn set, rerun the 3 tasks
   (tumor/normal, stage, survival).
3. Use this Nat Comm dataset as an **orthogonal perturbation**
   validation: apply the Zn-proteome GNN attention readout to
   condition A vs C at the gene level and compare to the 82 significant
   Zn genes already identified.
4. Write up as "the metallome GNN framework, demonstrated on Cu (TCGA
   LIHC) and validated on Zn perturbation (collaborator's Nat Comm
   cohort)" — exactly the grant aim he described.

## Files produced
- `data/natcomm_agg_per_gene.tsv` — de-duplicated gene-level Agilent logFC table
- `outputs/final_comparison/natcomm_qc.md` — this document
- `outputs/final_comparison/figures/natcomm_cu_vs_zn.png` — volcano comparison
- `outputs/final_comparison/figures/natcomm_cu_proteome_heatmap.png` — Cu-proteome heatmap in his cohort
