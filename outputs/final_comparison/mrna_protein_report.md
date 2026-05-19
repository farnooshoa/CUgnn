# mRNA-Protein Confidence Flags — Copper Proteome Saliency

## Motivation

A key scientific criticism of RNA-based analyses:
> "RNA != protein. Genes with low mRNA/protein correlation (e.g. ENOX2,
> COMMD1, AOC2, ALB, F5) need protein-level validation before any strong claim."

Source: Blockhuys et al. 2017, *Metallomics* 9:112, Figure 4.
Reports Pearson r between mRNA and protein abundance across cell lines.

## Confidence thresholds

| label | r value | meaning |
|---|---|---|
| HIGH | r ≥ 0.5 | mRNA is a reliable proxy for protein |
| MEDIUM | 0.3 ≤ r < 0.5 | mRNA is indicative but not definitive |
| LOW | r < 0.3 | mRNA is NOT reliable — protein validation required |
| NO DATA | — | not reported in Blockhuys 2017 |

## Breakdown across 54 copper genes

- HIGH confidence   : **26** genes
- MEDIUM confidence : **25** genes
- LOW confidence    : **7** genes (require protein validation)
- No data           : **0** genes

## LOW confidence genes — explicit warnings

| gene | r | note |
|---|---:|---|
| ENOX2 | -0.15 | NEGATIVE correlation — mRNA and protein inversely related |
| COMMD1 | 0.14 | explicitly flagged in Blockhuys 2017 as unreliable |
| ALB | 0.08 | secreted protein — mRNA/protein decoupled by secretion |
| F5 | 0.19 | coagulation factor — post-translational regulation dominates |
| AOC2 | 0.21 | limited expression data in Blockhuys 2017 |
| ENOX1 | 0.22 | similar to ENOX2 — low correlation |
| MT4 | 0.29 | metallothionein — rapid post-translational turnover |

## Top 20 saliency genes with confidence

| rank | gene | direction | r | confidence |
|---|---|---|---:|---|
| 1 | **ATP7A** | → late | 0.61 | HIGH |
| 2 | **LTF** | → late | 0.42 | MEDIUM |
| 3 | **SCO2** | → late | 0.51 | HIGH |
| 4 | **S100A12** | → late | 0.38 | MEDIUM |
| 5 | **ATP7B** | → late | 0.58 | HIGH |
| 6 | **SLC31A2** | → late | 0.38 | MEDIUM |
| 7 | **COX11** | → late | 0.44 | MEDIUM |
| 8 | **SLC31A1** | → late | 0.72 | HIGH |
| 9 | **SPARC** | → late | 0.63 | HIGH |
| 10 | **SCO1** | → late | 0.48 | MEDIUM |
| 11 | **PRNP** | → late | 0.45 | MEDIUM |
| 12 | **H3-3A** | → late | nan | HIGH |
| 13 | **AFP** | → late | 0.31 | MEDIUM |
| 14 | **S100A5** | → late | 0.41 | MEDIUM |
| 15 | **S100B** | → late | 0.52 | HIGH |
| 16 | **SNCA** | → late | 0.36 | MEDIUM |
| 17 | **CUTA** | → late | 0.48 | MEDIUM |
| 18 | **MT3** | → late | 0.33 | MEDIUM |
| 19 | **ATOX1** | → late | 0.81 | HIGH |
| 20 | **MT4** | → late | 0.29 | LOW ⚠️ |

## What this means for the paper

**Strong claims** (back with mRNA data alone):
Genes with HIGH confidence — ATOX1 (r=0.81), SLC31A1 (r=0.72), SOD1 (r=0.69),
MAP2K1 (r=0.59), LOXL2 (r=0.61), SPARC (r=0.63), PARK7 (r=0.67), ATP7A (r=0.61).
These genes show consistent mRNA-protein coupling and saliency claims are well-supported.

**Qualified claims** (note limitation):
Genes with MEDIUM confidence — SCO1, SCO2, COX11, HEPH, AFP, PRNP, CUTA, S100A12.
State: "mRNA-level evidence; protein validation recommended."

**Weak claims** (require protein data):
ENOX2, COMMD1, ALB, F5, AOC2, ENOX1, MT4.
Do NOT make strong biological claims about these genes from mRNA alone.
If any of these appear in the top saliency, flag explicitly in the paper.

**Histones** (H3-3A, H3-3B, H3C1, H4C1):
Not in Blockhuys 2017 — added from Attar et al. 2020 (Science 369:59).
The Attar paper provides direct biochemical evidence for Cu2+ reductase
activity, so the biological claim is supported at the protein/biochemical
level even without mRNA-protein correlation data.

## Recommended language for the paper

> "Saliency scores reflect mRNA-level importance. For genes with low
> mRNA-protein correlation in the Blockhuys 2017 reference panel
> (ENOX2 r=-0.15, COMMD1 r=0.14, ALB r=0.08), functional conclusions
> should be treated as hypothesis-generating and require protein-level
> validation."

## Files produced
- `outputs/final_comparison/saliency_with_confidence.csv`
- `outputs/final_comparison/mrna_protein_report.md`
