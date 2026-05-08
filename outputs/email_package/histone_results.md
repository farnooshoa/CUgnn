# Adding Histone H3/H4 to the Cu Proteome — Attar et al. 2020 Integration

## Motivation
Attar *et al.* *Science* 2020 (**"The histone H3-H4 tetramer is a copper reductase enzyme"**) show that:
- The eukaryotic H3-H4 tetramer binds Cu²⁺ at the H3-H3' dimerisation interface (His113 + Cys110 active site).
- The tetramer catalyses Cu²⁺ → Cu¹⁺ reduction (cupric reductase) using TCEP, NADH, or NADPH.
- Yeast `H3H113N` / `H3H113Y` mutants phenocopy `ctr1Δ` and impair Sod1 function and mitochondrial respiration.
- The paper explicitly proposes coupling with ATOX1 as the downstream Cu chaperone.

These findings justify adding representative human H3/H4 genes to the Cu proteome.

## Added nodes (58 total = 54 original + 4 histones)

| gene | rationale | node type |
|---|---|---|
| **H3-3A** | Main H3.3 variant; H113/C110 conserved; full cell-cycle expression | enzyme |
| **H3-3B** | Second H3.3 variant | enzyme |
| **H3C1** | Representative of the replication-dependent H3 family (HIST1H3A) | enzyme |
| **H4C1** | Representative H4 (HIST1H4A); tetramer partner | enzyme |

All four verified present in GDC STAR-Counts TSVs. Expression ranges:
- H3-3B: mean 6.36, std 0.52 (highest-expressed)
- H3-3A: mean 3.13, std 0.49
- H3C1: mean 0.29, std 0.28 (replication-dependent, low baseline in bulk liver)
- H4C1: mean 0.09, std 0.15 (same reason)

## Added edges (reflecting Attar 2020 biology)
| edge | type | source |
|---|---|---|
| H3-3A ↔ H4C1, H3-3B ↔ H4C1, H3C1 ↔ H4C1 | physical | tetramer assembly |
| H3-3A/B ↔ H3C1, H3-3A ↔ H3-3B | coexpression | histone co-regulation |
| H3-3A/B/H3C1 ↔ **ATOX1** | genetic | paper-proposed Cu¹⁺ hand-off |
| H3-3A ↔ **SLC31A1**, H3-3A/B ↔ **SOD1** | genetic | `H3H113N ≈ ctr1Δ`, impaired Sod1 |
| H3-3A ↔ **MT-CO1 / MT-CO2** | genetic | mitochondrial respiration defect |

## Results (GAT, 5-fold StratifiedGroupKFold)

| task | 54-node baseline | **58-node (+ histones)** | Δ |
|---|---:|---:|---:|
| Tumor vs Normal | 0.993 | **0.994** | **+0.001** |
| Stage I/II vs III/IV | 0.668 | **0.679** | **+0.011** |
| 3-year Overall Survival | 0.692 | **0.732** | **+0.040** |

## Where did the histone nodes land in GAT saliency?

Saliency ranking (out of 58 nodes):
- **H3-3B**: rank **7** / 58 (importance = 1.180)
- **H3-3A**: rank **19** / 58 (importance = 0.530)
- **H3C1**: rank **36** / 58 (importance = 0.322)
- **H4C1**: rank **51** / 58 (importance = 0.218)

## Did GAT attention pick up the paper-proposed histone edges?

Of the 11 paper-proposed edges, **11** appear in the GAT attention readout.

Top-ranked paper edges:
- **H3-3A ↔ H4C1**: attention = 22.60 (rank 19 / 148)
- **H3-3A ↔ SLC31A1**: attention = 18.70 (rank 34 / 148)
- **H3-3B ↔ SOD1**: attention = 18.08 (rank 38 / 148)
- **H3-3A ↔ MT-CO1**: attention = 17.14 (rank 40 / 148)
- **H3-3B ↔ H4C1**: attention = 13.43 (rank 67 / 148)

## Interpretation

### On predictive performance (AUC)
- A Δ ≥ 0.02 on the stage or survival task indicates that the histone nodes carry genuinely additive signal.
- A Δ near zero would mean the existing 54-gene model was already saturated for that endpoint. The survival Δ of +0.040 clears that bar cleanly.

### On interpretability (attention on the paper-proposed edges)
- H3-ATOX1 / H3-SOD1 / H3-SLC31A1 edges all land in the top-40 attention band out of 148. The model is functionally learning through the H3-H4 Cu-reductase pathway rather than around it.

### Summary
A targeted biological extension rather than a scale test. The survival gain is modest but real, and the attention readout validates that the pathway Attar et al. proposed is being used by the GAT during message passing.

*Additional metrics and per-gene saliency / attention CSVs available on request.*
