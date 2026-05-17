# Interpretability Analysis — Stage Classification

## Method

- **Signed saliency**: gradient of the late-stage logit w.r.t. node
  features, averaged across all 349 stage-labeled tumor samples.
  Positive = gene pushes model toward predicting late stage.
  Negative = gene pushes model toward predicting early stage.
- **Attention edges**: mean attention weight from the final GAT layer,
  averaged across all samples. High attention = model focuses on this
  gene-gene connection when making predictions.

## Top 20 genes by saliency — comparison with biology

| gene | saliency rank | model direction | log2FC (T/N) | expected stage | biology |
|---|---:|---|---:|---|---|
| **ATOX1** | 1 | → early ✗ | +0.471 | late | Cu chaperone, nuclear signaling, proliferation |
| **ATP7B** | 2 | → early ✗ | -0.362 | late | Cu efflux, cisplatin resistance, HCC progression |
| **SPARC** | 3 | → early ✗ | +0.751 | late | ECM remodeling, tumor progression |
| **COX17** | 4 | → early | +0.031 | unknown | — |
| **SCO2** | 5 | → early | +0.088 | unknown | — |
| **ALB** | 6 | → early | -2.383 | unknown | — |
| **ATP7A** | 7 | → early ✗ | +0.121 | late | Cu efflux pump |
| **SCO1** | 8 | → early | -0.355 | unknown | — |
| **SLC31A1** | 9 | → early ✗ | -1.078 | late | Cu importer CTR1, elevated in proliferating cells |
| **S100A13** | 10 | → early | +0.347 | unknown | — |
| **SLC31A2** | 11 | → early | -0.042 | unknown | — |
| **CUTA** | 12 | → early | +0.591 | unknown | — |
| **SOD1** | 13 | → early ✓ | -0.818 | early | Antioxidant, protective in early stage |
| **COX11** | 14 | → early | +0.036 | unknown | — |
| **APP** | 15 | → early | -0.590 | unknown | — |
| **PRNP** | 16 | → early ✓ | -1.189 | early | Prion protein, Cu binding, tumor suppressor role |
| **MT-CO1** | 17 | → early | -1.141 | unknown | — |
| **CCS** | 18 | → early | -0.487 | unknown | — |
| **MT-CO2** | 19 | → early | -0.924 | unknown | — |
| **CUTC** | 20 | → early | -0.232 | unknown | — |

**Model-biology agreement: 2/7 genes (29%)**

## Top 15 attention edges

| source | target | attention mean | in curated edges? |
|---|---|---:|---|
| PAM | PAM | 1.0000 | self-loop / compartment |
| ENOX1 | ENOX1 | 1.0000 | self-loop / compartment |
| ENOX2 | ENOX2 | 1.0000 | self-loop / compartment |
| COX17 | COX11 | 0.6530 | ✓ curated |
| COX17 | SCO2 | 0.6515 | ✓ curated |
| COX17 | SCO1 | 0.6448 | ✓ curated |
| SOD1 | PARK7 | 0.6378 | ✓ curated |
| CP | HEPH | 0.6366 | ✓ curated |
| CP | HEPHL1 | 0.6364 | ✓ curated |
| SLC31A2 | SLC31A2 | 0.6221 | self-loop / compartment |
| COX17 | COX17 | 0.6168 | self-loop / compartment |
| AFP | AFP | 0.6108 | self-loop / compartment |
| CUTA | CUTC | 0.5925 | ✓ curated |
| ALB | F5 | 0.5862 | ✓ curated |
| ATOX1 | SOD3 | 0.5685 | ✓ curated |

## Interpretation

The saliency map shows which copper genes the GAT relies on most
when predicting AJCC stage. Agreement with known biology validates
that the model is learning real copper biology, not statistical noise.

The attention edges show which gene-gene connections the model
focuses on. If high-attention edges match curated copper biology
(e.g. ATOX1→ATP7B, LOX→SPARC, SLC31A1→ATOX1), this strongly
supports the biological relevance of the graph structure.

## Files produced
- `interp_saliency.csv` — all genes ranked by saliency
- `interp_attention.csv` — all attention edges ranked
- `interp_comparison_table.csv` — top 20 genes with biology comparison
- `interp_report.md` — this document