# GNN interpretability summary

Best model by CV ROC-AUC: **gat** (AUC=0.995).

## Top-10 important genes (saliency)

```
gene_symbol  importance
         CP    1.768464
       SOD1    1.141387
     MT-CO1    1.085892
      ATOX1    1.052884
     MAP2K1    0.985290
        ALB    0.971084
      LOXL4    0.887514
     MT-CO2    0.826755
        LOX    0.808407
      ATP7B    0.708219
```


## Top-20 attention edges (GAT)

```
 source  target  attention_sum
  ENOX1   ENOX1      50.000000
  ENOX2   ENOX2      50.000000
    PAM     PAM      50.000000
  PARK7   PARK7      45.685221
    CCS SLC31A1      41.885603
  ATP7B     AFP      39.311686
   HEPH    HEPH      35.268848
   SOD3   ATOX1      34.491044
 HEPHL1  HEPHL1      34.313006
  PARK7    SOD1      33.998735
S100A12 S100A13      33.950987
S100A12   S100B      33.950987
S100A12  S100A5      33.950987
 MT-CO2    SCO2      31.186671
     F5      F5      29.682642
   SOD1     CCS      29.114729
  ATP7B  COMMD1      28.608684
   CUTC    CUTC      28.458271
  ATP7B   ATP7A      27.555302
  MOXD1  MOXD2P      27.270298
```

Notes: saliency is the mean absolute gradient of the tumor logit w.r.t.
node features, averaged across the first 50 graphs. Large values point
to Cu-proteome genes whose expression most strongly shifts the model's
tumor probability — candidates worth biological follow-up.