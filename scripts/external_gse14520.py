"""External validation on GSE14520 — Roessler et al. HCC cohort.

GSE14520:
  Platform: Affymetrix HG-U133A 2.0 (GPL3921)
  N: 488 samples (~225 tumor + ~220 adjacent-normal)
  Source: Chinese liver-cancer cohort, mostly HBV-related
  Paper: Roessler et al. Cancer Res. 2010

Pipeline:
  1. Download SOFT via GEOparse (cached under data/geo_cache)
  2. Build sample x probe expression matrix
  3. Map probes to HGNC symbols using the GSM platform table
  4. Aggregate probes to genes (max-variance probe per gene)
  5. log2-transform (already log2 on Affy, but ensure consistent scale)
  6. Subset to 54 Cu genes
  7. Harmonise with TCGA: Pearson-z each gene across samples
  8. Train GAT on TCGA-LIHC, test on GSE14520 (no CV — pure external)
  9. Report ROC-AUC on external + gene-importance overlap with TCGA-trained model

Outputs:
  data/gse14520/*                                   - cached SOFT + expression
  outputs/final_comparison/external_gse14520.md
  outputs/final_comparison/external_gse14520_metrics.csv
"""
from __future__ import annotations
import sys
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from torch_geometric.loader import DataLoader
from sklearn.metrics import (
    accuracy_score, balanced_accuracy_score, f1_score, roc_auc_score,
)

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")

from src.preprocessing import load_lihc_dataset, load_copper_genes
from src.graph_building import build_functional_graph
from src.gnn_models.dataset import build_graph_dataset, graph_to_edge_tensors
from src.gnn_models.models import GATGraphClassifier
from src.gnn_models.train import _class_weights, TrainConfig
from src.utils import RANDOM_SEED

OUT = ROOT / "outputs" / "final_comparison"
CACHE = ROOT / "data" / "geo_cache"
CACHE.mkdir(parents=True, exist_ok=True)
GSE_EXPR = ROOT / "data" / "gse14520_expression.tsv"
GSE_META = ROOT / "data" / "gse14520_metadata.tsv"


def download_and_prepare(force: bool = False) -> tuple[pd.DataFrame, pd.DataFrame]:
    if not force and GSE_EXPR.exists() and GSE_META.exists():
        print("[gse] cached files present")
        return (pd.read_csv(GSE_EXPR, sep="\t", index_col=0),
                pd.read_csv(GSE_META, sep="\t", index_col=0))

    import GEOparse
    print("[gse] fetching GSE14520 via GEOparse (this hits NCBI) ...")
    gse = GEOparse.get_GEO(geo="GSE14520", destdir=str(CACHE),
                            include_data=True, silent=True)
    gsms = list(gse.gsms.values())
    print(f"[gse] {len(gsms)} samples downloaded")

    # Collect expression per sample. GSE14520 has two platforms
    # GPL3921 (HT-HG-U133A) and GPL571 (HG-U133A 2.0). Use GPL3921 (most samples).
    target_plat = "GPL3921"
    sel = [g for g in gsms if g.metadata.get("platform_id", [""])[0] == target_plat]
    print(f"[gse] samples on {target_plat}: {len(sel)}")

    # Expression matrix
    cols = {}
    for gsm in sel:
        t = gsm.table
        if "VALUE" not in t.columns:
            continue
        t = t.set_index("ID_REF")["VALUE"]
        cols[gsm.name] = t
    expr_raw = pd.DataFrame(cols).astype(float)
    print(f"[gse] raw expression: {expr_raw.shape}")

    # Probe annotation from the platform table
    gpl = gse.gpls[target_plat]
    anno = gpl.table.copy()
    anno_sym_col = "Gene Symbol" if "Gene Symbol" in anno.columns else \
                   ("ORF" if "ORF" in anno.columns else None)
    if anno_sym_col is None:
        raise RuntimeError(f"Cannot find gene-symbol column in {target_plat} table")
    anno = anno[["ID", anno_sym_col]].rename(columns={"ID": "probe", anno_sym_col: "gene"})
    anno = anno[anno["gene"].notna() & (anno["gene"].astype(str) != "")]
    anno["gene"] = anno["gene"].astype(str).str.split(r"\s*///\s*", regex=True)
    anno = anno.explode("gene")
    anno["gene"] = anno["gene"].str.upper().str.strip()
    anno = anno[anno["gene"] != ""]

    # Restrict to Cu proteome
    copper_genes = set(load_copper_genes()["gene_symbol"])
    anno_cu = anno[anno["gene"].isin(copper_genes)]
    print(f"[gse] probes targeting Cu genes: {anno_cu['probe'].nunique()} "
          f"covering {anno_cu['gene'].nunique()} / 54 Cu genes")

    # Keep only those probes in expression
    expr_raw = expr_raw.loc[expr_raw.index.intersection(anno_cu["probe"])]
    # Aggregate probes to genes — take the max-variance probe per gene
    expr_raw["_var"] = expr_raw.var(axis=1)
    best_per_gene = (anno_cu.set_index("probe")
                            .join(expr_raw[["_var"]], how="inner")
                            .sort_values("_var", ascending=False)
                            .drop_duplicates(subset=["gene"], keep="first"))
    expr_raw = expr_raw.drop(columns=["_var"])
    per_gene = expr_raw.loc[best_per_gene.index]
    per_gene.index = best_per_gene["gene"]
    per_gene.index.name = "gene_symbol"

    # GSE14520 is already log2 (RMA-normalised Affy)
    # Reindex to all 54 Cu genes — fill missing with gene-median if any
    copper_list = load_copper_genes()["gene_symbol"].tolist()
    per_gene = per_gene.reindex(copper_list)
    missing = per_gene.index[per_gene.isna().all(axis=1)].tolist()
    print(f"[gse] final: {per_gene.shape}. Missing Cu genes: {missing}")

    # Save
    per_gene.to_csv(GSE_EXPR, sep="\t")

    # Metadata — sample_type from title / description
    rows = []
    for gsm in sel:
        name = gsm.name
        meta = gsm.metadata
        title = meta.get("title", [""])[0]
        src = meta.get("source_name_ch1", [""])[0].lower()
        # GSE14520 titles are like "Liver Tumor ..." or "Liver non-Tumor ..."
        is_tumor = ("tumor" in title.lower()) and ("non-tumor" not in title.lower()) \
                    and ("non tumor" not in title.lower())
        # Fall back to source name
        if not is_tumor and ("cancer" in src or "tumor" in src) and "non" not in src:
            is_tumor = True
        rows.append({
            "sample_id": name, "title": title,
            "source": src,
            "sample_type": "Tumor" if is_tumor else "Normal",
        })
    md = pd.DataFrame(rows).set_index("sample_id")
    md.to_csv(GSE_META, sep="\t")
    print(f"[gse] metadata: {md['sample_type'].value_counts().to_dict()}")
    return per_gene, md


def train_on_tcga(ds, graph):
    """Train a single GAT on all TCGA-LIHC; return the trained model + gene order."""
    bundle = build_graph_dataset(ds, graph, zscore_train_mask=np.ones(ds.n_samples, dtype=bool))
    y = np.array([int(d.y.item()) for d in bundle.data_list])
    cfg = TrainConfig(model="gat", epochs=80)
    class_w = _class_weights(y, cfg.device)
    model = GATGraphClassifier(bundle.in_dim, cfg.hidden, n_classes=2,
                                n_heads=4, n_layers=cfg.n_layers,
                                dropout=cfg.dropout).to(cfg.device)
    opt = torch.optim.Adam(model.parameters(), lr=cfg.lr, weight_decay=cfg.weight_decay)
    loss_fn = torch.nn.CrossEntropyLoss(weight=class_w)
    dl = DataLoader(bundle.data_list, batch_size=cfg.batch_size, shuffle=True)
    for epoch in range(cfg.epochs):
        model.train()
        for batch in dl:
            batch = batch.to(cfg.device)
            opt.zero_grad()
            logits = model(batch.x, batch.edge_index, batch.batch,
                           edge_weight=getattr(batch, "edge_weight", None))
            loss = loss_fn(logits, batch.y)
            loss.backward()
            opt.step()
    return model, bundle, ds.expression.index.tolist()


def apply_to_external(model, ds, graph, gene_order, gse_expr: pd.DataFrame,
                       gse_md: pd.DataFrame):
    """Feed GSE14520 samples through the trained model."""
    cfg = TrainConfig(model="gat")
    # Build per-gene mean/std using TCGA for fold-safe scaling
    tcga_mu = ds.expression.mean(axis=1)
    tcga_sd = ds.expression.std(axis=1).replace(0, 1)

    # For GSE14520, fill missing genes with TCGA mean (so z-score = 0)
    gse = gse_expr.copy()
    gse = gse.reindex(gene_order)
    for g in gse.index:
        if gse.loc[g].isna().all():
            gse.loc[g] = tcga_mu[g]  # all NaN row → set to TCGA mean
    gse = gse.fillna(gse.mean(axis=1).fillna(tcga_mu))

    # Cross-platform harmonisation: rank-quantile-map each sample to TCGA distribution
    # For simplicity here, just use TCGA z-score normalisation.
    z_gse = gse.sub(tcga_mu, axis=0).div(tcga_sd, axis=0)

    # Category one-hot
    cop = ds.copper_genes.set_index("gene_symbol")
    cats = ["transporter", "enzyme", "other_or_unknown"]
    cat_feats = np.zeros((len(gene_order), len(cats)), dtype=np.float32)
    for i, g in enumerate(gene_order):
        c = cop.loc[g, "functional_category"] if g in cop.index else "other_or_unknown"
        cat_feats[i, cats.index(c) if c in cats else 2] = 1.0

    # Edges
    edge_index, edge_weight = graph_to_edge_tensors(graph, gene_order)

    # Labels
    labels = (gse_md.loc[gse.columns, "sample_type"].str.lower() == "tumor").astype(int).to_numpy()

    # Build Data per sample
    from torch_geometric.data import Data
    data_list = []
    for i, sid in enumerate(gse.columns):
        x_cols = [
            gse[sid].to_numpy().astype(np.float32).reshape(-1, 1),
            z_gse[sid].to_numpy().astype(np.float32).reshape(-1, 1),
            cat_feats,
        ]
        x = np.concatenate(x_cols, axis=1)
        data = Data(
            x=torch.tensor(x, dtype=torch.float32),
            edge_index=edge_index, edge_weight=edge_weight,
            y=torch.tensor([int(labels[i])], dtype=torch.long),
        )
        data_list.append(data)

    # Predict
    dl = DataLoader(data_list, batch_size=cfg.batch_size, shuffle=False)
    model.eval()
    ys, ps, preds = [], [], []
    with torch.no_grad():
        for batch in dl:
            batch = batch.to(cfg.device)
            logits = model(batch.x, batch.edge_index, batch.batch,
                           edge_weight=getattr(batch, "edge_weight", None))
            prob = torch.softmax(logits, dim=-1)[:, 1].cpu().numpy()
            pr = logits.argmax(dim=-1).cpu().numpy()
            ys.extend(batch.y.cpu().numpy().tolist())
            ps.extend(prob.tolist())
            preds.extend(pr.tolist())
    ys = np.array(ys); ps = np.array(ps); preds = np.array(preds)
    return {
        "accuracy": accuracy_score(ys, preds),
        "balanced_accuracy": balanced_accuracy_score(ys, preds),
        "f1": f1_score(ys, preds, zero_division=0),
        "roc_auc": roc_auc_score(ys, ps) if len(set(ys)) > 1 else float("nan"),
        "n_samples": len(ys),
        "n_tumor": int((ys == 1).sum()),
        "n_normal": int((ys == 0).sum()),
    }


def main():
    ds = load_lihc_dataset(require_real=True)
    genes = ds.expression.index.tolist()
    graph = build_functional_graph(genes, ds.copper_genes)

    gse_expr, gse_md = download_and_prepare()

    print(f"\n[gse] training GAT on full TCGA-LIHC ...")
    model, bundle, gene_order = train_on_tcga(ds, graph)

    print(f"[gse] applying trained model to GSE14520 ...")
    ext = apply_to_external(model, ds, graph, gene_order, gse_expr, gse_md)
    print(f"[gse] external cohort result: {ext}")

    summary_rows = [
        {"cohort": "TCGA-LIHC (CV, internal)", "model": "GAT",
         "n": 424, "n_tumor": 374, "n_normal": 50,
         "roc_auc": 0.993, "balanced_accuracy": 0.913, "note": "5-fold CV (from leakage audit)"},
        {"cohort": "GSE14520 (external)", "model": "GAT (TCGA-trained)",
         "n": ext["n_samples"], "n_tumor": ext["n_tumor"], "n_normal": ext["n_normal"],
         "roc_auc": ext["roc_auc"], "balanced_accuracy": ext["balanced_accuracy"],
         "note": "no CV — trained once on TCGA, predicted once on external"},
    ]
    df = pd.DataFrame(summary_rows)
    df.to_csv(OUT / "external_gse14520_metrics.csv", index=False)

    missing_genes = gse_expr.index[gse_expr.isna().all(axis=1)].tolist()
    (OUT / "external_gse14520.md").write_text(f"""# External Validation — GSE14520 (Roessler HCC Cohort)

## Source
- **NCBI GEO accession**: GSE14520
- **Platform**: Affymetrix HT-HG-U133A (GPL3921)
- **Downloaded via**: GEOparse 2.0.4 (cached under `data/geo_cache/`)
- **Samples after filtering to GPL3921**: {len(gse_md)} ({int((gse_md['sample_type']=='Tumor').sum())} tumor, {int((gse_md['sample_type']=='Normal').sum())} normal)
- **Paper**: Roessler S et al. *Cancer Research* 2010, 70:10202 — Chinese HBV-related HCC cohort
- **Probe → HGNC**: platform annotation (GPL3921), aggregated by max-variance probe per gene
- **Cu-gene coverage**: {(~gse_expr.isna().all(axis=1)).sum()} / 54 Cu genes ({'missing: ' + ', '.join(missing_genes) if missing_genes else 'all present'})

## Setup (no CV on external)
1. Train a **single GAT** on all 424 TCGA-LIHC samples.
2. Apply it once to each GSE14520 sample (no retraining, no label access).
3. Node features for GSE14520 are z-scored using **TCGA's** per-gene mean/std
   (fold-safe; the external cohort has no influence on training statistics).
4. Missing Cu genes (if any) are set to the TCGA per-gene mean; their z-score
   becomes zero and they contribute only through graph neighbours.

## Results

| cohort | N | ROC-AUC | balanced acc |
|---|---:|---:|---:|
| TCGA-LIHC 5-fold CV (internal) | 424 | **0.993** | 0.913 |
| **GSE14520 external** | {ext['n_samples']} | **{ext['roc_auc']:.3f}** | {ext['balanced_accuracy']:.3f} |
| Δ (external − internal) | | **{ext['roc_auc'] - 0.993:+.3f}** | {ext['balanced_accuracy'] - 0.913:+.3f} |

## Interpretation

### Quick read
- External ROC-AUC {ext['roc_auc']:.3f} vs TCGA-internal 0.993.
- AUC drop of {abs(ext['roc_auc'] - 0.993):.3f} is {'small — the model generalises well across cohorts, platforms (Affymetrix microarray vs Illumina RNA-seq), and populations (US/Western vs Chinese HBV-related HCC).' if ext['roc_auc'] > 0.85 else 'larger than expected — likely driven by platform differences (Affy probes + RMA normalisation do not map cleanly onto TCGA FPKM-UQ). Future work: quantile-normalise the external cohort to the TCGA distribution per gene.'}
- Balanced accuracy is more conservative than AUC on imbalanced cohorts; compare both.

### What this tells us
- **Positive external AUC (> 0.85)** would confirm the 54-gene Cu proteome carries generalisable tumor-vs-normal signal across cohorts. The pilot would be ready to write up as a short methods paper.
- **Lower external AUC** would point to platform batch effects (Affy microarray vs Illumina RNA-seq) rather than biological breakdown. Mitigation: quantile normalisation (rank-match each GSE14520 sample to the TCGA distribution) before inference.
- Either way, this is the first time we are testing the model on data it has never seen during training — far more meaningful than CV on a single cohort.

### Platform caveat
- TCGA is Illumina RNA-seq, log2(FPKM-UQ + 1). Scale roughly 0 to 15.
- GSE14520 is Affymetrix HG-U133A, RMA-normalised log2 signal. Scale roughly 2 to 14.
- Different dynamic ranges; z-scoring using TCGA's mu/sigma puts GSE14520 on
  the TCGA scale only to first order. A proper cross-platform experiment
  would add quantile normalisation or ComBat batch correction.

## Files produced
- `data/gse14520_expression.tsv` — Cu-gene × sample expression matrix for GSE14520
- `data/gse14520_metadata.tsv` — sample IDs + tumor/normal labels
- `data/geo_cache/` — raw GEOparse downloads
- `outputs/final_comparison/external_gse14520.md` — this document
- `outputs/final_comparison/external_gse14520_metrics.csv` — internal vs external comparison
""")
    print(f"\n[gse] wrote {OUT/'external_gse14520.md'}")
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
