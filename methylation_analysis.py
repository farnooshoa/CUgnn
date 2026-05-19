"""Methylation analysis — multi-omic validation of interpretability results.

Three analyses:
  1. Expression-methylation correlation per gene
     (hypomethylated promoter → higher expression = expected)

  2. Stage-associated methylation changes
     (which copper genes show promoter methylation differences
      between early I/II and late III/IV tumors?)

  3. Multi-omic convergence
     (do GAT-important genes show BOTH expression AND methylation
      changes in the same direction? That is the strongest claim.)

  4. Add methylation as node feature → rerun stage GAT
     (does multi-omic input improve AUC?)

Outputs:
  outputs/final_comparison/methylation_correlations.csv
  outputs/final_comparison/methylation_stage_diff.csv
  outputs/final_comparison/methylation_convergence.csv
  outputs/final_comparison/methylation_report.md
"""
from __future__ import annotations
import sys
import random
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import torch
import torch.nn as nn
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.metrics import roc_auc_score
from torch_geometric.loader import DataLoader

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from src.preprocessing import load_lihc_dataset
from src.graph_building import build_functional_graph
from src.gnn_models.dataset import build_graph_dataset
from src.gnn_models.models import GATGraphClassifier
from src.gnn_models.train import TrainConfig, _class_weights

OUT  = ROOT / "outputs" / "final_comparison"
OUT.mkdir(parents=True, exist_ok=True)
METH = ROOT / "data" / "lihc_methylation_promoter.tsv"
SEEDS = list(range(1, 11))


def set_seeds(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)


def classify_stage(s):
    if not isinstance(s, str): return None
    s = s.upper().replace("STAGE ", "").strip()
    if s in {"I", "II", "IIA", "IIB"}: return 0
    if s.startswith("III") or s.startswith("IV"): return 1
    return None


def bh_correction(pvals):
    """Benjamini-Hochberg FDR correction."""
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = np.array(pvals)[order]
    fdr = ranked * n / (np.arange(n) + 1)
    fdr = np.minimum.accumulate(fdr[::-1])[::-1]
    result = np.empty(n)
    result[order] = np.clip(fdr, 0, 1)
    return result


# ═══════════════════════════════════════════════════════════════════════════ #
#  LOAD DATA                                                                  #
# ═══════════════════════════════════════════════════════════════════════════ #
def load_all(ds):
    meth = pd.read_csv(METH, sep="\t", index_col=0)
    expr = ds.expression   # rows=genes, cols=samples
    md   = ds.metadata.loc[ds.expression.columns].copy()

    # align samples — keep samples present in both
    shared = [s for s in expr.columns
              if s in meth.columns and s in md.index]
    print(f"  expression samples : {expr.shape[1]}")
    print(f"  methylation samples: {meth.shape[1]}")
    print(f"  shared samples     : {len(shared)}")
    print(f"  copper genes in meth: {meth.shape[0]}")

    expr_sh = expr[shared]
    meth_sh = meth[shared]
    md_sh   = md.loc[shared]
    return expr_sh, meth_sh, md_sh, shared


# ═══════════════════════════════════════════════════════════════════════════ #
#  1. EXPRESSION-METHYLATION CORRELATION                                      #
# ═══════════════════════════════════════════════════════════════════════════ #
def expr_meth_correlation(expr, meth, shared):
    print("\n[meth] === 1. Expression-Methylation Correlation ===")
    rows = []
    for gene in meth.index:
        if gene not in expr.index:
            continue
        e = expr.loc[gene, shared].to_numpy(dtype=float)
        m = meth.loc[gene, shared].to_numpy(dtype=float)
        mask = np.isfinite(e) & np.isfinite(m)
        if mask.sum() < 20:
            continue
        r, p = stats.pearsonr(e[mask], m[mask])
        rows.append({"gene": gene, "pearson_r": round(r, 4),
                     "p_value": p, "n": int(mask.sum())})

    df = pd.DataFrame(rows).sort_values("pearson_r")
    df["adj_p"] = bh_correction(df["p_value"].tolist())
    df["significant"] = df["adj_p"] < 0.05
    df["expected"]    = df["pearson_r"] < 0  # negative = expected (meth↑ expr↓)

    n_sig = int(df["significant"].sum())
    n_neg = int((df["pearson_r"] < 0).sum())
    print(f"  {n_sig}/{len(df)} genes significant (adj.p<0.05)")
    print(f"  {n_neg}/{len(df)} genes show negative correlation (expected direction)")

    print(f"\n  Top negatively correlated (hypomethylation → expression):")
    print(f"  {'gene':<12} {'r':>8} {'adj.p':>10} {'sig':>5}")
    print("  " + "-"*38)
    for _, r in df.head(10).iterrows():
        sig = "✓" if r["significant"] else ""
        print(f"  {r['gene']:<12} {r['pearson_r']:>8.4f} "
              f"{r['adj_p']:>10.2e} {sig:>5}")

    df.to_csv(OUT / "methylation_correlations.csv", index=False)
    return df


# ═══════════════════════════════════════════════════════════════════════════ #
#  2. STAGE-ASSOCIATED METHYLATION                                            #
# ═══════════════════════════════════════════════════════════════════════════ #
def stage_methylation(meth, md, shared):
    print("\n[meth] === 2. Stage-Associated Methylation ===")

    stage_bin = md["stage"].map(classify_stage)
    tumor_mask = md["sample_type"] == "Tumor"
    keep = tumor_mask & stage_bin.notna()
    md_t = md[keep].copy()
    md_t["stage_bin"] = stage_bin[keep]

    early_ids = [s for s in md_t[md_t["stage_bin"] == 0].index
                 if s in meth.columns]
    late_ids  = [s for s in md_t[md_t["stage_bin"] == 1].index
                 if s in meth.columns]

    print(f"  early samples with methylation: {len(early_ids)}")
    print(f"  late samples with methylation : {len(late_ids)}")

    rows = []
    for gene in meth.index:
        e_vals = meth.loc[gene, early_ids].dropna().to_numpy(dtype=float)
        l_vals = meth.loc[gene, late_ids].dropna().to_numpy(dtype=float)
        if len(e_vals) < 5 or len(l_vals) < 5:
            continue
        delta = l_vals.mean() - e_vals.mean()
        _, p  = stats.mannwhitneyu(e_vals, l_vals, alternative="two-sided")
        rows.append({
            "gene":         gene,
            "mean_early":   round(e_vals.mean(), 4),
            "mean_late":    round(l_vals.mean(), 4),
            "delta_meth":   round(delta, 4),
            "direction":    "hypermethylated in late" if delta > 0
                            else "hypomethylated in late",
            "p_value":      p,
        })

    df = pd.DataFrame(rows)
    df["adj_p"] = bh_correction(df["p_value"].tolist())
    df["significant"] = df["adj_p"] < 0.05
    df = df.sort_values("adj_p")

    n_sig = int(df["significant"].sum())
    print(f"  {n_sig}/{len(df)} genes significant (adj.p<0.05)")

    print(f"\n  Top stage-associated methylation changes:")
    print(f"  {'gene':<12} {'delta':>8} {'adj.p':>10} {'direction':<28} {'sig':>4}")
    print("  " + "-"*66)
    for _, r in df.head(15).iterrows():
        sig = "✓" if r["significant"] else ""
        print(f"  {r['gene']:<12} {r['delta_meth']:>8.4f} "
              f"{r['adj_p']:>10.2e} {r['direction']:<28} {sig:>4}")

    df.to_csv(OUT / "methylation_stage_diff.csv", index=False)
    return df, early_ids, late_ids


# ═══════════════════════════════════════════════════════════════════════════ #
#  3. MULTI-OMIC CONVERGENCE                                                  #
# ═══════════════════════════════════════════════════════════════════════════ #
def convergence_analysis(expr, meth, md, shared,
                          corr_df, stage_df, sal_top_genes):
    print("\n[meth] === 3. Multi-Omic Convergence ===")

    stage_bin  = md["stage"].map(classify_stage)
    tumor_mask = md["sample_type"] == "Tumor"
    keep       = tumor_mask & stage_bin.notna()
    md_t       = md[keep]

    early_ids = [s for s in md_t[md_t["stage_bin"] == 0].index
                 if s in shared]
    late_ids  = [s for s in md_t[md_t["stage_bin"] == 1].index
                 if s in shared]

    rows = []
    for gene in meth.index:
        if gene not in expr.index:
            continue

        # expression change: late vs early
        e_expr = expr.loc[gene, early_ids].to_numpy(dtype=float)
        l_expr = expr.loc[gene, late_ids].to_numpy(dtype=float)
        delta_expr = l_expr.mean() - e_expr.mean()

        # methylation change: late vs early
        e_meth_vals = meth.loc[gene, early_ids].dropna().to_numpy(dtype=float)
        l_meth_vals = meth.loc[gene, late_ids].dropna().to_numpy(dtype=float)
        delta_meth  = (l_meth_vals.mean() - e_meth_vals.mean()
                       if len(e_meth_vals) > 0 and len(l_meth_vals) > 0
                       else float("nan"))

        # convergence: expression up + methylation down = active in late stage
        #              expression down + methylation up = silenced in late stage
        if np.isfinite(delta_meth):
            converges = (delta_expr > 0 and delta_meth < 0) or \
                        (delta_expr < 0 and delta_meth > 0)
        else:
            converges = False

        # saliency info
        in_top_sal = gene in sal_top_genes
        corr_r = corr_df.set_index("gene").loc[gene, "pearson_r"] \
                 if gene in corr_df["gene"].values else float("nan")

        rows.append({
            "gene":          gene,
            "delta_expr":    round(delta_expr, 4),
            "delta_meth":    round(delta_meth, 4) if np.isfinite(delta_meth) else None,
            "expr_direction": "up in late" if delta_expr > 0 else "down in late",
            "meth_direction": ("hypo in late" if delta_meth < 0
                               else "hyper in late") if np.isfinite(delta_meth) else "—",
            "converges":     converges,
            "in_top_saliency": in_top_sal,
            "expr_meth_r":   round(corr_r, 4) if np.isfinite(corr_r) else None,
        })

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "methylation_convergence.csv", index=False)

    n_conv  = int(df["converges"].sum())
    n_total = len(df)
    n_sal_conv = int(df[df["in_top_saliency"]]["converges"].sum())
    n_sal_total = int(df["in_top_saliency"].sum())

    print(f"  All genes converging: {n_conv}/{n_total}")
    print(f"  Top saliency genes converging: {n_sal_conv}/{n_sal_total}")

    print(f"\n  Multi-omic convergence table (top saliency genes):")
    print(f"  {'gene':<12} {'Δexpr':>8} {'Δmeth':>8} "
          f"{'expr dir':<14} {'meth dir':<14} {'conv':>6}")
    print("  " + "-"*66)
    for _, r in df[df["in_top_saliency"]].iterrows():
        conv = "✓" if r["converges"] else "✗"
        dm   = f"{r['delta_meth']:>8.4f}" if r["delta_meth"] is not None else "     —  "
        print(f"  {r['gene']:<12} {r['delta_expr']:>8.4f} {dm} "
              f"{r['expr_direction']:<14} {r['meth_direction']:<14} {conv:>6}")

    return df


# ═══════════════════════════════════════════════════════════════════════════ #
#  4. MULTI-OMIC GAT (expression + methylation as node features)             #
# ═══════════════════════════════════════════════════════════════════════════ #
def build_multiomics_dataset(ds, meth, graph, shared, y_map):
    """Build PyG dataset with methylation added as extra node feature."""
    from src.gnn_models.dataset import graph_to_edge_tensors, _categorical_features
    from torch_geometric.data import Data

    gene_order = ds.expression.index.tolist()
    edge_index, edge_weight = graph_to_edge_tensors(graph, gene_order)
    expr = ds.expression[shared]

    # z-score expression
    mu = expr.mean(axis=1); sd = expr.std(axis=1).replace(0, 1)
    z_expr = expr.sub(mu, axis=0).div(sd, axis=0)

    # z-score methylation
    meth_shared = meth[shared]
    mu_m = meth_shared.mean(axis=1); sd_m = meth_shared.std(axis=1).replace(0, 1)
    z_meth = meth_shared.sub(mu_m, axis=0).div(sd_m, axis=0)

    cat_feats = _categorical_features(ds.copper_genes, gene_order)

    data_list = []
    for sid in shared:
        if sid not in y_map:
            continue
        e_col = z_expr[sid].to_numpy(dtype=np.float32).reshape(-1, 1)
        m_col = np.array([
            z_meth.loc[g, sid] if g in z_meth.index else 0.0
            for g in gene_order
        ], dtype=np.float32).reshape(-1, 1)
        x = np.concatenate([e_col, m_col, cat_feats], axis=1)
        data = Data(
            x=torch.tensor(x, dtype=torch.float32),
            edge_index=edge_index,
            edge_weight=edge_weight,
            y=torch.tensor([y_map[sid]], dtype=torch.long),
        )
        data.sample_id = sid
        data_list.append(data)

    in_dim = data_list[0].x.shape[1] if data_list else 0
    return data_list, in_dim


def run_multiomics_gat(ds, meth, graph, shared, md, seeds):
    print("\n[meth] === 4. Multi-Omic GAT (expression + methylation) ===")

    stage_bin  = md["stage"].map(classify_stage)
    tumor_mask = md["sample_type"] == "Tumor"
    keep       = tumor_mask & stage_bin.notna()
    md_t       = md[keep]

    # only samples with methylation data
    valid = [s for s in md_t.index if s in shared]
    y_map = {s: int(stage_bin[s]) for s in valid}
    y_arr = np.array([y_map[s] for s in valid])
    g_arr = md_t.loc[valid, "case_submitter_id"].to_numpy()

    data_list, in_dim = build_multiomics_dataset(
        ds, meth, graph, valid, y_map)

    print(f"  n={len(data_list)}  in_dim={in_dim}  "
          f"early={int((y_arr==0).sum())}  late={int((y_arr==1).sum())}")

    seed_aucs = []
    for seed in seeds:
        print(f"  seed {seed:2d}/{seeds[-1]} ...", end=" ", flush=True)
        set_seeds(seed)

        min_class = int(min((y_arr==0).sum(), (y_arr==1).sum()))
        n_splits  = max(2, min(5, min_class))
        cv  = StratifiedGroupKFold(n_splits=n_splits,
                                   shuffle=True, random_state=seed)
        cfg = TrainConfig(model="gat", epochs=80)

        counts = np.bincount(y_arr, minlength=2).astype(np.float32)
        w      = (counts.sum() / (2 * counts)) ** 1.5
        w      = w / w.mean()
        class_w = torch.tensor(w, dtype=torch.float32, device=cfg.device)

        fold_aucs = []
        for fold, (tr_idx, va_idx) in enumerate(
                cv.split(np.zeros(len(y_arr)), y_arr, g_arr)):
            set_seeds(seed + fold * 100)
            tr_list = [data_list[i] for i in tr_idx]
            va_list = [data_list[i] for i in va_idx]

            model = GATGraphClassifier(
                in_dim, 64, n_classes=2,
                n_heads=4, n_layers=2, dropout=0.4,
            ).to(cfg.device)
            opt     = torch.optim.Adam(model.parameters(),
                                       lr=cfg.lr, weight_decay=cfg.weight_decay)
            loss_fn = nn.CrossEntropyLoss(weight=class_w)
            tr_dl   = DataLoader(tr_list, batch_size=32, shuffle=True)
            va_dl   = DataLoader(va_list, batch_size=32, shuffle=False)

            best_auc, best_state = -1.0, None
            for _ in range(cfg.epochs):
                model.train()
                for batch in tr_dl:
                    batch = batch.to(cfg.device)
                    opt.zero_grad()
                    loss_fn(model(batch.x, batch.edge_index, batch.batch,
                                  edge_weight=getattr(batch, "edge_weight", None)),
                            batch.y).backward()
                    opt.step()
                model.eval()
                ys, ps = [], []
                with torch.no_grad():
                    for batch in va_dl:
                        batch = batch.to(cfg.device)
                        p = torch.softmax(
                            model(batch.x, batch.edge_index, batch.batch,
                                  edge_weight=getattr(batch, "edge_weight", None)),
                            dim=-1)[:, 1].cpu().numpy()
                        ys.extend(batch.y.cpu().numpy())
                        ps.extend(p)
                if len(set(ys)) > 1:
                    v = roc_auc_score(ys, ps)
                    if v > best_auc:
                        best_auc = v
                        best_state = {k: v2.clone()
                                      for k, v2 in model.state_dict().items()}

            fold_aucs.append(best_auc)

        m = float(np.mean(fold_aucs))
        print(f"AUC={m:.4f}")
        seed_aucs.append(m)

    return float(np.mean(seed_aucs)), float(np.std(seed_aucs))


# ═══════════════════════════════════════════════════════════════════════════ #
#  MAIN                                                                       #
# ═══════════════════════════════════════════════════════════════════════════ #
def main():
    print("[meth] loading dataset ...")
    ds    = load_lihc_dataset(require_real=True)
    graph = build_functional_graph(ds.expression.index.tolist(), ds.copper_genes)
    meth  = pd.read_csv(METH, sep="\t", index_col=0)

    expr, meth_sh, md, shared = load_all(ds)

    # top saliency genes from interpretability run
    sal_top = ["ATP7A","LTF","SCO2","S100A12","ATP7B",
               "SLC31A2","COX11","SLC31A1","SPARC","SCO1",
               "PRNP","H3-3A","AFP","S100A5","S100B",
               "SNCA","CUTA","MT3","ATOX1","MT4"]

    # ── analyses ──────────────────────────────────────────────────────────
    corr_df  = expr_meth_correlation(expr, meth_sh, shared)
    stage_df, early_ids, late_ids = stage_methylation(meth_sh, md, shared)
    conv_df  = convergence_analysis(expr, meth_sh, md, shared,
                                     corr_df, stage_df, sal_top)

    # expression-only GAT AUC for comparison (from previous run)
    expr_only_auc = 0.6667
    expr_only_std = 0.0174

    mo_auc, mo_std = run_multiomics_gat(ds, meth_sh, graph, shared, md, SEEDS)

    # ── final summary ──────────────────────────────────────────────────────
    print("\n" + "="*60)
    print("METHYLATION ANALYSIS SUMMARY")
    print("="*60)

    n_sig_corr  = int(corr_df["significant"].sum())
    n_neg_corr  = int((corr_df["pearson_r"] < 0).sum())
    n_sig_stage = int(stage_df["significant"].sum())
    n_conv      = int(conv_df["converges"].sum())
    n_sal_conv  = int(conv_df[conv_df["in_top_saliency"]]["converges"].sum())
    n_sal_total = int(conv_df["in_top_saliency"].sum())

    improvement = mo_auc - expr_only_auc
    print(f"\n1. Expression-methylation correlation")
    print(f"   {n_sig_corr}/54 genes significant  |  "
          f"{n_neg_corr}/54 negative (expected direction)")
    print(f"\n2. Stage-associated methylation")
    print(f"   {n_sig_stage}/54 genes differ between early and late stage")
    print(f"\n3. Multi-omic convergence (expr + meth both change)")
    print(f"   All genes: {n_conv}/54  |  "
          f"Top saliency genes: {n_sal_conv}/{n_sal_total}")
    print(f"\n4. Multi-omic GAT (expression + methylation features)")
    print(f"   Expression only : {expr_only_auc:.4f} ± {expr_only_std:.4f}")
    print(f"   + Methylation   : {mo_auc:.4f} ± {mo_std:.4f}")
    print(f"   Improvement     : {improvement:+.4f}")

    if improvement > 0.01:
        verdict = "Methylation improves prediction — multi-omic features add value"
    elif improvement > -0.01:
        verdict = "Methylation neutral — expression alone sufficient for prediction"
    else:
        verdict = "Methylation hurts — adds noise, expression features are cleaner"
    print(f"   Verdict: {verdict}")

    # write report
    report = f"""# Methylation Analysis — Multi-Omic Validation

## Data
- Promoter methylation beta values for all 54 copper genes
- {len(shared)} samples with both expression and methylation
- Beta values 0–1 (0 = unmethylated, 1 = fully methylated)

## 1. Expression-Methylation Correlation
- **{n_sig_corr}/54** genes show significant correlation (adj.p < 0.05)
- **{n_neg_corr}/54** genes show negative correlation (expected: hypomethylation → expression)
- Negative correlation = methylation regulates expression as expected

## 2. Stage-Associated Methylation Changes
- **{n_sig_stage}/54** genes show significant methylation differences between early and late stage tumors
- These genes are candidates for epigenetic regulation of stage progression

## 3. Multi-Omic Convergence
- A gene converges if expression and methylation change in opposite directions between early and late stage (active gene = high expression + low methylation)
- All copper genes: **{n_conv}/54** converge
- Top GAT saliency genes: **{n_sal_conv}/{n_sal_total}** converge
- Convergence in saliency genes = the model's important genes are also epigenetically regulated

## 4. Stage Classification: Expression Only vs Multi-Omic
| model | ROC-AUC | std |
|---|---:|---:|
| GAT (expression only) | {expr_only_auc:.4f} | {expr_only_std:.4f} |
| GAT (expression + methylation) | {mo_auc:.4f} | {mo_std:.4f} |
| Improvement | {improvement:+.4f} | — |

**{verdict}**

## Biological Interpretation
Genes that appear in all three analyses (important for the GAT,
show stage-associated methylation changes, AND show expression-methylation
correlation) are the strongest candidates for epigenetic regulation of
copper metabolism in HCC progression.
"""
    (OUT / "methylation_report.md").write_text(report, encoding="utf-8")
    print(f"\n[meth] wrote all outputs to outputs/final_comparison/")
    print("Done.")


if __name__ == "__main__":
    main()
