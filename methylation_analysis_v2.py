"""Methylation analysis v2 — clean rewrite with sample ID fix.

Sample ID mismatch:
  expression: 'TCGA-DD-AACI-01A'  (ends in letter)
  methylation: 'TCGA-DD-AACI-01'  (no trailing letter)
Fix: strip last char from expression IDs to get base ID.

Four analyses:
  1. Expression-methylation correlation per gene
  2. Stage-associated methylation changes (early vs late)
  3. Multi-omic convergence (GAT important genes)
  4. Multi-omic GAT (expression + methylation as node features)
"""
from __future__ import annotations
import sys, random
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import torch
import torch.nn as nn
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.metrics import roc_auc_score
from torch_geometric.loader import DataLoader
from torch_geometric.data import Data

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from src.preprocessing import load_lihc_dataset
from src.graph_building import build_functional_graph
from src.gnn_models.dataset import graph_to_edge_tensors, _categorical_features
from src.gnn_models.models import GATGraphClassifier
from src.gnn_models.train import TrainConfig

OUT  = ROOT / "outputs" / "final_comparison"
OUT.mkdir(parents=True, exist_ok=True)
METH_FILE = ROOT / "data" / "lihc_methylation_promoter.tsv"
SEEDS = list(range(1, 11))

# top saliency genes from interpretability_v2.py
SAL_TOP = ["ATP7A","LTF","SCO2","S100A12","ATP7B","SLC31A2","COX11",
           "SLC31A1","SPARC","SCO1","PRNP","H3-3A","AFP","S100A5",
           "S100B","SNCA","CUTA","MT3","ATOX1","MT4"]


def set_seeds(s):
    random.seed(s); np.random.seed(s); torch.manual_seed(s)

def classify_stage(s):
    if not isinstance(s, str): return None
    s = s.upper().replace("STAGE ","").strip()
    if s in {"I","II","IIA","IIB"}: return 0
    if s.startswith("III") or s.startswith("IV"): return 1
    return None

def bh_fdr(pvals):
    n = len(pvals); order = np.argsort(pvals)
    ranked = np.array(pvals)[order]
    fdr = ranked * n / (np.arange(n)+1)
    fdr = np.minimum.accumulate(fdr[::-1])[::-1]
    out = np.empty(n); out[order] = np.clip(fdr, 0, 1)
    return out


# ── load & align ──────────────────────────────────────────────────────────
def load_aligned(ds):
    meth = pd.read_csv(METH_FILE, sep="\t", index_col=0)
    expr = ds.expression          # rows=genes, cols=expr_ids ('TCGA-XX-XXXX-01A')
    md   = ds.metadata.loc[ds.expression.columns].copy()

    # strip last char from expr IDs: 'TCGA-DD-AACI-01A' → 'TCGA-DD-AACI-01'
    base = {eid: eid[:-1] for eid in expr.columns}
    meth_set = set(meth.columns)

    paired_expr = [eid for eid in expr.columns if base[eid] in meth_set]
    paired_meth = [base[eid] for eid in paired_expr]

    print(f"  expression samples : {expr.shape[1]}")
    print(f"  methylation samples: {meth.shape[1]}")
    print(f"  shared (paired)    : {len(paired_expr)}")

    if len(paired_expr) == 0:
        raise ValueError("No shared samples. Check ID format.")

    expr_al = expr[paired_expr].copy()
    expr_al.columns = paired_meth       # rename to meth IDs
    meth_al = meth[paired_meth].copy()
    md_al   = md.loc[paired_expr].copy()
    md_al.index = paired_meth           # rename index to meth IDs

    return expr_al, meth_al, md_al, paired_meth


# ── 1. expression-methylation correlation ─────────────────────────────────
def corr_analysis(expr, meth, shared):
    print("\n[meth] === 1. Expression-Methylation Correlation ===")
    rows = []
    for gene in meth.index:
        if gene not in expr.index: continue
        e = expr.loc[gene, shared].to_numpy(dtype=float)
        m = meth.loc[gene, shared].to_numpy(dtype=float)
        ok = np.isfinite(e) & np.isfinite(m)
        if ok.sum() < 20: continue
        r, p = stats.pearsonr(e[ok], m[ok])
        rows.append({"gene":gene,"pearson_r":round(r,4),"p_value":p,"n":int(ok.sum())})
    df = pd.DataFrame(rows).sort_values("pearson_r")
    df["adj_p"] = bh_fdr(df["p_value"].tolist())
    df["sig"]   = df["adj_p"] < 0.05
    n_sig = int(df["sig"].sum())
    n_neg = int((df["pearson_r"]<0).sum())
    print(f"  {n_sig}/54 significant | {n_neg}/54 negative correlation (expected)")
    print(f"  {'gene':<12} {'r':>8} {'adj.p':>10} {'sig':>5}")
    print("  "+"-"*38)
    for _,r in df.head(10).iterrows():
        print(f"  {r['gene']:<12} {r['pearson_r']:>8.4f} "
              f"{r['adj_p']:>10.2e} {'✓' if r['sig'] else '':>5}")
    df.to_csv(OUT/"methylation_correlations.csv", index=False)
    return df


# ── 2. stage methylation differences ──────────────────────────────────────
def stage_meth(meth, md, shared):
    print("\n[meth] === 2. Stage-Associated Methylation ===")
    sb  = md["stage"].map(classify_stage)
    tum = md["sample_type"] == "Tumor"
    keep = tum & sb.notna()
    mdt = md[keep].copy(); mdt["sb"] = sb[keep]
    early = [s for s in mdt[mdt["sb"]==0].index if s in meth.columns]
    late  = [s for s in mdt[mdt["sb"]==1].index if s in meth.columns]
    print(f"  early={len(early)}  late={len(late)}")
    rows = []
    for gene in meth.index:
        ev = meth.loc[gene,early].dropna().to_numpy(dtype=float)
        lv = meth.loc[gene,late].dropna().to_numpy(dtype=float)
        if len(ev)<5 or len(lv)<5: continue
        delta = lv.mean()-ev.mean()
        _,p = stats.mannwhitneyu(ev,lv,alternative="two-sided")
        rows.append({"gene":gene,"mean_early":round(ev.mean(),4),
                     "mean_late":round(lv.mean(),4),"delta":round(delta,4),
                     "direction":"hyper in late" if delta>0 else "hypo in late",
                     "p_value":p})
    df = pd.DataFrame(rows)
    df["adj_p"] = bh_fdr(df["p_value"].tolist())
    df["sig"]   = df["adj_p"] < 0.05
    df = df.sort_values("adj_p")
    n_sig = int(df["sig"].sum())
    print(f"  {n_sig}/54 genes differ between early and late (adj.p<0.05)")
    print(f"  {'gene':<12} {'delta':>8} {'adj.p':>10} {'direction':<20} {'sig':>4}")
    print("  "+"-"*58)
    for _,r in df.head(15).iterrows():
        print(f"  {r['gene']:<12} {r['delta']:>8.4f} {r['adj_p']:>10.2e} "
              f"{r['direction']:<20} {'✓' if r['sig'] else '':>4}")
    df.to_csv(OUT/"methylation_stage_diff.csv", index=False)
    return df, early, late


# ── 3. convergence ────────────────────────────────────────────────────────
def convergence(expr, meth, md, shared, corr_df):
    print("\n[meth] === 3. Multi-Omic Convergence ===")
    sb  = md["stage"].map(classify_stage)
    tum = md["sample_type"] == "Tumor"
    keep = tum & sb.notna()
    mdt = md[keep]
    early = [s for s in mdt[mdt["sb"] if "sb" in mdt.columns
             else sb[keep]==0].index if s in shared]
    # recompute cleanly
    mdt2 = md[tum & sb.notna()].copy(); mdt2["sb"] = sb[tum & sb.notna()]
    early = [s for s in mdt2[mdt2["sb"]==0].index if s in shared]
    late  = [s for s in mdt2[mdt2["sb"]==1].index if s in shared]

    corr_map = corr_df.set_index("gene")["pearson_r"].to_dict() \
               if "pearson_r" in corr_df.columns else {}
    rows = []
    for gene in meth.index:
        if gene not in expr.index: continue
        de = expr.loc[gene,late].mean() - expr.loc[gene,early].mean()
        ev = meth.loc[gene,early].dropna().to_numpy(dtype=float)
        lv = meth.loc[gene,late].dropna().to_numpy(dtype=float)
        dm = lv.mean()-ev.mean() if len(ev)>0 and len(lv)>0 else float("nan")
        conv = (np.isfinite(dm) and
                ((de>0 and dm<0) or (de<0 and dm>0)))
        rows.append({"gene":gene,
                     "delta_expr":round(de,4),
                     "delta_meth":round(dm,4) if np.isfinite(dm) else None,
                     "expr_dir":"up in late" if de>0 else "down in late",
                     "meth_dir":("hypo in late" if dm<0 else "hyper in late")
                                 if np.isfinite(dm) else "—",
                     "converges":conv,
                     "in_top_sal":gene in SAL_TOP,
                     "r":round(corr_map.get(gene,float("nan")),4)})
    df = pd.DataFrame(rows)
    df.to_csv(OUT/"methylation_convergence.csv", index=False)
    nc = int(df["converges"].sum())
    ns = int(df[df["in_top_sal"]]["converges"].sum())
    nt = int(df["in_top_sal"].sum())
    print(f"  All genes: {nc}/54 converge")
    print(f"  Top saliency genes: {ns}/{nt} converge")
    print(f"\n  {'gene':<12} {'Δexpr':>8} {'Δmeth':>8} "
          f"{'expr':>12} {'meth':>14} {'conv':>6}")
    print("  "+"-"*64)
    for _,r in df[df["in_top_sal"]].iterrows():
        dm = f"{r['delta_meth']:>8.4f}" if r["delta_meth"] is not None else "       —"
        print(f"  {r['gene']:<12} {r['delta_expr']:>8.4f} {dm} "
              f"{r['expr_dir']:>12} {r['meth_dir']:>14} "
              f"{'✓' if r['converges'] else '✗':>6}")
    return df, nc, ns, nt


# ── 4. multi-omic GAT ─────────────────────────────────────────────────────
def run_multiomics_gat(expr, meth, md, graph, gene_order, copper_genes, seeds):
    print("\n[meth] === 4. Multi-Omic GAT (expression + methylation) ===")
    sb  = md["stage"].map(classify_stage)
    tum = md["sample_type"] == "Tumor"
    keep = tum & sb.notna()
    mdt = md[keep].copy(); mdt["sb"] = sb[keep]

    mdt = mdt[~mdt.index.duplicated(keep="first")]
    valid = [s for s in mdt.index if s in meth.columns and s in expr.columns]
    y_map = {s: int(mdt.at[s, "sb"]) for s in valid}
    y_arr = np.array([y_map[s] for s in valid])
    g_arr = np.array([mdt.at[s, "case_submitter_id"] for s in valid])

    print(f"  n={len(valid)}  early={int((y_arr==0).sum())}  late={int((y_arr==1).sum())}")

    # build dataset
    ei, ew = graph_to_edge_tensors(graph, gene_order)
    cat    = _categorical_features(copper_genes, gene_order)

    # z-score both modalities
    # deduplicate both rows AND columns before z-scoring
    expr_v = expr[valid].T.drop_duplicates().T   # dedup sample cols
    meth_v = meth[valid].T.drop_duplicates().T
    valid  = [s for s in valid if s in expr_v.columns and s in meth_v.columns]
    expr_v = expr_v[valid]; meth_v = meth_v[valid]

    e_mu = expr_v.mean(axis=1); e_sd = expr_v.std(axis=1).replace(0,1)
    ze   = expr_v.sub(e_mu,axis=0).div(e_sd,axis=0)
    ze   = ze[~ze.index.duplicated(keep='first')]
    m_mu = meth_v.mean(axis=1); m_sd = meth_v.std(axis=1).replace(0,1)
    zm   = meth_v.sub(m_mu,axis=0).div(m_sd,axis=0)
    zm   = zm[~zm.index.duplicated(keep='first')]
    y_arr = np.array([y_map[s] for s in valid])
    g_arr = np.array([mdt.at[s,'case_submitter_id'] for s in valid])

    data_list = []
    for sid in valid:
        ec = ze[sid].to_numpy(dtype=np.float32).reshape(-1,1)
        def _get(g, sid):
            if g not in zm.index: return 0.0
            v = zm.loc[g, sid]
            return float(v.iloc[0]) if hasattr(v, 'iloc') else float(v)
        mc = np.array([_get(g, sid) for g in gene_order],
                      dtype=np.float32).reshape(-1,1)
        x  = np.concatenate([ec, mc, cat], axis=1)
        d  = Data(x=torch.tensor(x,dtype=torch.float32),
                  edge_index=ei, edge_weight=ew,
                  y=torch.tensor([y_map[sid]],dtype=torch.long))
        data_list.append(d)

    in_dim = data_list[0].x.shape[1]

    # class weights
    counts = np.bincount(y_arr, minlength=2).astype(np.float32)
    w = (counts.sum()/(2*counts))**1.5; w = w/w.mean()

    seed_aucs = []
    for seed in seeds:
        print(f"  seed {seed:2d}/{seeds[-1]} ...", end=" ", flush=True)
        set_seeds(seed)
        min_c = int(min((y_arr==0).sum(),(y_arr==1).sum()))
        cv = StratifiedGroupKFold(n_splits=max(2,min(5,min_c)),
                                  shuffle=True, random_state=seed)
        cfg = TrainConfig(model="gat", epochs=80)
        cw  = torch.tensor(w, dtype=torch.float32, device=cfg.device)
        fold_aucs = []
        for fold,(tr,va) in enumerate(cv.split(np.zeros(len(y_arr)),y_arr,g_arr)):
            set_seeds(seed+fold*100)
            model = GATGraphClassifier(in_dim,64,n_classes=2,
                                       n_heads=4,n_layers=2,dropout=0.4).to(cfg.device)
            opt     = torch.optim.Adam(model.parameters(),lr=cfg.lr,
                                       weight_decay=cfg.weight_decay)
            lf      = nn.CrossEntropyLoss(weight=cw)
            trdl    = DataLoader([data_list[i] for i in tr],batch_size=32,shuffle=True)
            vadl    = DataLoader([data_list[i] for i in va],batch_size=32,shuffle=False)
            best=-1.0
            for _ in range(cfg.epochs):
                model.train()
                for b in trdl:
                    b=b.to(cfg.device); opt.zero_grad()
                    lf(model(b.x,b.edge_index,b.batch,
                             edge_weight=getattr(b,"edge_weight",None)),b.y).backward()
                    opt.step()
                model.eval(); ys,ps=[],[]
                with torch.no_grad():
                    for b in vadl:
                        b=b.to(cfg.device)
                        p=torch.softmax(model(b.x,b.edge_index,b.batch,
                                             edge_weight=getattr(b,"edge_weight",None)),
                                        dim=-1)[:,1].cpu().numpy()
                        ys.extend(b.y.cpu().numpy()); ps.extend(p)
                if len(set(ys))>1:
                    v=roc_auc_score(ys,ps)
                    if v>best: best=v
            fold_aucs.append(best)
        m=float(np.mean(fold_aucs))
        print(f"AUC={m:.4f}")
        seed_aucs.append(m)
    return float(np.mean(seed_aucs)), float(np.std(seed_aucs))


# ── main ──────────────────────────────────────────────────────────────────
def main():
    print("[meth] loading dataset ...")
    ds    = load_lihc_dataset(require_real=True)
    graph = build_functional_graph(ds.expression.index.tolist(), ds.copper_genes)
    gene_order = ds.expression.index.tolist()

    expr, meth, md, shared = load_aligned(ds)

    corr_df = corr_analysis(expr, meth, shared)
    stage_df, early, late = stage_meth(meth, md, shared)
    conv_df, nc, ns, nt   = convergence(expr, meth, md, shared, corr_df)

    expr_only_auc = 0.6667; expr_only_std = 0.0174
    mo_auc, mo_std = run_multiomics_gat(
        expr, meth, md, graph, gene_order, ds.copper_genes, SEEDS)

    imp = mo_auc - expr_only_auc
    verdict = ("Methylation improves prediction — multi-omic adds value" if imp>0.01
               else "Methylation neutral — expression alone sufficient" if imp>-0.01
               else "Methylation hurts — adds noise")

    print("\n"+"="*60)
    print("METHYLATION SUMMARY")
    print("="*60)
    n_sig_c = int(corr_df["sig"].sum())
    n_neg_c = int((corr_df["pearson_r"]<0).sum())
    n_sig_s = int(stage_df["sig"].sum())
    print(f"  Expr-meth correlation : {n_sig_c}/54 sig | {n_neg_c}/54 negative")
    print(f"  Stage methylation     : {n_sig_s}/54 genes differ early vs late")
    print(f"  Convergence (all)     : {nc}/54")
    print(f"  Convergence (top sal) : {ns}/{nt}")
    print(f"  Expression only AUC   : {expr_only_auc:.4f} ± {expr_only_std:.4f}")
    print(f"  + Methylation AUC     : {mo_auc:.4f} ± {mo_std:.4f}")
    print(f"  Improvement           : {imp:+.4f}")
    print(f"  Verdict               : {verdict}")

    report = f"""# Methylation Analysis

## Sample alignment
- Expression IDs end in letter (e.g. TCGA-DD-AACI-01A)
- Methylation IDs do not (e.g. TCGA-DD-AACI-01)
- Aligned by stripping last character: {len(shared)} shared samples

## 1. Expression-Methylation Correlation
- {n_sig_c}/54 genes significant (adj.p<0.05)
- {n_neg_c}/54 genes negative correlation (hypomethylation → expression, expected)

## 2. Stage-Associated Methylation
- {n_sig_s}/54 genes differ between early and late stage tumors

## 3. Multi-Omic Convergence
- All copper genes: {nc}/54 converge (expression + methylation both change)
- Top GAT saliency genes: {ns}/{nt} converge

## 4. Multi-Omic GAT
| model | AUC | std |
|---|---:|---:|
| Expression only | {expr_only_auc:.4f} | {expr_only_std:.4f} |
| + Methylation | {mo_auc:.4f} | {mo_std:.4f} |
| Improvement | {imp:+.4f} | — |

**{verdict}**
"""
    (OUT/"methylation_report.md").write_text(report, encoding="utf-8")
    print("\n[meth] all outputs written. Done.")


if __name__ == "__main__":
    main()
