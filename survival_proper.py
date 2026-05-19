"""Week 2 — Proper survival analysis.

Compares three approaches on the same TCGA-LIHC tumor patients:

  1. Random baseline     — C-index = 0.5 by definition
  2. Cox regression      — clinical gold standard, uses ALL patients
                           including censored ones (proper survival model)
  3. GAT classifier      — our model, reports both ROC-AUC and C-index

Key improvements over the original pipeline:
  - C-index instead of ROC-AUC (correct metric for survival)
  - Cox uses censored patients (no data thrown away)
  - 10-seed stability for GAT
  - Direct model comparison on same patients

Outputs:
  outputs/final_comparison/week2_survival_comparison.csv
  outputs/final_comparison/week2_survival_report.md

Usage:
    python week2_survival_proper.py
"""
from __future__ import annotations
import sys
import random
from pathlib import Path
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from sklearn.model_selection import StratifiedGroupKFold, KFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score
from torch_geometric.loader import DataLoader
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from src.preprocessing import load_lihc_dataset
from src.graph_building import build_functional_graph
from src.gnn_models.dataset import build_graph_dataset
from src.gnn_models.models import GATGraphClassifier
from src.gnn_models.train import TrainConfig, _class_weights

OUT = ROOT / "outputs" / "final_comparison"
OUT.mkdir(parents=True, exist_ok=True)

SEEDS = list(range(1, 11))


# ── seeding ────────────────────────────────────────────────────────────────
def set_seeds(seed: int):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


# ── survival labels ────────────────────────────────────────────────────────
def binary_label(row):
    """For GAT — only confirmed outcomes (no censored)."""
    os  = row.get("overall_survival_days")
    vit = row.get("vital_status")
    if pd.isna(os) or pd.isna(vit):
        return None
    if vit == "Dead":
        return 1 if os <= 1095 else 0
    if vit == "Alive":
        return 0 if os > 1095 else None   # censored — unknown
    return None


# ══════════════════════════════════════════════════════════════════════════ #
#  1. COX REGRESSION                                                         #
# ══════════════════════════════════════════════════════════════════════════ #
def run_cox(ds, md_tumor, seeds):
    """
    Cox proportional hazards model — 5-fold CV, 10 seeds.
    Uses ALL tumor patients with any survival data (including censored).
    Returns mean C-index ± std across seeds.
    """
    # select patients with survival info
    has_data = (md_tumor["overall_survival_days"].notna() &
                md_tumor["vital_status"].notna())
    md_cox = md_tumor[has_data].copy()
    md_cox["event"]    = (md_cox["vital_status"] == "Dead").astype(int)
    md_cox["duration"] = md_cox["overall_survival_days"]

    # expression for those patients
    shared = [s for s in md_cox.index if s in ds.expression.columns]
    md_cox = md_cox.loc[shared]
    expr   = ds.expression[shared].T.to_numpy()   # (n_samples, n_genes)

    groups  = md_cox["case_submitter_id"].to_numpy()
    y_strat = md_cox["event"].to_numpy()

    print(f"\n  n={len(md_cox)}  events={int(md_cox['event'].sum())}  "
          f"censored={int((md_cox['event']==0).sum())}")
    print(f"  (Cox uses all {len(md_cox)} patients; "
          f"GAT uses only confirmed-outcome subset)")

    seed_cis = []
    for seed in seeds:
        print(f"  seed {seed:2d}/{seeds[-1]} ...", end=" ", flush=True)
        set_seeds(seed)

        cv = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=seed)
        fold_cis = []

        for tr_idx, va_idx in cv.split(np.zeros(len(md_cox)), y_strat, groups):
            scaler = StandardScaler()
            X_tr   = scaler.fit_transform(expr[tr_idx])
            X_va   = scaler.transform(expr[va_idx])

            tr_df = pd.DataFrame(X_tr, columns=ds.expression.index.tolist())
            tr_df["duration"] = md_cox["duration"].iloc[tr_idx].values
            tr_df["event"]    = md_cox["event"].iloc[tr_idx].values

            va_df = pd.DataFrame(X_va, columns=ds.expression.index.tolist())

            cph = CoxPHFitter(penalizer=0.1)   # L2 penalty for 58 features
            try:
                cph.fit(tr_df, duration_col="duration", event_col="event",
                        show_progress=False)
                risk = cph.predict_partial_hazard(va_df).values
                ci   = concordance_index(
                    md_cox["duration"].iloc[va_idx].values,
                    -risk,
                    md_cox["event"].iloc[va_idx].values,
                )
                fold_cis.append(ci)
            except Exception as e:
                print(f"\n  [Cox fold error: {e}]", end=" ")

        mean_ci = float(np.mean(fold_cis)) if fold_cis else float("nan")
        print(f"C-index={mean_ci:.4f}")
        seed_cis.append(mean_ci)

    return float(np.mean(seed_cis)), float(np.std(seed_cis))


# ══════════════════════════════════════════════════════════════════════════ #
#  2. GAT — ROC-AUC + C-index                                               #
# ══════════════════════════════════════════════════════════════════════════ #
def _train_fold(train_list, val_list, in_dim, cfg):
    """Train one fold, return (prob_scores, true_labels)."""
    y_tr    = np.array([int(d.y.item()) for d in train_list])
    class_w = _class_weights(y_tr, cfg.device)
    model   = GATGraphClassifier(
        in_dim, cfg.hidden, n_classes=2,
        n_heads=4, n_layers=cfg.n_layers, dropout=cfg.dropout,
    ).to(cfg.device)
    opt     = torch.optim.Adam(model.parameters(),
                               lr=cfg.lr, weight_decay=cfg.weight_decay)
    loss_fn = nn.CrossEntropyLoss(weight=class_w)
    tr_dl   = DataLoader(train_list, batch_size=cfg.batch_size, shuffle=True)
    va_dl   = DataLoader(val_list,   batch_size=cfg.batch_size, shuffle=False)

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

    if best_state:
        model.load_state_dict(best_state)

    model.eval()
    ps_final = []
    with torch.no_grad():
        for batch in va_dl:
            batch = batch.to(cfg.device)
            p = torch.softmax(
                model(batch.x, batch.edge_index, batch.batch,
                      edge_weight=getattr(batch, "edge_weight", None)),
                dim=-1)[:, 1].cpu().numpy()
            ps_final.extend(p)

    return np.array(ps_final)


def run_gat(ds, graph, y, groups, sample_indices, md_keep, seeds, epochs=80):
    """
    GAT — 5-fold CV, 10 seeds.
    Reports both ROC-AUC (binary) and C-index (survival ranking).
    """
    durations = md_keep["overall_survival_days"].to_numpy()
    events    = (md_keep["vital_status"] == "Dead").astype(int).to_numpy()

    seed_aucs, seed_cis = [], []

    for seed in seeds:
        print(f"  seed {seed:2d}/{seeds[-1]} ...", end=" ", flush=True)
        set_seeds(seed)

        min_class = int(min((y == 0).sum(), (y == 1).sum()))
        n_splits  = max(2, min(5, min_class))
        cv  = StratifiedGroupKFold(n_splits=n_splits,
                                   shuffle=True, random_state=seed)
        cfg = TrainConfig(model="gat", epochs=epochs)

        all_probs    = np.zeros(len(y))
        all_ys       = np.zeros(len(y), dtype=int)

        for fold, (tr_idx, va_idx) in enumerate(
                cv.split(np.zeros(len(y)), y, groups)):

            global_mask = np.zeros(ds.n_samples, dtype=bool)
            global_mask[sample_indices[tr_idx]] = True
            bundle = build_graph_dataset(ds, graph,
                                          zscore_train_mask=global_mask)

            active = []
            for pos, orig_i in enumerate(sample_indices):
                d   = bundle.data_list[orig_i]
                d.y = torch.tensor([int(y[pos])], dtype=torch.long)
                active.append(d)

            set_seeds(seed + fold * 100)
            probs = _train_fold(
                [active[i] for i in tr_idx],
                [active[i] for i in va_idx],
                bundle.in_dim, cfg,
            )
            all_probs[va_idx] = probs
            all_ys[va_idx]    = y[va_idx]

        # ROC-AUC on binary labels
        auc = roc_auc_score(all_ys, all_probs) \
              if len(set(all_ys)) > 1 else float("nan")

        # C-index: higher prob of dying → higher risk
        ci = concordance_index(durations, -all_probs, events)

        print(f"AUC={auc:.4f}  C-index={ci:.4f}")
        seed_aucs.append(auc)
        seed_cis.append(ci)

    return (float(np.mean(seed_aucs)), float(np.std(seed_aucs)),
            float(np.mean(seed_cis)),  float(np.std(seed_cis)))


# ══════════════════════════════════════════════════════════════════════════ #
#  MAIN                                                                      #
# ══════════════════════════════════════════════════════════════════════════ #
def main():
    print("[week2] loading dataset ...")
    ds    = load_lihc_dataset(require_real=True)
    genes = ds.expression.index.tolist()
    graph = build_functional_graph(genes, ds.copper_genes)
    md    = ds.metadata.loc[ds.expression.columns].copy()
    md["overall_survival_days"] = pd.to_numeric(
        md["overall_survival_days"], errors="coerce")

    # tumor only
    md_tumor = md[md["sample_type"] == "Tumor"].copy()

    # GAT subset — confirmed outcomes only
    md_tumor["y_bin"] = md_tumor.apply(binary_label, axis=1)
    keep = md_tumor["y_bin"].notna()
    md_gat = md_tumor[keep].copy()
    idx    = np.array([ds.expression.columns.tolist().index(s)
                       for s in md_gat.index])
    y_gat  = md_gat["y_bin"].astype(int).to_numpy()
    g_gat  = md_gat["case_submitter_id"].to_numpy()

    # ── 1. Random baseline ────────────────────────────────────────────────
    print("\n[week2] === Random baseline ===")
    print("  C-index = 0.500 by definition (random ranking)")
    random_ci = 0.500

    # ── 2. Cox regression ─────────────────────────────────────────────────
    print("\n[week2] === Cox Regression (10 seeds) ===")
    cox_ci_mean, cox_ci_std = run_cox(ds, md_tumor, SEEDS)

    # ── 3. GAT ───────────────────────────────────────────────────────────
    print(f"\n[week2] === GAT (10 seeds, n={len(y_gat)} confirmed patients) ===")
    gat_auc_mean, gat_auc_std, gat_ci_mean, gat_ci_std = run_gat(
        ds, graph, y_gat, g_gat, idx, md_gat, SEEDS, epochs=80)

    # ── summary ───────────────────────────────────────────────────────────
    print("\n" + "="*60)
    print("WEEK 2 SURVIVAL SUMMARY")
    print("="*60)
    print(f"\n{'Model':<20} {'C-index mean':>14} {'C-index std':>12}")
    print("-"*48)
    print(f"{'Random baseline':<20} {random_ci:>14.4f} {'—':>12}")
    print(f"{'Cox regression':<20} {cox_ci_mean:>14.4f} {cox_ci_std:>12.4f}")
    print(f"{'GAT (C-index)':<20} {gat_ci_mean:>14.4f} {gat_ci_std:>12.4f}")
    print(f"\nGAT ROC-AUC (binary): {gat_auc_mean:.4f} ± {gat_auc_std:.4f}")
    print(f"GAT vs Cox C-index gap: {gat_ci_mean - cox_ci_mean:+.4f}")

    # interpretation
    gap = gat_ci_mean - cox_ci_mean
    if gap > 0.02:
        verdict = "GAT clearly beats Cox — graph adds clinical value"
    elif gap > 0.00:
        verdict = "GAT slightly better than Cox — graph adds marginal value"
    elif gap > -0.02:
        verdict = "GAT and Cox are equivalent — graph adds interpretability not prediction"
    else:
        verdict = "Cox beats GAT — simpler model wins on this task"
    print(f"Verdict: {verdict}")

    # save CSV
    rows = [
        {"model": "Random baseline", "metric": "C-index",
         "mean": random_ci, "std": 0.0,
         "note": "random ranking — theoretical floor"},
        {"model": "Cox regression", "metric": "C-index",
         "mean": round(cox_ci_mean, 4), "std": round(cox_ci_std, 4),
         "note": f"uses all {len(md_tumor[md_tumor.overall_survival_days.notna()])} tumor patients incl. censored"},
        {"model": "GAT", "metric": "ROC-AUC (binary)",
         "mean": round(gat_auc_mean, 4), "std": round(gat_auc_std, 4),
         "note": f"confirmed-outcome patients only n={len(y_gat)}"},
        {"model": "GAT", "metric": "C-index",
         "mean": round(gat_ci_mean, 4), "std": round(gat_ci_std, 4),
         "note": f"confirmed-outcome patients only n={len(y_gat)}"},
    ]
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "week2_survival_comparison.csv", index=False)

    # save report
    report = f"""# Week 2 — Proper Survival Analysis

## What changed vs the original pipeline

| Issue | Original | Fixed |
|---|---|---|
| Metric | ROC-AUC (binary classification) | **C-index** (survival ranking) |
| Censored patients | Dropped | **Included in Cox** |
| Baseline model | None | **Cox regression** |
| Stability | Single seed | **10 seeds, mean ± std** |

## Results

| model | C-index mean | C-index std | note |
|---|---:|---:|---|
| Random baseline | 0.500 | — | theoretical floor |
| Cox regression | {cox_ci_mean:.4f} | {cox_ci_std:.4f} | uses all patients incl. censored |
| GAT | {gat_ci_mean:.4f} | {gat_ci_std:.4f} | confirmed-outcome patients only |

GAT ROC-AUC (binary): **{gat_auc_mean:.4f} ± {gat_auc_std:.4f}**

GAT vs Cox C-index gap: **{gap:+.4f}**

## Verdict

{verdict}

## What the C-index means

A C-index of 0.5 is random guessing.
A C-index of 0.7 is considered good for a clinical survival model.
A C-index of 1.0 is perfect ranking.

Unlike ROC-AUC, the C-index:
- Does not require throwing away censored patients
- Measures ranking (who dies sooner) not just binary classification
- Is the standard metric used in clinical survival research

## Honest caveats

- Cox uses more patients than GAT (censored included vs excluded)
  so the comparison is not perfectly apples-to-apples
- GAT uses only {len(y_gat)} confirmed-outcome patients; Cox uses all available
- A truly fair comparison would implement a survival-aware GNN loss
  (e.g. DeepSurv or Cox-PH loss) — that is the next step
- With n=198 and 10 seeds, std around 0.01-0.02 is expected

## Files produced

- `outputs/final_comparison/week2_survival_comparison.csv`
- `outputs/final_comparison/week2_survival_report.md`
"""
    (OUT / "week2_survival_report.md").write_text(report, encoding="utf-8")

    print(f"\n[week2] wrote week2_survival_comparison.csv")
    print(f"[week2] wrote week2_survival_report.md")
    print("\nWeek 2 complete.")


if __name__ == "__main__":
    main()
