"""Leakage audit — re-run classical + GNN with:
  1. sklearn Pipeline (no global StandardScaler leakage)
  2. StratifiedGroupKFold grouped by case_submitter_id (patient-level splits)
  3. permutation-label sanity check (labels shuffled, expect AUC ~ 0.5)

Outputs:
  outputs/final_comparison/leakage_audit.md
  outputs/final_comparison/leakage_audit_metrics.csv
"""
from __future__ import annotations
import sys
import json
import time
from pathlib import Path
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.model_selection import (
    StratifiedKFold, StratifiedGroupKFold, cross_validate,
)
from sklearn.metrics import roc_auc_score, balanced_accuracy_score, f1_score
from torch_geometric.loader import DataLoader

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.preprocessing import load_lihc_dataset
from src.graph_building import build_functional_graph
from src.gnn_models.dataset import build_graph_dataset
from src.gnn_models.train import (
    TrainConfig, _class_weights, _train_one_fold, _make_model,
)
from src.utils import RANDOM_SEED

OUT = ROOT / "outputs" / "final_comparison"
OUT.mkdir(parents=True, exist_ok=True)


# ------------------------------------------------------------------ #
# Classical — Pipeline + GroupKFold                                   #
# ------------------------------------------------------------------ #
def classical_cv(X, y, groups, cv, tag: str) -> pd.DataFrame:
    models = {
        "logreg": Pipeline([
            ("scaler", StandardScaler()),
            ("clf", LogisticRegression(max_iter=2000, class_weight="balanced",
                                        random_state=RANDOM_SEED)),
        ]),
        "random_forest": Pipeline([
            ("clf", RandomForestClassifier(n_estimators=300, class_weight="balanced",
                                           random_state=RANDOM_SEED, n_jobs=-1)),
        ]),
        "svm_rbf": Pipeline([
            ("scaler", StandardScaler()),
            ("clf", SVC(kernel="rbf", class_weight="balanced", probability=True,
                        random_state=RANDOM_SEED)),
        ]),
    }
    scoring = {"accuracy": "accuracy", "balanced_accuracy": "balanced_accuracy",
               "f1": "f1", "roc_auc": "roc_auc"}
    rows = []
    for name, model in models.items():
        res = cross_validate(model, X, y, cv=cv, groups=groups,
                              scoring=scoring, n_jobs=1)
        rows.append({
            "model": name, "setting": tag,
            "accuracy": res["test_accuracy"].mean(),
            "balanced_accuracy": res["test_balanced_accuracy"].mean(),
            "f1": res["test_f1"].mean(),
            "roc_auc": res["test_roc_auc"].mean(),
            "roc_auc_std": res["test_roc_auc"].std(),
        })
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ #
# GNN — per-fold bundle + GroupKFold                                  #
# ------------------------------------------------------------------ #
def gnn_cv(ds, graph, y, groups, cv, tag: str, model_name: str = "gat") -> pd.DataFrame:
    cfg = TrainConfig(model=model_name, epochs=60)
    all_fold = []
    rng_state = torch.get_rng_state()
    torch.manual_seed(RANDOM_SEED)

    fold_iter = cv.split(np.zeros(len(y)), y, groups) if groups is not None \
                else cv.split(np.zeros(len(y)), y)

    for fold, (tr_idx, va_idx) in enumerate(fold_iter):
        train_mask = np.zeros(len(y), dtype=bool)
        train_mask[tr_idx] = True
        bundle = build_graph_dataset(ds, graph, zscore_train_mask=train_mask)
        # overwrite labels with the (possibly-permuted) y so permutation
        # tests actually shuffle what the model trains on.
        for i, data in enumerate(bundle.data_list):
            data.y = torch.tensor([int(y[i])], dtype=torch.long)
        class_w = _class_weights(y, cfg.device)

        train_dl = DataLoader([bundle.data_list[i] for i in tr_idx],
                               batch_size=cfg.batch_size, shuffle=True)
        val_dl = DataLoader([bundle.data_list[i] for i in va_idx],
                             batch_size=cfg.batch_size, shuffle=False)
        model = _make_model(cfg, bundle.in_dim).to(cfg.device)
        model, _ = _train_one_fold(model, train_dl, val_dl, cfg, class_w)

        model.eval()
        ys, ps, preds = [], [], []
        with torch.no_grad():
            for batch in val_dl:
                batch = batch.to(cfg.device)
                logits = model(batch.x, batch.edge_index, batch.batch,
                               edge_weight=getattr(batch, "edge_weight", None))
                prob = torch.softmax(logits, dim=-1)[:, 1].cpu().numpy()
                pr = logits.argmax(dim=-1).cpu().numpy()
                ys.extend(batch.y.cpu().numpy().tolist())
                ps.extend(prob.tolist())
                preds.extend(pr.tolist())
        ys = np.array(ys); ps = np.array(ps); preds = np.array(preds)
        all_fold.append({
            "model": model_name, "setting": tag, "fold": fold,
            "balanced_accuracy": balanced_accuracy_score(ys, preds),
            "f1": f1_score(ys, preds, zero_division=0),
            "roc_auc": roc_auc_score(ys, ps) if len(set(ys)) > 1 else float("nan"),
            "n_val": len(ys),
        })
    torch.set_rng_state(rng_state)
    df = pd.DataFrame(all_fold)
    agg = (df.drop(columns=["fold", "n_val"]).groupby(["model", "setting"])
             .agg(["mean", "std"]).reset_index())
    agg.columns = ["_".join(c).strip("_") for c in agg.columns]
    return agg, df


# ------------------------------------------------------------------ #
# Main audit                                                          #
# ------------------------------------------------------------------ #
def main():
    ds = load_lihc_dataset(require_real=True)
    print(f"[audit] loaded {ds.n_samples} samples, {ds.n_genes} Cu genes")

    y = (ds.metadata.loc[ds.expression.columns, "sample_type"]
         .str.lower().eq("tumor").astype(int).to_numpy())
    groups = ds.metadata.loc[ds.expression.columns, "case_submitter_id"].to_numpy()
    X = ds.expression.T.to_numpy()  # samples x genes
    graph = build_functional_graph(ds.expression.index.tolist(), ds.copper_genes)

    print(f"[audit] groups = {len(np.unique(groups))} unique patients")
    print(f"[audit] class balance = {int(y.sum())} tumor / {int((1-y).sum())} normal")

    results = []
    gnn_per_fold_records = []

    # --- CLASSICAL ------------------------------------------------------- #
    print("\n=== Classical — legacy (StratifiedKFold, global scaler) ===")
    legacy_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=RANDOM_SEED)
    from sklearn.preprocessing import StandardScaler as _SS
    X_legacy = _SS().fit_transform(X)
    from src.baseline_models.run_baselines import run_classical_baselines as _legacy
    # the legacy function takes the LIHCDataset directly; we re-implement minimally
    legacy_models = {
        "logreg": LogisticRegression(max_iter=2000, class_weight="balanced",
                                     random_state=RANDOM_SEED),
        "random_forest": RandomForestClassifier(n_estimators=300, class_weight="balanced",
                                                random_state=RANDOM_SEED, n_jobs=-1),
        "svm_rbf": SVC(kernel="rbf", class_weight="balanced", probability=True,
                        random_state=RANDOM_SEED),
    }
    scoring = {"accuracy": "accuracy", "balanced_accuracy": "balanced_accuracy",
               "f1": "f1", "roc_auc": "roc_auc"}
    for name, m in legacy_models.items():
        res = cross_validate(m, X_legacy, y, cv=legacy_cv, scoring=scoring)
        results.append({
            "model": name, "setting": "classical_legacy",
            "accuracy": res["test_accuracy"].mean(),
            "balanced_accuracy": res["test_balanced_accuracy"].mean(),
            "f1": res["test_f1"].mean(),
            "roc_auc": res["test_roc_auc"].mean(),
            "roc_auc_std": res["test_roc_auc"].std(),
        })

    print("\n=== Classical — Pipeline + StratifiedGroupKFold ===")
    sgkf = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=RANDOM_SEED)
    df_fixed = classical_cv(X, y, groups, sgkf, "classical_pipeline_groupkfold")
    results.extend(df_fixed.to_dict("records"))
    print(df_fixed.round(4).to_string(index=False))

    print("\n=== Classical — PERMUTATION (shuffled labels) ===")
    rng = np.random.default_rng(RANDOM_SEED)
    perm_rows = []
    for seed in range(5):
        y_perm = rng.permutation(y)
        df_perm = classical_cv(X, y_perm, groups, sgkf,
                                f"classical_permutation_seed{seed}")
        perm_rows.append(df_perm)
    perm_df = pd.concat(perm_rows)
    perm_summary = perm_df.groupby("model")[["roc_auc", "balanced_accuracy"]].agg(["mean", "std"]).round(4)
    print(perm_summary)

    # --- GNN ------------------------------------------------------------- #
    print("\n=== GNN — legacy (StratifiedKFold, global z-score) ===")
    agg_legacy, df_legacy = gnn_cv(ds, graph, y, None,
                                    StratifiedKFold(n_splits=5, shuffle=True,
                                                     random_state=RANDOM_SEED),
                                    "gnn_legacy", model_name="gat")
    gnn_per_fold_records.append(df_legacy)
    print(agg_legacy)

    print("\n=== GNN — per-fold z-score + StratifiedGroupKFold ===")
    agg_fixed, df_fixed = gnn_cv(ds, graph, y, groups, sgkf,
                                  "gnn_groupkfold", model_name="gat")
    gnn_per_fold_records.append(df_fixed)
    print(agg_fixed)

    print("\n=== GNN — PERMUTATION (shuffled labels, 3 seeds to save time) ===")
    gnn_perm_rows = []
    for seed in range(3):
        y_perm = np.random.default_rng(RANDOM_SEED + seed).permutation(y)
        agg_perm, df_perm = gnn_cv(ds, graph, y_perm, groups, sgkf,
                                    f"gnn_permutation_seed{seed}", model_name="gat")
        gnn_per_fold_records.append(df_perm)
        gnn_perm_rows.append(agg_perm)
    gnn_perm = pd.concat(gnn_perm_rows)
    print(gnn_perm)

    # Save tidy tables
    classical_df = pd.DataFrame(results)
    classical_df.to_csv(OUT / "leakage_audit_metrics.csv", index=False)

    all_gnn = pd.concat(gnn_per_fold_records, ignore_index=True)
    all_gnn.to_csv(OUT / "leakage_audit_gnn_folds.csv", index=False)

    # Markdown report
    real_rows = classical_df[classical_df["setting"].isin(
        ["classical_legacy", "classical_pipeline_groupkfold"])]
    perm_rows_df = classical_df[classical_df["setting"].str.startswith("classical_permutation")]
    perm_by_model = perm_rows_df.groupby("model")[["roc_auc", "balanced_accuracy"]].mean().round(3)

    gnn_real = all_gnn[all_gnn["setting"].isin(["gnn_legacy", "gnn_groupkfold"])]
    gnn_real_agg = gnn_real.groupby("setting")[["roc_auc", "balanced_accuracy", "f1"]].mean().round(3)
    gnn_perm_fold = all_gnn[all_gnn["setting"].str.startswith("gnn_permutation")]
    gnn_perm_agg = gnn_perm_fold[["roc_auc", "balanced_accuracy"]].mean().round(3)

    (OUT / "leakage_audit.md").write_text(f"""# Leakage Audit — TCGA-LIHC Copper Proteome

## Fixes applied

1. **sklearn Pipeline** — StandardScaler is fit *inside* each CV fold instead of globally, so no test-fold statistics leak into training-fold scaling.
2. **StratifiedGroupKFold (groups = case_submitter_id)** — all samples from the same patient stay in the same fold. 50 matched normals each come from a different patient that also contributes a tumor sample; without grouping, the previous CV could place a patient's tumor in train and their normal in test.
3. **Per-fold z-score for GNN node features** — `build_graph_dataset(..., zscore_train_mask=...)` now computes per-gene mean/std from training samples only.

## Permutation-label sanity check
Labels shuffled, CV repeated. A model with no leakage should drop to ROC-AUC ≈ 0.5.

## Classical models

| model | setting | ROC-AUC | balanced acc |
|---|---|---:|---:|
""")

    lines = []
    for _, r in classical_df[
        classical_df["setting"].isin(["classical_legacy", "classical_pipeline_groupkfold"])
    ].iterrows():
        lines.append(f"| {r['model']} | {r['setting']} | {r['roc_auc']:.3f} | {r['balanced_accuracy']:.3f} |")
    lines.append(f"| **permutation (shuffled labels, mean over 5 seeds)** | | | |")
    for m, row in perm_by_model.iterrows():
        lines.append(f"| {m} | classical_permutation | {row['roc_auc']:.3f} | {row['balanced_accuracy']:.3f} |")

    gnn_block = ["", "## GNN (GAT)", "",
                 "| setting | ROC-AUC | balanced acc | F1 |",
                 "|---|---:|---:|---:|"]
    for setting, row in gnn_real_agg.iterrows():
        gnn_block.append(f"| {setting} | {row['roc_auc']:.3f} | {row['balanced_accuracy']:.3f} | {row['f1']:.3f} |")
    gnn_block.append(f"| **gnn_permutation (mean over 3 seeds)** | {gnn_perm_agg['roc_auc']:.3f} | {gnn_perm_agg['balanced_accuracy']:.3f} | - |")

    interp = f"""

## Interpretation

- If the *pipeline_groupkfold* ROC-AUC is close to the *legacy* ROC-AUC, leakage was small and the task is intrinsically easy.
- If the permutation-label ROC-AUC is near 0.5, the model is genuinely using the labels, not a data artefact.
- Large gaps (legacy − groupkfold > 0.05) would indicate meaningful patient-level leakage.

**Read the Classical table row-by-row and the GNN block before drawing conclusions.**
"""
    with open(OUT / "leakage_audit.md", "a") as f:
        f.write("\n".join(lines) + "\n")
        f.write("\n".join(gnn_block) + "\n")
        f.write(interp)

    print(f"\n[audit] wrote {OUT/'leakage_audit.md'}")
    print(f"[audit] wrote {OUT/'leakage_audit_metrics.csv'}")
    print(f"[audit] wrote {OUT/'leakage_audit_gnn_folds.csv'}")


if __name__ == "__main__":
    main()
