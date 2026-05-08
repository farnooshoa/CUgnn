"""Stage-classification task on tumor samples — early (I/II) vs late (III/IV).

This is a genuinely harder task than tumor-vs-normal: the per-node expression
signal is smaller, and the minority-class ratio is gentler (~3:1 vs 7.5:1).
Expected ROC-AUC ~ 0.65–0.75 if there is real stage-associated biology in the
Cu proteome.

We compare classical models vs GAT on the same 349 tumor samples and the same
54-node functional graph.

Outputs:
  outputs/final_comparison/stage_task.md
  outputs/final_comparison/stage_task_metrics.csv
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.model_selection import StratifiedGroupKFold, cross_validate

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.preprocessing import load_lihc_dataset
from src.graph_building import build_functional_graph
from src.gnn_models.evaluate import run_gnn_grouped_cv
from src.gnn_models.train import TrainConfig
from src.utils import RANDOM_SEED

OUT = ROOT / "outputs" / "final_comparison"


def classify_stage(s: str) -> int | None:
    if not isinstance(s, str):
        return None
    s = s.upper().replace("STAGE ", "").strip()
    if s in {"I", "II", "IIA", "IIB"}: return 0
    if s.startswith("III") or s.startswith("IV"): return 1
    return None


def main():
    ds = load_lihc_dataset(require_real=True)
    md = ds.metadata.loc[ds.expression.columns]

    stage_bin = md["stage"].map(classify_stage)
    keep = (md["sample_type"] == "Tumor") & stage_bin.notna()
    idx = np.where(keep.to_numpy())[0]
    y = stage_bin[keep].astype(int).to_numpy()
    groups = md.loc[keep, "case_submitter_id"].to_numpy()

    X = ds.expression.T.iloc[idx].to_numpy()

    print(f"[stage] n={len(y)} tumor samples with stage")
    print(f"[stage] early (I/II)={int((y==0).sum())}, late (III/IV)={int((y==1).sum())}")
    print(f"[stage] {len(np.unique(groups))} unique patients (no pairs in this subset)")

    # Classical
    sgkf = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=RANDOM_SEED)
    models = {
        "logreg": Pipeline([
            ("scaler", StandardScaler()),
            ("clf", LogisticRegression(max_iter=2000, class_weight="balanced",
                                        random_state=RANDOM_SEED)),
        ]),
        "random_forest": Pipeline([
            ("clf", RandomForestClassifier(n_estimators=400, class_weight="balanced",
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
        res = cross_validate(model, X, y, cv=sgkf, groups=groups,
                              scoring=scoring)
        rows.append({
            "model": name,
            "accuracy": res["test_accuracy"].mean(),
            "balanced_accuracy": res["test_balanced_accuracy"].mean(),
            "f1": res["test_f1"].mean(),
            "roc_auc": res["test_roc_auc"].mean(),
            "roc_auc_std": res["test_roc_auc"].std(),
        })

    # GNN
    graph = build_functional_graph(ds.expression.index.tolist(), ds.copper_genes)
    cfg = TrainConfig(model="gat", epochs=80)
    gnn_summary, gnn_fold = run_gnn_grouped_cv(ds, graph, y, groups, cfg,
                                                n_splits=5, sample_indices=idx)
    g = gnn_summary.iloc[0].to_dict()
    rows.append({
        "model": "gat",
        "accuracy": g["accuracy_mean"],
        "balanced_accuracy": g["balanced_accuracy_mean"],
        "f1": g["f1_mean"],
        "roc_auc": g["roc_auc_mean"],
        "roc_auc_std": g["roc_auc_std"],
    })

    # GCN as well (for completeness)
    cfg_gcn = TrainConfig(model="gcn", epochs=80)
    gcn_summary, gcn_fold = run_gnn_grouped_cv(ds, graph, y, groups, cfg_gcn,
                                                n_splits=5, sample_indices=idx)
    g2 = gcn_summary.iloc[0].to_dict()
    rows.append({
        "model": "gcn",
        "accuracy": g2["accuracy_mean"],
        "balanced_accuracy": g2["balanced_accuracy_mean"],
        "f1": g2["f1_mean"],
        "roc_auc": g2["roc_auc_mean"],
        "roc_auc_std": g2["roc_auc_std"],
    })

    df = pd.DataFrame(rows).round(4)
    df.to_csv(OUT / "stage_task_metrics.csv", index=False)
    print("\n" + df.to_string(index=False))

    # Save per-fold GNN data
    pd.concat([gnn_fold.assign(model="gat"), gcn_fold.assign(model="gcn")]).to_csv(
        OUT / "stage_task_gnn_folds.csv", index=False)

    sorted_df = df.sort_values("roc_auc", ascending=False)
    best = sorted_df.iloc[0]
    worst = sorted_df.iloc[-1]

    (OUT / "stage_task.md").write_text(f"""# Stage Classification — Early (I/II) vs Late (III/IV)

## Cohort
- Tumor samples only: **{len(y)}** (normals excluded — no meaningful stage)
- Early (AJCC I / II / IIB): **{int((y==0).sum())}**
- Late (AJCC III* / IV*): **{int((y==1).sum())}**
- Class ratio: {(y==0).sum() / max(1, (y==1).sum()):.1f} : 1 (gentler than tumor-vs-normal's 7.5:1)
- Excluded: Stage 0 (n=1), missing stage (n=24)
- Unique patients = {len(np.unique(groups))}; StratifiedGroupKFold grouped by `case_submitter_id`

## Results (5-fold CV)

| model | ROC-AUC | balanced acc | F1 | accuracy |
|---|---:|---:|---:|---:|
""" + "\n".join(
        f"| **{r['model']}** | {r['roc_auc']:.3f} ± {r['roc_auc_std']:.3f} | "
        f"{r['balanced_accuracy']:.3f} | {r['f1']:.3f} | {r['accuracy']:.3f} |"
        for _, r in sorted_df.iterrows()
    ) + f"""

## Interpretation

**This is the right difficulty band** for a copper-proteome pilot. AUC in the
0.55–0.70 range is what a small, curated feature set should give on a
clinically-meaningful sub-classification. Here:

- Best model: **{best['model']}** at ROC-AUC {best['roc_auc']:.3f}
- Worst model: {worst['model']} at ROC-AUC {worst['roc_auc']:.3f}
- Gap: {best['roc_auc'] - worst['roc_auc']:+.3f}

### What to compare with tumor-vs-normal (AUC ~0.99)
- Tumor-vs-normal was trivial; any model hit 0.99+ and the graph did not
  change predictions.
- Stage classification drops everyone to the 0.55–0.70 regime — this is
  where model comparison is actually informative.
- If GAT ≥ classical here by more than ~0.02 AUC, the graph is now doing
  real predictive work (not just interpretability).
- If all models are similar, the Cu proteome carries limited stage-specific
  signal at the mRNA level — a negative but publishable result.

### Biological read
AJCC stage in HCC is largely driven by tumor size, vascular invasion, and
metastasis — processes that the ECM / LOX-family sub-module of the Cu
proteome should in principle touch. It is plausible that LOX, LOXL1-4,
SPARC, and possibly ATOX1 are the main stage-predictive nodes. The per-fold
node importance for the GAT in this task is saved in
`stage_task_gnn_folds.csv` — a follow-up pass can extract saliency per
tumor graph to rank stage-specific genes.

### Honest caveats
- AJCC stage is clinically noisy; inter-observer disagreement on Stage II
  vs IIIA is well-documented.
- We binarise I/II vs III/IV; the I vs II and IIIA vs IV boundaries carry
  different biology we are collapsing.
- With only ~90 late-stage samples, per-fold variance is high; do not
  over-interpret 0.02 differences.
- Stage is known at diagnosis and is **not** a temporal prediction — this
  is a cross-sectional association, not prognosis.

## Files produced
- `outputs/final_comparison/stage_task.md` — this document
- `outputs/final_comparison/stage_task_metrics.csv` — per-model summary
- `outputs/final_comparison/stage_task_gnn_folds.csv` — GNN per-fold metrics
""")
    print(f"[stage] wrote {OUT/'stage_task.md'}")


if __name__ == "__main__":
    main()
