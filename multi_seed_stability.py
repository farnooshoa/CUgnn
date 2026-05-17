"""Multi-seed stability test — properly reseeds every component.

Because RANDOM_SEED is imported as a value (not a reference) in evaluate.py
and train.py, we cannot patch it at runtime. Instead this script re-implements
the survival CV loop directly, seeding torch / numpy / sklearn at the start
of every seed iteration.

Usage:
    python multi_seed_stability.py          # survival only (~30 min)
    python multi_seed_stability.py --all    # all three tasks (~90 min)

Outputs:
    outputs/final_comparison/stability_results.csv
    outputs/final_comparison/stability_report.md
"""
from __future__ import annotations
import sys
import argparse
import random
from pathlib import Path
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.metrics import roc_auc_score, balanced_accuracy_score
from torch_geometric.loader import DataLoader

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from src.preprocessing import load_lihc_dataset
from src.graph_building import build_functional_graph
from src.gnn_models.dataset import build_graph_dataset
from src.gnn_models.models import GATGraphClassifier
from src.gnn_models.train import TrainConfig, _class_weights

OUT = ROOT / "outputs" / "final_comparison"
OUT.mkdir(parents=True, exist_ok=True)

SEEDS = list(range(1, 11))   # 1 through 10


def set_all_seeds(seed: int):
    """Reseed everything that matters."""
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def make_gat(in_dim: int, cfg: TrainConfig) -> GATGraphClassifier:
    return GATGraphClassifier(
        in_dim, cfg.hidden, n_classes=2,
        n_heads=4, n_layers=cfg.n_layers, dropout=cfg.dropout,
    ).to(cfg.device)


def train_and_eval(train_list, val_list, in_dim, cfg):
    """Train one fold and return (roc_auc, balanced_accuracy)."""
    y_train = np.array([int(d.y.item()) for d in train_list])
    class_w  = _class_weights(y_train, cfg.device)
    model    = make_gat(in_dim, cfg)
    opt      = torch.optim.Adam(model.parameters(),
                                lr=cfg.lr, weight_decay=cfg.weight_decay)
    loss_fn  = nn.CrossEntropyLoss(weight=class_w)

    tr_dl = DataLoader(train_list, batch_size=cfg.batch_size, shuffle=True)
    va_dl = DataLoader(val_list,   batch_size=cfg.batch_size, shuffle=False)

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
                prob = torch.softmax(
                    model(batch.x, batch.edge_index, batch.batch,
                          edge_weight=getattr(batch, "edge_weight", None)),
                    dim=-1)[:, 1].cpu().numpy()
                ys.extend(batch.y.cpu().numpy())
                ps.extend(prob)
        if len(set(ys)) > 1:
            v = roc_auc_score(ys, ps)
            if v > best_auc:
                best_auc = v
                best_state = {k: v2.clone() for k, v2 in model.state_dict().items()}

    if best_state:
        model.load_state_dict(best_state)

    model.eval()
    ys, ps, preds = [], [], []
    with torch.no_grad():
        for batch in va_dl:
            batch = batch.to(cfg.device)
            logits = model(batch.x, batch.edge_index, batch.batch,
                           edge_weight=getattr(batch, "edge_weight", None))
            prob = torch.softmax(logits, dim=-1)[:, 1].cpu().numpy()
            pr   = logits.argmax(dim=-1).cpu().numpy()
            ys.extend(batch.y.cpu().numpy())
            ps.extend(prob)
            preds.extend(pr)

    ys = np.array(ys); ps = np.array(ps); preds = np.array(preds)
    auc = roc_auc_score(ys, ps) if len(set(ys)) > 1 else float("nan")
    bal = balanced_accuracy_score(ys, preds)
    return auc, bal


def run_task_multiseed(task_name, ds, graph, y, groups,
                       sample_indices, seeds, epochs):
    rows = []
    for seed in seeds:
        print(f"  seed {seed:2d}/{seeds[-1]} ...", end=" ", flush=True)
        set_all_seeds(seed)

        # build bundle with this seed's train-fold z-score
        min_class = int(min((y == 0).sum(), (y == 1).sum()))
        n_splits  = max(2, min(5, min_class))
        cv = StratifiedGroupKFold(n_splits=n_splits,
                                  shuffle=True, random_state=seed)

        cfg = TrainConfig(model="gat", epochs=epochs)

        fold_aucs, fold_bals = [], []
        for fold, (tr_idx, va_idx) in enumerate(
                cv.split(np.zeros(len(y)), y, groups)):

            global_train_mask = np.zeros(ds.n_samples, dtype=bool)
            global_train_mask[sample_indices[tr_idx]] = True
            bundle = build_graph_dataset(ds, graph,
                                          zscore_train_mask=global_train_mask)

            # overwrite labels with task-specific y
            active_data = []
            for pos, orig_i in enumerate(sample_indices):
                d = bundle.data_list[orig_i]
                d.y = torch.tensor([int(y[pos])], dtype=torch.long)
                active_data.append(d)

            set_all_seeds(seed + fold * 100)   # fold-level re-seed
            auc, bal = train_and_eval(
                [active_data[i] for i in tr_idx],
                [active_data[i] for i in va_idx],
                bundle.in_dim, cfg,
            )
            fold_aucs.append(auc)
            fold_bals.append(bal)

        mean_auc = float(np.nanmean(fold_aucs))
        mean_bal = float(np.nanmean(fold_bals))
        print(f"AUC={mean_auc:.4f}")
        rows.append({"task": task_name, "seed": seed,
                     "roc_auc": mean_auc, "balanced_accuracy": mean_bal})

    return pd.DataFrame(rows)


# ── label helpers ──────────────────────────────────────────────────────────
def label_survival(row):
    os  = row.get("overall_survival_days")
    vit = row.get("vital_status")
    if pd.isna(os) or pd.isna(vit): return None
    if vit == "Dead":  return 1 if os <= 1095 else 0
    if vit == "Alive": return 0 if os > 1095  else None
    return None

def classify_stage(s):
    if not isinstance(s, str): return None
    s = s.upper().replace("STAGE ", "").strip()
    if s in {"I","II","IIA","IIB"}: return 0
    if s.startswith("III") or s.startswith("IV"): return 1
    return None


# ── main ───────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--all", action="store_true")
    args = ap.parse_args()

    print("[stability] loading dataset ...")
    ds    = load_lihc_dataset(require_real=True)
    genes = ds.expression.index.tolist()
    graph = build_functional_graph(genes, ds.copper_genes)
    md    = ds.metadata.loc[ds.expression.columns].copy()
    md["overall_survival_days"] = pd.to_numeric(
        md["overall_survival_days"], errors="coerce")

    all_results = []

    # ── survival ───────────────────────────────────────────────────────────
    print("\n[stability] === 3-year Overall Survival ===")
    md["y_3y"] = md.apply(label_survival, axis=1)
    keep = (md["sample_type"] == "Tumor") & md["y_3y"].notna()
    idx  = np.where(keep.to_numpy())[0]
    y    = md.loc[keep, "y_3y"].astype(int).to_numpy()
    g    = md.loc[keep, "case_submitter_id"].to_numpy()
    print(f"  n={len(y)}  dead={int((y==1).sum())}  alive={int((y==0).sum())}")
    all_results.append(
        run_task_multiseed("survival_3yr", ds, graph, y, g, idx, SEEDS, 80))

    if args.all:
        # ── tumor vs normal ────────────────────────────────────────────────
        print("\n[stability] === Tumor vs Normal ===")
        y_tn = (md["sample_type"].str.lower().eq("tumor").astype(int).to_numpy())
        g_tn = md["case_submitter_id"].to_numpy()
        all_results.append(run_task_multiseed(
            "tumor_vs_normal", ds, graph, y_tn, g_tn,
            np.arange(len(y_tn)), SEEDS, 60))

        # ── stage ──────────────────────────────────────────────────────────
        print("\n[stability] === Stage I/II vs III/IV ===")
        sb   = md["stage"].map(classify_stage)
        keep = (md["sample_type"] == "Tumor") & sb.notna()
        idx  = np.where(keep.to_numpy())[0]
        y_st = sb[keep].astype(int).to_numpy()
        g_st = md.loc[keep, "case_submitter_id"].to_numpy()
        all_results.append(run_task_multiseed(
            "stage_early_vs_late", ds, graph, y_st, g_st, idx, SEEDS, 80))

    # ── summary ────────────────────────────────────────────────────────────
    results = pd.concat(all_results, ignore_index=True)
    results.to_csv(OUT / "stability_results.csv", index=False)

    print("\n" + "="*55)
    print("STABILITY SUMMARY")
    print("="*55)

    report_lines = [
        "# Multi-Seed Stability Report",
        "",
        f"- torch           : {torch.__version__}",
        f"- seeds           : {SEEDS}",
        f"- model           : GAT, 5-fold StratifiedGroupKFold",
        "",
        "## Results",
        "",
        "| task | AUC mean | AUC std | AUC min | AUC max |",
        "|---|---:|---:|---:|---:|",
    ]

    for task, grp in results.groupby("task"):
        m = grp["roc_auc"].mean()
        s = grp["roc_auc"].std()
        lo = grp["roc_auc"].min()
        hi = grp["roc_auc"].max()
        print(f"\n{task}")
        print(f"  ROC-AUC : {m:.4f} ± {s:.4f}  [{lo:.4f} – {hi:.4f}]")
        verdict = ("STABLE" if s < 0.02
                   else "MODERATE" if s < 0.04
                   else "UNSTABLE")
        print(f"  verdict : {verdict}")
        report_lines.append(
            f"| {task} | **{m:.4f}** | {s:.4f} | {lo:.4f} | {hi:.4f} |")

    report_lines += [
        "",
        "## Verdict",
        "- std < 0.02 → **STABLE** — report mean ± std",
        "- std 0.02–0.04 → **MODERATE** — report mean ± std, note limitation",
        "- std > 0.04 → **UNSTABLE** — do not quote a single number",
        "",
        "## Per-seed detail",
        "",
    ]
    report_lines.append(results.to_string(index=False))

    report_path = OUT / "stability_report.md"
    report_path.write_text("\n".join(report_lines), encoding="utf-8")

    print(f"\n[stability] wrote stability_results.csv")
    print(f"[stability] wrote stability_report.md")


if __name__ == "__main__":
    main()
