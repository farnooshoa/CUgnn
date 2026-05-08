"""Reusable GNN cross-validation with patient-grouped splits and per-fold z-score.

Import this for any downstream experiment (graph ablation, stage classification,
alternative edge sources, multi-omics features). It re-uses the already-fixed
training loop from ``train.py`` but exposes a clean interface:

    from src.gnn_models.evaluate import run_gnn_grouped_cv
    summary_df, per_fold_df = run_gnn_grouped_cv(ds, graph, y, groups, cfg)

All of this code reflects the post-audit methodology from leakage_audit.py:
  * StratifiedGroupKFold on ``groups``
  * per-fold z-score fit on training samples only (zscore_train_mask)
  * labels explicitly written into each Data.y so permutation / external labels
    correctly flow into training.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
import pandas as pd
import torch
from sklearn.model_selection import StratifiedGroupKFold, StratifiedKFold
from sklearn.metrics import (
    accuracy_score, balanced_accuracy_score, f1_score, roc_auc_score,
)
from torch_geometric.loader import DataLoader

from src.gnn_models.dataset import build_graph_dataset
from src.gnn_models.train import (
    TrainConfig, _class_weights, _train_one_fold, _make_model,
)
from src.utils import RANDOM_SEED


def run_gnn_grouped_cv(ds, graph, y: np.ndarray, groups: np.ndarray | None,
                       cfg: TrainConfig, n_splits: int = 5,
                       sample_indices: np.ndarray | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run stratified-group CV for a GAT/GCN graph-classification task.

    Parameters
    ----------
    ds : LIHCDataset
    graph : networkx.Graph over the 54 Cu genes
    y : array-like, one label per sample (same ordering as ds.expression.columns)
    groups : array-like of patient ids (same ordering) or None for plain stratified
    cfg : TrainConfig
    sample_indices : optional positional indices into ds.expression.columns
        for when the task uses a subset (e.g. tumor-only for stage classification)

    Returns
    -------
    summary_df : one-row-per-model aggregation (mean, std across folds)
    per_fold_df : per-fold metrics
    """
    if sample_indices is None:
        sample_indices = np.arange(len(y))

    min_class = int(min((y == 0).sum(), (y == 1).sum()))
    n_splits = max(2, min(n_splits, min_class))

    if groups is not None:
        cv = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=RANDOM_SEED)
        fold_iter = cv.split(np.zeros(len(y)), y, groups)
    else:
        cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=RANDOM_SEED)
        fold_iter = cv.split(np.zeros(len(y)), y)

    torch.manual_seed(RANDOM_SEED)
    np.random.seed(RANDOM_SEED)

    per_fold = []
    for fold, (tr_idx, va_idx) in enumerate(fold_iter):
        # Build a full bundle but restrict labels / active samples
        global_train_mask = np.zeros(ds.n_samples, dtype=bool)
        active_train = sample_indices[tr_idx]
        global_train_mask[active_train] = True
        bundle = build_graph_dataset(ds, graph, zscore_train_mask=global_train_mask)

        # Pick out the active subset of graphs and overwrite labels
        active_all = sample_indices
        sub_data = []
        for pos, orig_i in enumerate(active_all):
            d = bundle.data_list[orig_i]
            d.y = torch.tensor([int(y[pos])], dtype=torch.long)
            sub_data.append(d)
        train_list = [sub_data[i] for i in tr_idx]
        val_list = [sub_data[i] for i in va_idx]

        class_w = _class_weights(y, cfg.device)
        train_dl = DataLoader(train_list, batch_size=cfg.batch_size, shuffle=True)
        val_dl = DataLoader(val_list, batch_size=cfg.batch_size, shuffle=False)
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
        per_fold.append({
            "model": cfg.model, "fold": fold,
            "accuracy": accuracy_score(ys, preds),
            "balanced_accuracy": balanced_accuracy_score(ys, preds),
            "f1": f1_score(ys, preds, zero_division=0),
            "roc_auc": roc_auc_score(ys, ps) if len(set(ys)) > 1 else float("nan"),
            "n_val": len(ys),
        })

    per_fold_df = pd.DataFrame(per_fold)
    summary = (per_fold_df.drop(columns=["fold", "n_val"])
                 .groupby("model").agg(["mean", "std"]).round(4).reset_index())
    summary.columns = ["_".join(c).strip("_") for c in summary.columns]
    return summary, per_fold_df
