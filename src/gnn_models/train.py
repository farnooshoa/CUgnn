"""Train GCN / GAT on graph-classification (tumor vs normal) with CV."""
from __future__ import annotations
from pathlib import Path
from dataclasses import dataclass
import json
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from torch_geometric.loader import DataLoader
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import (
    accuracy_score, balanced_accuracy_score, f1_score, roc_auc_score,
)

from src.preprocessing import LIHCDataset, load_lihc_dataset
from src.graph_building import build_ppi_graph, build_functional_graph
from src.gnn_models.models import GCNGraphClassifier, GATGraphClassifier
from src.gnn_models.dataset import build_graph_dataset, GraphDataBundle
from src.utils import GNN_DIR, RANDOM_SEED


@dataclass
class TrainConfig:
    model: str = "gcn"                    # 'gcn' or 'gat'
    hidden: int = 64
    n_layers: int = 2
    dropout: float = 0.4
    epochs: int = 80
    lr: float = 5e-3
    weight_decay: float = 5e-4
    batch_size: int = 32
    n_splits: int = 5
    device: str = "cuda" if torch.cuda.is_available() else "cpu"


def _make_model(cfg: TrainConfig, in_dim: int) -> nn.Module:
    if cfg.model == "gcn":
        return GCNGraphClassifier(in_dim, cfg.hidden, n_classes=2,
                                  n_layers=cfg.n_layers, dropout=cfg.dropout)
    if cfg.model == "gat":
        return GATGraphClassifier(in_dim, cfg.hidden, n_classes=2,
                                  n_layers=cfg.n_layers, dropout=cfg.dropout)
    raise ValueError(cfg.model)


def _class_weights(y: np.ndarray, device: str) -> torch.Tensor:
    counts = np.bincount(y, minlength=2).astype(np.float32)
    counts = np.where(counts == 0, 1.0, counts)
    w = counts.sum() / (2 * counts)
    return torch.tensor(w, dtype=torch.float32, device=device)


def _train_one_fold(model, train_loader, val_loader, cfg, class_w):
    opt = torch.optim.Adam(model.parameters(), lr=cfg.lr, weight_decay=cfg.weight_decay)
    loss_fn = nn.CrossEntropyLoss(weight=class_w)
    best_val = -1.0
    best_state = None
    for ep in range(cfg.epochs):
        model.train()
        for batch in train_loader:
            batch = batch.to(cfg.device)
            opt.zero_grad()
            logits = model(batch.x, batch.edge_index, batch.batch,
                           edge_weight=getattr(batch, "edge_weight", None))
            loss = loss_fn(logits, batch.y)
            loss.backward()
            opt.step()
        # validation
        model.eval()
        ys, ps = [], []
        with torch.no_grad():
            for batch in val_loader:
                batch = batch.to(cfg.device)
                logits = model(batch.x, batch.edge_index, batch.batch,
                               edge_weight=getattr(batch, "edge_weight", None))
                prob = torch.softmax(logits, dim=-1)[:, 1].cpu().numpy()
                ys.extend(batch.y.cpu().numpy().tolist())
                ps.extend(prob.tolist())
        try:
            val = roc_auc_score(ys, ps) if len(set(ys)) > 1 else 0.0
        except Exception:
            val = 0.0
        if val > best_val:
            best_val = val
            best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}
    if best_state is not None:
        model.load_state_dict(best_state)
    return model, best_val


def cross_validate(bundle: GraphDataBundle, cfg: TrainConfig) -> pd.DataFrame:
    torch.manual_seed(RANDOM_SEED)
    np.random.seed(RANDOM_SEED)
    data_list = bundle.data_list
    y = np.array([int(d.y.item()) for d in data_list])

    min_class = int(min(y.sum(), (1 - y).sum()))
    n_splits = max(2, min(cfg.n_splits, min_class))
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=RANDOM_SEED)
    class_w = _class_weights(y, cfg.device)

    fold_rows = []
    for fold, (tr_idx, va_idx) in enumerate(skf.split(np.zeros(len(data_list)), y)):
        train_dl = DataLoader([data_list[i] for i in tr_idx],
                              batch_size=cfg.batch_size, shuffle=True)
        val_dl = DataLoader([data_list[i] for i in va_idx],
                            batch_size=cfg.batch_size, shuffle=False)
        model = _make_model(cfg, bundle.in_dim).to(cfg.device)
        model, best_val = _train_one_fold(model, train_dl, val_dl, cfg, class_w)
        # final fold metrics
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
        ys = np.array(ys); preds = np.array(preds); ps = np.array(ps)
        fold_rows.append({
            "model": cfg.model, "fold": fold,
            "accuracy": accuracy_score(ys, preds),
            "balanced_accuracy": balanced_accuracy_score(ys, preds),
            "f1": f1_score(ys, preds, zero_division=0),
            "roc_auc": roc_auc_score(ys, ps) if len(set(ys)) > 1 else float("nan"),
            "n_val": len(ys),
        })
    return pd.DataFrame(fold_rows)


def aggregate(df: pd.DataFrame) -> pd.DataFrame:
    return (df.drop(columns=["fold"])
              .groupby("model")
              .agg(["mean", "std"])
              .round(4))


def train_full_then_embed(bundle: GraphDataBundle, cfg: TrainConfig):
    """Train on everything, then return embeddings + labels for viz."""
    torch.manual_seed(RANDOM_SEED)
    y = np.array([int(d.y.item()) for d in bundle.data_list])
    class_w = _class_weights(y, cfg.device)
    model = _make_model(cfg, bundle.in_dim).to(cfg.device)
    opt = torch.optim.Adam(model.parameters(), lr=cfg.lr, weight_decay=cfg.weight_decay)
    loss_fn = nn.CrossEntropyLoss(weight=class_w)
    dl = DataLoader(bundle.data_list, batch_size=cfg.batch_size, shuffle=True)
    for _ in range(cfg.epochs):
        model.train()
        for batch in dl:
            batch = batch.to(cfg.device)
            opt.zero_grad()
            logits = model(batch.x, batch.edge_index, batch.batch,
                           edge_weight=getattr(batch, "edge_weight", None))
            loss = loss_fn(logits, batch.y)
            loss.backward()
            opt.step()

    model.eval()
    ordered_dl = DataLoader(bundle.data_list, batch_size=cfg.batch_size, shuffle=False)
    embs = []
    ys = []
    with torch.no_grad():
        for batch in ordered_dl:
            batch = batch.to(cfg.device)
            _, g_emb, _ = model(batch.x, batch.edge_index, batch.batch,
                                 edge_weight=getattr(batch, "edge_weight", None),
                                 return_embedding=True)
            embs.append(g_emb.cpu().numpy())
            ys.extend(batch.y.cpu().numpy().tolist())
    return model, np.concatenate(embs, axis=0), np.array(ys)


def embedding_scatter(emb: np.ndarray, labels: np.ndarray, out_path: Path,
                       title: str):
    try:
        import umap
        coords = umap.UMAP(n_components=2, random_state=RANDOM_SEED).fit_transform(emb)
        method = "UMAP"
    except Exception:
        from sklearn.manifold import TSNE
        perp = min(30, max(5, len(emb) // 5))
        coords = TSNE(n_components=2, random_state=RANDOM_SEED,
                       perplexity=perp, init="pca").fit_transform(emb)
        method = "t-SNE"
    fig, ax = plt.subplots(figsize=(6, 5))
    for cls, color, name in [(0, "#4682B4", "Normal"), (1, "#B22222", "Tumor")]:
        mask = labels == cls
        ax.scatter(coords[mask, 0], coords[mask, 1], s=25, alpha=0.75,
                   label=name, c=color)
    ax.set_title(f"{method} of GNN graph embeddings — {title}")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def node_importance_from_gradient(model, bundle: GraphDataBundle, cfg: TrainConfig,
                                   n_graphs: int = 20) -> pd.DataFrame:
    """Saliency-style node importance: |d logit_tumor / d x| pooled per gene."""
    model.eval()
    loader = DataLoader(bundle.data_list[:n_graphs], batch_size=1, shuffle=False)
    per_gene = np.zeros(len(bundle.gene_order), dtype=np.float64)
    count = 0
    for batch in loader:
        batch = batch.to(cfg.device)
        batch.x.requires_grad_(True)
        logits = model(batch.x, batch.edge_index, batch.batch,
                       edge_weight=getattr(batch, "edge_weight", None))
        tumor_logit = logits[:, 1].sum()
        tumor_logit.backward()
        grad = batch.x.grad.detach().cpu().numpy()
        per_gene += np.abs(grad).sum(axis=1)
        count += 1
    per_gene /= max(count, 1)
    df = (pd.DataFrame({"gene_symbol": bundle.gene_order, "importance": per_gene})
            .sort_values("importance", ascending=False))
    return df


def attention_summary(model: GATGraphClassifier, bundle: GraphDataBundle,
                       cfg: TrainConfig, n_graphs: int = 20) -> pd.DataFrame:
    model.eval()
    loader = DataLoader(bundle.data_list[:n_graphs], batch_size=1, shuffle=False)
    agg: dict[tuple[int, int], float] = {}
    with torch.no_grad():
        for batch in loader:
            batch = batch.to(cfg.device)
            _ = model(batch.x, batch.edge_index, batch.batch,
                      edge_weight=getattr(batch, "edge_weight", None),
                      return_attention=True)
            attn = model.get_last_attention()
            if not attn:
                continue
            ei, alpha = attn[-1]
            ei = ei.cpu().numpy()
            alpha = alpha.mean(dim=1).cpu().numpy() if alpha.ndim > 1 else alpha.cpu().numpy()
            for (s, t), a in zip(ei.T, alpha):
                key = (int(s), int(t))
                agg[key] = agg.get(key, 0.0) + float(a)
    rows = []
    n_nodes = len(bundle.gene_order)
    for (s, t), a in agg.items():
        if s < n_nodes and t < n_nodes:
            rows.append({"source": bundle.gene_order[s % n_nodes],
                         "target": bundle.gene_order[t % n_nodes],
                         "attention_sum": a})
    return (pd.DataFrame(rows)
              .groupby(["source", "target"], as_index=False)["attention_sum"].sum()
              .sort_values("attention_sum", ascending=False))


def run_all(ds: LIHCDataset | None = None, out_dir: Path = GNN_DIR) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    if ds is None:
        ds = load_lihc_dataset()
    genes = ds.expression.index.tolist()
    graph = build_functional_graph(genes, ds.copper_genes)
    bundle = build_graph_dataset(ds, graph)

    all_metrics = []
    best_auc = -np.inf
    best_model_name = "gcn"
    summaries = {}
    for model_name in ["gcn", "gat"]:
        cfg = TrainConfig(model=model_name)
        df = cross_validate(bundle, cfg)
        df["model"] = model_name
        all_metrics.append(df)
        summaries[model_name] = {
            "mean_roc_auc": float(df["roc_auc"].mean()),
            "mean_balanced_accuracy": float(df["balanced_accuracy"].mean()),
            "mean_f1": float(df["f1"].mean()),
            "mean_accuracy": float(df["accuracy"].mean()),
            "n_folds": len(df),
        }
        if df["roc_auc"].mean() > best_auc:
            best_auc = df["roc_auc"].mean()
            best_model_name = model_name

    metrics_df = pd.concat(all_metrics, ignore_index=True)
    metrics_df.to_csv(out_dir / "model_metrics.csv", index=False)
    (out_dir / "model_metrics_summary.json").write_text(
        json.dumps(summaries, indent=2))

    # Final model + embeddings + importance
    cfg_best = TrainConfig(model=best_model_name)
    model, emb, labels = train_full_then_embed(bundle, cfg_best)
    embedding_scatter(emb, labels, out_dir / "graph_embedding_umap.png",
                       title=f"best={best_model_name}")

    imp = node_importance_from_gradient(model, bundle, cfg_best, n_graphs=min(50, len(bundle.data_list)))
    imp.to_csv(out_dir / "node_importance.csv", index=False)

    attn_text = ""
    if best_model_name == "gat":
        attn_df = attention_summary(model, bundle, cfg_best,
                                     n_graphs=min(50, len(bundle.data_list)))
        attn_df.head(50).to_csv(out_dir / "top_attention_edges.csv", index=False)
        top_rows = attn_df.head(20).to_string(index=False)
        attn_text = f"\n\n## Top-20 attention edges (GAT)\n\n```\n{top_rows}\n```\n"

    top_imp = imp.head(10).to_string(index=False)
    (out_dir / "top_subgraph_or_attention_summary.md").write_text(
        f"# GNN interpretability summary\n\n"
        f"Best model by CV ROC-AUC: **{best_model_name}** (AUC={best_auc:.3f}).\n\n"
        f"## Top-10 important genes (saliency)\n\n```\n{top_imp}\n```\n"
        f"{attn_text}\n"
        f"Notes: saliency is the mean absolute gradient of the tumor logit w.r.t.\n"
        f"node features, averaged across the first 50 graphs. Large values point\n"
        f"to Cu-proteome genes whose expression most strongly shifts the model's\n"
        f"tumor probability — candidates worth biological follow-up."
    )
    return {"metrics": metrics_df, "summaries": summaries,
            "best_model": best_model_name, "importance": imp}


if __name__ == "__main__":
    run_all()
