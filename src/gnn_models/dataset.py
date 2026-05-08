"""Build per-sample PyG Data objects for graph classification.

Each patient sample -> one ``torch_geometric.data.Data``:
  - node features: gene expression (+ optional z-score + category one-hot)
  - fixed edge_index / edge_weight from the chosen graph variant
  - y = 0 (normal) or 1 (tumor)
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
import pandas as pd
import torch
from torch_geometric.data import Data
import networkx as nx

from src.preprocessing import LIHCDataset


@dataclass
class GraphDataBundle:
    data_list: list[Data]
    gene_order: list[str]
    edge_index: torch.Tensor
    edge_weight: torch.Tensor
    in_dim: int


def _categorical_features(copper: pd.DataFrame, genes: list[str]) -> np.ndarray:
    cat_map = copper.set_index("gene_symbol")["functional_category"].to_dict()
    cats = ["transporter", "enzyme", "other_or_unknown"]
    feats = np.zeros((len(genes), len(cats)), dtype=np.float32)
    for i, g in enumerate(genes):
        c = cat_map.get(g, "other_or_unknown")
        feats[i, cats.index(c) if c in cats else 2] = 1.0
    return feats


def graph_to_edge_tensors(G: nx.Graph, gene_order: list[str]) -> tuple[torch.Tensor, torch.Tensor]:
    index = {g: i for i, g in enumerate(gene_order)}
    edges = []
    weights = []
    for u, v, d in G.edges(data=True):
        if u in index and v in index:
            edges.append((index[u], index[v]))
            edges.append((index[v], index[u]))
            w = float(d.get("weight", 1.0))
            weights.extend([w, w])
    if not edges:
        return (torch.zeros((2, 0), dtype=torch.long),
                torch.zeros((0,), dtype=torch.float32))
    ei = torch.tensor(edges, dtype=torch.long).t().contiguous()
    ew = torch.tensor(weights, dtype=torch.float32)
    return ei, ew


def build_graph_dataset(ds: LIHCDataset, graph: nx.Graph,
                        include_zscore: bool = True,
                        include_category: bool = True,
                        zscore_train_mask: np.ndarray | None = None) -> GraphDataBundle:
    """Build per-sample PyG Data objects.

    ``zscore_train_mask`` (length = n_samples) selects which columns are used
    for computing per-gene mean/std when include_zscore=True. If None, use all
    samples (legacy behaviour; has mild leakage inside CV).
    """
    gene_order = ds.expression.index.tolist()
    edge_index, edge_weight = graph_to_edge_tensors(graph, gene_order)

    expr = ds.expression
    if zscore_train_mask is not None:
        train_cols = expr.columns[zscore_train_mask]
        mu = expr[train_cols].mean(axis=1)
        sd = expr[train_cols].std(axis=1).replace(0, 1)
    else:
        mu = expr.mean(axis=1)
        sd = expr.std(axis=1).replace(0, 1)
    z = expr.sub(mu, axis=0).div(sd, axis=0)
    cat_feats = _categorical_features(ds.copper_genes, gene_order) if include_category else None

    labels = (ds.metadata.loc[expr.columns, "sample_type"]
              .str.lower().eq("tumor").astype(int).to_numpy())

    data_list: list[Data] = []
    for i, sample_id in enumerate(expr.columns):
        cols = [expr[sample_id].to_numpy().astype(np.float32).reshape(-1, 1)]
        if include_zscore:
            cols.append(z[sample_id].to_numpy().astype(np.float32).reshape(-1, 1))
        if cat_feats is not None:
            cols.append(cat_feats)
        x = np.concatenate(cols, axis=1)
        data = Data(
            x=torch.tensor(x, dtype=torch.float32),
            edge_index=edge_index,
            edge_weight=edge_weight,
            y=torch.tensor([int(labels[i])], dtype=torch.long),
        )
        data.sample_id = sample_id
        data_list.append(data)

    return GraphDataBundle(
        data_list=data_list, gene_order=gene_order,
        edge_index=edge_index, edge_weight=edge_weight,
        in_dim=data_list[0].x.shape[1],
    )
