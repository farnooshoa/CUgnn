"""GCN and GAT graph-classification models for the Cu proteome."""
from __future__ import annotations
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GCNConv, GATConv, global_mean_pool, global_add_pool


class GCNGraphClassifier(nn.Module):
    def __init__(self, in_dim: int, hidden_dim: int = 64, n_classes: int = 2,
                 n_layers: int = 2, dropout: float = 0.4, pool: str = "mean"):
        super().__init__()
        self.convs = nn.ModuleList()
        self.bns = nn.ModuleList()
        prev = in_dim
        for _ in range(n_layers):
            self.convs.append(GCNConv(prev, hidden_dim))
            self.bns.append(nn.BatchNorm1d(hidden_dim))
            prev = hidden_dim
        self.dropout = dropout
        self.pool = global_mean_pool if pool == "mean" else global_add_pool
        self.classifier = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim), nn.ReLU(),
            nn.Dropout(dropout), nn.Linear(hidden_dim, n_classes),
        )

    def forward(self, x, edge_index, batch, edge_weight=None, return_embedding: bool = False):
        for conv, bn in zip(self.convs, self.bns):
            x = conv(x, edge_index, edge_weight=edge_weight)
            x = bn(x)
            x = F.relu(x)
            x = F.dropout(x, p=self.dropout, training=self.training)
        graph_emb = self.pool(x, batch)
        out = self.classifier(graph_emb)
        if return_embedding:
            return out, graph_emb, x
        return out


class GATGraphClassifier(nn.Module):
    def __init__(self, in_dim: int, hidden_dim: int = 64, n_classes: int = 2,
                 n_heads: int = 4, n_layers: int = 2, dropout: float = 0.4,
                 pool: str = "mean"):
        super().__init__()
        self.convs = nn.ModuleList()
        self.bns = nn.ModuleList()
        prev = in_dim
        for i in range(n_layers):
            heads = n_heads if i < n_layers - 1 else 1
            out_dim = hidden_dim
            self.convs.append(GATConv(prev, out_dim, heads=heads,
                                      concat=(i < n_layers - 1),
                                      dropout=dropout, add_self_loops=True))
            multi = heads if i < n_layers - 1 else 1
            self.bns.append(nn.BatchNorm1d(out_dim * multi))
            prev = out_dim * multi
        self.dropout = dropout
        self.pool = global_mean_pool if pool == "mean" else global_add_pool
        self.classifier = nn.Sequential(
            nn.Linear(prev, prev), nn.ReLU(),
            nn.Dropout(dropout), nn.Linear(prev, n_classes),
        )
        self._last_attn: list[tuple[torch.Tensor, torch.Tensor]] = []

    def forward(self, x, edge_index, batch, edge_weight=None,
                return_embedding: bool = False, return_attention: bool = False):
        self._last_attn = []
        for i, (conv, bn) in enumerate(zip(self.convs, self.bns)):
            if return_attention:
                x, (ei, alpha) = conv(x, edge_index, return_attention_weights=True)
                self._last_attn.append((ei.detach(), alpha.detach()))
            else:
                x = conv(x, edge_index)
            x = bn(x)
            x = F.elu(x)
            x = F.dropout(x, p=self.dropout, training=self.training)
        graph_emb = self.pool(x, batch)
        out = self.classifier(graph_emb)
        if return_embedding:
            return out, graph_emb, x
        return out

    def get_last_attention(self):
        return self._last_attn
