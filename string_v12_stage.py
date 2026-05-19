"""Stage classification using STRING v12 edges instead of hand-curated edges.

This directly addresses the scientific criticism that the curated edge list
is not a substitute for a proper versioned database query.

STRING v12 edges are in data/string_v12_copper_edges.tsv:
  columns: source, target, score (combined score 0-1)
  77 edges, all copper proteome gene pairs

We compare three graph variants on stage classification:
  1. Hand-curated edges (our previous result, baseline)
  2. STRING v12 edges only
  3. STRING v12 + curated (union graph)

All use the same 349 tumor patients, same 5-fold StratifiedGroupKFold, 10 seeds.

Outputs:
  outputs/final_comparison/string_v12_stage_comparison.csv
  outputs/final_comparison/string_v12_stage_report.md
"""
from __future__ import annotations
import sys
import random
from pathlib import Path
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import networkx as nx
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.metrics import roc_auc_score
from torch_geometric.loader import DataLoader

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from src.preprocessing import load_lihc_dataset
from src.graph_building import build_functional_graph
from src.graph_building.build_graph import CURATED_CU_EDGES
from src.gnn_models.dataset import build_graph_dataset
from src.gnn_models.models import GATGraphClassifier
from src.gnn_models.train import TrainConfig, _class_weights

OUT = ROOT / "outputs" / "final_comparison"
OUT.mkdir(parents=True, exist_ok=True)
STRING_FILE = ROOT / "data" / "string_v12_copper_edges.tsv"
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


# ── graph builders ─────────────────────────────────────────────────────────
def build_string_graph(genes: list[str],
                       score_threshold: float = 0.4) -> nx.Graph:
    """Build graph from STRING v12 edges only.
    
    score_threshold: minimum combined score to include edge (0-1).
    Default 0.4 = medium confidence, standard in the field.
    """
    G = nx.Graph()
    G.add_nodes_from(genes)
    genes_set = set(genes)

    string_df = pd.read_csv(STRING_FILE, sep="\t")
    string_df.columns = [c.lower().strip() for c in string_df.columns]

    n_total = 0
    n_added = 0
    for _, row in string_df.iterrows():
        s = str(row["source"]).upper().strip()
        t = str(row["target"]).upper().strip()
        w = float(row["score"])
        n_total += 1
        if s in genes_set and t in genes_set and s != t and w >= score_threshold:
            G.add_edge(s, t, weight=w, edge_type="string_ppi")
            n_added += 1

    print(f"  STRING graph: {n_added}/{n_total} edges above score {score_threshold}")
    print(f"  Nodes: {G.number_of_nodes()}  Edges: {G.number_of_edges()}")
    isolates = [n for n in G.nodes() if G.degree(n) == 0]
    print(f"  Isolated nodes: {len(isolates)}")
    return G


def build_union_graph(genes: list[str],
                      copper: pd.DataFrame,
                      score_threshold: float = 0.4) -> nx.Graph:
    """Union of STRING v12 + hand-curated edges."""
    # start with functional graph (curated + compartment)
    G = build_functional_graph(genes, copper)

    genes_set = set(genes)
    string_df = pd.read_csv(STRING_FILE, sep="\t")
    string_df.columns = [c.lower().strip() for c in string_df.columns]

    n_added = 0
    for _, row in string_df.iterrows():
        s = str(row["source"]).upper().strip()
        t = str(row["target"]).upper().strip()
        w = float(row["score"])
        if s in genes_set and t in genes_set and s != t and w >= score_threshold:
            if not G.has_edge(s, t):
                G.add_edge(s, t, weight=w, edge_type="string_ppi")
                n_added += 1

    print(f"  Union graph: {n_added} STRING edges added on top of curated")
    print(f"  Nodes: {G.number_of_nodes()}  Edges: {G.number_of_edges()}")
    return G


# ── training ───────────────────────────────────────────────────────────────
def _train_eval_fold(train_list, val_list, in_dim, cfg):
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
    auc = roc_auc_score(ys, ps) if len(set(ys)) > 1 else float("nan")
    return auc


def run_gat(graph_name, graph, ds, y, groups, sample_indices, seeds, epochs=80):
    """Run GAT with given graph across seeds."""
    seed_aucs = []
    for seed in seeds:
        print(f"  seed {seed:2d}/{seeds[-1]} ...", end=" ", flush=True)
        set_seeds(seed)

        min_class = int(min((y == 0).sum(), (y == 1).sum()))
        n_splits  = max(2, min(5, min_class))
        cv  = StratifiedGroupKFold(n_splits=n_splits,
                                   shuffle=True, random_state=seed)
        cfg = TrainConfig(model="gat", epochs=epochs)

        fold_aucs = []
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
            auc = _train_eval_fold(
                [active[i] for i in tr_idx],
                [active[i] for i in va_idx],
                bundle.in_dim, cfg,
            )
            fold_aucs.append(auc)

        m = float(np.nanmean(fold_aucs))
        print(f"AUC={m:.4f}")
        seed_aucs.append(m)

    return float(np.mean(seed_aucs)), float(np.std(seed_aucs))


# ── main ───────────────────────────────────────────────────────────────────
def main():
    print("[string] loading dataset ...")
    ds    = load_lihc_dataset(require_real=True)
    genes = ds.expression.index.tolist()
    md    = ds.metadata.loc[ds.expression.columns].copy()

    # stage subset
    stage_bin = md["stage"].map(classify_stage)
    keep      = (md["sample_type"] == "Tumor") & stage_bin.notna()
    md_stage  = md[keep].copy()
    y         = stage_bin[keep].astype(int).to_numpy()
    groups    = md_stage["case_submitter_id"].to_numpy()
    idx       = np.array([ds.expression.columns.tolist().index(s)
                          for s in md_stage.index])

    n_early = int((y == 0).sum())
    n_late  = int((y == 1).sum())
    print(f"  n={len(y)}  early={n_early}  late={n_late}")

    rows = []

    # ── 1. hand-curated (previous result — no rerun needed) ────────────────
    print("\n[string] === Graph 1: Hand-curated edges (previous result) ===")
    print("  Skipping rerun — using validated result from week2_stage_comparison.py")
    rows.append({
        "graph":        "hand_curated",
        "n_edges":      74,
        "source":       "CURATED_CU_EDGES + shared_compartment",
        "auc_mean":     0.6667,
        "auc_std":      0.0174,
        "note":         "from week2_stage_comparison.py, 10 seeds"
    })

    # ── 2. STRING v12 only ─────────────────────────────────────────────────
    print("\n[string] === Graph 2: STRING v12 edges only (score >= 0.4) ===")
    string_graph = build_string_graph(genes, score_threshold=0.4)
    auc_m, auc_s = run_gat("string_v12", string_graph,
                             ds, y, groups, idx, SEEDS)
    rows.append({
        "graph":    "string_v12",
        "n_edges":  string_graph.number_of_edges(),
        "source":   "STRING v12 combined score >= 0.4",
        "auc_mean": auc_m,
        "auc_std":  auc_s,
        "note":     "10 seeds"
    })

    # ── 3. union graph ─────────────────────────────────────────────────────
    print("\n[string] === Graph 3: STRING v12 + hand-curated (union) ===")
    union_graph = build_union_graph(genes, ds.copper_genes, score_threshold=0.4)
    auc_m2, auc_s2 = run_gat("union", union_graph,
                               ds, y, groups, idx, SEEDS)
    rows.append({
        "graph":    "union_string_curated",
        "n_edges":  union_graph.number_of_edges(),
        "source":   "STRING v12 + CURATED_CU_EDGES + shared_compartment",
        "auc_mean": auc_m2,
        "auc_std":  auc_s2,
        "note":     "10 seeds"
    })

    # ── summary ────────────────────────────────────────────────────────────
    df = pd.DataFrame(rows).sort_values("auc_mean", ascending=False)
    df.to_csv(OUT / "string_v12_stage_comparison.csv", index=False)

    print("\n" + "="*65)
    print("STRING v12 GRAPH COMPARISON — STAGE CLASSIFICATION")
    print("="*65)
    print(f"\n{'graph':<28} {'edges':>6} {'AUC':>8} {'± std':>7}")
    print("-"*52)
    for _, r in df.iterrows():
        print(f"{r['graph']:<28} {r['n_edges']:>6} "
              f"{r['auc_mean']:>8.4f} {r['auc_std']:>7.4f}")

    best = df.iloc[0]
    curated = df[df["graph"] == "hand_curated"].iloc[0]
    string  = df[df["graph"] == "string_v12"].iloc[0]
    union   = df[df["graph"] == "union_string_curated"].iloc[0]

    gap_string = string["auc_mean"] - curated["auc_mean"]
    gap_union  = union["auc_mean"]  - curated["auc_mean"]

    print(f"\nSTRING vs curated : {gap_string:+.4f}")
    print(f"Union  vs curated : {gap_union:+.4f}")

    if string["auc_mean"] >= curated["auc_mean"] - 0.01:
        verdict_string = "STRING edges reproduce the result — biological interpretation is robust to edge source"
    else:
        verdict_string = "STRING edges score lower — curated edges carry domain knowledge STRING misses"

    if union["auc_mean"] > curated["auc_mean"] + 0.01:
        verdict_union = "Union graph improves on both — combining sources adds value"
    elif union["auc_mean"] >= curated["auc_mean"] - 0.01:
        verdict_union = "Union graph is equivalent — additional STRING edges do not hurt"
    else:
        verdict_union = "Union graph is worse — additional edges add noise"

    print(f"\nVerdict STRING : {verdict_string}")
    print(f"Verdict Union  : {verdict_union}")

    # ── report ─────────────────────────────────────────────────────────────
    report = f"""# STRING v12 Graph Comparison — Stage Classification

## Motivation

A key scientific criticism of the original pipeline is that the hand-curated
edge list (`CURATED_CU_EDGES`) is not a substitute for a proper versioned
database query. This analysis addresses that criticism by rebuilding the
copper proteome graph using STRING v12 edges
(`data/string_v12_copper_edges.tsv`) and rerunning the stage classification.

## Graph variants compared

| graph | edges | source |
|---|---:|---|
| hand_curated | 74 | CURATED_CU_EDGES + shared_compartment |
| string_v12 | {string_graph.number_of_edges()} | STRING v12 combined score ≥ 0.4 |
| union | {union_graph.number_of_edges()} | STRING v12 + hand-curated |

STRING score threshold 0.4 = medium confidence (standard in the field).

## Results (GAT, 10 seeds, 5-fold StratifiedGroupKFold)

| graph | ROC-AUC | std |
|---|---:|---:|
| hand_curated | {curated['auc_mean']:.4f} | {curated['auc_std']:.4f} |
| string_v12 | {string['auc_mean']:.4f} | {string['auc_std']:.4f} |
| union | {union['auc_mean']:.4f} | {union['auc_std']:.4f} |

## Interpretation

**STRING vs curated: {gap_string:+.4f} AUC**
{verdict_string}

**Union vs curated: {gap_union:+.4f} AUC**
{verdict_union}

## What this means for the paper

The result {'holds' if string['auc_mean'] >= curated['auc_mean'] - 0.01 else 'changes'}
when using a versioned database source instead of hand-curated edges.
{'This demonstrates that the biological interpretation is robust to the choice of edge source — a key requirement for a publishable graph-based analysis.' if string['auc_mean'] >= curated['auc_mean'] - 0.01 else 'This suggests the hand-curated edges carry domain knowledge that STRING combined scores do not fully capture at this threshold.'}

The STRING v12 edge file is pinned at version 12 and is fully reproducible
from `data/string_v12_copper_edges.tsv`. All future analyses should use this
as the primary edge source, with hand-curated edges as a biological prior
for edges not captured in STRING (e.g. the histone Cu-reductase edges from
Attar et al. 2020).

## Recommended edge strategy going forward

1. Primary: STRING v12 combined score ≥ 0.4 (reproducible, versioned)
2. Supplement: hand-curated edges for biology not in STRING
   (histone H3-H4 Cu-reductase, Au-compound interactions)
3. Report both as a sensitivity analysis

## Files produced
- `outputs/final_comparison/string_v12_stage_comparison.csv`
- `outputs/final_comparison/string_v12_stage_report.md`
"""
    (OUT / "string_v12_stage_report.md").write_text(report, encoding="utf-8")
    print(f"\n[string] wrote string_v12_stage_comparison.csv")
    print(f"[string] wrote string_v12_stage_report.md")
    print("\nDone.")


if __name__ == "__main__":
    main()
