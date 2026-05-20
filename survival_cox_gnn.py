"""Survival GNN with Cox proportional hazards loss.

Replaces binary classification with a proper survival model:
  - Output: single risk score per patient (not dead/alive)
  - Loss: Cox partial likelihood (standard survival analysis loss)
  - Uses ALL 369 tumor patients including censored ones
  - Metric: C-index (concordance index)

This is a fair comparison with Cox regression because:
  1. Same loss function (Cox partial likelihood)
  2. Same patients (all 369, including censored)
  3. Same metric (C-index)
  
The only difference is the GNN adds graph structure on top.

Expected improvement: from C-index 0.586 (binary GAT) toward 0.638 (Cox regression)
If GAT with Cox loss beats Cox regression -> graph adds survival value.
If they are equal -> graph adds interpretability only.

Outputs:
  outputs/final_comparison/survival_cox_gnn_results.csv
  outputs/final_comparison/survival_cox_gnn_report.md
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
from lifelines.utils import concordance_index
from torch_geometric.loader import DataLoader
from torch_geometric.data import Data

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from src.preprocessing import load_lihc_dataset
from src.graph_building import build_functional_graph
from src.gnn_models.dataset import graph_to_edge_tensors, _categorical_features
from src.gnn_models.models import GATGraphClassifier

OUT = ROOT / "outputs" / "final_comparison"
OUT.mkdir(parents=True, exist_ok=True)

SEEDS = list(range(1, 11))


def set_seeds(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


# ── Cox GNN model ──────────────────────────────────────────────────────────
class CoxGAT(nn.Module):
    """GAT with a single risk score output instead of classification head.
    
    Everything is identical to GATGraphClassifier except:
    - Output is 1 scalar (risk score) not 2 logits
    - No softmax — raw risk score fed to Cox loss
    """
    def __init__(self, in_dim: int, hidden_dim: int = 64,
                 n_heads: int = 4, n_layers: int = 2,
                 dropout: float = 0.4):
        super().__init__()
        from torch_geometric.nn import GATConv, global_mean_pool
        self.convs = nn.ModuleList()
        self.bns   = nn.ModuleList()
        prev = in_dim
        for i in range(n_layers):
            heads   = n_heads if i < n_layers - 1 else 1
            out_dim = hidden_dim
            self.convs.append(GATConv(prev, out_dim, heads=heads,
                                      concat=(i < n_layers - 1),
                                      dropout=dropout, add_self_loops=True))
            multi = heads if i < n_layers - 1 else 1
            self.bns.append(nn.BatchNorm1d(out_dim * multi))
            prev = out_dim * multi
        self.dropout = dropout
        self.pool    = global_mean_pool
        # single risk score output
        self.risk_head = nn.Sequential(
            nn.Linear(prev, prev // 2),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(prev // 2, 1),
        )

    def forward(self, x, edge_index, batch, edge_weight=None):
        from torch_geometric.nn import global_mean_pool
        import torch.nn.functional as F
        for conv, bn in zip(self.convs, self.bns):
            x = conv(x, edge_index)
            x = bn(x)
            x = F.elu(x)
            x = F.dropout(x, p=self.dropout, training=self.training)
        graph_emb = self.pool(x, batch)
        risk = self.risk_head(graph_emb).squeeze(-1)  # (batch_size,)
        return risk


# ── Cox partial likelihood loss ────────────────────────────────────────────
def cox_loss(risk_scores: torch.Tensor,
             durations:   torch.Tensor,
             events:      torch.Tensor) -> torch.Tensor:
    """
    Cox partial likelihood loss.

    For each patient i who experienced the event (died):
      loss_i = risk_i - log( sum of exp(risk_j) for all j at risk at time i )

    'At risk at time i' means duration_j >= duration_i.

    Parameters
    ----------
    risk_scores : (n,) higher = higher risk of dying sooner
    durations   : (n,) survival time in days
    events      : (n,) 1 = died, 0 = censored

    Returns
    -------
    negative mean partial log-likelihood (scalar, minimise this)
    """
    # sort by duration descending so risk set is a prefix
    order       = torch.argsort(durations, descending=True)
    risk_scores = risk_scores[order]
    events      = events[order]

    # log-sum-exp of risk scores for patients still at risk
    # (cumulative from the highest duration downward)
    log_cumsum  = torch.logcumsumexp(risk_scores, dim=0)

    # partial log-likelihood for events only
    event_mask  = events.bool()
    if event_mask.sum() == 0:
        return torch.tensor(0.0, requires_grad=True)

    pll = (risk_scores - log_cumsum)[event_mask]
    return -pll.mean()


# ── dataset builder ────────────────────────────────────────────────────────
def build_survival_dataset(ds, graph, md_tumor):
    """Build one Data object per patient with survival labels.
    
    Uses ALL tumor patients with any survival data — including censored.
    Labels stored as (duration, event) not binary dead/alive.
    """
    gene_order = ds.expression.index.tolist()
    ei, ew     = graph_to_edge_tensors(graph, gene_order)
    cat        = _categorical_features(ds.copper_genes, gene_order)

    # z-score expression across all samples
    expr = ds.expression
    mu   = expr.mean(axis=1); sd = expr.std(axis=1).replace(0, 1)
    ze   = expr.sub(mu, axis=0).div(sd, axis=0)

    data_list  = []
    durations  = []
    events     = []
    sample_ids = []
    groups     = []

    for sid in md_tumor.index:
        if sid not in expr.columns:
            continue
        dur = md_tumor.loc[sid, "overall_survival_days"]
        vit = md_tumor.loc[sid, "vital_status"]
        if pd.isna(dur) or pd.isna(vit):
            continue

        ev  = 1 if vit == "Dead" else 0
        ec  = ze[sid].to_numpy(dtype=np.float32).reshape(-1, 1)
        x   = np.concatenate([ec, cat], axis=1)
        d   = Data(x=torch.tensor(x, dtype=torch.float32),
                   edge_index=ei, edge_weight=ew,
                   y=torch.tensor([ev], dtype=torch.long))
        d.sample_id = sid
        d.duration  = float(dur)
        d.event     = ev

        data_list.append(d)
        durations.append(float(dur))
        events.append(ev)
        sample_ids.append(sid)
        groups.append(md_tumor.loc[sid, "case_submitter_id"])

    return (data_list,
            np.array(durations),
            np.array(events),
            np.array(groups),
            sample_ids)


# ── training loop ──────────────────────────────────────────────────────────
def train_eval_fold(train_list, val_list, in_dim, cfg_dict, seed):
    set_seeds(seed)
    device = "cuda" if torch.cuda.is_available() else "cpu"

    model = CoxGAT(
        in_dim,
        hidden_dim=cfg_dict["hidden"],
        n_heads=4,
        n_layers=2,
        dropout=cfg_dict["dropout"],
    ).to(device)

    opt = torch.optim.Adam(model.parameters(),
                           lr=cfg_dict["lr"],
                           weight_decay=cfg_dict["wd"])

    tr_dl = DataLoader(train_list, batch_size=cfg_dict["batch"], shuffle=True)
    va_dl = DataLoader(val_list,   batch_size=cfg_dict["batch"], shuffle=False)

    best_ci    = -1.0
    best_state = None

    for ep in range(cfg_dict["epochs"]):
        # ── train ──────────────────────────────────────────────────────────
        model.train()
        for batch in tr_dl:
            batch = batch.to(device)
            risk  = model(batch.x, batch.edge_index, batch.batch,
                          edge_weight=getattr(batch, "edge_weight", None))
            dur   = torch.tensor([d.duration for d in train_list
                                  if d.sample_id in
                                  [train_list[i].sample_id
                                   for i in range(len(train_list))]],
                                 dtype=torch.float32).to(device)
            # collect durations and events for this batch
            batch_dur = torch.tensor(
                [train_list[i].duration
                 for i in range(len(batch.y))],
                dtype=torch.float32).to(device)
            batch_ev  = batch.y.float().squeeze()

            loss = cox_loss(risk, batch_dur, batch_ev)
            if torch.isnan(loss) or torch.isinf(loss):
                continue
            opt.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            opt.step()

        # ── validate ───────────────────────────────────────────────────────
        model.eval()
        val_risks = []
        val_durs  = []
        val_evs   = []
        with torch.no_grad():
            for batch in va_dl:
                batch = batch.to(device)
                risk  = model(batch.x, batch.edge_index, batch.batch,
                              edge_weight=getattr(batch, "edge_weight", None))
                val_risks.extend(risk.cpu().numpy().tolist())
                val_durs.extend([d.duration for d in val_list])
                val_evs.extend(batch.y.cpu().numpy().tolist())

        val_risks = np.array(val_risks)
        val_durs  = np.array(val_durs)
        val_evs   = np.array(val_evs)

        if val_evs.sum() > 0:
            ci = concordance_index(val_durs, -val_risks, val_evs)
            if ci > best_ci:
                best_ci    = ci
                best_state = {k: v.clone()
                              for k, v in model.state_dict().items()}

    if best_state:
        model.load_state_dict(best_state)
    return model, best_ci, val_risks, val_durs, val_evs


# ── main ───────────────────────────────────────────────────────────────────
def main():
    print("[cox-gnn] loading dataset ...")
    ds    = load_lihc_dataset(require_real=True)
    graph = build_functional_graph(ds.expression.index.tolist(),
                                   ds.copper_genes)
    md    = ds.metadata.loc[ds.expression.columns].copy()
    md["overall_survival_days"] = pd.to_numeric(
        md["overall_survival_days"], errors="coerce")

    md_tumor = md[md["sample_type"] == "Tumor"].copy()

    print("[cox-gnn] building survival dataset (ALL patients incl. censored)...")
    data_list, durations, events, groups, sample_ids = \
        build_survival_dataset(ds, graph, md_tumor)

    n_total  = len(data_list)
    n_events = int(events.sum())
    n_cens   = n_total - n_events
    print(f"  n={n_total}  events={n_events}  censored={n_cens}")
    print(f"  (previous binary GAT used only {n_events + int((durations[events==0]>1095).sum())} confirmed patients)")

    # stratify by event for CV
    min_class = int(min(events.sum(), (1-events).sum()))
    n_splits  = max(2, min(5, min_class))

    cfg = {
        "hidden":  64,
        "dropout": 0.4,
        "lr":      1e-3,
        "wd":      5e-4,
        "batch":   32,
        "epochs":  100,
    }

    seed_cis = []

    for seed in SEEDS:
        print(f"  seed {seed:2d}/{SEEDS[-1]} ...", end=" ", flush=True)
        set_seeds(seed)

        cv = StratifiedGroupKFold(n_splits=n_splits,
                                  shuffle=True, random_state=seed)

        fold_cis = []
        for fold, (tr_idx, va_idx) in enumerate(
                cv.split(np.zeros(len(data_list)), events, groups)):

            set_seeds(seed + fold * 100)
            tr_list = [data_list[i] for i in tr_idx]
            va_list = [data_list[i] for i in va_idx]

            _, ci, _, _, _ = train_eval_fold(
                tr_list, va_list,
                data_list[0].x.shape[1],
                cfg, seed + fold * 100,
            )
            fold_cis.append(ci)

        mean_ci = float(np.mean(fold_cis))
        print(f"C-index={mean_ci:.4f}")
        seed_cis.append(mean_ci)

    cox_gnn_mean = float(np.mean(seed_cis))
    cox_gnn_std  = float(np.std(seed_cis))

    # ── comparison table ───────────────────────────────────────────────────
    print("\n" + "="*60)
    print("SURVIVAL COMPARISON — FINAL")
    print("="*60)
    print(f"  Random baseline      : 0.500")
    print(f"  Binary GAT (old)     : 0.586 ± 0.020  (198 patients, binary loss)")
    print(f"  Cox regression       : 0.638 ± 0.010  (369 patients, Cox loss)")
    print(f"  Cox GNN (new)        : {cox_gnn_mean:.3f} ± {cox_gnn_std:.3f}  "
          f"({n_total} patients, Cox loss + graph)")

    gap = cox_gnn_mean - 0.638
    if gap > 0.02:
        verdict = "Cox GNN beats Cox regression — graph adds survival value"
    elif gap > -0.01:
        verdict = "Cox GNN matches Cox regression — graph adds interpretability"
    else:
        verdict = "Cox regression still wins — larger dataset or better features needed"

    print(f"\n  Cox GNN vs Cox regression: {gap:+.3f}")
    print(f"  Verdict: {verdict}")

    # save
    rows = [
        {"model": "Random baseline",    "n_patients": n_total,
         "loss": "—",           "ci_mean": 0.500, "ci_std": 0.000},
        {"model": "Binary GAT (old)",   "n_patients": 198,
         "loss": "CrossEntropy", "ci_mean": 0.586, "ci_std": 0.020},
        {"model": "Cox regression",     "n_patients": 369,
         "loss": "Cox PL",       "ci_mean": 0.638, "ci_std": 0.010},
        {"model": "Cox GNN (new)",      "n_patients": n_total,
         "loss": "Cox PL",       "ci_mean": round(cox_gnn_mean, 4),
         "ci_std": round(cox_gnn_std, 4)},
    ]
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "survival_cox_gnn_results.csv", index=False)

    report = f"""# Survival GNN with Cox Proportional Hazards Loss

## What changed vs original pipeline

| aspect | binary GAT (old) | Cox GNN (new) |
|---|---|---|
| loss function | CrossEntropy (dead/alive) | Cox partial likelihood |
| patients used | 198 (confirmed outcome only) | **{n_total} (all incl. censored)** |
| metric | ROC-AUC | **C-index** |
| fair vs Cox? | No | **Yes** |

## Results (10 seeds, 5-fold StratifiedGroupKFold)

| model | n patients | loss | C-index | std |
|---|---:|---|---:|---:|
| Random baseline | {n_total} | — | 0.500 | — |
| Binary GAT (old) | 198 | CrossEntropy | 0.586 | 0.020 |
| Cox regression | 369 | Cox PL | 0.638 | 0.010 |
| **Cox GNN (new)** | **{n_total}** | **Cox PL** | **{cox_gnn_mean:.3f}** | **{cox_gnn_std:.3f}** |

## Verdict

Cox GNN vs Cox regression gap: **{gap:+.3f}**

{verdict}

## Why this comparison is now fair

The original binary GAT threw away {369 - 198} censored patients and used
a classification loss — fundamentally different from Cox regression.
This Cox GNN uses the same loss function (partial likelihood) and the
same patients (all {n_total} with survival data). The only difference
is the GNN adds copper proteome graph structure on top of Cox.

## Biological interpretation

If Cox GNN > Cox regression: the copper proteome graph edges carry
survival-relevant information beyond what individual gene expression
captures. The attention weights on the survival model would show
which copper biology edges are most predictive of survival.

If Cox GNN ≈ Cox regression: the graph structure does not add
survival prediction value but adds interpretability — the attention
edges still recover canonical copper biology.

## Files produced
- `outputs/final_comparison/survival_cox_gnn_results.csv`
- `outputs/final_comparison/survival_cox_gnn_report.md`
"""
    (OUT / "survival_cox_gnn_report.md").write_text(report, encoding="utf-8")
    print(f"\n[cox-gnn] wrote survival_cox_gnn_results.csv")
    print(f"[cox-gnn] wrote survival_cox_gnn_report.md")
    print("\nDone.")


if __name__ == "__main__":
    main()
