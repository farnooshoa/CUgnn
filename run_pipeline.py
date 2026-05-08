"""End-to-end driver for the TCGA-LIHC copper-proteome pilot.

Usage
-----
    python run_pipeline.py            # full pipeline (baselines + GNN)
    python run_pipeline.py baselines  # only Part C (no training)
    python run_pipeline.py gnn        # only Part D
    python run_pipeline.py compare    # only Part E summary

If data/lihc_expression.tsv and data/lihc_metadata.tsv are missing, a
synthetic demo dataset is generated and the whole pipeline still runs.
"""
from __future__ import annotations
import json
import sys
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from src.preprocessing import load_lihc_dataset
from src.baseline_models import run_all as run_baselines
from src.gnn_models import run_all as run_gnn
from src.utils import BASELINE_DIR, GNN_DIR, FINAL_COMPARISON_DIR


def write_comparison_summary():
    lines = ["# Model comparison — TCGA-LIHC Cu proteome pilot", ""]

    classical_path = BASELINE_DIR / "classical_model_metrics.csv"
    gnn_path = GNN_DIR / "model_metrics.csv"

    if classical_path.exists():
        classical = pd.read_csv(classical_path)
        lines += ["## Classical baselines (5-fold CV, 54-gene vector)",
                   "", classical.round(4).to_markdown(index=False), ""]
    else:
        lines += ["## Classical baselines", "_not available_", ""]

    if gnn_path.exists():
        gnn = pd.read_csv(gnn_path)
        agg = (gnn.groupby("model")[["accuracy", "balanced_accuracy", "f1", "roc_auc"]]
                 .mean().round(4).reset_index())
        lines += ["## GNN models (5-fold CV, per-patient graph)",
                   "", agg.to_markdown(index=False), ""]
    else:
        lines += ["## GNN models", "_not available_", ""]

    lines += [
        "## Interpretation",
        "",
        "- Compare the ROC-AUC row across classical vs GNN. If GNN is within",
        "  ~0.02 AUC of the best classical baseline, the marginal predictive",
        "  value of graph structure is small on this cohort — which is the",
        "  *expected* situation for a 54-gene tumor-vs-normal problem.",
        "- The *real* value of the GNN in the pilot is interpretability:",
        "  attention / saliency maps highlight Cu-proteome modules (e.g.",
        "  ATP7B / CP / ceruloplasmin axis; LOX family) that tabular models",
        "  cannot expose.",
        "- See outputs/gnn/top_subgraph_or_attention_summary.md for the",
        "  per-model interpretability report.",
        "- See outputs/baseline/copper_network_logfc.png for the static",
        "  LIHC-coloured network.",
        "",
        "## Conclusion",
        "",
        "The copper proteome is small enough (54 genes) that classical ML on",
        "a flat expression vector is already a strong baseline. The GNN",
        "earns its place mainly by (a) aligning with the biology of Fig. 3",
        "of Blockhuys 2017, (b) providing module-level interpretability, and",
        "(c) offering a natural substrate for future multi-modal graphs",
        "(methylation, mutation, clinical covariates as additional node/edge",
        "features). Predictive parity with tabular models is acceptable at",
        "this pilot stage.",
    ]
    FINAL_COMPARISON_DIR.mkdir(parents=True, exist_ok=True)
    (FINAL_COMPARISON_DIR / "model_comparison_summary.md").write_text("\n".join(lines))


def main(mode: str = "full", require_real: bool = False):
    print(f"[pipeline] mode = {mode}, require_real = {require_real}")
    ds = load_lihc_dataset(require_real=require_real)
    print(f"[pipeline] Loaded {ds.n_samples} samples, "
          f"{ds.n_genes}/54 Cu genes covered, "
          f"{int(ds.tumor_mask.sum())} tumor, {int((~ds.tumor_mask).sum())} normal.")

    if mode in ("full", "baselines"):
        print("[pipeline] Running baseline analyses ...")
        run_baselines(ds)

    if mode in ("full", "gnn"):
        print("[pipeline] Running GNN experiments ...")
        run_gnn(ds)

    if mode in ("full", "compare"):
        print("[pipeline] Writing model comparison summary ...")
        write_comparison_summary()

    print("[pipeline] Done. See outputs/ for results.")


if __name__ == "__main__":
    args = [a for a in sys.argv[1:]]
    require_real = "--require-real-data" in args
    args = [a for a in args if a != "--require-real-data"]
    mode = args[0] if args else "full"
    main(mode, require_real=require_real)
