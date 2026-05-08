"""Baseline analyses — run *before* any GNN experiment.

Steps
-----
C1. Gene coverage report
C2. Differential expression (tumor vs normal) on the Cu proteome
C3. Z-scored heat map of Cu-gene expression
C4. PCA (+ UMAP/t-SNE fallback)
C5. Static network visualisation with node color = logFC, size = centrality
C6. Community / module detection
Classical ML tumor-vs-normal classifiers (LR, RF, SVM) for later GNN comparison.

All outputs are written to ``outputs/baseline/``.
"""
from __future__ import annotations
from pathlib import Path
import json
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.model_selection import StratifiedKFold, cross_validate
from sklearn.metrics import (
    accuracy_score, balanced_accuracy_score, f1_score, roc_auc_score,
)

from src.preprocessing import LIHCDataset, load_lihc_dataset
from src.graph_building import (
    build_ppi_graph, build_functional_graph, build_coexpression_graph,
    save_graph, summarize_graph,
)
from src.utils import BASELINE_DIR, RANDOM_SEED

warnings.filterwarnings("ignore", category=RuntimeWarning)


# --------------------------------------------------------------------------- #
# C1. Gene coverage                                                            #
# --------------------------------------------------------------------------- #
def write_coverage_report(ds: LIHCDataset, out_dir: Path = BASELINE_DIR) -> None:
    present = ds.expression.index.tolist()
    missing = ds.missing_genes
    n_total = len(ds.copper_genes)
    pct = len(present) / n_total * 100 if n_total else 0.0
    lines = [
        "# Gene Coverage Report — Cu Proteome vs TCGA-LIHC",
        "",
        f"- Copper proteome size (Blockhuys 2017): **{n_total}**",
        f"- Genes found in LIHC expression matrix: **{len(present)}** ({pct:.1f}%)",
        f"- Genes missing: **{len(missing)}**",
        "",
        "## Missing genes",
    ]
    if missing:
        for g in missing:
            lines.append(f"- {g}")
    else:
        lines.append("- _none_")
    lines += ["", "## Present genes"]
    lines += [f"- {g}" for g in present]
    (out_dir / "gene_coverage_report.md").write_text("\n".join(lines))


# --------------------------------------------------------------------------- #
# C2. Differential expression                                                  #
# --------------------------------------------------------------------------- #
def differential_expression(ds: LIHCDataset, out_dir: Path = BASELINE_DIR) -> pd.DataFrame:
    """t-test per gene + Benjamini-Hochberg FDR."""
    tumor_mask = ds.tumor_mask
    normal_mask = ~tumor_mask
    expr = ds.expression.to_numpy()
    tumor = expr[:, tumor_mask]
    normal = expr[:, normal_mask]

    log2fc = tumor.mean(axis=1) - normal.mean(axis=1)
    tstat, pval = stats.ttest_ind(tumor, normal, axis=1, equal_var=False)
    pval = np.nan_to_num(pval, nan=1.0)

    order = np.argsort(pval)
    ranked = pval[order]
    m = len(pval)
    fdr_ranked = ranked * m / (np.arange(m) + 1)
    fdr_ranked = np.minimum.accumulate(fdr_ranked[::-1])[::-1]
    fdr = np.empty_like(fdr_ranked)
    fdr[order] = np.clip(fdr_ranked, 0, 1)

    df = pd.DataFrame({
        "gene_symbol": ds.expression.index,
        "mean_tumor": tumor.mean(axis=1),
        "mean_normal": normal.mean(axis=1),
        "log2FC": log2fc,
        "t_stat": tstat,
        "p_value": pval,
        "adj_p_BH": fdr,
        "significant_0.05_BH": fdr < 0.05,
    }).sort_values("adj_p_BH")
    df.to_csv(out_dir / "copper_de_results.csv", index=False)
    return df


# --------------------------------------------------------------------------- #
# C3. Heat map                                                                 #
# --------------------------------------------------------------------------- #
def plot_heatmap(ds: LIHCDataset, out_dir: Path = BASELINE_DIR) -> None:
    expr = ds.expression
    z = expr.sub(expr.mean(axis=1), axis=0).div(expr.std(axis=1).replace(0, 1), axis=0)
    sample_types = ds.metadata.loc[expr.columns, "sample_type"]
    col_order = (sample_types.sort_values(kind="stable").index.tolist())
    z = z[col_order]
    annotation = sample_types.loc[col_order]
    col_colors = annotation.map({"Tumor": "#B22222", "Normal": "#4682B4"}).values

    g = sns.clustermap(
        z, cmap="RdBu_r", center=0, vmin=-3, vmax=3,
        col_cluster=False, row_cluster=True,
        col_colors=col_colors, figsize=(12, 10),
        xticklabels=False, yticklabels=True, cbar_pos=(0.02, 0.8, 0.02, 0.15),
    )
    g.ax_heatmap.set_xlabel("Samples (blue=Normal, red=Tumor)")
    g.ax_heatmap.set_ylabel("Cu-proteome gene")
    g.fig.suptitle("TCGA-LIHC — Cu proteome expression (z-score)", y=1.02)
    g.savefig(out_dir / "copper_heatmap.png", dpi=150, bbox_inches="tight")
    plt.close("all")


# --------------------------------------------------------------------------- #
# C4. PCA / UMAP / t-SNE                                                       #
# --------------------------------------------------------------------------- #
def _scatter(emb: np.ndarray, labels: np.ndarray, path: Path, title: str) -> None:
    fig, ax = plt.subplots(figsize=(6, 5))
    for lab, color in zip(["Normal", "Tumor"], ["#4682B4", "#B22222"]):
        mask = labels == lab
        ax.scatter(emb[mask, 0], emb[mask, 1], s=25, alpha=0.75, label=lab, c=color)
    ax.set_title(title)
    ax.set_xlabel("component 1")
    ax.set_ylabel("component 2")
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_pca_and_umap(ds: LIHCDataset, out_dir: Path = BASELINE_DIR) -> None:
    X = ds.expression.T.to_numpy()
    X_s = StandardScaler().fit_transform(X)
    labels = ds.metadata.loc[ds.expression.columns, "sample_type"].to_numpy()

    pca = PCA(n_components=2, random_state=RANDOM_SEED).fit_transform(X_s)
    _scatter(pca, labels, out_dir / "pca_scatter.png",
             "PCA of TCGA-LIHC samples (Cu proteome features)")

    try:
        import umap  # noqa
        embed = umap.UMAP(n_components=2, random_state=RANDOM_SEED,
                          n_neighbors=min(15, max(2, X.shape[0] - 1))).fit_transform(X_s)
        method = "UMAP"
    except Exception:
        from sklearn.manifold import TSNE
        perp = min(30, max(5, X.shape[0] // 5))
        embed = TSNE(n_components=2, random_state=RANDOM_SEED,
                     perplexity=perp, init="pca").fit_transform(X_s)
        method = "t-SNE"
    _scatter(embed, labels, out_dir / "umap_scatter.png",
             f"{method} of TCGA-LIHC samples (Cu proteome features)")


# --------------------------------------------------------------------------- #
# C5. Network visualisation with logFC coloring                                #
# --------------------------------------------------------------------------- #
def plot_network_logfc(ds: LIHCDataset, de_df: pd.DataFrame,
                        out_dir: Path = BASELINE_DIR) -> nx.Graph:
    genes = ds.expression.index.tolist()
    G = build_functional_graph(genes, ds.copper_genes)

    logfc = dict(zip(de_df["gene_symbol"], de_df["log2FC"]))
    vmax = max(abs(min(logfc.values(), default=0)), abs(max(logfc.values(), default=0)), 1.0)
    centrality = nx.degree_centrality(G)

    pos = nx.spring_layout(G, seed=RANDOM_SEED, k=0.7)

    fig, ax = plt.subplots(figsize=(12, 10))
    node_colors = [logfc.get(n, 0.0) for n in G.nodes()]
    node_sizes = [300 + 2500 * centrality.get(n, 0.01) for n in G.nodes()]

    edge_styles = {"physical": ("red", 1.5),
                   "coexpression": ("purple", 1.0),
                   "genetic": ("green", 1.0),
                   "shared_compartment": ("lightgray", 0.6)}
    for etype, (color, width) in edge_styles.items():
        elist = [(u, v) for u, v, d in G.edges(data=True)
                 if d.get("edge_type") == etype]
        if elist:
            nx.draw_networkx_edges(G, pos, edgelist=elist, edge_color=color,
                                   width=width, alpha=0.7, ax=ax)

    nx.draw_networkx_nodes(G, pos, node_color=node_colors, cmap="RdBu_r",
                           vmin=-vmax, vmax=vmax, node_size=node_sizes,
                           edgecolors="black", linewidths=0.5, ax=ax)
    nx.draw_networkx_labels(G, pos, font_size=8, ax=ax)

    sm = plt.cm.ScalarMappable(cmap="RdBu_r",
                                norm=plt.Normalize(vmin=-vmax, vmax=vmax))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.6, pad=0.02)
    cbar.set_label("log2FC (tumor / normal)")
    ax.set_title("Cu-proteome functional network, colored by LIHC tumor-vs-normal log2FC")
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(out_dir / "copper_network_logfc.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    return G


# --------------------------------------------------------------------------- #
# C6. Community detection                                                      #
# --------------------------------------------------------------------------- #
def module_detection(G: nx.Graph, out_dir: Path = BASELINE_DIR) -> pd.DataFrame:
    try:
        communities = nx.community.greedy_modularity_communities(G)
    except Exception:
        communities = [set(G.nodes())]
    rows = []
    for cid, comm in enumerate(communities):
        for g in sorted(comm):
            rows.append({"gene_symbol": g, "module_id": cid, "module_size": len(comm)})
    df = pd.DataFrame(rows)
    df.to_csv(out_dir / "copper_modules.csv", index=False)

    summary = [
        "# Cu-proteome module analysis (LIHC graph)",
        "",
        f"- Graph: functional (curated Cu-edges + shared-compartment)",
        f"- Nodes: {G.number_of_nodes()}  Edges: {G.number_of_edges()}",
        f"- Modules detected: **{len(communities)}**",
        "",
    ]
    for cid, comm in enumerate(communities):
        summary.append(f"## Module {cid} ({len(comm)} genes)")
        summary.append(", ".join(sorted(comm)))
        summary.append("")
    (out_dir / "module_summary.md").write_text("\n".join(summary))
    return df


# --------------------------------------------------------------------------- #
# Classical ML baselines                                                       #
# --------------------------------------------------------------------------- #
def run_classical_baselines(ds: LIHCDataset, out_dir: Path = BASELINE_DIR) -> pd.DataFrame:
    X = ds.expression.T.to_numpy()
    y = (ds.metadata.loc[ds.expression.columns, "sample_type"]
         .str.lower().eq("tumor").astype(int).to_numpy())
    X = StandardScaler().fit_transform(X)

    min_class = int(min(y.sum(), (1 - y).sum()))
    n_splits = max(2, min(5, min_class))
    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=RANDOM_SEED)

    models = {
        "logreg": LogisticRegression(max_iter=2000, class_weight="balanced",
                                     random_state=RANDOM_SEED),
        "random_forest": RandomForestClassifier(n_estimators=300, class_weight="balanced",
                                                random_state=RANDOM_SEED, n_jobs=-1),
        "svm_rbf": SVC(kernel="rbf", class_weight="balanced", probability=True,
                        random_state=RANDOM_SEED),
    }
    scoring = {"accuracy": "accuracy", "balanced_accuracy": "balanced_accuracy",
               "f1": "f1", "roc_auc": "roc_auc"}

    rows = []
    for name, model in models.items():
        res = cross_validate(model, X, y, cv=cv, scoring=scoring, n_jobs=1)
        rows.append({
            "model": name,
            "cv_folds": n_splits,
            "accuracy": res["test_accuracy"].mean(),
            "balanced_accuracy": res["test_balanced_accuracy"].mean(),
            "f1": res["test_f1"].mean(),
            "roc_auc": res["test_roc_auc"].mean(),
            "accuracy_std": res["test_accuracy"].std(),
        })
    df = pd.DataFrame(rows)
    df.to_csv(out_dir / "classical_model_metrics.csv", index=False)
    return df


# --------------------------------------------------------------------------- #
# Main driver                                                                  #
# --------------------------------------------------------------------------- #
def run_all(ds: LIHCDataset | None = None) -> dict:
    if ds is None:
        ds = load_lihc_dataset()

    write_coverage_report(ds)
    de = differential_expression(ds)
    plot_heatmap(ds)
    plot_pca_and_umap(ds)

    genes = ds.expression.index.tolist()
    ppi = build_ppi_graph(genes)
    functional = build_functional_graph(genes, ds.copper_genes)
    coexpr = build_coexpression_graph(ds.expression)

    graph_dir = BASELINE_DIR / "graphs"
    save_graph(ppi, graph_dir, "ppi_graph", genes)
    save_graph(functional, graph_dir, "functional_graph", genes)
    save_graph(coexpr, graph_dir, "coexpression_graph", genes)

    graph_summary = {
        "ppi_graph": summarize_graph(ppi),
        "functional_graph": summarize_graph(functional),
        "coexpression_graph": summarize_graph(coexpr),
    }
    (BASELINE_DIR / "graphs" / "graph_summary.json").write_text(
        json.dumps(graph_summary, indent=2, default=list))

    plot_network_logfc(ds, de)
    modules = module_detection(functional)
    classical = run_classical_baselines(ds)

    return {
        "de_table": de, "graph_summary": graph_summary,
        "modules": modules, "classical_metrics": classical,
    }


if __name__ == "__main__":
    run_all()
