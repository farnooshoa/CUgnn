"""Construct fixed-topology graphs over the 54-gene Cu proteome.

Three graph variants are supported:
  1. ppi_graph        - STRING / BioGRID-style physical interactions
  2. functional_graph - union of physical + curated pathway + shared compartment
  3. coexpression_graph - exploratory; built from TCGA-LIHC expression only

The primary graph is the PPI graph. If no external interaction table is
available, a biologically informed fallback is produced from hand-curated
Cu-homeostasis edges + shared-localisation edges + shared-category edges, so
the pipeline remains runnable offline.
"""
from __future__ import annotations
from pathlib import Path
import itertools
import numpy as np
import pandas as pd
import networkx as nx

from src.utils import RANDOM_SEED

CURATED_CU_EDGES: list[tuple[str, str, str]] = [
    ("SLC31A1", "ATOX1", "physical"),
    ("SLC31A1", "CCS", "physical"),
    ("SLC31A2", "ATOX1", "physical"),
    ("ATOX1", "ATP7A", "physical"),
    ("ATOX1", "ATP7B", "physical"),
    ("ATOX1", "SOD3", "genetic"),
    ("CCS", "SOD1", "physical"),
    ("COMMD1", "ATP7B", "physical"),
    ("COX17", "COX11", "physical"),
    ("COX17", "SCO1", "physical"),
    ("COX17", "SCO2", "physical"),
    ("SCO1", "MT-CO1", "physical"),
    ("SCO2", "MT-CO2", "physical"),
    ("MT-CO1", "MT-CO2", "physical"),
    ("COX11", "MT-CO1", "physical"),
    ("ATP7A", "ATP7B", "physical"),
    ("ATP7B", "CP", "physical"),
    ("CP", "HEPH", "physical"),
    ("CP", "HEPHL1", "physical"),
    ("LOX", "LOXL1", "coexpression"),
    ("LOX", "LOXL2", "coexpression"),
    ("LOX", "LOXL3", "coexpression"),
    ("LOX", "LOXL4", "coexpression"),
    ("LOXL1", "LOXL2", "coexpression"),
    ("LOXL2", "LOXL3", "coexpression"),
    ("LOXL3", "LOXL4", "coexpression"),
    ("LOX", "SPARC", "physical"),
    ("SPARC", "LOXL1", "physical"),
    ("AOC1", "AOC2", "coexpression"),
    ("AOC2", "AOC3", "coexpression"),
    ("AOC1", "AOC3", "coexpression"),
    ("SOD1", "SOD3", "coexpression"),
    ("S100A5", "S100A12", "coexpression"),
    ("S100A5", "S100A13", "coexpression"),
    ("S100A5", "S100B", "coexpression"),
    ("S100A12", "S100A13", "coexpression"),
    ("S100A12", "S100B", "coexpression"),
    ("S100A13", "S100B", "coexpression"),
    ("MT3", "MT4", "coexpression"),
    ("MOXD1", "MOXD2P", "coexpression"),
    ("TYR", "TYRP1", "physical"),
    ("TYR", "DBH", "coexpression"),
    ("TYRP1", "DBH", "coexpression"),
    ("MAP2K1", "MEMO1", "genetic"),
    ("MEMO1", "SPARC", "genetic"),
    ("PARK7", "SOD1", "genetic"),
    ("APP", "PRNP", "physical"),
    ("APP", "SNCA", "physical"),
    ("PRNP", "SNCA", "physical"),
    ("ATP7B", "AFP", "coexpression"),
    ("LTF", "CP", "coexpression"),
    ("LTF", "ALB", "coexpression"),
    ("ALB", "CP", "coexpression"),
    ("F5", "ALB", "coexpression"),
    ("CUTA", "CUTC", "coexpression"),
    ("CUTA", "COMMD1", "physical"),
    ("GPC1", "SPARC", "physical"),
    # Histone H3-H4 tetramer Cu2+ reductase (Attar et al. Science 2020)
    ("H3-3A", "H4C1", "physical"),
    ("H3-3B", "H4C1", "physical"),
    ("H3C1", "H4C1", "physical"),
    ("H3-3A", "H3-3B", "coexpression"),
    ("H3-3A", "H3C1", "coexpression"),
    ("H3-3B", "H3C1", "coexpression"),
    # Paper-proposed coupling between histone Cu-reductase and the cytosolic Cu system
    ("H3-3A", "ATOX1", "genetic"),
    ("H3-3B", "ATOX1", "genetic"),
    ("H3C1", "ATOX1", "genetic"),
    # H3H113N phenocopies ctr1Δ and impairs Sod1 (paper Fig. 5)
    ("H3-3A", "SLC31A1", "genetic"),
    ("H3-3A", "SOD1", "genetic"),
    ("H3-3B", "SOD1", "genetic"),
    # Mitochondrial respiration defect in H3H113N (paper Fig. 5A-C)
    ("H3-3A", "MT-CO1", "genetic"),
    ("H3-3A", "MT-CO2", "genetic"),
]


def _empty_graph(genes: list[str]) -> nx.Graph:
    G = nx.Graph()
    G.add_nodes_from(genes)
    return G


def build_ppi_graph(genes: list[str],
                    external_edges: pd.DataFrame | None = None) -> nx.Graph:
    """PPI graph restricted to the Cu proteome.

    ``external_edges`` (optional) should have columns ['source', 'target',
    'score'] and be pre-filtered to the Cu proteome.
    """
    G = _empty_graph(genes)
    genes_set = set(genes)
    if external_edges is not None and not external_edges.empty:
        for row in external_edges.itertuples(index=False):
            s = str(row.source).upper().strip()
            t = str(row.target).upper().strip()
            if s in genes_set and t in genes_set and s != t:
                w = float(getattr(row, "score", 1.0))
                G.add_edge(s, t, weight=w, edge_type="physical")
    else:
        for s, t, et in CURATED_CU_EDGES:
            if s in genes_set and t in genes_set and et == "physical":
                G.add_edge(s, t, weight=1.0, edge_type="physical")
    return G


def build_functional_graph(genes: list[str], copper: pd.DataFrame) -> nx.Graph:
    """PPI edges + shared-compartment edges + shared-category edges."""
    G = _empty_graph(genes)
    genes_set = set(genes)

    for s, t, et in CURATED_CU_EDGES:
        if s in genes_set and t in genes_set:
            G.add_edge(s, t, weight=1.0, edge_type=et)

    cop = copper.set_index("gene_symbol")
    cop = cop.loc[[g for g in genes if g in cop.index]]
    loc = cop["subcellular_localization"].fillna("").str.split("|")

    compartment_members: dict[str, list[str]] = {}
    for gene, comps in loc.items():
        for c in comps:
            c = c.strip()
            if not c:
                continue
            compartment_members.setdefault(c, []).append(gene)

    for comp, members in compartment_members.items():
        if 2 <= len(members) <= 8:
            for a, b in itertools.combinations(members, 2):
                if not G.has_edge(a, b):
                    G.add_edge(a, b, weight=0.5, edge_type="shared_compartment")
    return G


def build_coexpression_graph(expr: pd.DataFrame, top_k: int = 4,
                              min_corr: float = 0.35) -> nx.Graph:
    """kNN graph on Pearson correlation between gene rows.

    ``expr`` rows must be the Cu proteome only. Edges are symmetrised.
    Warning: sample-set-dependent; secondary / exploratory only.
    """
    genes = list(expr.index)
    G = _empty_graph(genes)
    if len(genes) < 2:
        return G
    corr = expr.T.corr().abs()
    np.fill_diagonal(corr.values, 0.0)
    for gene in genes:
        partners = corr.loc[gene].sort_values(ascending=False).head(top_k)
        for p, c in partners.items():
            if c >= min_corr and gene != p:
                w = float(c)
                if G.has_edge(gene, p):
                    G[gene][p]["weight"] = max(G[gene][p]["weight"], w)
                else:
                    G.add_edge(gene, p, weight=w, edge_type="coexpression")
    return G


def graph_to_edge_list(G: nx.Graph) -> pd.DataFrame:
    return pd.DataFrame(
        [(u, v, d.get("weight", 1.0), d.get("edge_type", "unknown"))
         for u, v, d in G.edges(data=True)],
        columns=["source", "target", "weight", "edge_type"],
    )


def graph_to_adjacency(G: nx.Graph, genes: list[str]) -> pd.DataFrame:
    adj = nx.to_pandas_adjacency(G, nodelist=genes, weight="weight")
    return adj


def save_graph(G: nx.Graph, out_dir: Path, name: str, genes: list[str]) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    edge_df = graph_to_edge_list(G)
    adj_df = graph_to_adjacency(G, genes)
    edge_df.to_csv(out_dir / f"{name}_edges.tsv", sep="\t", index=False)
    adj_df.to_csv(out_dir / f"{name}_adjacency.tsv", sep="\t")


def summarize_graph(G: nx.Graph) -> dict:
    degs = dict(G.degree())
    n_nodes = G.number_of_nodes()
    n_edges = G.number_of_edges()
    isolates = [n for n, d in degs.items() if d == 0]
    density = nx.density(G)
    try:
        n_components = nx.number_connected_components(G)
    except Exception:
        n_components = None
    return {
        "n_nodes": n_nodes,
        "n_edges": n_edges,
        "density": density,
        "isolated_nodes": isolates,
        "n_isolates": len(isolates),
        "n_connected_components": n_components,
    }
