"""Generate a cartoon-cell version of the interactive Cu-proteome network.

The regular interactive plot is a force-directed graph. This version keeps the
same node/edge encodings but places genes over a simple cell drawing according
to their annotated subcellular localization.

Produces:
  outputs/visualizations/copper_cell_cartoon_network.html
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import networkx as nx
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.graph_building import build_functional_graph, graph_to_edge_list

OUT = ROOT / "outputs" / "visualizations"
OUT.mkdir(parents=True, exist_ok=True)

COPPER_CSV = ROOT / "outputs/paper_2017_extraction/copper_gene_list.csv"
DE_CSV = ROOT / "outputs/baseline/copper_de_results.csv"
IMP_CSV = ROOT / "outputs/gnn/node_importance.csv"
HISTONE_IMP_CSV = ROOT / "outputs/final_comparison/histone_node_importance.csv"
ATTN_CSV = ROOT / "outputs/gnn/top_attention_edges.csv"
HISTONE_ATTN_CSV = ROOT / "outputs/final_comparison/histone_top_attention.csv"
MODULES_CSV = ROOT / "outputs/baseline/copper_modules.csv"

AU_NODE = {
    "id": "AuCompounds",
    "label": "Au\ncompounds",
    "position": (1160, 400),
    "summary": "Extracellular Au(I)/Au(III) drug panel (auranofin, Aubipyc, Au(I)-phosphine-thioamide, Au(III) coumarin-phosphine). Shown with connecting lines to biomolecules with reported Au interaction.",
}

AU_EDGES = [
    # (target, confidence, mechanism, citation)
    ("ATOX1", "lit", "Au(I) forms 1:1/1:2 adduct at the Cu(I) CXXC site (Cys12-Cys15); displaces Cu(I)",
     "Gabbiani 2012 Chem Commun; Marzo 2015 Biometals"),
    ("ATP7A", "lit", "Cu-pathway efflux pump; Au(III) co-transported with Cu, implicated in Au efflux",
     "Spreckelmeyer 2018 Front Chem"),
    ("ATP7B", "lit", "Cu efflux pump; same efflux route in cisplatin-resistant A2780 cells",
     "Spreckelmeyer 2018 Front Chem"),
    ("LOX", "lit", "Au(I)-phosphine-thioamide inhibits lipoxygenase (IC50 1-39 uM). NOTE: 2011 paper tested lipoxygenase (ALOX class). Our graph LOX = lysyl oxidase (different Cu enzyme). Edge drawn provisionally pending clarification.",
     "Georgiou 2011 J Enzyme Inhib Med Chem"),
    ("H3-3A", "unpublished", "Histone complex - Au interaction suggested by collaborator's unpublished data",
     "Schwartz-Duval (unpublished)"),
    ("H3-3B", "unpublished", "Histone complex - Au interaction suggested by collaborator's unpublished data",
     "Schwartz-Duval (unpublished)"),
    ("H3C1", "unpublished", "Histone complex - Au interaction suggested by collaborator's unpublished data",
     "Schwartz-Duval (unpublished)"),
    ("H4C1", "unpublished", "Histone complex - Au interaction suggested by collaborator's unpublished data",
     "Schwartz-Duval (unpublished)"),
    ("PRNP", "unpublished", "Transcriptomics-derived candidate (significant FC in 2020 Nat Comm cohort); Au interaction pending validation",
     "Schwartz-Duval (unpublished) + Nat Comm 2020"),
    ("ENOX2", "unpublished", "Transcriptomics-derived candidate (significant FC in 2020 Nat Comm cohort); Au interaction pending validation",
     "Schwartz-Duval (unpublished) + Nat Comm 2020"),
]

CANONICAL_EDGES = {
    frozenset(("CCS", "SLC31A1")),
    frozenset(("SOD1", "CCS")),
    frozenset(("SOD3", "ATOX1")),
    frozenset(("ATP7B", "COMMD1")),
    frozenset(("MT-CO2", "SCO2")),
    frozenset(("ATP7B", "ATP7A")),
    frozenset(("ATOX1", "ATP7A")),
    frozenset(("ATOX1", "ATP7B")),
    frozenset(("COX17", "SCO1")),
    frozenset(("COX17", "SCO2")),
    frozenset(("COX17", "COX11")),
    frozenset(("H3-3A", "H4C1")),
    frozenset(("H3-3B", "H4C1")),
    frozenset(("H3C1", "H4C1")),
    frozenset(("H3-3A", "ATOX1")),
    frozenset(("H3-3B", "ATOX1")),
    frozenset(("H3C1", "ATOX1")),
    frozenset(("H3-3A", "SLC31A1")),
    frozenset(("H3-3A", "SOD1")),
    frozenset(("H3-3B", "SOD1")),
    frozenset(("H3-3A", "MT-CO1")),
    frozenset(("H3-3A", "MT-CO2")),
}

MANUAL_POS = {
    # Cell membrane / copper import + surface enzymes
    "SLC31A1": (620, 92),
    "SLC31A2": (735, 98),
    "AOC2": (940, 205),
    "AOC3": (975, 260),
    "HEPH": (1010, 155),
    "HEPHL1": (1085, 180),
    "APP": (1045, 575),
    "GPC1": (965, 515),
    "ATP7A": (820, 218),
    "ENOX1": (1075, 470),
    "ENOX2": (1075, 520),
    "PRNP": (855, 120),
    # Cytosolic copper-handling axis
    "ATOX1": (610, 265),
    "CCS": (565, 485),
    "SOD1": (680, 505),
    "SOD3": (735, 430),
    "COMMD1": (535, 345),
    "CUTC": (520, 415),
    "PARK7": (760, 540),
    "MAP2K1": (455, 315),
    "MEMO1": (405, 385),
    "SNCA": (625, 555),
    "MT3": (470, 445),
    # S100 family cluster (cytoplasm)
    "S100A5": (465, 475),
    "S100A12": (500, 450),
    "S100A13": (500, 485),
    "S100B": (530, 465),
    # Golgi / ER / vesicles
    "ATP7B": (405, 255),
    "MOXD1": (360, 235),
    "MOXD2P": (365, 300),
    "TYR": (430, 340),
    "TYRP1": (455, 335),
    # Nucleus — histone tetramer cluster (upper-right of nucleus, clear of nucleolus at (272,355) r=54)
    "H3-3A": (340, 240),
    "H3-3B": (370, 270),
    "H3C1":  (345, 295),
    "H4C1":  (375, 285),
    # Mitochondrial COX module
    "COX17": (800, 585),
    "COX11": (875, 605),
    "SCO1": (860, 665),
    "SCO2": (945, 650),
    "MT-CO1": (1015, 610),
    "MT-CO2": (1030, 690),
    "CUTA": (930, 550),
    # Secreted / ECM / plasma proteins outside the cell
    "ALB": (1000, 62),
    "CP": (1080, 82),
    "AFP": (930, 92),
    "F5": (1118, 110),
    "LTF": (1055, 130),
    "LOX": (110, 75),
    "LOXL1": (175, 72),
    "LOXL2": (240, 78),
    "LOXL3": (305, 82),
    "LOXL4": (370, 82),
    "MT4": (440, 58),
    "SPARC": (450, 110),
    "AOC1": (525, 90),
    "DBH": (585, 64),
    "PAM": (880, 74),
}

COMPARTMENT_CENTERS = {
    "cell_membrane": (900, 225),
    "golgi_apparatus": (390, 275),
    "endoplasmic_reticulum": (325, 260),
    "intracellular_vesicles": (250, 335),
    "mitochondrion": (915, 625),
    "nucleus": (245, 530),
    "cytoplasm": (575, 395),
    "cytoskeleton": (430, 505),
    "extracellular_space": (565, 88),
}


def first_localization(localization: str) -> str:
    locs = [x.strip() for x in str(localization).split("|") if x.strip()]
    priority = [
        "cell_membrane",
        "golgi_apparatus",
        "endoplasmic_reticulum",
        "intracellular_vesicles",
        "mitochondrion",
        "nucleus",
        "cytoskeleton",
        "cytoplasm",
        "extracellular_space",
    ]
    for loc in priority:
        if loc in locs:
            return loc
    return locs[0] if locs else "cytoplasm"


def spread_positions(genes: list[str], center: tuple[float, float], radius: float) -> dict[str, tuple[float, float]]:
    if not genes:
        return {}
    out = {}
    for i, gene in enumerate(sorted(genes)):
        angle = (2 * math.pi * i / max(len(genes), 1)) - math.pi / 2
        ring = radius * (0.55 + 0.45 * ((i % 3) / 2))
        out[gene] = (center[0] + ring * math.cos(angle), center[1] + ring * math.sin(angle))
    return out


def load_table(path: Path, index: str) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path).set_index(index)


def load_attention(path: Path) -> dict[frozenset[str], float]:
    if not path.exists():
        return {}
    df = pd.read_csv(path)
    if not {"source", "target", "attention_sum"}.issubset(df.columns):
        return {}
    df = df[df["source"] != df["target"]]
    out: dict[frozenset[str], float] = {}
    for row in df.itertuples(index=False):
        key = frozenset((str(row.source), str(row.target)))
        out[key] = max(out.get(key, 0.0), float(row.attention_sum))
    return out


def build_positions(copper: pd.DataFrame) -> dict[str, tuple[float, float]]:
    positions = dict(MANUAL_POS)
    used = set(positions)
    by_compartment: dict[str, list[str]] = {}
    for row in copper.itertuples(index=False):
        gene = str(row.gene_symbol)
        if gene in used:
            continue
        comp = first_localization(row.subcellular_localization)
        by_compartment.setdefault(comp, []).append(gene)
    for comp, genes in by_compartment.items():
        center = COMPARTMENT_CENTERS.get(comp, COMPARTMENT_CENTERS["cytoplasm"])
        radius = 65 if comp != "extracellular_space" else 120
        positions.update(spread_positions(genes, center, radius))
    return positions


def load_payload() -> dict:
    copper = pd.read_csv(COPPER_CSV)
    genes = copper["gene_symbol"].tolist()
    de = load_table(DE_CSV, "gene_symbol")
    imp_path = HISTONE_IMP_CSV if HISTONE_IMP_CSV.exists() else IMP_CSV
    imp = load_table(imp_path, "gene_symbol")
    modules = load_table(MODULES_CSV, "gene_symbol")
    attn_map = load_attention(HISTONE_ATTN_CSV if HISTONE_ATTN_CSV.exists() else ATTN_CSV)
    top_attn = set(sorted(attn_map, key=attn_map.get, reverse=True)[:18])

    graph = build_functional_graph(genes, copper)
    edges = graph_to_edge_list(graph)
    positions = build_positions(copper)

    imp_max = float(imp["importance"].max()) if "importance" in imp and len(imp) else 1.0
    top_imp = set(imp.sort_values("importance", ascending=False).head(12).index.tolist()) if len(imp) else set()

    nodes = []
    for row in copper.itertuples(index=False):
        gene = str(row.gene_symbol)
        de_row = de.loc[gene] if gene in de.index else None
        imp_val = float(imp.loc[gene, "importance"]) if gene in imp.index else 0.0
        mod_id = int(modules.loc[gene, "module_id"]) if gene in modules.index else -1
        x, y = positions[gene]
        nodes.append({
            "data": {
                "id": gene,
                "label": gene,
                "category": row.functional_category,
                "localization": row.subcellular_localization,
                "compartment": first_localization(row.subcellular_localization),
                "notes": row.notes_from_paper,
                "log2FC": float(de_row["log2FC"]) if de_row is not None else 0.0,
                "adj_p": float(de_row["adj_p_BH"]) if de_row is not None else 1.0,
                "significant": bool(de_row["significant_0.05_BH"]) if de_row is not None else False,
                "importance": imp_val,
                "importance_norm": imp_val / imp_max if imp_max else 0.0,
                "module": mod_id,
                "is_top_importance": gene in top_imp,
            },
            "position": {"x": x, "y": y},
        })

    edge_list = []
    for row in edges.itertuples(index=False):
        key = frozenset((row.source, row.target))
        edge_list.append({
            "data": {
                "id": f"{row.source}__{row.target}",
                "source": row.source,
                "target": row.target,
                "weight": float(row.weight),
                "edge_type": row.edge_type,
                "attention": float(attn_map.get(key, 0.0)),
                "is_top_attention": key in top_attn,
                "is_canonical": key in CANONICAL_EDGES,
            }
        })

    # Append the Au compounds node and its connecting edges. The Au node is
    # conceptually separate from the expression-level gene nodes: it is a drug
    # perturbation, not an mRNA/protein; it carries no logFC or saliency.
    nodes.append({
        "data": {
            "id": AU_NODE["id"],
            "label": AU_NODE["label"],
            "category": "au_compound",
            "localization": "extracellular_space",
            "compartment": "extracellular_space",
            "notes": AU_NODE["summary"],
            "log2FC": 0.0,
            "adj_p": 1.0,
            "significant": False,
            "importance": 0.0,
            "importance_norm": 0.0,
            "module": -1,
            "is_top_importance": False,
        },
        "position": {"x": AU_NODE["position"][0], "y": AU_NODE["position"][1]},
    })

    present_genes = {n["data"]["id"] for n in nodes}
    au_edges_added = 0
    for target, confidence, mechanism, citation in AU_EDGES:
        if target not in present_genes:
            continue
        edge_list.append({
            "data": {
                "id": f"{AU_NODE['id']}__{target}",
                "source": AU_NODE["id"],
                "target": target,
                "weight": 1.0,
                "edge_type": "au_interaction",
                "attention": 0.0,
                "is_top_attention": False,
                "is_canonical": False,
                "au_confidence": confidence,
                "au_mechanism": mechanism,
                "au_citation": citation,
            }
        })
        au_edges_added += 1

    attn_values = [e["data"]["attention"] for e in edge_list if e["data"]["attention"] > 0]
    return {
        "nodes": nodes,
        "edges": edge_list,
        "meta": {
            "attn_max": max(attn_values) if attn_values else 1.0,
            "imp_max": imp_max,
            "n_nodes": len(nodes),
            "n_edges": len(edge_list),
            "n_au_edges": au_edges_added,
        },
    }


HTML = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>TCGA-LIHC Copper Proteome Cartoon Cell Network</title>
<script src="https://unpkg.com/cytoscape@3.28.1/dist/cytoscape.min.js"></script>
<style>
  html, body { margin: 0; height: 100%; font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; color: #1f2933; background: #f6f8fb; }
  #app { display: grid; grid-template-columns: 288px 1fr; grid-template-rows: 54px 1fr; height: 100vh; }
  #header { grid-column: 1 / 3; display: flex; align-items: center; justify-content: space-between; padding: 0 18px; background: #243b53; color: #fff; box-shadow: 0 1px 5px rgba(0,0,0,.18); }
  #header h1 { margin: 0; font-size: 16px; letter-spacing: 0; font-weight: 650; }
  #header .subtitle { font-size: 12px; opacity: .84; margin-top: 2px; }
  #sidebar { background: #fff; border-right: 1px solid #d9e2ec; overflow: auto; padding: 14px; }
  #stage { position: relative; overflow: hidden; background: #f8fbff; }
  #cy { position: absolute; inset: 0; z-index: 2; background: transparent; }
  #cell-bg { position: absolute; inset: 0; z-index: 1; pointer-events: none; }
  #tooltip { position: absolute; z-index: 5; display: none; pointer-events: none; max-width: 330px; padding: 10px 12px; border-radius: 7px; background: rgba(19, 32, 44, .96); color: white; font-size: 12px; line-height: 1.42; box-shadow: 0 8px 24px rgba(16,24,40,.28); }
  .tip-title { font-size: 13px; font-weight: 700; margin-bottom: 5px; }
  .tip-row { display: flex; justify-content: space-between; gap: 12px; }
  .tip-label { color: #b8c5d3; }
  .tip-value { font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace; color: #fff; text-align: right; }
  .up { color: #ffb4a8; }
  .down { color: #a8d8ff; }
  .panel { margin-bottom: 15px; }
  .panel h3 { margin: 0 0 8px; padding-bottom: 5px; border-bottom: 1px solid #e6edf5; font-size: 12px; text-transform: uppercase; letter-spacing: .06em; color: #52606d; }
  .control-row { display: flex; align-items: center; gap: 8px; margin: 7px 0; font-size: 13px; }
  .control-row input[type=search] { width: 100%; border: 1px solid #cbd5e1; border-radius: 5px; padding: 7px 8px; font-size: 13px; }
  button { border: 1px solid #bcccdc; border-radius: 5px; background: #f0f4f8; padding: 6px 9px; cursor: pointer; font-size: 12px; color: #243b53; }
  button:hover { background: #e6edf5; }
  label { cursor: pointer; }
  .legend-item { display: flex; align-items: center; gap: 8px; margin: 6px 0; font-size: 12px; }
  .swatch { width: 18px; height: 18px; border: 1px solid #9fb3c8; background: #fff; flex: 0 0 auto; }
  .round { border-radius: 50%; }
  .edge-swatch { width: 34px; height: 3px; flex: 0 0 auto; }
  .small { font-size: 11px; color: #627d98; line-height: 1.45; }
</style>
</head>
<body>
<div id="app">
  <div id="header">
    <div>
      <h1>TCGA-LIHC Copper Proteome — Cartoon Cell Network</h1>
      <div class="subtitle">same graph data, arranged by cellular compartment · zoom/pan/drag freely, click Reset to restore layout</div>
    </div>
  </div>

  <aside id="sidebar">
    <div class="panel">
      <h3>Search</h3>
      <div class="control-row"><input id="search" type="search" placeholder="gene symbol, e.g. ATOX1"></div>
      <div class="control-row"><button id="fit">Fit</button><button id="reset">Reset</button><button id="lock">Lock nodes</button></div>
      <div class="control-row">
        <button id="save_layout" title="Download current node positions as a JSON file">Save layout</button>
        <button id="load_layout" title="Load node positions from a previously saved JSON file">Load layout</button>
        <input type="file" id="load_file" accept="application/json,.json" style="display:none">
      </div>
      <div class="small">Drag nodes to tidy up, then Save layout to keep them. Load layout restores a saved arrangement.</div>
    </div>

    <div class="panel">
      <h3>Map</h3>
      <div class="legend-item"><span class="swatch round" style="background:#dbeafe"></span>cell membrane / surface</div>
      <div class="legend-item"><span class="swatch round" style="background:#efe0ff"></span>nucleus / histones</div>
      <div class="legend-item"><span class="swatch round" style="background:#d9f99d"></span>mitochondrion</div>
      <div class="legend-item"><span class="swatch round" style="background:#bfdbfe"></span>ER / Golgi / vesicles</div>
      <div class="legend-item"><span class="swatch round" style="background:#fef3c7"></span>extracellular / secreted</div>
    </div>

    <div class="panel">
      <h3>Node Encoding</h3>
      <div class="legend-item"><span class="swatch round" style="background:#b22222"></span>red = up in tumor</div>
      <div class="legend-item"><span class="swatch round" style="background:#f5f5f5"></span>white = unchanged / missing</div>
      <div class="legend-item"><span class="swatch round" style="background:#1f4e79"></span>blue = down in tumor</div>
      <div class="small">Size = GNN saliency. Black outline = top-ranked node.</div>
    </div>

    <div class="panel">
      <h3>Edges</h3>
      <div class="control-row"><input type="checkbox" id="show_physical" checked><label for="show_physical">physical</label><span class="edge-swatch" style="background:#cf3d35"></span></div>
      <div class="control-row"><input type="checkbox" id="show_coexpression" checked><label for="show_coexpression">co-expression</label><span class="edge-swatch" style="background:#8e44ad"></span></div>
      <div class="control-row"><input type="checkbox" id="show_genetic" checked><label for="show_genetic">genetic</label><span class="edge-swatch" style="background:#2f855a"></span></div>
      <div class="control-row"><input type="checkbox" id="show_shared_compartment" checked><label for="show_shared_compartment">shared compartment</label><span class="edge-swatch" style="background:#94a3b8"></span></div>
      <div class="control-row"><input type="checkbox" id="show_au_interaction" checked><label for="show_au_interaction">Au interaction</label><span class="edge-swatch" style="background:#8b6914"></span></div>
      <div class="small">Width = GAT attention when available. Canonical copper-biology edges are drawn darker. Au edges: solid = published lit, dashed = collaborator's unpublished data.</div>
    </div>

    <div class="panel">
      <h3>Au compounds</h3>
      <div class="legend-item"><span class="swatch" style="background:#e0b324;border-color:#8b6914;transform:rotate(45deg);width:14px;height:14px"></span>extracellular Au drug panel</div>
      <div class="small">Targets drawn with solid lines are supported by the 2011 / 2012 / 2015 / 2018 references. Dashed lines (histones, PRNP, ENOX2) are provisional pending collaborator's unpublished data.</div>
    </div>

    <div class="panel">
      <h3>Stats</h3>
      <div id="stats" class="small"></div>
    </div>
  </aside>

  <main id="stage">
    <svg id="cell-bg" viewBox="0 0 1200 760" preserveAspectRatio="xMidYMid meet" aria-hidden="true">
      <defs>
        <filter id="soft"><feGaussianBlur stdDeviation="1.1"/></filter>
      </defs>
      <rect x="0" y="0" width="1200" height="760" fill="#f8fbff"/>
      <text x="575" y="36" text-anchor="middle" fill="#6b7280" font-size="13" font-weight="650">extracellular / secreted copper proteins</text>
      <path d="M69,151 C143,58 335,94 467,75 C593,57 776,53 967,80 C1114,101 1160,193 1148,337 C1135,503 1058,702 840,711 C652,719 550,662 381,662 C229,662 57,622 45,448 C35,304 25,205 69,151 Z"
            fill="#fffef7" stroke="#d8be45" stroke-width="15" opacity=".72"/>
      <path d="M84,162 C158,75 342,109 471,91 C594,74 773,70 952,94 C1091,113 1139,199 1129,334 C1116,488 1041,681 839,688 C655,694 555,642 383,642 C236,642 79,607 67,443 C57,309 42,214 84,162 Z"
            fill="none" stroke="#b59b2f" stroke-width="3" opacity=".48"/>
      <path d="M146,335 C151,255 207,210 284,218 C367,226 400,280 387,359 C374,440 315,481 236,463 C176,449 139,405 146,335 Z"
            fill="#efe0ff" stroke="#ad7dd1" stroke-width="10" opacity=".72"/>
      <circle cx="272" cy="355" r="54" fill="#854bbd" opacity=".72"/>
      <path d="M408,232 C339,213 336,286 421,285 C477,284 474,222 408,232 Z M397,300 C347,296 338,342 389,354 C451,369 466,305 397,300 Z"
            fill="#bfdbfe" stroke="#77a9d6" stroke-width="5" opacity=".68"/>
      <text x="300" y="206" text-anchor="middle" fill="#6d48a1" font-size="13" font-weight="650">nucleus</text>
      <text x="398" y="213" text-anchor="middle" fill="#426b9a" font-size="12" font-weight="650">ER / Golgi</text>
      <path d="M844,503 C903,451 1002,463 1050,527 C1092,584 1052,668 970,691 C895,711 804,681 779,619 C762,577 782,539 844,503 Z"
            fill="#d9f99d" stroke="#3f7a35" stroke-width="9" opacity=".8"/>
      <path d="M837,584 C878,549 907,592 934,557 C956,528 981,599 1010,565 C1033,537 1054,592 1019,617 C983,643 965,601 941,633 C911,672 873,623 837,655 C808,681 794,621 837,584 Z"
            fill="none" stroke="#315c2d" stroke-width="5" opacity=".85"/>
      <text x="931" y="490" text-anchor="middle" fill="#315c2d" font-size="13" font-weight="650">mitochondrion</text>
      <text x="602" y="400" text-anchor="middle" fill="#64748b" font-size="13" font-weight="650">cytoplasm</text>
      <text x="1000" y="231" text-anchor="middle" fill="#64748b" font-size="12" font-weight="650">membrane</text>
    </svg>
    <div id="cy"></div>
    <div id="tooltip"></div>
  </main>
</div>

<script>
const DATA = __DATA_JSON__;
const META = __META_JSON__;
const EDGE_COLOR = {
  physical: "#cf3d35",
  coexpression: "#8e44ad",
  genetic: "#2f855a",
  shared_compartment: "#94a3b8",
};
const NODE_SHAPE = { transporter: "ellipse", enzyme: "round-rectangle", other_or_unknown: "triangle", au_compound: "diamond" };
const AU_GOLD = "#e0b324";
const AU_GOLD_DARK = "#8b6914";

function lerp(a,b,t){ return a + (b-a)*t; }
function lerpColor(a,b,t){ return [Math.round(lerp(a[0],b[0],t)), Math.round(lerp(a[1],b[1],t)), Math.round(lerp(a[2],b[2],t))]; }
function logfcColor(fc) {
  const blue=[31,78,121], white=[245,245,245], red=[178,34,34];
  const clipped = Math.max(-3, Math.min(3, fc || 0));
  const rgb = clipped >= 0 ? lerpColor(white, red, clipped/3) : lerpColor(white, blue, -clipped/3);
  return `rgb(${rgb[0]},${rgb[1]},${rgb[2]})`;
}
function nodeSize(v){ return 20 + 42 * (v || 0); }
function edgeWidth(attn){ return attn > 0 ? 1.6 + 6.5 * (attn / META.attn_max) : 1.35; }

DATA.nodes.forEach(n => {
  const d = n.data;
  if (d.category === "au_compound") {
    d.color = AU_GOLD;
    d.size = 64;
    d.shape = "diamond";
    d.border_width = 3.2;
    d.border_color = AU_GOLD_DARK;
  } else {
    d.color = logfcColor(d.log2FC);
    d.size = nodeSize(d.importance_norm);
    d.shape = NODE_SHAPE[d.category] || "ellipse";
    d.border_width = d.is_top_importance ? 3.4 : 1.1;
    d.border_color = d.is_top_importance ? "#111827" : "#64748b";
  }
});
DATA.edges.forEach(e => {
  const d = e.data;
  if (d.edge_type === "au_interaction") {
    d.color = AU_GOLD_DARK;
    d.width = d.au_confidence === "lit" ? 3.2 : 2.2;
    d.opacity = 0.85;
  } else {
    d.color = d.is_canonical ? "#7f1d1d" : (EDGE_COLOR[d.edge_type] || "#94a3b8");
    d.width = d.is_top_attention ? Math.max(3.2, edgeWidth(d.attention)) : edgeWidth(d.attention);
    d.opacity = d.edge_type === "shared_compartment" ? .42 : .72;
  }
});

const HOME_POS = Object.fromEntries(DATA.nodes.map(n => [n.data.id, {x: n.position.x, y: n.position.y}]));

const cy = cytoscape({
  container: document.getElementById("cy"),
  elements: [...DATA.nodes, ...DATA.edges],
  layout: { name: "preset", fit: true, padding: 35 },
  wheelSensitivity: 0.22,
  minZoom: 0.45,
  maxZoom: 3.2,
  style: [
    { selector: "node", style: {
      "background-color": "data(color)",
      "shape": "data(shape)",
      "width": "data(size)",
      "height": "data(size)",
      "border-color": "data(border_color)",
      "border-width": "data(border_width)",
      "label": "data(label)",
      "font-size": 10,
      "font-weight": 650,
      "text-valign": "center",
      "text-halign": "center",
      "color": "#111827",
      "text-outline-color": "#ffffff",
      "text-outline-width": 1.7,
    }},
    { selector: "edge", style: {
      "width": "data(width)",
      "line-color": "data(color)",
      "opacity": "data(opacity)",
      "curve-style": "bezier",
    }},
    { selector: "edge[edge_type = 'shared_compartment']", style: { "line-style": "dashed" }},
    { selector: "edge[edge_type = 'au_interaction'][au_confidence = 'unpublished']", style: { "line-style": "dashed" }},
    { selector: "node[category = 'au_compound']", style: {
      "font-size": 13,
      "font-weight": 800,
      "text-wrap": "wrap",
      "color": "#3d2b00",
      "text-outline-color": "#fff8dc",
      "text-outline-width": 2.2,
    }},
    { selector: ".faded", style: { "opacity": .12 }},
    { selector: "node.selected", style: { "border-width": 5, "border-color": "#0f172a" }},
    { selector: "edge.selected", style: { "line-color": "#0f172a", "opacity": 1 }},
  ],
});

const tip = document.getElementById("tooltip");
let pinned = false;
let locked = false;

// Sync SVG cartoon background with cytoscape pan/zoom so the cartoon scales with the graph.
const cellBg = document.getElementById("cell-bg");
cellBg.style.transformOrigin = "0 0";
cellBg.style.willChange = "transform";
let baseZoom = null, basePan = {x: 0, y: 0};
function syncSvg() {
  if (baseZoom === null) return;
  const s = cy.zoom() / baseZoom;
  const p = cy.pan();
  const tx = p.x - basePan.x * s;
  const ty = p.y - basePan.y * s;
  cellBg.style.transform = `translate(${tx}px, ${ty}px) scale(${s})`;
}
cy.on("viewport", syncSvg);

function fmt(x, digits=3) {
  if (x === null || x === undefined || Number.isNaN(x)) return "-";
  if (Math.abs(x) < 0.001 && x !== 0) return Number(x).toExponential(2);
  return Number(x).toFixed(digits);
}
function nodeTip(n) {
  const d = n.data();
  if (d.category === "au_compound") {
    return `<div class="tip-title">Au compounds (extracellular)</div>
      <div style="margin-top:4px;color:#ffe6a0;font-size:11px">Drug panel: auranofin, Aubipyc, Au(I)-phosphine-thioamide, Au(III) coumarin-phosphine.</div>
      <div style="margin-top:6px;color:#d8e2ee;font-size:11px">Solid lines = supported by published literature. Dashed lines = collaborator's unpublished data (pending shared references).</div>`;
  }
  const fcClass = d.log2FC > 0 ? "up" : (d.log2FC < 0 ? "down" : "");
  return `<div class="tip-title">${d.label}</div>
    <div class="tip-row"><span class="tip-label">compartment</span><span class="tip-value">${d.compartment}</span></div>
    <div class="tip-row"><span class="tip-label">log2FC</span><span class="tip-value ${fcClass}">${fmt(d.log2FC)}</span></div>
    <div class="tip-row"><span class="tip-label">adj p</span><span class="tip-value">${fmt(d.adj_p)}</span></div>
    <div class="tip-row"><span class="tip-label">GNN importance</span><span class="tip-value">${fmt(d.importance)}</span></div>
    <div class="tip-row"><span class="tip-label">category</span><span class="tip-value">${d.category}</span></div>
    <div style="margin-top:6px;color:#d8e2ee;font-size:11px">${d.notes || ""}</div>`;
}
function edgeTip(e) {
  const d = e.data();
  if (d.edge_type === "au_interaction") {
    const confLabel = d.au_confidence === "lit" ? "published literature" : "unpublished (collaborator data)";
    return `<div class="tip-title">Au compounds - ${d.target}</div>
      <div class="tip-row"><span class="tip-label">evidence</span><span class="tip-value">${confLabel}</span></div>
      <div class="tip-row"><span class="tip-label">citation</span><span class="tip-value" style="white-space:normal;text-align:right">${d.au_citation || "-"}</span></div>
      <div style="margin-top:6px;color:#ffe6a0;font-size:11px">${d.au_mechanism || ""}</div>`;
  }
  const canon = d.is_canonical ? " · canonical Cu edge" : "";
  return `<div class="tip-title">${d.source} - ${d.target}${canon}</div>
    <div class="tip-row"><span class="tip-label">type</span><span class="tip-value">${String(d.edge_type).replaceAll("_"," ")}</span></div>
    <div class="tip-row"><span class="tip-label">attention</span><span class="tip-value">${fmt(d.attention)}</span></div>
    <div class="tip-row"><span class="tip-label">prior weight</span><span class="tip-value">${fmt(d.weight)}</span></div>`;
}
function showTip(content, evt) {
  tip.innerHTML = content;
  tip.style.left = (evt.clientX + 14) + "px";
  tip.style.top = (evt.clientY + 14) + "px";
  tip.style.display = "block";
}
function hideTip(){ if (!pinned) tip.style.display = "none"; }

cy.on("mouseover", "node", e => { if (!pinned) showTip(nodeTip(e.target), e.originalEvent); });
cy.on("mouseover", "edge", e => { if (!pinned) showTip(edgeTip(e.target), e.originalEvent); });
cy.on("mousemove", "node, edge", e => { if (!pinned) { tip.style.left = (e.originalEvent.clientX + 14) + "px"; tip.style.top = (e.originalEvent.clientY + 14) + "px"; }});
cy.on("mouseout", "node, edge", hideTip);
cy.on("tap", "node", e => {
  pinned = true;
  cy.elements().removeClass("selected faded");
  const hood = e.target.closedNeighborhood();
  cy.elements().not(hood).addClass("faded");
  e.target.addClass("selected");
  e.target.connectedEdges().addClass("selected");
  showTip(nodeTip(e.target), e.originalEvent);
});
cy.on("tap", "edge", e => {
  pinned = true;
  cy.elements().removeClass("selected faded");
  const keep = e.target.connectedNodes().union(e.target);
  cy.elements().not(keep).addClass("faded");
  e.target.addClass("selected");
  showTip(edgeTip(e.target), e.originalEvent);
});
cy.on("tap", e => {
  if (e.target === cy) {
    pinned = false;
    tip.style.display = "none";
    cy.elements().removeClass("selected faded");
  }
});

document.getElementById("search").addEventListener("input", e => {
  const q = e.target.value.trim().toUpperCase();
  cy.elements().removeClass("faded selected");
  if (!q) return;
  cy.nodes().forEach(n => {
    const hit = n.id().toUpperCase().includes(q);
    if (!hit) n.addClass("faded");
    else n.addClass("selected");
  });
});
document.getElementById("fit").addEventListener("click", () => cy.fit(undefined, 35));
document.getElementById("reset").addEventListener("click", () => {
  pinned = false;
  tip.style.display = "none";
  document.getElementById("search").value = "";
  cy.elements().removeClass("selected faded");
  // restore cartoon positions even if user has dragged nodes around
  cy.nodes().forEach(n => {
    const home = HOME_POS[n.id()];
    if (home) n.position(home);
  });
  cy.fit(undefined, 35);
});
document.getElementById("lock").addEventListener("click", e => {
  locked = !locked;
  cy.nodes().forEach(n => locked ? n.lock() : n.unlock());
  e.target.textContent = locked ? "Unlock nodes" : "Lock nodes";
});
["physical", "coexpression", "genetic", "shared_compartment", "au_interaction"].forEach(t => {
  document.getElementById("show_" + t).addEventListener("change", e => {
    cy.edges().forEach(edge => {
      if (edge.data("edge_type") === t) edge.style("display", e.target.checked ? "element" : "none");
    });
  });
});

// Save layout: write the current node positions to a JSON file the user downloads.
document.getElementById("save_layout").addEventListener("click", () => {
  const positions = {};
  cy.nodes().forEach(n => {
    const p = n.position();
    positions[n.id()] = { x: Math.round(p.x * 100) / 100, y: Math.round(p.y * 100) / 100 };
  });
  const payload = {
    version: 1,
    saved_at: new Date().toISOString(),
    n_nodes: cy.nodes().length,
    positions: positions,
  };
  const blob = new Blob([JSON.stringify(payload, null, 2)], { type: "application/json" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = `cartoon_layout_${new Date().toISOString().replace(/[:.]/g, "-").slice(0, 19)}.json`;
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  URL.revokeObjectURL(url);
});

// Load layout: read a JSON file and apply its positions to the current graph.
document.getElementById("load_layout").addEventListener("click", () => {
  document.getElementById("load_file").click();
});
document.getElementById("load_file").addEventListener("change", e => {
  const file = e.target.files[0];
  if (!file) return;
  const reader = new FileReader();
  reader.onload = () => {
    try {
      const parsed = JSON.parse(reader.result);
      const pos = parsed && parsed.positions;
      if (!pos || typeof pos !== "object") {
        alert("Layout file is missing a 'positions' field.");
        return;
      }
      let applied = 0;
      cy.nodes().forEach(n => {
        const p = pos[n.id()];
        if (p && typeof p.x === "number" && typeof p.y === "number") {
          n.position({ x: p.x, y: p.y });
          applied += 1;
        }
      });
      cy.fit(undefined, 35);
      requestAnimationFrame(() => {
        baseZoom = cy.zoom();
        basePan = { ...cy.pan() };
        syncSvg();
      });
      alert(`Loaded ${applied} / ${cy.nodes().length} node positions.`);
    } catch (err) {
      alert("Could not parse layout file: " + err.message);
    }
  };
  reader.readAsText(file);
  e.target.value = "";
});

document.getElementById("stats").innerHTML =
  `nodes: ${META.n_nodes}<br>edges: ${META.n_edges}<br>` +
  `largest saliency: ${fmt(META.imp_max)}<br>` +
  `layout: cartoon (drag freely, Reset to restore)`;

cy.ready(() => {
  cy.fit(undefined, 35);
  // Record cy viewport AFTER fit so we can scale/pan the SVG cartoon to match.
  requestAnimationFrame(() => {
    baseZoom = cy.zoom();
    basePan = { ...cy.pan() };
    syncSvg();
  });
});
window.addEventListener("resize", () => {
  // On resize, refit and re-baseline so SVG stays aligned.
  cy.fit(undefined, 35);
  requestAnimationFrame(() => {
    baseZoom = cy.zoom();
    basePan = { ...cy.pan() };
    cellBg.style.transform = "none";
  });
});
</script>
</body>
</html>
"""


def main() -> None:
    print("[cell-cartoon] loading data")
    payload = load_payload()
    html = (
        HTML
        .replace("__DATA_JSON__", json.dumps({"nodes": payload["nodes"], "edges": payload["edges"]}))
        .replace("__META_JSON__", json.dumps(payload["meta"]))
    )
    out_html = OUT / "copper_cell_cartoon_network.html"
    out_html.write_text(html)
    print(f"[cell-cartoon] wrote {out_html} ({out_html.stat().st_size / 1024:.1f} KB)")

    notes = OUT / "copper_cell_cartoon_network_notes.md"
    notes.write_text("""# Cartoon Cell Copper Network

Open `copper_cell_cartoon_network.html` in a browser.

Compartment-mapped version of the Cu-proteome graph:

- Cell-membrane importer / surface proteins are placed along the outer membrane.
- ATOX1 / CCS / SOD1 / COMMD1 / CUTC sit in the cytosol.
- ATP7A / ATP7B and ER/Golgi/vesicle proteins sit near the ER/Golgi drawing.
- COX / SCO / MT-CO genes sit over the mitochondrion.
- H3 / H4 histone nodes sit in the nucleus (upper-right, clear of the nucleolus).
- Secreted / ECM / plasma proteins sit outside the cell.
- An extracellular gold-diamond **Au compounds** node is connected to published
  Au targets (ATOX1, ATP7A, ATP7B, LOX) with solid lines and to provisional
  collaborator-data targets (histones, PRNP, ENOX2) with dashed lines.

Visual encodings:

- Node color = LIHC tumor-vs-normal log2FC (red up, blue down).
- Node size = GNN saliency / importance.
- Edge color = edge type; Au edges are gold.
- Edge width = GAT attention when available (Au edges: solid = literature,
  dashed = unpublished collaborator data).

Controls:

- The cartoon SVG scales and pans with the graph (mouse wheel to zoom).
- Nodes are draggable — move them around to tidy overlapping neighbourhoods.
- `Save layout` downloads a JSON file of the current node positions.
- `Load layout` restores a previously saved JSON file.
- `Reset` restores the built-in cartoon positions and re-fits the view.
- `Lock nodes` freezes positions if you want to export a figure panel.
""")
    print(f"[cell-cartoon] wrote {notes}")


if __name__ == "__main__":
    main()
