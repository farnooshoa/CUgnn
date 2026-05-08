"""Generate an interactive draggable HTML network (Cytoscape.js) for the LIHC Cu graph.

Produces:
  outputs/visualizations/copper_interactive_network.html
"""
from __future__ import annotations
import json
import sys
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

OUT = ROOT / "outputs" / "visualizations"
OUT.mkdir(parents=True, exist_ok=True)

COPPER_CSV = ROOT / "outputs/paper_2017_extraction/copper_gene_list.csv"
DE_CSV = ROOT / "outputs/baseline/copper_de_results.csv"
IMP_CSV = ROOT / "outputs/gnn/node_importance.csv"
ATTN_CSV = ROOT / "outputs/gnn/top_attention_edges.csv"
EDGES_TSV = ROOT / "outputs/baseline/graphs/functional_graph_edges.tsv"
MODULES_CSV = ROOT / "outputs/baseline/copper_modules.csv"

BIO_EDGES = {
    frozenset(("CCS", "SLC31A1")),
    frozenset(("SOD3", "ATOX1")),
    frozenset(("ATP7B", "COMMD1")),
    frozenset(("MT-CO2", "SCO2")),
    frozenset(("SOD1", "CCS")),
    frozenset(("ATP7B", "ATP7A")),
    frozenset(("ATOX1", "ATP7A")),
    frozenset(("ATOX1", "ATP7B")),
    frozenset(("COX17", "SCO1")),
    frozenset(("COX17", "SCO2")),
    frozenset(("COX17", "COX11")),
}


def load_all() -> dict:
    copper = pd.read_csv(COPPER_CSV)
    de = pd.read_csv(DE_CSV).set_index("gene_symbol")
    imp = pd.read_csv(IMP_CSV).set_index("gene_symbol")
    modules = pd.read_csv(MODULES_CSV).set_index("gene_symbol")
    edges = pd.read_csv(EDGES_TSV, sep="\t")

    attn = pd.DataFrame(columns=["source", "target", "attention_sum"])
    if ATTN_CSV.exists():
        attn = pd.read_csv(ATTN_CSV)
        attn = attn[attn["source"] != attn["target"]]
    attn_map = {}
    for r in attn.itertuples(index=False):
        key = frozenset((r.source, r.target))
        attn_map[key] = max(attn_map.get(key, 0.0), float(r.attention_sum))

    top_imp = set(imp.sort_values("importance", ascending=False).head(10).index.tolist())
    top_attn_edges = set()
    if attn_map:
        sorted_keys = sorted(attn_map, key=attn_map.get, reverse=True)
        top_attn_edges = set(sorted_keys[:15])

    nodes = []
    imp_max = float(imp["importance"].max()) if len(imp) else 1.0
    for row in copper.itertuples(index=False):
        g = row.gene_symbol
        de_row = de.loc[g] if g in de.index else None
        imp_val = float(imp.loc[g, "importance"]) if g in imp.index else 0.0
        mod_id = int(modules.loc[g, "module_id"]) if g in modules.index else -1
        node = {
            "data": {
                "id": g,
                "label": g,
                "category": row.functional_category,
                "localization": row.subcellular_localization,
                "notes": row.notes_from_paper,
                "log2FC": float(de_row["log2FC"]) if de_row is not None else 0.0,
                "adj_p": float(de_row["adj_p_BH"]) if de_row is not None else 1.0,
                "p_value": float(de_row["p_value"]) if de_row is not None else 1.0,
                "significant": bool(de_row["significant_0.05_BH"]) if de_row is not None else False,
                "importance": imp_val,
                "importance_norm": imp_val / imp_max if imp_max else 0.0,
                "module": mod_id,
                "is_top_importance": g in top_imp,
            }
        }
        nodes.append(node)

    edge_list = []
    seen = set()
    for row in edges.itertuples(index=False):
        key = frozenset((row.source, row.target))
        if key in seen:
            continue
        seen.add(key)
        attn_val = float(attn_map.get(key, 0.0))
        is_bio_canonical = key in BIO_EDGES
        is_top_attn = key in top_attn_edges
        edge_list.append({
            "data": {
                "id": f"{row.source}__{row.target}",
                "source": row.source,
                "target": row.target,
                "weight": float(row.weight) if hasattr(row, "weight") else 1.0,
                "edge_type": row.edge_type,
                "attention": attn_val,
                "is_bio_canonical": is_bio_canonical,
                "is_top_attention": is_top_attn,
            }
        })

    if attn_map:
        existing = seen
        for key, val in sorted(attn_map.items(), key=lambda kv: kv[1], reverse=True)[:15]:
            if key in existing:
                continue
            a, b = tuple(key)
            edge_list.append({
                "data": {
                    "id": f"{a}__{b}",
                    "source": a, "target": b,
                    "weight": 0.0,
                    "edge_type": "attention_only",
                    "attention": float(val),
                    "is_bio_canonical": key in BIO_EDGES,
                    "is_top_attention": True,
                }
            })

    attn_values = [e["data"]["attention"] for e in edge_list if e["data"]["attention"] > 0]
    attn_max = max(attn_values) if attn_values else 1.0

    return {
        "nodes": nodes,
        "edges": edge_list,
        "meta": {
            "attn_max": attn_max,
            "imp_max": imp_max,
            "n_top_importance": len(top_imp),
            "n_top_attention": len(top_attn_edges),
        },
    }


HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>TCGA-LIHC Copper Proteome Interactive Network</title>
<script src="https://unpkg.com/cytoscape@3.28.1/dist/cytoscape.min.js"></script>
<script src="https://unpkg.com/layout-base@2.0.1/layout-base.js"></script>
<script src="https://unpkg.com/cose-base@2.2.0/cose-base.js"></script>
<script src="https://unpkg.com/cytoscape-fcose@2.2.0/cytoscape-fcose.js"></script>
<script src="https://unpkg.com/@popperjs/core@2.11.8/dist/umd/popper.min.js"></script>
<script src="https://unpkg.com/cytoscape-popper@2.0.0/cytoscape-popper.js"></script>
<style>
  html, body { margin: 0; padding: 0; font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; background: #fafafa; }
  #app { display: grid; grid-template-columns: 260px 1fr; grid-template-rows: 48px 1fr; height: 100vh; }
  #header { grid-column: 1 / 3; display: flex; align-items: center; justify-content: space-between;
             padding: 0 16px; background: #1a237e; color: white; box-shadow: 0 1px 4px rgba(0,0,0,0.1); }
  #header h1 { margin: 0; font-size: 16px; font-weight: 500; }
  #header .subtitle { font-size: 12px; opacity: 0.8; }
  #sidebar { border-right: 1px solid #e0e0e0; overflow-y: auto; padding: 12px; background: white; }
  #cy { width: 100%; height: 100%; background: #fafafa; position: relative; }

  .panel { margin-bottom: 16px; }
  .panel h3 { font-size: 12px; text-transform: uppercase; letter-spacing: 0.5px;
              color: #555; margin: 0 0 8px; border-bottom: 1px solid #eee; padding-bottom: 4px; }
  .control-row { display: flex; align-items: center; gap: 8px; margin-bottom: 6px; font-size: 13px; }
  .control-row input[type=search] { flex: 1; padding: 4px 6px; border: 1px solid #ddd; border-radius: 3px; font-size: 12px; }
  .control-row button { padding: 4px 8px; border: 1px solid #ccc; background: #f5f5f5; border-radius: 3px; cursor: pointer; font-size: 12px; }
  .control-row button:hover { background: #eee; }
  .control-row label { flex: 1; font-size: 12px; color: #333; cursor: pointer; }
  .control-row input[type=checkbox] { cursor: pointer; }

  .legend-item { display: flex; align-items: center; gap: 8px; margin: 4px 0; font-size: 12px; }
  .swatch { width: 18px; height: 18px; border-radius: 3px; flex-shrink: 0; border: 1px solid #aaa; }
  .swatch.round { border-radius: 50%; }
  .edge-swatch { width: 28px; height: 3px; flex-shrink: 0; }

  #tooltip { position: absolute; background: rgba(30, 30, 40, 0.96); color: white;
             padding: 10px 12px; border-radius: 6px; font-size: 12px; line-height: 1.4;
             pointer-events: none; max-width: 320px; box-shadow: 0 4px 16px rgba(0,0,0,0.3);
             z-index: 1000; display: none; }
  #tooltip .tip-title { font-weight: 600; font-size: 13px; margin-bottom: 4px; }
  #tooltip .tip-row { display: flex; justify-content: space-between; gap: 10px; }
  #tooltip .tip-label { color: #bbb; }
  #tooltip .tip-value { color: white; font-family: ui-monospace, monospace; }
  #tooltip .up { color: #ff8a80; }
  #tooltip .down { color: #80d8ff; }

  #info { position: absolute; top: 8px; right: 8px; background: white; border: 1px solid #e0e0e0;
          padding: 6px 10px; border-radius: 4px; font-size: 11px; color: #666; box-shadow: 0 1px 3px rgba(0,0,0,0.08); }
</style>
</head>
<body>
<div id="app">
  <div id="header">
    <div>
      <h1>TCGA-LIHC Copper Proteome — Interactive Network</h1>
      <div class="subtitle">54 Cu genes · node size = GNN importance · node colour = LIHC log2FC · edge width = GAT attention</div>
    </div>
  </div>

  <div id="sidebar">
    <div class="panel">
      <h3>Search</h3>
      <div class="control-row">
        <input id="search" type="search" placeholder="gene symbol, e.g. ATOX1">
      </div>
      <div class="control-row">
        <button id="fit">Fit view</button>
        <button id="reset">Reset highlight</button>
      </div>
    </div>

    <div class="panel">
      <h3>Highlight</h3>
      <div class="control-row">
        <input type="checkbox" id="hl_top_nodes" checked>
        <label for="hl_top_nodes">Top-10 important genes (thick border)</label>
      </div>
      <div class="control-row">
        <input type="checkbox" id="hl_top_edges" checked>
        <label for="hl_top_edges">Top-15 attention edges (bold)</label>
      </div>
      <div class="control-row">
        <input type="checkbox" id="hl_bio" checked>
        <label for="hl_bio">Canonical Cu-biology edges (red glow)</label>
      </div>
    </div>

    <div class="panel">
      <h3>Edge types — show/hide</h3>
      <div class="control-row">
        <input type="checkbox" id="show_physical" checked>
        <label for="show_physical">physical</label>
        <span class="edge-swatch" style="background:#d32f2f"></span>
      </div>
      <div class="control-row">
        <input type="checkbox" id="show_coexpression" checked>
        <label for="show_coexpression">co-expression</label>
        <span class="edge-swatch" style="background:#7b1fa2"></span>
      </div>
      <div class="control-row">
        <input type="checkbox" id="show_genetic" checked>
        <label for="show_genetic">genetic</label>
        <span class="edge-swatch" style="background:#2e7d32"></span>
      </div>
      <div class="control-row">
        <input type="checkbox" id="show_shared_compartment" checked>
        <label for="show_shared_compartment">shared compartment</label>
        <span class="edge-swatch" style="background:#9e9e9e"></span>
      </div>
      <div class="control-row">
        <input type="checkbox" id="show_attention_only" checked>
        <label for="show_attention_only">attention-only edges (not in graph)</label>
        <span class="edge-swatch" style="background:#1565c0;border-bottom:1px dashed #1565c0"></span>
      </div>
    </div>

    <div class="panel">
      <h3>Node colour — log2FC</h3>
      <div class="legend-item"><span class="swatch round" style="background:#b22222"></span>up in tumor</div>
      <div class="legend-item"><span class="swatch round" style="background:#f5f5f5"></span>unchanged</div>
      <div class="legend-item"><span class="swatch round" style="background:#1f4e79"></span>down in tumor</div>
    </div>

    <div class="panel">
      <h3>Node shape — category</h3>
      <div class="legend-item"><span class="swatch" style="border-radius:50%;background:#f5f5f5"></span>transporter (circle)</div>
      <div class="legend-item"><span class="swatch" style="border-radius:4px;background:#f5f5f5"></span>enzyme (rounded square)</div>
      <div class="legend-item"><span class="swatch" style="background:#f5f5f5;clip-path:polygon(50% 0,100% 100%,0 100%)"></span>other (triangle)</div>
    </div>

    <div class="panel">
      <h3>Interact</h3>
      <div style="font-size:11px;color:#666;line-height:1.5">
        • drag nodes to rearrange<br>
        • scroll to zoom<br>
        • click-drag background to pan<br>
        • hover over a node or edge for details<br>
        • click a node to pin its tooltip
      </div>
    </div>

    <div class="panel">
      <h3>Stats</h3>
      <div style="font-size:11px;color:#666" id="stats"></div>
    </div>
  </div>

  <div id="cy">
    <div id="info">Right-click / hover for tooltip. Scroll to zoom.</div>
    <div id="tooltip"></div>
  </div>
</div>

<script>
const DATA = __DATA_JSON__;
const META = __META_JSON__;

const EDGE_COLOR = {
  physical: "#d32f2f",
  coexpression: "#7b1fa2",
  genetic: "#2e7d32",
  shared_compartment: "#9e9e9e",
  attention_only: "#1565c0",
};
const NODE_SHAPE = {
  transporter: "ellipse",
  enzyme: "round-rectangle",
  other_or_unknown: "triangle",
};

const FC_MIN = -3, FC_MAX = 3;
function lerp(a, b, t) { return a + (b - a) * t; }
function lerpColor(rgbA, rgbB, t) {
  return [Math.round(lerp(rgbA[0], rgbB[0], t)),
          Math.round(lerp(rgbA[1], rgbB[1], t)),
          Math.round(lerp(rgbA[2], rgbB[2], t))];
}
function logfcColor(fc) {
  const blue = [31, 78, 121], white = [245, 245, 245], red = [178, 34, 34];
  const clipped = Math.max(FC_MIN, Math.min(FC_MAX, fc));
  if (clipped >= 0) {
    const t = clipped / FC_MAX;
    const [r,g,b] = lerpColor(white, red, t);
    return `rgb(${r},${g},${b})`;
  } else {
    const t = -clipped / (-FC_MIN);
    const [r,g,b] = lerpColor(white, blue, t);
    return `rgb(${r},${g},${b})`;
  }
}
function nodeSize(imp_norm) { return 18 + 42 * imp_norm; }
function edgeWidth(attn, has_attn) {
  if (!has_attn || attn <= 0) return 1.5;
  return 1.5 + 7.0 * (attn / META.attn_max);
}

const elements = [];
DATA.nodes.forEach(n => {
  const d = n.data;
  d.color = logfcColor(d.log2FC);
  d.size = nodeSize(d.importance_norm);
  d.shape = NODE_SHAPE[d.category] || "ellipse";
  d.border_color = d.is_top_importance ? "#111" : "#777";
  d.border_width = d.is_top_importance ? 3 : 1;
  elements.push(n);
});
DATA.edges.forEach(e => {
  const d = e.data;
  d.color = EDGE_COLOR[d.edge_type] || "#aaa";
  d.line_style = d.edge_type === "attention_only" ? "dashed" : "solid";
  d.width = edgeWidth(d.attention, d.attention > 0);
  d.glow = d.is_bio_canonical ? "#ff1744" : "transparent";
  elements.push(e);
});

const cy = cytoscape({
  container: document.getElementById("cy"),
  elements: elements,
  layout: {
    name: "fcose",
    animate: false,
    quality: "proof",
    idealEdgeLength: 90,
    nodeSeparation: 85,
    nodeRepulsion: 6500,
    gravity: 0.25,
    padding: 40,
  },
  style: [
    { selector: "node",
      style: {
        "background-color": "data(color)",
        "label": "data(label)",
        "text-valign": "center",
        "text-halign": "center",
        "font-size": 10,
        "font-weight": 500,
        "color": "#111",
        "text-outline-color": "#fff",
        "text-outline-width": 1.5,
        "width": "data(size)",
        "height": "data(size)",
        "shape": "data(shape)",
        "border-color": "data(border_color)",
        "border-width": "data(border_width)",
      }
    },
    { selector: "edge",
      style: {
        "width": "data(width)",
        "line-color": "data(color)",
        "line-style": "data(line_style)",
        "curve-style": "bezier",
        "opacity": 0.7,
      }
    },
    { selector: "edge[?is_bio_canonical]",
      style: { "line-color": "#c62828", "opacity": 0.95, "width": "data(width)" }
    },
    { selector: "edge[?is_top_attention]",
      style: { "target-arrow-color": "#333" }
    },
    { selector: ".faded",
      style: { "opacity": 0.12 }
    },
    { selector: ".selected-edge",
      style: { "line-color": "#000", "opacity": 1.0 }
    },
    { selector: "node.selected-node",
      style: { "border-color": "#000", "border-width": 4 }
    },
  ],
  wheelSensitivity: 0.25,
  minZoom: 0.3,
  maxZoom: 4.0,
});

// Stats
const sig = DATA.nodes.filter(n => n.data.significant).length;
document.getElementById("stats").innerHTML =
  `nodes: ${DATA.nodes.length} (${sig} sig. BH<0.05)<br>` +
  `edges: ${DATA.edges.length}<br>` +
  `top-10 important: <b>${META.n_top_importance}</b><br>` +
  `top attention edges: ${META.n_top_attention}`;

// Tooltip
const tip = document.getElementById("tooltip");
let tipPinned = false;

function fmt(x, d = 3) {
  if (x === null || x === undefined || Number.isNaN(x)) return "-";
  if (typeof x !== "number") return x;
  if (Math.abs(x) < 1e-3 && x !== 0) return x.toExponential(2);
  return x.toFixed(d);
}

function nodeTip(n) {
  const d = n.data();
  const fcClass = d.log2FC > 0 ? "up" : (d.log2FC < 0 ? "down" : "");
  const sigBadge = d.significant ? "<span style='color:#69f0ae'>✓ sig</span>" : "<span style='color:#bbb'>ns</span>";
  return `
    <div class="tip-title">${d.label}   ${sigBadge}</div>
    <div class="tip-row"><span class="tip-label">log2FC</span><span class="tip-value ${fcClass}">${fmt(d.log2FC)}</span></div>
    <div class="tip-row"><span class="tip-label">adj p (BH)</span><span class="tip-value">${fmt(d.adj_p, 3)}</span></div>
    <div class="tip-row"><span class="tip-label">GNN importance</span><span class="tip-value">${fmt(d.importance)}</span></div>
    <div class="tip-row"><span class="tip-label">category</span><span class="tip-value">${d.category}</span></div>
    <div class="tip-row"><span class="tip-label">module</span><span class="tip-value">#${d.module}</span></div>
    <div class="tip-row"><span class="tip-label">localization</span><span class="tip-value" style="text-align:right;max-width:180px">${(d.localization || '').replaceAll('|',', ')}</span></div>
    <div style="margin-top:6px;color:#ccc;font-size:11px;font-style:italic;max-width:300px">${d.notes || ''}</div>`;
}

function edgeTip(e) {
  const d = e.data();
  let label = d.edge_type.replaceAll("_", " ");
  const badge = d.is_bio_canonical ? " <span style='color:#ff6e6e'>★ canonical Cu biology</span>" : "";
  return `
    <div class="tip-title">${d.source} — ${d.target}${badge}</div>
    <div class="tip-row"><span class="tip-label">type</span><span class="tip-value">${label}</span></div>
    <div class="tip-row"><span class="tip-label">GAT attention</span><span class="tip-value">${fmt(d.attention)}</span></div>
    <div class="tip-row"><span class="tip-label">prior weight</span><span class="tip-value">${fmt(d.weight)}</span></div>
    <div class="tip-row"><span class="tip-label">top attention?</span><span class="tip-value">${d.is_top_attention ? "yes" : "no"}</span></div>`;
}

function showTipAt(content, x, y) {
  tip.innerHTML = content;
  tip.style.left = (x + 14) + "px";
  tip.style.top = (y + 14) + "px";
  tip.style.display = "block";
}
function hideTip() { if (!tipPinned) tip.style.display = "none"; }

cy.on("mouseover", "node", evt => { if (!tipPinned) showTipAt(nodeTip(evt.target), evt.originalEvent.clientX, evt.originalEvent.clientY); });
cy.on("mouseover", "edge", evt => { if (!tipPinned) showTipAt(edgeTip(evt.target), evt.originalEvent.clientX, evt.originalEvent.clientY); });
cy.on("mousemove", "node, edge", evt => { if (!tipPinned) { tip.style.left = (evt.originalEvent.clientX + 14) + "px"; tip.style.top = (evt.originalEvent.clientY + 14) + "px"; }});
cy.on("mouseout", "node, edge", () => hideTip());

cy.on("tap", "node", evt => {
  tipPinned = true;
  cy.elements().removeClass("selected-node selected-edge faded");
  const n = evt.target;
  const hood = n.closedNeighborhood();
  cy.elements().not(hood).addClass("faded");
  n.addClass("selected-node");
  n.connectedEdges().addClass("selected-edge");
  showTipAt(nodeTip(n), evt.originalEvent.clientX, evt.originalEvent.clientY);
});
cy.on("tap", "edge", evt => {
  tipPinned = true;
  cy.elements().removeClass("selected-node selected-edge faded");
  const e = evt.target;
  cy.elements().not(e.connectedNodes().union(e)).addClass("faded");
  e.addClass("selected-edge");
  showTipAt(edgeTip(e), evt.originalEvent.clientX, evt.originalEvent.clientY);
});
cy.on("tap", evt => {
  if (evt.target === cy) { tipPinned = false; hideTip(); cy.elements().removeClass("selected-node selected-edge faded"); }
});

// Search
document.getElementById("search").addEventListener("input", evt => {
  const q = evt.target.value.trim().toUpperCase();
  cy.elements().removeClass("faded");
  if (!q) return;
  cy.nodes().forEach(n => {
    if (!n.data("label").toUpperCase().includes(q)) n.addClass("faded");
  });
});
document.getElementById("fit").addEventListener("click", () => cy.fit(undefined, 40));
document.getElementById("reset").addEventListener("click", () => {
  tipPinned = false; hideTip();
  cy.elements().removeClass("selected-node selected-edge faded");
  document.getElementById("search").value = "";
});

function applyEdgeVisibility() {
  const types = ["physical", "coexpression", "genetic", "shared_compartment", "attention_only"];
  types.forEach(t => {
    const on = document.getElementById("show_" + t).checked;
    cy.edges().forEach(e => {
      if (e.data("edge_type") === t) e.style("display", on ? "element" : "none");
    });
  });
}
["physical", "coexpression", "genetic", "shared_compartment", "attention_only"].forEach(t => {
  document.getElementById("show_" + t).addEventListener("change", applyEdgeVisibility);
});

document.getElementById("hl_top_nodes").addEventListener("change", evt => {
  cy.nodes().forEach(n => {
    if (n.data("is_top_importance")) {
      n.style("border-width", evt.target.checked ? 3 : 1);
      n.style("border-color", evt.target.checked ? "#111" : "#777");
    }
  });
});
document.getElementById("hl_top_edges").addEventListener("change", evt => {
  cy.edges().forEach(e => {
    if (e.data("is_top_attention")) {
      e.style("opacity", evt.target.checked ? 1.0 : 0.5);
      e.style("width", evt.target.checked ? Math.max(3, e.data("width")) : e.data("width"));
    }
  });
});
document.getElementById("hl_bio").addEventListener("change", evt => {
  cy.edges().forEach(e => {
    if (e.data("is_bio_canonical")) {
      e.style("line-color", evt.target.checked ? "#c62828" : e.data("color"));
      e.style("width", evt.target.checked ? Math.max(4, e.data("width")) : e.data("width"));
    }
  });
});

cy.ready(() => cy.fit(undefined, 50));
</script>
</body>
</html>
"""


def main():
    print("[interactive] loading data")
    payload = load_all()
    print(f"[interactive] {len(payload['nodes'])} nodes, {len(payload['edges'])} edges")
    html = (HTML_TEMPLATE
            .replace("__DATA_JSON__", json.dumps({"nodes": payload["nodes"], "edges": payload["edges"]}))
            .replace("__META_JSON__", json.dumps(payload["meta"])))
    out_html = OUT / "copper_interactive_network.html"
    out_html.write_text(html)
    print(f"[interactive] wrote {out_html} ({out_html.stat().st_size/1024:.1f} KB)")

    notes = OUT / "copper_interactive_network_notes.md"
    notes.write_text(f"""# Interactive Copper Network — Legend & Notes

Open [`copper_interactive_network.html`](copper_interactive_network.html) in any modern browser. No server needed — it is a single self-contained HTML file that pulls Cytoscape.js from CDN.

## Layout

The initial layout is computed with **fCoSE** — a force-directed algorithm tuned for biological networks. After load you can drag any node to rearrange it, scroll to zoom, and click-drag the background to pan.

## Nodes — 54 Cu proteome genes

| visual | meaning |
|---|---|
| **size** | GNN saliency (`outputs/gnn/node_importance.csv`). Bigger = more influential on the tumor/normal decision. Top-10 genes have a thick black border. |
| **colour** | TCGA-LIHC log2FC from `copper_de_results.csv`. Red = up in tumor, blue = down, near-white = no change. Colour scale is clipped at ±3. |
| **shape** | functional category: circle = transporter, rounded square = enzyme, triangle = other/unknown. |
| **tooltip (hover)** | gene name, log2FC, adj p (BH), GNN importance, category, module id, subcellular localization, and the paper-derived notes field. |

## Edges

| visual | meaning |
|---|---|
| **colour** | interaction type — red = physical / PPI, purple = co-expression, green = genetic, grey = shared subcellular compartment, blue dashed = attention-only (GAT top edge that is not in the fixed topology) |
| **width** | GAT attention sum (`outputs/gnn/top_attention_edges.csv`). Thicker = higher attention. Edges with no attention data get default thin width. |
| **bold red** | canonical Cu-handling pair (CCS–SLC31A1, SOD3–ATOX1, ATP7B–COMMD1, MT-CO2–SCO2, …). Toggleable via the "Canonical Cu-biology edges" checkbox. |
| **tooltip (hover)** | source–target pair, edge type, GAT attention, prior weight, top-attention flag. |

## Controls

- **Search box** — type a gene symbol or prefix to fade everything else.
- **Fit view / Reset highlight** — self-explanatory.
- **Highlight toggles** — turn on/off the top-10 important nodes, top-15 attention edges, canonical Cu-biology edges.
- **Edge type toggles** — show/hide edges by interaction type (useful to see the attention-only blue edges on their own).
- **Click a node** — pins its tooltip and fades everything outside its 1-hop neighbourhood.
- **Click an edge** — pins its tooltip and fades everything except its two endpoints.
- **Click the background** — unpins.

## What to look for

1. The **blue cluster** in the lower-right (ALB, CP, DBH, MT-CO1/CO2, SOD1) is the hepatocyte-secretome + mitochondrial-COX axis, strongly downregulated.
2. The **red cluster** (LOX, LOXL2, SPARC, AFP) is the ECM-remodelling + HCC-marker axis, upregulated.
3. Switch OFF all edge types except "attention-only" to see **purely what the GAT learned** — the canonical Cu-handling pairs (CCS↔SLC31A1, SOD3↔ATOX1, ATP7B↔COMMD1, MT-CO2↔SCO2) should dominate.
4. Click **ATOX1** or **ATP7B** — both have modest log2FC but high GNN importance because their 1-hop neighbourhood is central to Cu homeostasis.
5. Toggle off the grey "shared compartment" edges to declutter the graph and focus on biology-derived interactions.

## Stats
- Nodes: {len(payload['nodes'])} (of which {sum(1 for n in payload['nodes'] if n['data']['significant'])} significant at BH<0.05)
- Edges: {len(payload['edges'])} (including {sum(1 for e in payload['edges'] if e['data']['edge_type']=='attention_only')} attention-only)
- Top-10 important genes: highlighted with thick black border
- Canonical Cu-biology edges: red glow when the "Canonical Cu-biology edges" toggle is on
""")
    print(f"[interactive] wrote {notes}")


if __name__ == "__main__":
    main()
