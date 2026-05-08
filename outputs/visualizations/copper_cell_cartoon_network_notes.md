# Cartoon Cell Copper Network

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
