"""Project-wide configuration: paths, constants, seeds."""
from __future__ import annotations
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = PROJECT_ROOT / "data"
SRC_DIR = PROJECT_ROOT / "src"
OUTPUTS_DIR = PROJECT_ROOT / "outputs"

PAPER_EXTRACTION_DIR = OUTPUTS_DIR / "paper_2017_extraction"
BASELINE_DIR = OUTPUTS_DIR / "baseline"
GNN_DIR = OUTPUTS_DIR / "gnn"
FINAL_COMPARISON_DIR = OUTPUTS_DIR / "final_comparison"

EXPRESSION_FILE = DATA_DIR / "lihc_expression.tsv"
METADATA_FILE = DATA_DIR / "lihc_metadata.tsv"
COPPER_GENE_FILE = PAPER_EXTRACTION_DIR / "copper_gene_list.csv"

RANDOM_SEED = 42

for d in (DATA_DIR, BASELINE_DIR, GNN_DIR, FINAL_COMPARISON_DIR):
    d.mkdir(parents=True, exist_ok=True)
