"""Load TCGA-LIHC expression + metadata and subset to the Cu proteome.

Expected file formats
---------------------
- data/lihc_expression.tsv
    TSV with rows=genes, columns=samples. First column: gene_symbol (HGNC).
- data/lihc_metadata.tsv
    TSV with one row per sample and columns: sample_id, sample_type (values
    'Tumor' or 'Normal'), plus any optional clinical columns
    (stage, grade, overall_survival_days, vital_status).

If neither file is present, a small synthetic dataset is generated so the
whole pipeline can still be exercised end-to-end. See ``build_synthetic``.
"""
from __future__ import annotations
from pathlib import Path
from dataclasses import dataclass
import numpy as np
import pandas as pd

from src.utils import (
    EXPRESSION_FILE, METADATA_FILE, COPPER_GENE_FILE, RANDOM_SEED,
)


@dataclass
class LIHCDataset:
    expression: pd.DataFrame       # rows=genes (Cu proteome subset), cols=sample_id
    metadata: pd.DataFrame         # indexed by sample_id, has sample_type column
    copper_genes: pd.DataFrame     # the 54-gene annotation table
    missing_genes: list[str]       # Cu genes not found in the expression matrix

    @property
    def n_samples(self) -> int:
        return self.expression.shape[1]

    @property
    def n_genes(self) -> int:
        return self.expression.shape[0]

    @property
    def tumor_mask(self) -> np.ndarray:
        return (self.metadata.loc[self.expression.columns, "sample_type"]
                .str.lower().eq("tumor").to_numpy())


def load_copper_genes(path: Path = COPPER_GENE_FILE) -> pd.DataFrame:
    df = pd.read_csv(path)
    df["gene_symbol"] = df["gene_symbol"].str.upper().str.strip()
    return df


def _standardize_expression(expr: pd.DataFrame) -> pd.DataFrame:
    expr = expr.copy()
    expr.columns = [str(c).strip() for c in expr.columns]
    if expr.index.name is None:
        expr.index.name = "gene_symbol"
    expr.index = expr.index.astype(str).str.upper().str.strip()
    expr = expr[~expr.index.duplicated(keep="first")]
    return expr


def load_expression(path: Path = EXPRESSION_FILE) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Expression file not found: {path}")
    df = pd.read_csv(path, sep="\t", index_col=0)
    return _standardize_expression(df)


def load_metadata(path: Path = METADATA_FILE) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Metadata file not found: {path}")
    md = pd.read_csv(path, sep="\t")
    md.columns = [c.strip() for c in md.columns]
    md["sample_id"] = md["sample_id"].astype(str).str.strip()
    md = md.set_index("sample_id")
    if "sample_type" not in md.columns:
        raise ValueError("metadata must contain a 'sample_type' column "
                         "with values 'Tumor' or 'Normal'")
    md["sample_type"] = md["sample_type"].str.strip().str.title()
    return md


def subset_to_copper(expr: pd.DataFrame, copper: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    wanted = copper["gene_symbol"].tolist()
    present = [g for g in wanted if g in expr.index]
    missing = [g for g in wanted if g not in expr.index]
    return expr.loc[present], missing


def load_lihc_dataset(require_real: bool = False) -> LIHCDataset:
    """Main entry point. Falls back to synthetic data when files are absent.

    If ``require_real=True``, raises FileNotFoundError instead of synthesising —
    used by Phase 2 real-data runs where a silent fallback is dangerous.
    """
    copper = load_copper_genes()
    if not EXPRESSION_FILE.exists() or not METADATA_FILE.exists():
        if require_real:
            raise FileNotFoundError(
                "Real TCGA-LIHC data required but not found. Expected:\n"
                f"  {EXPRESSION_FILE}\n  {METADATA_FILE}\n"
                "See data/EXPECTED_INPUT_FORMAT.md for the required format."
            )
        print("[load_data] Real TCGA-LIHC files not found -> building synthetic demo dataset.")
        print(f"[load_data] Provide {EXPRESSION_FILE} and {METADATA_FILE} to run on real data.")
        expr, md = build_synthetic(copper)
    else:
        expr = load_expression()
        md = load_metadata()
        shared = [s for s in expr.columns if s in md.index]
        expr = expr[shared]
        md = md.loc[shared]

    expr_cu, missing = subset_to_copper(expr, copper)
    return LIHCDataset(
        expression=expr_cu, metadata=md, copper_genes=copper, missing_genes=missing,
    )


def build_synthetic(copper: pd.DataFrame, n_tumor: int = 120, n_normal: int = 30,
                    seed: int = RANDOM_SEED) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Generate plausible log-scale expression so downstream code can run.

    Tumor samples get a shift on a random subset of genes to create a learnable
    tumor-vs-normal signal. Do NOT interpret synthetic results biologically.
    """
    rng = np.random.default_rng(seed)
    genes = copper["gene_symbol"].tolist()
    n_genes = len(genes)
    tumor_ids = [f"LIHC_T{i:03d}" for i in range(n_tumor)]
    normal_ids = [f"LIHC_N{i:03d}" for i in range(n_normal)]
    sample_ids = tumor_ids + normal_ids

    baseline = rng.normal(loc=6.0, scale=1.5, size=(n_genes, 1))
    noise = rng.normal(loc=0.0, scale=0.8, size=(n_genes, len(sample_ids)))
    data = baseline + noise

    up_idx = rng.choice(n_genes, size=max(1, n_genes // 3), replace=False)
    down_idx = rng.choice([i for i in range(n_genes) if i not in up_idx],
                          size=max(1, n_genes // 4), replace=False)
    data[up_idx, :n_tumor] += rng.uniform(0.8, 2.2, size=(len(up_idx), n_tumor))
    data[down_idx, :n_tumor] -= rng.uniform(0.6, 1.8, size=(len(down_idx), n_tumor))

    expr = pd.DataFrame(data, index=genes, columns=sample_ids)
    md = pd.DataFrame({
        "sample_id": sample_ids,
        "sample_type": ["Tumor"] * n_tumor + ["Normal"] * n_normal,
        "stage": rng.choice(["I", "II", "III", "IV"], size=len(sample_ids),
                             p=[0.35, 0.35, 0.2, 0.1]),
    }).set_index("sample_id")
    md.loc[normal_ids, "stage"] = np.nan
    return _standardize_expression(expr), md
