"""Download TCGA-LIHC RNA-seq (STAR-Counts) from the GDC API.

Equivalent to the TCGAbiolinks R snippet:

    query <- GDCquery(project="TCGA-LIHC",
                      data.category="Transcriptome Profiling",
                      data.type="Gene Expression Quantification",
                      workflow.type="STAR - Counts")
    GDCdownload(query)
    data <- GDCprepare(query)

Produces:
  data/lihc_expression.tsv   rows=HGNC symbol, cols=sample_id,
                             values=log2(fpkm_uq_unstranded + 1)
  data/lihc_metadata.tsv     sample_id, sample_type (Tumor/Normal),
                             plus extra TCGA-level annotations
  data/gdc_raw/              per-sample TSVs (kept; can be deleted later)
  data/manifest.tsv          the GDC manifest used

Re-running is idempotent — files already downloaded are skipped.
"""
from __future__ import annotations
import argparse
import gzip
import io
import json
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import requests
import pandas as pd

GDC_API = "https://api.gdc.cancer.gov"
ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
RAW_DIR = DATA_DIR / "gdc_raw"
RAW_DIR.mkdir(parents=True, exist_ok=True)

FIELDS = [
    "file_id", "file_name", "file_size", "md5sum",
    "cases.submitter_id",
    "cases.samples.submitter_id",
    "cases.samples.sample_type",
    "cases.samples.sample_type_id",
    "cases.samples.tissue_type",
    "cases.diagnoses.ajcc_pathologic_stage",
    "cases.diagnoses.tumor_grade",
    "cases.diagnoses.days_to_last_follow_up",
    "cases.demographic.vital_status",
    "cases.demographic.days_to_death",
    "cases.demographic.gender",
    "cases.demographic.race",
    "cases.demographic.age_at_index",
]


def get_file_manifest() -> pd.DataFrame:
    filters = {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": ["TCGA-LIHC"]}},
            {"op": "in", "content": {"field": "data_category", "value": ["Transcriptome Profiling"]}},
            {"op": "in", "content": {"field": "data_type", "value": ["Gene Expression Quantification"]}},
            {"op": "in", "content": {"field": "analysis.workflow_type", "value": ["STAR - Counts"]}},
            {"op": "in", "content": {"field": "access", "value": ["open"]}},
        ],
    }
    payload = {
        "filters": filters, "size": 10000, "format": "JSON",
        "fields": ",".join(FIELDS),
    }
    r = requests.post(f"{GDC_API}/files", json=payload, timeout=60)
    r.raise_for_status()
    hits = r.json()["data"]["hits"]
    print(f"[download] GDC returned {len(hits)} files for TCGA-LIHC STAR-Counts")

    rows = []
    for h in hits:
        case = h["cases"][0] if h.get("cases") else {}
        sample = case.get("samples", [{}])[0]
        diag = (case.get("diagnoses") or [{}])[0]
        demo = case.get("demographic", {}) or {}
        rows.append({
            "file_id": h["file_id"],
            "file_name": h["file_name"],
            "file_size": h["file_size"],
            "md5sum": h.get("md5sum"),
            "case_submitter_id": case.get("submitter_id"),
            "sample_submitter_id": sample.get("submitter_id"),
            "sample_type": sample.get("sample_type"),
            "sample_type_id": sample.get("sample_type_id"),
            "tissue_type": sample.get("tissue_type"),
            "ajcc_pathologic_stage": diag.get("ajcc_pathologic_stage"),
            "tumor_grade": diag.get("tumor_grade"),
            "days_to_last_follow_up": diag.get("days_to_last_follow_up"),
            "vital_status": demo.get("vital_status"),
            "days_to_death": demo.get("days_to_death"),
            "gender": demo.get("gender"),
            "race": demo.get("race"),
            "age_at_index": demo.get("age_at_index"),
        })
    mf = pd.DataFrame(rows)
    mf.to_csv(DATA_DIR / "manifest.tsv", sep="\t", index=False)
    return mf


def download_one(file_id: str, file_name: str, dest: Path,
                  session: requests.Session, max_retries: int = 3) -> bool:
    if dest.exists() and dest.stat().st_size > 0:
        return True
    url = f"{GDC_API}/data/{file_id}"
    for attempt in range(max_retries):
        try:
            with session.get(url, stream=True, timeout=180) as r:
                r.raise_for_status()
                tmp = dest.with_suffix(dest.suffix + ".part")
                with open(tmp, "wb") as f:
                    for chunk in r.iter_content(chunk_size=1 << 20):
                        if chunk:
                            f.write(chunk)
                tmp.replace(dest)
            return True
        except Exception as e:
            print(f"[download] retry {attempt+1}/{max_retries} for {file_id}: {e}")
            time.sleep(2 ** attempt)
    return False


def download_all(mf: pd.DataFrame, max_workers: int = 8) -> list[str]:
    session = requests.Session()
    session.headers.update({"Accept": "application/json"})

    failures = []
    total = len(mf)
    done = 0
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        futures = {
            ex.submit(download_one, row.file_id, row.file_name,
                      RAW_DIR / row.file_name, session): row.file_id
            for row in mf.itertuples(index=False)
        }
        for fut in as_completed(futures):
            done += 1
            fid = futures[fut]
            ok = False
            try:
                ok = fut.result()
            except Exception as e:
                print(f"[download] {fid} failed: {e}")
            if not ok:
                failures.append(fid)
            if done % 20 == 0 or done == total:
                dt = time.time() - t0
                print(f"[download] {done}/{total} files "
                      f"({dt:.0f}s, {len(failures)} failed)")
    session.close()
    return failures


def parse_one(path: Path) -> pd.Series | None:
    """Extract fpkm_uq_unstranded keyed by gene_name, averaging duplicates."""
    try:
        df = pd.read_csv(path, sep="\t", comment="#", skiprows=1, dtype=str)
    except Exception as e:
        print(f"[parse] read failed for {path.name}: {e}")
        return None
    required = {"gene_id", "gene_name", "fpkm_uq_unstranded"}
    if not required.issubset(df.columns):
        print(f"[parse] missing columns in {path.name}; got {df.columns.tolist()[:8]}")
        return None
    df = df[~df["gene_id"].str.startswith(("N_", "__"), na=False)]
    df = df[df["gene_name"].notna() & (df["gene_name"].astype(str).str.strip() != "")]
    df["gene_name"] = df["gene_name"].str.upper().str.strip()
    df["fpkm_uq_unstranded"] = pd.to_numeric(df["fpkm_uq_unstranded"], errors="coerce")
    df = df.dropna(subset=["fpkm_uq_unstranded"])
    series = (df.groupby("gene_name")["fpkm_uq_unstranded"].max())
    return series


def build_expression_matrix(mf: pd.DataFrame, copper_genes: list[str]) -> pd.DataFrame:
    """Build Cu-gene only expression matrix to keep memory small."""
    copper_set = set(copper_genes)
    cols = {}
    for row in mf.itertuples(index=False):
        path = RAW_DIR / row.file_name
        if not path.exists():
            continue
        series = parse_one(path)
        if series is None:
            continue
        series = series[series.index.isin(copper_set)]
        cols[row.sample_submitter_id] = series
    if not cols:
        raise RuntimeError("no files parsed")
    matrix = pd.concat(cols, axis=1)
    matrix.index.name = "gene_symbol"
    matrix = matrix.reindex(copper_genes)
    import numpy as np
    matrix = np.log2(matrix + 1.0)
    return matrix


def derive_sample_type(sample_type_raw: str | None) -> str:
    if not sample_type_raw:
        return "Unknown"
    s = str(sample_type_raw).lower()
    if "normal" in s:
        return "Normal"
    if "tumor" in s or "tumour" in s or "metasta" in s or "neoplasm" in s or "primary" in s:
        return "Tumor"
    return "Unknown"


def build_metadata(mf: pd.DataFrame, expression_columns: list[str]) -> pd.DataFrame:
    md = mf.copy()
    md["sample_id"] = md["sample_submitter_id"]
    md["sample_type"] = md["sample_type"].apply(derive_sample_type)

    overall_survival = md["days_to_death"].fillna(md["days_to_last_follow_up"])
    md["overall_survival_days"] = pd.to_numeric(overall_survival, errors="coerce")

    md = md[[
        "sample_id", "sample_type", "case_submitter_id",
        "ajcc_pathologic_stage", "tumor_grade",
        "overall_survival_days", "vital_status",
        "gender", "race", "age_at_index",
    ]]
    md = md.rename(columns={
        "ajcc_pathologic_stage": "stage",
        "tumor_grade": "grade",
    })
    md = md[md["sample_id"].isin(expression_columns)]
    md = md.drop_duplicates(subset="sample_id", keep="first")
    return md


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-workers", type=int, default=8)
    ap.add_argument("--skip-download", action="store_true",
                    help="use files already in data/gdc_raw")
    args = ap.parse_args()

    print(f"[download] project = TCGA-LIHC, workflow = STAR - Counts")
    mf = get_file_manifest()
    print(f"[download] manifest saved to {DATA_DIR/'manifest.tsv'}")
    total_gb = mf["file_size"].sum() / 1e9
    print(f"[download] total payload: {total_gb:.2f} GB across {len(mf)} files")

    if not args.skip_download:
        fails = download_all(mf, max_workers=args.max_workers)
        if fails:
            print(f"[download] {len(fails)} files failed after retries — see stderr")
    else:
        print("[download] --skip-download set; assuming files exist in gdc_raw")

    from src.preprocessing import load_copper_genes
    copper = load_copper_genes()
    gene_list = copper["gene_symbol"].tolist()
    print(f"[prepare] building expression matrix for {len(gene_list)} Cu genes")
    expr = build_expression_matrix(mf, gene_list)
    print(f"[prepare] expression matrix: {expr.shape[0]} genes x {expr.shape[1]} samples")

    md = build_metadata(mf, expr.columns.tolist())
    print(f"[prepare] metadata: {md.shape[0]} rows; "
          f"tumor={int((md['sample_type']=='Tumor').sum())}, "
          f"normal={int((md['sample_type']=='Normal').sum())}")

    expr.to_csv(DATA_DIR / "lihc_expression.tsv", sep="\t")
    md.to_csv(DATA_DIR / "lihc_metadata.tsv", sep="\t", index=False)
    print(f"[prepare] wrote {DATA_DIR/'lihc_expression.tsv'}")
    print(f"[prepare] wrote {DATA_DIR/'lihc_metadata.tsv'}")


if __name__ == "__main__":
    sys.path.insert(0, str(ROOT))
    main()
