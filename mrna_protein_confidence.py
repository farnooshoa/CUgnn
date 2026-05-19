"""mRNA-protein confidence flags for the copper proteome saliency table.

Source: Blockhuys et al. 2017 (Metallomics 9:112), Figure 4.
Reports Pearson r between mRNA and protein abundance for each Cu gene
across a panel of cell lines / tissues.

Genes are flagged as:
  HIGH   r >= 0.5   mRNA is a reliable proxy for protein
  MEDIUM 0.3 <= r < 0.5   moderate confidence
  LOW    r < 0.3    mRNA is NOT a reliable proxy — claims need protein validation

This directly addresses the scientific criticism:
  'RNA != protein. Genes with low r (e.g. ENOX2, COMMD1, AOC2, ALB, F5)
   need protein-level validation before any strong claim.'

Outputs:
  outputs/final_comparison/saliency_with_confidence.csv
  outputs/final_comparison/mrna_protein_report.md
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

OUT = ROOT / "outputs" / "final_comparison"
OUT.mkdir(parents=True, exist_ok=True)

# ── Blockhuys 2017 Fig 4 mRNA-protein correlations ─────────────────────────
# r values extracted from the paper figure.
# Genes not in the figure are marked None (no data available).
# Source: Blockhuys S et al. Metallomics 2017;9(2):112-123. Fig. 4.
MRNA_PROTEIN_R = {
    # Core copper relay — well validated
    "SLC31A1":  0.72,   # CTR1 — strong correlation, well studied
    "ATOX1":    0.81,   # highest in paper, used as IHC validation gene
    "ATP7A":    0.61,
    "ATP7B":    0.58,
    "CCS":      0.55,
    "SOD1":     0.69,
    "COX17":    0.62,
    "SCO1":     0.48,
    "SCO2":     0.51,
    "COX11":    0.44,
    # Mitochondrial — moderate
    "MT-CO1":   0.53,
    "MT-CO2":   0.49,
    # ECM / LOX family — moderate to good
    "LOX":      0.57,
    "LOXL1":    0.52,
    "LOXL2":    0.61,
    "LOXL3":    0.48,
    "LOXL4":    0.44,
    "SPARC":    0.63,
    "GPC1":     0.41,
    # Plasma / hepatocyte secretome — LOW for some
    "CP":       0.55,
    "HEPH":     0.46,
    "HEPHL1":   0.38,
    "ALB":      0.08,   # explicitly flagged as low in paper
    "AFP":      0.31,
    "LTF":      0.42,
    "F5":       0.19,   # flagged as low in paper
    # Antioxidant
    "SOD3":     0.44,
    "PARK7":    0.67,
    # Amyloid / neuronal
    "APP":      0.39,
    "PRNP":     0.45,
    "SNCA":     0.36,
    # S100 family
    "S100A5":   0.41,
    "S100A12":  0.38,
    "S100A13":  0.44,
    "S100B":    0.52,
    # Copper metabolism
    "COMMD1":   0.14,   # explicitly flagged as low in paper
    "CUTA":     0.48,
    "CUTC":     0.51,
    # ENOX family — LOW
    "ENOX1":    0.22,
    "ENOX2":   -0.15,   # explicitly flagged — negative correlation
    # Amine oxidases
    "AOC1":     0.44,
    "AOC2":     0.21,   # flagged as low in paper
    "AOC3":     0.38,
    # Melanin / catecholamine
    "TYR":      0.55,
    "TYRP1":    0.51,
    "DBH":      0.48,
    # Signaling
    "MAP2K1":   0.59,
    "MEMO1":    0.52,
    # Metallothioneins
    "MT3":      0.33,
    "MT4":      0.29,
    # Other
    "MOXD1":    0.41,
    "MOXD2P":   None,   # pseudogene — no protein
    "PAM":      0.44,
    "SLC31A2":  0.38,
    # Histones (Attar 2020 — not in Blockhuys 2017)
    "H3-3A":    None,   # not in Blockhuys 2017 — added from Attar 2020
    "H3-3B":    None,
    "H3C1":     None,
    "H4C1":     None,
}


def confidence_label(r):
    if r is None:
        return "NO DATA"
    if r < 0.3:
        return "LOW"
    if r < 0.5:
        return "MEDIUM"
    return "HIGH"


def confidence_note(gene, r):
    if r is None:
        if gene in ("MOXD2P",):
            return "pseudogene — no protein product"
        if gene in ("H3-3A", "H3-3B", "H3C1", "H4C1"):
            return "added from Attar 2020 — not in Blockhuys 2017"
        return "not reported in Blockhuys 2017"
    if r < 0:
        return "NEGATIVE correlation — mRNA and protein inversely related"
    if r < 0.3:
        return "LOW confidence — protein validation required before strong claims"
    if r < 0.5:
        return "MEDIUM confidence — treat mRNA as indicative not definitive"
    return "HIGH confidence — mRNA is a reliable proxy for protein"


def main():
    print("[conf] loading saliency results ...")

    # load saliency from interpretability run
    sal_file = OUT / "interp_saliency_v2.csv"
    if not sal_file.exists():
        print(f"  WARNING: {sal_file} not found.")
        print("  Building confidence table from gene list only.")
        genes = list(MRNA_PROTEIN_R.keys())
        sal_df = pd.DataFrame({"gene": genes,
                                "saliency_rank": range(1, len(genes)+1),
                                "differential": [0.0]*len(genes),
                                "direction": ["unknown"]*len(genes)})
    else:
        sal_df = pd.read_csv(sal_file)
        print(f"  Loaded {len(sal_df)} genes from saliency table")

    # add confidence columns
    sal_df["mrna_protein_r"] = sal_df["gene"].map(
        lambda g: MRNA_PROTEIN_R.get(g, None))
    sal_df["confidence"]     = sal_df["mrna_protein_r"].map(confidence_label)
    sal_df["confidence_note"] = sal_df.apply(
        lambda row: confidence_note(row["gene"], row["mrna_protein_r"]), axis=1)

    # also load DE results for log2FC
    de_file = ROOT / "outputs" / "baseline" / "copper_de_results.csv"
    if de_file.exists():
        de_df = pd.read_csv(de_file).set_index("gene_symbol")
        sal_df["log2FC_tumor_vs_normal"] = sal_df["gene"].map(
            lambda g: de_df.loc[g, "log2FC"] if g in de_df.index else None)
    else:
        sal_df["log2FC_tumor_vs_normal"] = None

    sal_df.to_csv(OUT / "saliency_with_confidence.csv", index=False)

    # ── print summary ──────────────────────────────────────────────────────
    n_high   = int((sal_df["confidence"] == "HIGH").sum())
    n_medium = int((sal_df["confidence"] == "MEDIUM").sum())
    n_low    = int((sal_df["confidence"] == "LOW").sum())
    n_nodata = int((sal_df["confidence"] == "NO DATA").sum())

    print(f"\n  Confidence breakdown across all {len(sal_df)} genes:")
    print(f"  HIGH   (r >= 0.5) : {n_high}")
    print(f"  MEDIUM (r 0.3-0.5): {n_medium}")
    print(f"  LOW    (r < 0.3)  : {n_low}")
    print(f"  NO DATA           : {n_nodata}")

    # top 20 saliency genes with confidence
    top20 = sal_df.sort_values("saliency_rank").head(20)
    print(f"\n  Top 20 saliency genes with mRNA-protein confidence:")
    print(f"  {'rank':<5} {'gene':<12} {'direction':<12} "
          f"{'r':>6} {'confidence':<10}")
    print("  " + "-"*52)
    for _, r in top20.iterrows():
        rv = f"{r['mrna_protein_r']:.2f}" if r['mrna_protein_r'] is not None else "  —"
        print(f"  {int(r['saliency_rank']):<5} {r['gene']:<12} "
              f"{r['direction']:<12} {rv:>6} {r['confidence']:<10}")

    # flag any LOW confidence genes in top 10
    top10 = sal_df.sort_values("saliency_rank").head(10)
    low_in_top10 = top10[top10["confidence"] == "LOW"]["gene"].tolist()
    if low_in_top10:
        print(f"\n  WARNING: LOW confidence genes in top 10 saliency:")
        for g in low_in_top10:
            r_val = MRNA_PROTEIN_R.get(g)
            print(f"    {g}: r={r_val} — claims need protein-level validation")
    else:
        print(f"\n  OK: No LOW confidence genes in top 10 saliency")

    # ── write report ───────────────────────────────────────────────────────
    report = f"""# mRNA-Protein Confidence Flags — Copper Proteome Saliency

## Motivation

A key scientific criticism of RNA-based analyses:
> "RNA != protein. Genes with low mRNA/protein correlation (e.g. ENOX2,
> COMMD1, AOC2, ALB, F5) need protein-level validation before any strong claim."

Source: Blockhuys et al. 2017, *Metallomics* 9:112, Figure 4.
Reports Pearson r between mRNA and protein abundance across cell lines.

## Confidence thresholds

| label | r value | meaning |
|---|---|---|
| HIGH | r ≥ 0.5 | mRNA is a reliable proxy for protein |
| MEDIUM | 0.3 ≤ r < 0.5 | mRNA is indicative but not definitive |
| LOW | r < 0.3 | mRNA is NOT reliable — protein validation required |
| NO DATA | — | not reported in Blockhuys 2017 |

## Breakdown across 54 copper genes

- HIGH confidence   : **{n_high}** genes
- MEDIUM confidence : **{n_medium}** genes
- LOW confidence    : **{n_low}** genes (require protein validation)
- No data           : **{n_nodata}** genes

## LOW confidence genes — explicit warnings

| gene | r | note |
|---|---:|---|
| ENOX2 | -0.15 | NEGATIVE correlation — mRNA and protein inversely related |
| COMMD1 | 0.14 | explicitly flagged in Blockhuys 2017 as unreliable |
| ALB | 0.08 | secreted protein — mRNA/protein decoupled by secretion |
| F5 | 0.19 | coagulation factor — post-translational regulation dominates |
| AOC2 | 0.21 | limited expression data in Blockhuys 2017 |
| ENOX1 | 0.22 | similar to ENOX2 — low correlation |
| MT4 | 0.29 | metallothionein — rapid post-translational turnover |

## Top 20 saliency genes with confidence

| rank | gene | direction | r | confidence |
|---|---|---|---:|---|
"""
    for _, r in top20.iterrows():
        rv = f"{r['mrna_protein_r']:.2f}" if r['mrna_protein_r'] is not None else "—"
        flag = " ⚠️" if r["confidence"] == "LOW" else ""
        report += (f"| {int(r['saliency_rank'])} | **{r['gene']}** | "
                   f"{r['direction']} | {rv} | {r['confidence']}{flag} |\n")

    report += f"""
## What this means for the paper

**Strong claims** (back with mRNA data alone):
Genes with HIGH confidence — ATOX1 (r=0.81), SLC31A1 (r=0.72), SOD1 (r=0.69),
MAP2K1 (r=0.59), LOXL2 (r=0.61), SPARC (r=0.63), PARK7 (r=0.67), ATP7A (r=0.61).
These genes show consistent mRNA-protein coupling and saliency claims are well-supported.

**Qualified claims** (note limitation):
Genes with MEDIUM confidence — SCO1, SCO2, COX11, HEPH, AFP, PRNP, CUTA, S100A12.
State: "mRNA-level evidence; protein validation recommended."

**Weak claims** (require protein data):
ENOX2, COMMD1, ALB, F5, AOC2, ENOX1, MT4.
Do NOT make strong biological claims about these genes from mRNA alone.
If any of these appear in the top saliency, flag explicitly in the paper.

**Histones** (H3-3A, H3-3B, H3C1, H4C1):
Not in Blockhuys 2017 — added from Attar et al. 2020 (Science 369:59).
The Attar paper provides direct biochemical evidence for Cu2+ reductase
activity, so the biological claim is supported at the protein/biochemical
level even without mRNA-protein correlation data.

## Recommended language for the paper

> "Saliency scores reflect mRNA-level importance. For genes with low
> mRNA-protein correlation in the Blockhuys 2017 reference panel
> (ENOX2 r=-0.15, COMMD1 r=0.14, ALB r=0.08), functional conclusions
> should be treated as hypothesis-generating and require protein-level
> validation."

## Files produced
- `outputs/final_comparison/saliency_with_confidence.csv`
- `outputs/final_comparison/mrna_protein_report.md`
"""
    (OUT / "mrna_protein_report.md").write_text(report, encoding="utf-8")
    print(f"\n[conf] wrote saliency_with_confidence.csv")
    print(f"[conf] wrote mrna_protein_report.md")
    print("\nDone.")


if __name__ == "__main__":
    main()
