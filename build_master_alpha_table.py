#!/usr/bin/env python3
"""
build_master_alpha_table.py
============================
Generates: DISSERTATION_master_alpha_table.csv

This script compiles the master alpha table from three sources:

  1. TCGA HiChIP tumours — per-cancer-type summary (mean ± SD)
     Source: TCGA_HiChIP_hic/TCGA_powerlaw_results.csv
     (individual sample α values fitted via process_tcga_hichip.py)

  2. GSE185069 — pancreatic cell lines (BxPC3, PANC1) and normal (HPDE6C7)
     Source: GSE185069_analysis/cancer_powerlaw_results.csv
     (per-chromosome α values from Du et al. 2022 pipeline, averaged here)

  3. ENCODE / 4DN cell lines — 8 reference cell lines
     Source: published α values from Rao et al. 2014 (GM12878, K562, etc.)
     These were extracted from the original publications and GSE63525 data.

NOTE ON DISSERTATION_master_alpha_table.csv CONTENTS:
  The table has 14 TCGA tumour rows (per cancer type, with N = number of
  individual patients), 3 GSE185069 rows, and 8 ENCODE cell line rows.
  Alpha convention: NEGATIVE (e.g. -1.85 = steeper decay).
  Std column: SD of α across individual samples (TCGA) or chromosomes (GSE185069).

  The cellline_alpha_results.csv and gbm_cellline_alpha_results.csv contain
  the newer 29+4 cell line samples processed by the standardised cooltools /
  hicstraw pipeline. These are NOT in DISSERTATION_master_alpha_table.csv —
  they are used separately in the figure scripts.

Requirements:
  pip install numpy pandas scipy

Data locations:
  TCGA_HiChIP_hic/TCGA_powerlaw_results.csv  (run process_tcga_hichip.py first)
  GSE185069_analysis/cancer_powerlaw_results.csv

Usage:
  cd /mnt/c/Users/jaipa/OneDrive/Documents/Dissertation
  python build_master_alpha_table.py
"""

import os
import platform
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

# ── Paths ──────────────────────────────────────────────────────────────────────
if platform.system() == "Linux":
    ROOT = "/mnt/c/Users/jaipa/OneDrive/Documents/Dissertation"
else:
    ROOT = r"C:\Users\jaipa\OneDrive\Documents\Dissertation"

def P(*parts):
    return os.path.join(ROOT, *parts)


# ── Section 1: TCGA tumours ────────────────────────────────────────────────────
def load_tcga_summary():
    """
    Load individual TCGA sample α values and summarise per cancer type.
    Input file produced by process_tcga_hichip.py.
    Columns expected: Cancer_Type, Sample, Alpha (negative), Assay, N
    """
    path = P("TCGA_HiChIP_hic", "TCGA_powerlaw_results.csv")
    if not os.path.exists(path):
        print(f"  MISSING: {path}")
        print("  Run process_tcga_hichip.py first.")
        return pd.DataFrame()

    df = pd.read_csv(path)
    print(f"  TCGA: {len(df)} individual samples")

    # Standardise column names
    df.columns = [c.strip() for c in df.columns]
    col_map = {c: c for c in df.columns}
    for c in df.columns:
        if c.lower() in ("cancer_type","cancertype","cancer"):
            col_map[c] = "Cancer_Type"
        elif c.lower() == "alpha":
            col_map[c] = "Alpha"
    df = df.rename(columns=col_map)

    # Ensure Alpha is negative convention
    if df["Alpha"].mean() > 0:
        df["Alpha"] = -df["Alpha"]

    # Summarise per cancer type
    summary = df.groupby("Cancer_Type")["Alpha"].agg(
        Alpha_mean="mean",
        Std="std",
        N="count",
    ).reset_index()

    rows = []
    for _, row in summary.iterrows():
        rows.append({
            "Dataset": "TCGA",
            "Sample":  row["Cancer_Type"],
            "Assay":   "HiChIP (H3K27ac)",
            "Alpha":   row["Alpha_mean"],
            "Std":     row["Std"],
            "N":       int(row["N"]),
            "Type":    "Tumor",
        })

    print(f"  TCGA summary: {len(rows)} cancer types")
    return pd.DataFrame(rows)


# ── Section 2: GSE185069 (pancreatic) ─────────────────────────────────────────
def load_gse185069():
    """
    Per-chromosome α values from Du et al. 2022 (GSE185069).
    File: cancer_powerlaw_results.csv
    Columns: sample, chromosome, alpha, r_squared  (alpha negative)
    Averaged across chromosomes per cell line.
    """
    path = P("GSE185069_analysis", "cancer_powerlaw_results.csv")
    if not os.path.exists(path):
        print(f"  MISSING: {path}")
        return pd.DataFrame()

    df = pd.read_csv(path)
    print(f"  GSE185069: {len(df)} chromosome-level entries")

    # Standardise
    df.columns = [c.strip() for c in df.columns]
    rename = {}
    for c in df.columns:
        if c.lower() in ("sample","cell_line","name"):
            rename[c] = "sample"
        elif c.lower() in ("alpha","exponent","slope"):
            rename[c] = "alpha"
    df = df.rename(columns=rename)

    if df["alpha"].mean() > 0:
        df["alpha"] = -df["alpha"]

    agg = df.groupby("sample")["alpha"].agg(
        Alpha_mean="mean",
        Std="std",
        N="count",
    ).reset_index()

    type_map = {
        "BxPC3":   "Pancreatic Cancer",
        "PANC1":   "Pancreatic Cancer",
        "HPDE6C7": "Normal",
    }

    rows = []
    for _, row in agg.iterrows():
        rows.append({
            "Dataset": "GSE185069",
            "Sample":  row["sample"],
            "Assay":   "Hi-C",
            "Alpha":   row["Alpha_mean"],
            "Std":     row["Std"],
            "N":       int(row["N"]),
            "Type":    type_map.get(row["sample"], "Cell Line"),
        })

    print(f"  GSE185069 summary: {len(rows)} cell lines")
    return pd.DataFrame(rows)


# ── Section 3: ENCODE / 4DN reference cell lines ───────────────────────────────
def encode_reference_lines():
    """
    Published α values for 8 ENCODE/4DN reference cell lines.
    Source: Rao et al. 2014 (Cell) and 4DN data portal.
    Dataset GSE63525. Values are mean ± SD across replicates / chromosomes.
    These are the same values present in DISSERTATION_master_alpha_table.csv.
    Alpha convention: negative (e.g. -0.929 for IMR90).
    """
    encode_data = [
        # (Sample, Alpha_mean, Std, N_replicates)
        ("IMR90",   -0.9294131512158312,  0.0842292760228974, 3),
        ("HCT116",  -0.8613007561792900,  0.0,                1),
        ("GM12878", -0.8250899165801640,  0.1215530226584147, 5),
        ("K562",    -0.7939134276231676,  0.0773603361807752, 4),
        ("HepG2",   -0.7874495575141220,  0.0974634588094820, 5),
        ("HUVEC",   -0.7598507592719798,  0.0673618843322495, 3),
        ("HMEC",    -0.6630133683055421,  0.0079745511280623, 2),
        ("NHEK",    -0.6365191826070863,  0.0763127611274711, 5),
    ]

    rows = []
    for sample, alpha, std, n in encode_data:
        rows.append({
            "Dataset": "GSE63525/ENCODE/4DN",
            "Sample":  sample,
            "Assay":   "Hi-C",
            "Alpha":   alpha,
            "Std":     std,
            "N":       n,
            "Type":    "Cell Line",
        })

    print(f"  ENCODE/4DN: {len(rows)} reference cell lines (published values)")
    return pd.DataFrame(rows)


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    parts = []

    print("\nLoading TCGA tumour data...")
    tcga = load_tcga_summary()
    if not tcga.empty:
        parts.append(tcga)

    print("\nLoading GSE185069 pancreatic data...")
    gse185 = load_gse185069()
    if not gse185.empty:
        parts.append(gse185)

    print("\nLoading ENCODE reference cell lines...")
    encode = encode_reference_lines()
    parts.append(encode)

    if not parts:
        print("\nERROR: No data loaded.")
        return

    master = pd.concat(parts, ignore_index=True)

    print(f"\n{'='*60}")
    print(f"Master table: {len(master)} rows")
    print(f"  TCGA tumours:   {len(master[master['Type']=='Tumor'])}")
    print(f"  Cancer lines:   {len(master[master['Type']=='Pancreatic Cancer'])}")
    print(f"  Normal lines:   {len(master[master['Type']=='Normal'])}")
    print(f"  Cell lines:     {len(master[master['Type']=='Cell Line'])}")
    print()
    print(master[["Dataset","Sample","Alpha","Std","N","Type"]].to_string(index=False))

    out = P("DISSERTATION_master_alpha_table.csv")
    master.to_csv(out, index=False)
    print(f"\nSaved → {out}")
    print("\nDone.")


if __name__ == "__main__":
    main()
