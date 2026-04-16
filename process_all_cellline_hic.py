#!/usr/bin/env python3
"""
process_all_cellline_hic.py
============================
Generates: cellline_alpha_results.csv

Processes cancer cell line Hi-C datasets across 8 matched cancer types.
Computes P(s) using cooltools expected_cis, log-bins the curve, fits a
power law in 100 kb – 10 Mb by OLS in log-log space, and stores
α = -slope (positive).

Output columns:
  cancer_type, sample, source, genome, resolution, alpha, r_squared, n_points

Verified output (29 samples):
  BLCA   HT1376  α=1.1577  SCABER α=1.4119  RT4 α=1.2056  SW780 α=1.3640
  BRCA   T47D    α=1.3299  ZR7530 α=1.3140  HCC1954 α=1.3546 HCC70 α=1.2395 BT549 α=1.2068
  LUAD   A549_rep1–4  α≈1.28–1.34
  LIHC   Hep3B–SNU449 α≈1.22–1.27
  PRAD   22Rv1   α=1.2082
  SKCM   LC676–LC821  α≈1.01–1.32   (GSE253407, 10 kb .cool)
  GBM    GBM_P4CC3    α=0.8463      (GSE253407, 10 kb .cool)
  COAD   SW480   α=1.3168

Requirements:
  pip install cooler cooltools hic2cool numpy pandas scipy

Data locations (all relative to DATA_ROOT):
  BLCA:  GSE148079/GSE148079_{HT1376,RT4,SCABER,SW780}_HiC.mcool
  BRCA:  GSE167154/{T47D,ZR7530,HCC1954,HCC70,BT549}.hic
  LUAD:  GSE92804/GSE92804_ENCFF{089KBG,378RZT,401ZAN,939ARM}_chromatin_interactions.hic
  LIHC:  GSE226216/GSE226216_{Hep3B,Huh1,Huh7,SNU449}_Inter_Intra_Res10_20_40_100_500kb.hic
  PRAD:  GSE118629_22Rv1_HiC_40k.normalized.matrix.txt  (directly in DATA_ROOT)
  SKCM/GBM: GSE253407_analysis/{LC676,LC677,LC499,LC500,LC731,LC648,LC801,LC821,LC789,GBM_P4CC3}.cool
  COAD:  GSE133928/SW480_allValidPairs.txt

Usage:
  cd /mnt/c/Users/jaipa/OneDrive/Documents/Dissertation
  python process_all_cellline_hic.py
"""

import os
import gzip
import platform
import subprocess
import warnings
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy import stats as sp_stats
import cooler
import cooltools

warnings.filterwarnings("ignore")

# ── Paths ──────────────────────────────────────────────────────────────────────
if platform.system() == "Linux":
    DATA_ROOT = "/mnt/c/Users/jaipa/OneDrive/Documents/Dissertation"
else:
    DATA_ROOT = r"C:\Users\jaipa\OneDrive\Documents\Dissertation"

def P(*parts):
    return os.path.join(DATA_ROOT, *parts)

# ── Parameters ─────────────────────────────────────────────────────────────────
FIT_MIN    = 100_000     # 100 kb
FIT_MAX    = 10_000_000  # 10 Mb
N_BINS     = 200
MIN_POINTS = 5

# ── Dataset registry ───────────────────────────────────────────────────────────
# (cancer_type, sample, filepath, format, genome, resolution, GSE_source)
DATASETS = [
    # ── BLCA ─ GSE148079 ─ .mcool ─ hg19 ─ 40 kb ──────────────────────────────
    ("BLCA","HT1376", P("GSE148079","GSE148079_HT1376_HiC.mcool"), "mcool","hg19",40_000,"GSE148079"),
    ("BLCA","RT4",    P("GSE148079","GSE148079_RT4_HiC.mcool"),    "mcool","hg19",40_000,"GSE148079"),
    ("BLCA","SCABER", P("GSE148079","GSE148079_SCABER_HiC.mcool"), "mcool","hg19",40_000,"GSE148079"),
    ("BLCA","SW780",  P("GSE148079","GSE148079_SW780_HiC.mcool"),  "mcool","hg19",40_000,"GSE148079"),

    # ── BRCA ─ GSE167154 ─ .hic ─ hg38 ─ 40 kb ───────────────────────────────
    ("BRCA","T47D",    P("GSE167154","T47D.hic"),    "hic","hg38",40_000,"GSE167154"),
    ("BRCA","ZR7530",  P("GSE167154","ZR7530.hic"),  "hic","hg38",40_000,"GSE167154"),
    ("BRCA","HCC1954", P("GSE167154","HCC1954.hic"), "hic","hg38",40_000,"GSE167154"),
    ("BRCA","HCC70",   P("GSE167154","HCC70.hic"),   "hic","hg38",40_000,"GSE167154"),
    ("BRCA","BT549",   P("GSE167154","BT549.hic"),   "hic","hg38",40_000,"GSE167154"),

    # ── LUAD ─ GSE92804 ─ .hic ─ hg38 ─ 50 kb (native ENCODE resolution) ────
    ("LUAD","A549_rep1",P("GSE92804","GSE92804_ENCFF089KBG_chromatin_interactions.hic"),"hic","hg38",50_000,"GSE92804"),
    ("LUAD","A549_rep2",P("GSE92804","GSE92804_ENCFF378RZT_chromatin_interactions.hic"),"hic","hg38",50_000,"GSE92804"),
    ("LUAD","A549_rep3",P("GSE92804","GSE92804_ENCFF401ZAN_chromatin_interactions.hic"),"hic","hg38",50_000,"GSE92804"),
    ("LUAD","A549_rep4",P("GSE92804","GSE92804_ENCFF939ARM_chromatin_interactions.hic"),"hic","hg38",50_000,"GSE92804"),

    # ── LIHC ─ GSE226216 ─ .hic ─ hg38 ─ 40 kb ───────────────────────────────
    ("LIHC","Hep3B", P("GSE226216","GSE226216_Hep3B_Inter_Intra_Res10_20_40_100_500kb.hic"), "hic","hg38",40_000,"GSE226216"),
    ("LIHC","Huh1",  P("GSE226216","GSE226216_Huh1_Inter_Intra_Res10_20_40_100_500kb.hic"),  "hic","hg38",40_000,"GSE226216"),
    ("LIHC","Huh7",  P("GSE226216","GSE226216_Huh7_Inter_Intra_Res10_20_40_100_500kb.hic"),  "hic","hg38",40_000,"GSE226216"),
    ("LIHC","SNU449",P("GSE226216","GSE226216_SNU449_Inter_Intra_Res10_20_40_100_500kb.hic"),"hic","hg38",40_000,"GSE226216"),

    # ── PRAD ─ GSE118629 ─ ICE-normalised sparse .txt ─ hg19 ─ 40 kb ─────────
    ("PRAD","22Rv1",P("GSE118629_22Rv1_HiC_40k.normalized.matrix.txt"),"txt_ice","hg19",40_000,"GSE118629"),

    # ── SKCM ─ GSE253407 ─ .cool ─ hg38 ─ 10 kb ──────────────────────────────
    ("SKCM","LC676",P("GSE253407_analysis","LC676.cool"),"cool","hg38",10_000,"GSE253407"),
    ("SKCM","LC677",P("GSE253407_analysis","LC677.cool"),"cool","hg38",10_000,"GSE253407"),
    ("SKCM","LC499",P("GSE253407_analysis","LC499.cool"),"cool","hg38",10_000,"GSE253407"),
    ("SKCM","LC500",P("GSE253407_analysis","LC500.cool"),"cool","hg38",10_000,"GSE253407"),
    ("SKCM","LC731",P("GSE253407_analysis","LC731.cool"),"cool","hg38",10_000,"GSE253407"),
    ("SKCM","LC648",P("GSE253407_analysis","LC648.cool"),"cool","hg38",10_000,"GSE253407"),
    ("SKCM","LC801",P("GSE253407_analysis","LC801.cool"),"cool","hg38",10_000,"GSE253407"),
    ("SKCM","LC821",P("GSE253407_analysis","LC821.cool"),"cool","hg38",10_000,"GSE253407"),
    ("SKCM","LC789",P("GSE253407_analysis","LC789.cool"),"cool","hg38",10_000,"GSE253407"),

    # ── GBM ─ GSE253407 ─ .cool ─ hg38 ─ 10 kb ───────────────────────────────
    ("GBM","GBM_P4CC3",P("GSE253407_analysis","GBM_P4CC3.cool"),"cool","hg38",10_000,"GSE253407"),

    # ── COAD ─ GSE133928 ─ allValidPairs ─ hg19 ─ 40 kb ──────────────────────
    ("COAD","SW480",P("GSE133928","SW480_allValidPairs.txt"),"allvalidpairs","hg19",40_000,"GSE133928"),
]


# ── Core functions ─────────────────────────────────────────────────────────────

def autosome_view(clr):
    cs = clr.chromsizes
    nums = {str(i) for i in range(1, 23)}
    chroms = [c for c in cs.index if c.lstrip("chr") in nums]
    return pd.DataFrame({
        "chrom": chroms,
        "start": 0,
        "end":   [cs[c] for c in chroms],
        "name":  chroms,
    })


def run_cooltools(clr):
    view_df = autosome_view(clr)
    weight  = "weight" if "weight" in clr.bins().columns else None
    exp = cooltools.expected_cis(
        clr=clr,
        view_df=view_df,
        clr_weight_name=weight,
        smooth=True,
        aggregate_smoothed=True,
        nproc=1,
        ignore_diags=2,
    )
    col = next(
        (c for c in ["count.avg.smoothed.agg","count.avg.smoothed","count.avg"]
         if c in exp.columns), None
    )
    if col is None:
        raise ValueError(f"No P(s) column. Available: {list(exp.columns)}")
    agg = exp.groupby("dist")[col].mean().reset_index()
    return agg["dist"].values * clr.binsize, agg[col].values


def log_bin_fit(dist_bp, contacts):
    mask = (dist_bp > 0) & (contacts > 0) & np.isfinite(contacts)
    d, c = dist_bp[mask], contacts[mask]
    if len(d) < MIN_POINTS:
        return None, None, None, 0

    log_d = np.log10(d)
    edges = np.linspace(log_d.min(), log_d.max(), N_BINS + 1)
    bc, bm = [], []
    for i in range(len(edges) - 1):
        sel = (log_d >= edges[i]) & (log_d < edges[i + 1])
        if sel.sum():
            bc.append(10 ** ((edges[i] + edges[i + 1]) / 2))
            bm.append(c[sel].mean())

    bc, bm = np.array(bc), np.array(bm)
    ps_df  = pd.DataFrame({"genomic_distance": bc, "contact_freq": bm})

    fit = (bc >= FIT_MIN) & (bc <= FIT_MAX) & (bm > 0)
    n   = int(fit.sum())
    if n < MIN_POINTS:
        return ps_df, None, None, n

    slope, _, r, _, _ = sp_stats.linregress(np.log10(bc[fit]), np.log10(bm[fit]))
    return ps_df, -slope, r**2, n  # alpha stored positive


# ── Format handlers ────────────────────────────────────────────────────────────

def handle_mcool(fp, res):
    clr = cooler.Cooler(f"{fp}::/resolutions/{res}")
    return run_cooltools(clr)


def handle_cool(fp):
    clr = cooler.Cooler(fp)
    return run_cooltools(clr)


def handle_hic(fp, res):
    cp = fp.replace(".hic", f"_{res // 1000}kb.cool")
    if not os.path.exists(cp):
        print(f"    hic2cool → {res//1000} kb ...", end=" ", flush=True)
        r = subprocess.run(
            ["hic2cool","convert", fp, cp, "-r", str(res)],
            capture_output=True, text=True
        )
        if r.returncode != 0:
            raise RuntimeError(f"hic2cool failed:\n{r.stderr[:400]}")
        print("done")
    else:
        print("    Cached .cool found")
    return handle_cool(cp)


def handle_txt_ice(fp):
    opener = gzip.open if fp.endswith(".gz") else open
    rows, cols, vals = [], [], []
    with opener(fp, "rt") as fh:
        for line in fh:
            p = line.strip().split()
            if len(p) < 3:
                continue
            try:
                rows.append(float(p[0])); cols.append(float(p[1])); vals.append(float(p[2]))
            except ValueError:
                continue
    rows, cols, vals = map(np.array, [rows, cols, vals])
    coords = np.unique(np.concatenate([rows, cols]))
    diffs  = np.diff(np.sort(coords))
    res    = int(diffs[diffs > 0].min()) if (diffs > 0).any() else 40_000
    dist   = np.abs(cols - rows) * res
    ok     = (dist > 0) & (vals > 0)
    return dist[ok], vals[ok]


def handle_allvalidpairs(fp, res):
    opener = gzip.open if fp.endswith(".gz") else open
    autos  = {str(i) for i in range(1,23)} | {f"chr{i}" for i in range(1,23)}
    bins   = defaultdict(float)
    n = 0
    with opener(fp, "rt") as fh:
        for line in fh:
            p = line.strip().split()
            if len(p) < 7:
                continue
            try:
                c1, x1 = p[1], int(p[2])
                c2, x2 = p[5], int(p[6])
            except (ValueError, IndexError):
                continue
            if c1 != c2 or c1 not in autos:
                continue
            d = abs(x2 - x1)
            if d == 0:
                continue
            bins[(d // res) * res] += 1
            n += 1
    print(f"    {n:,} cis pairs")
    if not bins:
        raise ValueError("No cis contacts found")
    dists  = np.array(sorted(bins), dtype=float)
    counts = np.array([bins[int(d)] for d in dists], dtype=float)
    return dists, counts


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    results, all_ps = [], []
    n_ok = n_fail = 0

    for cancer, sample, fp, fmt, genome, res, source in DATASETS:
        print(f"\n{'─'*60}")
        print(f"  {cancer} | {sample} | {source} | {fmt} | {res//1000} kb | {genome}")

        if not os.path.exists(fp):
            print(f"  NOT FOUND — update path in DATASETS: {fp}")
            n_fail += 1
            continue

        try:
            if   fmt == "mcool":        dist_bp, contacts = handle_mcool(fp, res)
            elif fmt == "cool":         dist_bp, contacts = handle_cool(fp)
            elif fmt == "hic":          dist_bp, contacts = handle_hic(fp, res)
            elif fmt == "txt_ice":      dist_bp, contacts = handle_txt_ice(fp)
            elif fmt == "allvalidpairs":dist_bp, contacts = handle_allvalidpairs(fp, res)
            else:
                raise ValueError(f"Unknown format '{fmt}'")

            ps_df, alpha, r2, n_pts = log_bin_fit(
                np.asarray(dist_bp, dtype=float),
                np.asarray(contacts, dtype=float),
            )

            if alpha is None:
                print(f"  Fit failed ({n_pts} bins in fit range)")
                n_fail += 1
                continue

        except Exception as e:
            print(f"  ERROR: {e}")
            n_fail += 1
            continue

        print(f"  α = {alpha:.10f}   R² = {r2:.10f}   n_bins = {n_pts}")
        n_ok += 1

        results.append({
            "cancer_type": cancer,
            "sample":      sample,
            "source":      source,
            "genome":      genome,
            "resolution":  res,
            "alpha":       alpha,
            "r_squared":   r2,
            "n_points":    n_pts,
        })

        if ps_df is not None:
            ps_df["cancer_type"] = cancer
            ps_df["sample"]      = sample
            all_ps.append(ps_df)

    print(f"\n{'='*60}")
    print(f"Done: {n_ok} succeeded, {n_fail} failed / {len(DATASETS)} total")

    if results:
        df  = pd.DataFrame(results)
        out = P("cellline_alpha_results.csv")
        df.to_csv(out, index=False)
        print(f"\nSaved → {out}")
        print(df[["cancer_type","sample","alpha","r_squared","n_points"]].to_string(index=False))

    if all_ps:
        pd.concat(all_ps, ignore_index=True).to_csv(P("cellline_ps_curves.csv"), index=False)
        print(f"Saved → {P('cellline_ps_curves.csv')}")


if __name__ == "__main__":
    main()
