#!/usr/bin/env python3
"""
Figure 1: Violin plot — TCGA tumours (n=48) vs cell lines (n=44)
UPDATED: includes 4 new GSC cell lines from GSE280831

Reads: TCGA_powerlaw_results.csv, cellline_alpha_results.csv,
       DISSERTATION_master_alpha_table.csv, gbm_cellline_alpha_results.csv
Output: fig1_violin.png, fig1_violin.pdf

Usage:
    cd /mnt/c/Users/jaipa/OneDrive/Documents/Dissertation
    python fig1_violin.py
"""

import os, platform
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import stats

# === Paths ===
if platform.system() == "Linux":
    ROOT = "/mnt/c/Users/jaipa/OneDrive/Documents/Dissertation"
else:
    ROOT = r"C:\Users\jaipa\OneDrive\Documents\Dissertation"

TCGA_FILE   = os.path.join(ROOT, "TCGA_HiChIP_hic", "TCGA_powerlaw_results.csv")
CL_NEW_FILE = os.path.join(ROOT, "cellline_alpha_results.csv")
MASTER_FILE = os.path.join(ROOT, "DISSERTATION_master_alpha_table.csv")
GBM_CL_FILE = os.path.join(ROOT, "gbm_cellline_alpha_results.csv")

# === Load ===
tcga_raw  = pd.read_csv(TCGA_FILE)
cl_new    = pd.read_csv(CL_NEW_FILE)
master    = pd.read_csv(MASTER_FILE)
gbm_cl    = pd.read_csv(GBM_CL_FILE)

# TCGA: negate (stored positive)
tcga_alpha = -tcga_raw["Alpha"].values

# Old cell lines (11): already negative
old_cl = master[master["Dataset"] != "TCGA"]
old_alpha = old_cl["Alpha"].values
old_sources = old_cl["Dataset"].values

# New cell lines (29): negate
new_alpha = -cl_new["alpha"].values
new_sources = cl_new["source"].values

# GBM cell lines (4): negate
gbm_alpha = -gbm_cl["alpha"].values
gbm_sources = np.array(["GSE280831"] * len(gbm_cl))

# Combined
cl_alpha   = np.concatenate([old_alpha, new_alpha, gbm_alpha])
cl_sources = np.concatenate([old_sources, new_sources, gbm_sources])

# Source categories for colouring
def src_cat(s):
    if "ENCODE" in s or "GSE63525" in s or "4DN" in s:
        return "ENCODE"
    elif "GSE253407" in s:
        return "GSE253407"
    elif "GSE185069" in s:
        return "GSE185069"
    elif "GSE280831" in s:
        return "GSE280831"
    else:
        return "Other GEO"

cl_cats = [src_cat(s) for s in cl_sources]

# === Stats ===
mwu_stat, mwu_p = stats.mannwhitneyu(tcga_alpha, cl_alpha, alternative="two-sided")
welch_t, welch_p = stats.ttest_ind(tcga_alpha, cl_alpha, equal_var=False)
pooled_sd = np.sqrt((np.std(tcga_alpha, ddof=1)**2 + np.std(cl_alpha, ddof=1)**2) / 2)
cohens_d = (np.mean(tcga_alpha) - np.mean(cl_alpha)) / pooled_sd

print(f"TCGA: n={len(tcga_alpha)}, mean={np.mean(tcga_alpha):.4f}")
print(f"CL:   n={len(cl_alpha)}, mean={np.mean(cl_alpha):.4f}")
print(f"MW-U p={mwu_p:.4e}, Welch t={welch_t:.2f}, Cohen d={cohens_d:.2f}")

# === Colours ===
TCGA_COL = "#E07B54"
CL_COL   = "#6BA3D6"
CAT_COLS = {"ENCODE": "#C0392B", "GSE253407": "#E67E22", "GSE185069": "#8E44AD",
            "GSE280831": "#2980B9", "Other GEO": "#C0392B"}

# === Plot ===
plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 11,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 300, "savefig.dpi": 300, "savefig.bbox": "tight",
})

fig, ax = plt.subplots(figsize=(6, 7))

# Violins
vp1 = ax.violinplot([tcga_alpha], positions=[0], showmedians=False, showextrema=False)
vp2 = ax.violinplot([cl_alpha], positions=[1], showmedians=False, showextrema=False)
for b in vp1["bodies"]:
    b.set_facecolor(TCGA_COL); b.set_alpha(0.3); b.set_edgecolor(TCGA_COL); b.set_linewidth(1.5)
for b in vp2["bodies"]:
    b.set_facecolor(CL_COL); b.set_alpha(0.3); b.set_edgecolor(CL_COL); b.set_linewidth(1.5)

# Medians
ax.hlines(np.median(tcga_alpha), -0.35, 0.35, colors=TCGA_COL, linewidth=2.5)
ax.hlines(np.median(cl_alpha), 0.65, 1.35, colors=CL_COL, linewidth=2.5)

# Jittered dots
rng = np.random.default_rng(42)
ax.scatter(rng.uniform(-0.25, 0.25, len(tcga_alpha)), tcga_alpha,
           color=TCGA_COL, alpha=0.55, s=30, edgecolors="white", linewidth=0.5, zorder=5)

for i, (a, cat) in enumerate(zip(cl_alpha, cl_cats)):
    ax.scatter(1 + rng.uniform(-0.25, 0.25), a,
               color=CAT_COLS[cat], alpha=0.65, s=35, edgecolors="white", linewidth=0.5, zorder=5)

# Significance bracket
sig = "***" if mwu_p < 0.001 else ("**" if mwu_p < 0.01 else ("*" if mwu_p < 0.05 else "ns"))
ybar = max(tcga_alpha.max(), cl_alpha.max()) + 0.05
ax.plot([0, 0, 1, 1], [ybar, ybar + 0.03, ybar + 0.03, ybar], color="k", lw=1.2)
ax.text(0.5, ybar + 0.04, f"{sig}  p = {mwu_p:.1e}", ha="center", fontsize=11, fontweight="bold")

# Stats inset
ax.text(0.98, 0.95,
        f"Cohen's d = {cohens_d:.2f}\nWelch's t = {welch_t:.2f}\nMW-U p = {mwu_p:.2e}",
        transform=ax.transAxes, fontsize=9, va="top", ha="right",
        bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="gray", alpha=0.85))

# Source legend
cat_order = ["ENCODE", "GSE253407", "GSE185069", "GSE280831", "Other GEO"]
cat_counts = {c: sum(1 for x in cl_cats if x == c) for c in cat_order}
legend_el = [
    Line2D([0], [0], marker="o", color="w", markerfacecolor=CAT_COLS[c], markersize=7,
           label=f"{c} (n={cat_counts[c]})")
    for c in cat_order if cat_counts[c] > 0
]
ax.legend(handles=legend_el, loc="lower left", fontsize=8, framealpha=0.9)

# Labels
ax.set_xticks([0, 1])
ax.set_xticklabels([f"TCGA\nTumours\nn={len(tcga_alpha)}", f"Cell\nLines\nn={len(cl_alpha)}"],
                    fontsize=12, fontweight="bold")
ax.set_ylabel(r"Power-law exponent  $\alpha$", fontsize=13)
ax.set_title("Tumours vs cell lines", fontsize=15, fontweight="bold", pad=15)
ax.set_xlim(-0.7, 1.7)

fig.tight_layout()
for ext in ("png", "pdf"):
    out = os.path.join(ROOT, f"fig1_violin.{ext}")
    fig.savefig(out)
    print(f"Saved: {out}")
plt.close()
