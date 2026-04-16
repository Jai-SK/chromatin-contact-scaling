#!/usr/bin/env python3
"""
Figure 2: Box plot of cell line alpha by cancer type (cancer lines only)
UPDATED: GBM now n=5 (GBM_P4CC3 + 4 GSC lines from GSE280831)

Reads: cellline_alpha_results.csv, DISSERTATION_master_alpha_table.csv,
       gbm_cellline_alpha_results.csv
Output: fig2_boxplot_by_type.png / .pdf

Usage:
    cd /mnt/c/Users/jaipa/OneDrive/Documents/Dissertation
    python fig2_boxplot_by_type.py
"""

import os, platform
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# === Paths ===
if platform.system() == "Linux":
    ROOT = "/mnt/c/Users/jaipa/OneDrive/Documents/Dissertation"
else:
    ROOT = r"C:\Users\jaipa\OneDrive\Documents\Dissertation"

CL_NEW_FILE = os.path.join(ROOT, "cellline_alpha_results.csv")
MASTER_FILE = os.path.join(ROOT, "DISSERTATION_master_alpha_table.csv")
GBM_CL_FILE = os.path.join(ROOT, "gbm_cellline_alpha_results.csv")

# === Load ===
cl_new = pd.read_csv(CL_NEW_FILE)
master = pd.read_csv(MASTER_FILE)
old_cl = master[master["Dataset"] != "TCGA"].copy()
gbm_cl = pd.read_csv(GBM_CL_FILE)

# === Build combined — CANCER LINES ONLY ===
old_cancer = {
    "K562": "CML", "HepG2": "HCC", "HCT116": "CRC",
    "BxPC3": "PDAC", "PANC1": "PDAC",
}
old_rows = []
for _, r in old_cl.iterrows():
    if r["Sample"] in old_cancer:
        old_rows.append({
            "sample": r["Sample"],
            "alpha_neg": r["Alpha"],
            "cancer_type": old_cancer[r["Sample"]],
        })
old_df = pd.DataFrame(old_rows)

new_df = pd.DataFrame({
    "sample": cl_new["sample"].values,
    "alpha_neg": -cl_new["alpha"].values,
    "cancer_type": cl_new["cancer_type"].values,
})

gbm_df = pd.DataFrame({
    "sample": gbm_cl["sample"].values,
    "alpha_neg": -gbm_cl["alpha"].values,
    "cancer_type": ["GBM"] * len(gbm_cl),
})

cl_all = pd.concat([old_df, new_df, gbm_df], ignore_index=True)

print(f"Total cancer cell lines: {len(cl_all)}")
print(cl_all.groupby("cancer_type").size().sort_values(ascending=False).to_string())

# === Sort by median ===
type_medians = cl_all.groupby("cancer_type")["alpha_neg"].median().sort_values()
ordered = list(type_medians.index)

# === Colours ===
TYPE_COLS = {
    "BLCA": "#C62828", "BRCA": "#AD1457", "COAD": "#6A1B9A", "GBM": "#1565C0",
    "LIHC": "#2E7D32", "LUAD": "#00838F", "PRAD": "#4E342E", "SKCM": "#E65100",
    "PDAC": "#8E24AA", "CML": "#546E7A", "HCC": "#388E3C", "CRC": "#7B1FA2",
}

# === Plot ===
plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 11,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 300, "savefig.dpi": 300, "savefig.bbox": "tight",
})

fig, ax = plt.subplots(figsize=(10, 6))
rng = np.random.default_rng(42)

for pos, ct in enumerate(ordered):
    vals = cl_all.loc[cl_all["cancer_type"] == ct, "alpha_neg"].values
    col = TYPE_COLS.get(ct, "#607D8B")
    n = len(vals)

    if n >= 2:
        bp = ax.boxplot([vals], positions=[pos], widths=0.55, patch_artist=True,
                        showfliers=False,
                        medianprops=dict(color="black", linewidth=2),
                        whiskerprops=dict(color=col, linewidth=1.2),
                        capprops=dict(color=col, linewidth=1.2))
        bp["boxes"][0].set_facecolor(col)
        bp["boxes"][0].set_alpha(0.5)
        bp["boxes"][0].set_edgecolor(col)
        jit = rng.uniform(-0.15, 0.15, n)
        ax.scatter(pos + jit, vals, color=col, s=35, alpha=0.85,
                   edgecolors="white", linewidth=0.5, zorder=5)
    else:
        ax.scatter(pos, vals[0], color=col, s=80, alpha=0.9,
                   edgecolors="black", linewidth=1.2, zorder=5, marker="D")

    y_top = max(vals) + 0.025
    ax.text(pos, y_top, f"n={n}", ha="center", fontsize=8.5, fontweight="bold", color="#333")

ax.axhline(-1.0, color="gray", ls="--", lw=1, alpha=0.5)
ax.text(len(ordered) - 0.5, -0.98, r"$\alpha = -1$", fontsize=9, color="gray", ha="right")

ax.set_xticks(range(len(ordered)))
ax.set_xticklabels(ordered, rotation=40, ha="right", fontsize=10, fontweight="bold")
ax.set_ylabel(r"Power-law exponent  $\alpha$", fontsize=13)
ax.set_title("Cell line scaling exponents by cancer type", fontsize=14, fontweight="bold")
ax.set_xlim(-0.6, len(ordered) - 0.4)

fig.tight_layout()
for ext in ("png", "pdf"):
    out = os.path.join(ROOT, f"fig2_boxplot_by_type.{ext}")
    fig.savefig(out)
    print(f"Saved: {out}")
plt.close()
