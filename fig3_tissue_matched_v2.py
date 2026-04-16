#!/usr/bin/env python3
"""
Figure 3 UPDATED: Tissue-matched alpha + P(s) curves
- Triangles for cell line means, squares for TCGA
- n=* labels on each column
- GBM now has 5 cell lines (GBM_P4CC3 + 4 GSCs)
Reads: TCGA_powerlaw_results.csv, cellline_alpha_results.csv, 
       gbm_cellline_alpha_results.csv, cellline_ps_curves.csv, gbm_cellline_ps_curves.csv
Output: fig3_tissue_matched_v2.png/.pdf
CHANGES: Panel B y-axis label -> Contact probability, P(s); removed grey axvspan
"""
import os
import platform
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

R = "/mnt/c/Users/jaipa/OneDrive/Documents/Dissertation" if platform.system() == "Linux" else r"C:\Users\jaipa\OneDrive\Documents\Dissertation"

tcga_raw = pd.read_csv(os.path.join(R, "TCGA_HiChIP_hic", "TCGA_powerlaw_results.csv"))

# Load 30 GBM patient tumour samples from GSE229962 and merge with TCGA tumour data
gbm_patient_file = os.path.join(R, "gbm_patient_alpha_results.csv")
if os.path.exists(gbm_patient_file):
    gbm_patient = pd.read_csv(gbm_patient_file)
    # alpha column is stored positive (same convention as TCGA Alpha column)
    gbm_extra = pd.DataFrame({
        "Alpha":       gbm_patient["alpha"].values,
        "Cancer_Type": ["GBMx"] * len(gbm_patient),
    })
    tcga_raw = pd.concat([tcga_raw, gbm_extra], ignore_index=True)
    print(f"Added {len(gbm_patient)} GBM patient samples (GSE229962) to tumour cohort")
    gbm_all = tcga_raw[tcga_raw["Cancer_Type"] == "GBMx"]
    print(f"Total GBM tumour samples now: {len(gbm_all)}, mean alpha = {-gbm_all['Alpha'].mean():.3f}")
else:
    print(f"WARNING: {gbm_patient_file} not found - GBM tumour n will be TCGA only (n=3)")
cl_raw   = pd.read_csv(os.path.join(R, "cellline_alpha_results.csv"))
gbm_cl   = pd.read_csv(os.path.join(R, "gbm_cellline_alpha_results.csv"))
ps_cl    = pd.read_csv(os.path.join(R, "cellline_ps_curves.csv"))
ps_gbm   = pd.read_csv(os.path.join(R, "gbm_cellline_ps_curves.csv"))

tcga_raw["alpha"]    = -tcga_raw["Alpha"]
cl_raw["alpha_neg"]  = -cl_raw["alpha"]
cl_raw["ct_match"]   = cl_raw["cancer_type"].replace({"GBM": "GBMx"})

gbm_add = pd.DataFrame({
    "sample":      gbm_cl["sample"],
    "alpha_neg":   -gbm_cl["alpha"],
    "cancer_type": ["GBM"] * len(gbm_cl),
    "ct_match":    ["GBMx"] * len(gbm_cl),
    "source":      ["GSE280831"] * len(gbm_cl),
})
cl_combined = pd.concat([cl_raw, gbm_add], ignore_index=True)
ps_combined = pd.concat([ps_cl, ps_gbm], ignore_index=True)

matched = sorted(set(tcga_raw["Cancer_Type"]) & set(cl_combined["ct_match"]))

TYPE_COLS = {
    "BLCA": "#C62828", "BRCA": "#AD1457", "COAD": "#6A1B9A", "GBMx": "#1565C0",
    "LIHC": "#2E7D32", "LUAD": "#00838F", "PRAD": "#4E342E", "SKCM": "#E65100",
}
DISPLAY = {
    "BLCA": "Bladder\n(BLCA)", "BRCA": "Breast\n(BRCA)", "COAD": "Colon\n(COAD)",
    "GBMx": "Brain\n(GBM)",    "LIHC": "Liver\n(LIHC)", "LUAD": "Lung Adeno\n(LUAD)",
    "PRAD": "Prostate\n(PRAD)", "SKCM": "Skin\n(SKCM)",
}

per_type = []
for ct in matched:
    t = tcga_raw.loc[tcga_raw["Cancer_Type"] == ct, "alpha"].values
    c = cl_combined.loc[cl_combined["ct_match"] == ct, "alpha_neg"].values
    per_type.append({
        "type":   ct,
        "t_vals": t, "c_vals": c,
        "t_mean": np.mean(t),
        "t_std":  np.std(t, ddof=1) if len(t) > 1 else 0,
        "c_mean": np.mean(c),
        "c_std":  np.std(c, ddof=1) if len(c) > 1 else 0,
        "delta":  np.mean(c) - np.mean(t),
        "t_n":    len(t), "c_n": len(c),
    })

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 11,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 300, "savefig.dpi": 300, "savefig.bbox": "tight",
})

fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(16, 7),
                                  gridspec_kw={"width_ratios": [1, 1.1]})

# ── Panel A — UNCHANGED ──────────────────────────────────────────────────────
rng = np.random.default_rng(123)
w   = 0.35
for i, s in enumerate(per_type):
    col = TYPE_COLS.get(s["type"], "gray")
    jt  = rng.uniform(-0.08, 0.08, len(s["t_vals"]))
    ax_a.scatter(i - w/2 + jt, s["t_vals"], color=col, alpha=0.35, s=25, zorder=4)
    ax_a.errorbar(i - w/2, s["t_mean"], yerr=s["t_std"], fmt="s", color=col,
                  ms=10, mec="white", mew=1.5, capsize=4, capthick=1.5, lw=1.5, zorder=5)
    jc = rng.uniform(-0.08, 0.08, len(s["c_vals"]))
    ax_a.scatter(i + w/2 + jc, s["c_vals"], color=col, alpha=0.35, s=25, marker="D", zorder=4)
    ax_a.errorbar(i + w/2, s["c_mean"], yerr=s["c_std"], fmt="^", color=col,
                  ms=10, mec="white", mew=1.5, capsize=4, capthick=1.5, lw=1.5, zorder=5)
    ax_a.plot([i - w/2, i + w/2], [s["t_mean"], s["c_mean"]], color=col, lw=1, alpha=0.5, zorder=3)
    yl = min(s["t_mean"], s["c_mean"]) - 0.08
    ax_a.text(i, yl, r"$\Delta\alpha$ = " + f"{s['delta']:+.2f}",
              ha="center", fontsize=8, fontweight="bold", color=col)
    yt_t = max(s["t_vals"]) + 0.04
    yt_c = max(s["c_vals"]) + 0.04
    ax_a.text(i - w/2, yt_t, f"n={s['t_n']}", ha="center", fontsize=7, color="gray")
    ax_a.text(i + w/2, yt_c, f"n={s['c_n']}", ha="center", fontsize=7, color="gray")

ax_a.axhline(-1.0, color="gray", ls="--", lw=1, alpha=0.5)
ax_a.set_xticks(range(len(per_type)))
ax_a.set_xticklabels([DISPLAY.get(s["type"], s["type"]) for s in per_type],
                      fontsize=9, fontweight="bold")
ax_a.set_ylabel(r"Power-law exponent  $\alpha$", fontsize=13)
ax_a.set_title(r"A   Tissue-matched $\alpha$ comparison", fontsize=14,
               fontweight="bold", loc="left")
ax_a.legend(handles=[
    Line2D([0], [0], marker="s", color="gray", ls="None", ms=9,
           markerfacecolor="gray", label="TCGA tumour"),
    Line2D([0], [0], marker="^", color="gray", ls="None", ms=9,
           markerfacecolor="gray", label="Cell line"),
], loc="upper right", fontsize=10, framealpha=0.9)

# ── Panel B — y-axis label updated, grey background removed ──────────────────
for ct in matched:
    ct_cl = "GBM" if ct == "GBMx" else ct
    col   = TYPE_COLS.get(ct, "gray")
    samps = ps_combined.loc[ps_combined["cancer_type"] == ct_cl, "sample"].unique()

    for samp in samps:
        sub = ps_combined[
            (ps_combined["cancer_type"] == ct_cl) &
            (ps_combined["sample"] == samp)
        ].sort_values("genomic_distance")
        ax_b.loglog(sub["genomic_distance"] / 1e6, sub["contact_freq"],
                    color=col, alpha=max(0.12, 0.4 / len(samps)), lw=0.7)

    if len(samps) > 0:
        type_ps = ps_combined[ps_combined["cancer_type"] == ct_cl].copy()
        grid    = np.logspace(
            np.log10(type_ps["genomic_distance"].min()),
            np.log10(type_ps["genomic_distance"].max()), 120
        )
        interps = []
        for samp in samps:
            sub = type_ps[type_ps["sample"] == samp].sort_values("genomic_distance")
            interps.append(np.interp(
                np.log10(grid),
                np.log10(sub["genomic_distance"].values),
                np.log10(sub["contact_freq"].values),
            ))
        ax_b.loglog(grid / 1e6, 10 ** np.mean(interps, axis=0),
                    color=col, lw=2.5, alpha=0.9, label=ct)

xr = np.logspace(-1, 2, 50)
ax_b.loglog(xr, 0.3 * xr ** (-1.0), "k--", alpha=0.25, lw=1)
ax_b.text(25, 0.3 * 25 ** (-1.0) * 0.65, r"$\alpha = -1$",
          fontsize=9, color="gray", rotation=-28)

# NO axvspan grey background
ax_b.set_facecolor("white")
fig.patch.set_facecolor("white")

ax_b.set_xlabel("Genomic distance (Mb)", fontsize=13)
ax_b.set_ylabel("Contact probability, P(s)", fontsize=13)   # UPDATED
ax_b.set_title("B   P(s) curves - cell lines by cancer type",
               fontsize=14, fontweight="bold", loc="left")
ax_b.legend(fontsize=9, loc="lower left", framealpha=0.9, ncol=2)
ax_b.set_xlim(0.01, 100)

fig.tight_layout()
for e in ("png", "pdf"):
    fig.savefig(os.path.join(R, "fig3_tissue_matched_v2." + e), dpi=300, facecolor="white")
    print("Saved fig3_tissue_matched_v2." + e)
plt.close()
print("Done.")
