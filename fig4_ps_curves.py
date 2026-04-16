#!/usr/bin/env python3
"""
Figure 4 UPDATED: Cell line P(s) curves — all individual + mean per cancer type
Now includes 4 new GBM cell lines from GSE280831
Reads: cellline_alpha_results.csv, cellline_ps_curves.csv,
       gbm_cellline_alpha_results.csv, gbm_cellline_ps_curves.csv
Output: fig4_ps_curves.png/.pdf
CHANGES: Both panel y-axis labels -> Contact probability, P(s)
"""
import os
import platform
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = "/mnt/c/Users/jaipa/OneDrive/Documents/Dissertation" if platform.system() == "Linux" else r"C:\Users\jaipa\OneDrive\Documents\Dissertation"

cl_raw = pd.read_csv(os.path.join(R, "cellline_alpha_results.csv"))
ps_cl  = pd.read_csv(os.path.join(R, "cellline_ps_curves.csv"))
gbm_cl = pd.read_csv(os.path.join(R, "gbm_cellline_alpha_results.csv"))
ps_gbm = pd.read_csv(os.path.join(R, "gbm_cellline_ps_curves.csv"))

ps_all  = pd.concat([ps_cl, ps_gbm], ignore_index=True)
n_total = cl_raw["sample"].nunique() + gbm_cl["sample"].nunique()

TYPE_COLS = {
    "BLCA": "#C62828", "BRCA": "#AD1457", "COAD": "#6A1B9A", "GBM": "#1565C0",
    "LIHC": "#2E7D32", "LUAD": "#00838F", "PRAD": "#4E342E", "SKCM": "#E65100",
}

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 11,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 300, "savefig.dpi": 300, "savefig.bbox": "tight",
})

fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(14, 6))
fig.patch.set_facecolor("white")

# ── Panel A — all individual curves ──────────────────────────────────────────
for ct in sorted(ps_all["cancer_type"].unique()):
    col   = TYPE_COLS.get(ct, "gray")
    samps = ps_all.loc[ps_all["cancer_type"] == ct, "sample"].unique()
    for j, samp in enumerate(samps):
        sub = ps_all[
            (ps_all["cancer_type"] == ct) &
            (ps_all["sample"] == samp)
        ].sort_values("genomic_distance")
        ax_a.loglog(
            sub["genomic_distance"] / 1e6,
            sub["contact_freq"],
            color=col, alpha=0.4, lw=0.8,
            label=ct if j == 0 else None,
        )

xr = np.logspace(-1, 2, 50)
ax_a.loglog(xr, 0.3 * xr ** (-1.0), "k--", alpha=0.25, lw=1)
ax_a.text(20, 0.3 * 20 ** (-1.0) * 0.6, r"$\alpha = -1.0$",
          fontsize=9, color="gray", rotation=-28)
ax_a.set_xlabel("Genomic distance (Mb)", fontsize=13)
ax_a.set_ylabel("Contact probability, P(s)", fontsize=13)   # UPDATED
ax_a.set_title(f"A   Cell line P(s) curves (N={n_total})",
               fontsize=14, fontweight="bold", loc="left")
ax_a.legend(fontsize=8, loc="lower left", framealpha=0.9, ncol=2)
ax_a.set_xlim(0.01, 100)
ax_a.set_facecolor("white")

# ── Panel B — mean per cancer type with SD band ───────────────────────────────
for ct in sorted(ps_all["cancer_type"].unique()):
    col     = TYPE_COLS.get(ct, "gray")
    type_ps = ps_all[ps_all["cancer_type"] == ct].copy()
    samps   = type_ps["sample"].unique()
    grid    = np.logspace(
        np.log10(type_ps["genomic_distance"].min()),
        np.log10(type_ps["genomic_distance"].max()),
        120,
    )
    interps = []
    for samp in samps:
        sub = type_ps[type_ps["sample"] == samp].sort_values("genomic_distance")
        interps.append(np.interp(
            np.log10(grid),
            np.log10(sub["genomic_distance"].values),
            np.log10(sub["contact_freq"].values),
        ))
    interps  = np.array(interps)
    mean_log = np.mean(interps, axis=0)
    lbl      = f"{ct} (n={len(samps)})" if len(samps) > 1 else ct
    ax_b.loglog(grid / 1e6, 10 ** mean_log,
                color=col, lw=2.5, alpha=0.9, label=lbl)
    if len(samps) > 1:
        sd = np.std(interps, axis=0, ddof=1)
        ax_b.fill_between(grid / 1e6,
                          10 ** (mean_log - sd),
                          10 ** (mean_log + sd),
                          color=col, alpha=0.15)

ax_b.loglog(xr, 0.3 * xr ** (-1.0), "k--", alpha=0.25, lw=1)
ax_b.text(20, 0.3 * 20 ** (-1.0) * 0.6, r"$\alpha = -1.0$",
          fontsize=9, color="gray", rotation=-28)
ax_b.set_xlabel("Genomic distance (Mb)", fontsize=13)
ax_b.set_ylabel("Contact probability, P(s)", fontsize=13)   # UPDATED
ax_b.set_title("B   Mean P(s) by cancer type",
               fontsize=14, fontweight="bold", loc="left")
ax_b.legend(fontsize=8, loc="lower left", framealpha=0.9, ncol=2)
ax_b.set_xlim(0.01, 100)
ax_b.set_facecolor("white")

fig.tight_layout()
for e in ("png", "pdf"):
    fig.savefig(os.path.join(R, "fig4_ps_curves." + e), dpi=300, facecolor="white")
    print("Saved fig4_ps_curves." + e)
plt.close()
print("Done.")
