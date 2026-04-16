#!/usr/bin/env python3
"""
Figure 5: 8-panel P(s) — REAL TCGA tumour P(s) vs REAL cell line P(s)
UPDATED: Y-axis normalised to contact probability (0–1)
UPDATED: Alpha values annotated on each panel (per Bingxin)

Reads: TCGA contact txt files, cellline_ps_curves.csv, gbm_cellline_ps_curves.csv,
       gbm_patient_ps_curves.csv, TCGA_powerlaw_results.csv, cellline_alpha_results.csv,
       gbm_patient_alpha_results.csv, gbm_cellline_alpha_results.csv,
       DISSERTATION_master_alpha_table.csv
Output: fig5_8panel_ps.png / .pdf

Usage:
    cd /mnt/c/Users/jaipa/OneDrive/Documents/Dissertation
    python fig5_8panel_ps.py
"""
import os, platform, glob, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = "/mnt/c/Users/jaipa/OneDrive/Documents/Dissertation" if platform.system()=="Linux" else r"C:\Users\jaipa\OneDrive\Documents\Dissertation"
TCGA_DIR = os.path.join(R, "TCGA_HiChIP_hic")
RESOLUTION = 100000  # 100kb for TCGA contacts


def normalise_ps(df):
    """Normalise contact_freq to sum to 1 (contact probability)."""
    total = df["contact_freq"].sum()
    if total > 0:
        df = df.copy()
        df["contact_freq"] = df["contact_freq"] / total
    return df


# === Step 0: Load alpha values for annotation ===
tcga_alpha_df = pd.read_csv(os.path.join(R, "TCGA_HiChIP_hic", "TCGA_powerlaw_results.csv"))
tcga_alpha_df["alpha_neg"] = -tcga_alpha_df["Alpha"]

cl_alpha_df = pd.read_csv(os.path.join(R, "cellline_alpha_results.csv"))
cl_alpha_df["alpha_neg"] = -cl_alpha_df["alpha"]

master_df = pd.read_csv(os.path.join(R, "DISSERTATION_master_alpha_table.csv"))
old_cl = master_df[master_df["Dataset"] != "TCGA"].copy()

# GBM patient alphas
gbm_pat_alpha = pd.read_csv(os.path.join(R, "gbm_patient_alpha_results.csv"))
gbm_pat_alpha["alpha_neg"] = -gbm_pat_alpha["alpha"]

# GBM cell line alphas
gbm_cl_alpha = pd.read_csv(os.path.join(R, "gbm_cellline_alpha_results.csv"))
gbm_cl_alpha["alpha_neg"] = -gbm_cl_alpha["alpha"]

# Build per-cancer-type mean alpha for TCGA tumours
tcga_type_alpha = {}
for ct in ["BLCA", "BRCA", "COAD", "GBMx", "LIHC", "LUAD", "PRAD", "SKCM"]:
    vals = tcga_alpha_df.loc[tcga_alpha_df["Cancer_Type"] == ct, "alpha_neg"].values
    if ct == "GBMx":
        # Merge TCGA GBMx + GSE229962 patients
        gbm_vals = gbm_pat_alpha["alpha_neg"].values
        vals = np.concatenate([vals, gbm_vals])
    if len(vals) > 0:
        tcga_type_alpha[ct] = np.mean(vals)

# Build per-cancer-type mean alpha for cell lines
# Map old cell lines to cancer types
old_cancer_map = {"HCT116": "COAD", "HepG2": "LIHC"}
cl_type_alpha = {}
for ct_code in ["BLCA", "BRCA", "COAD", "GBM", "LIHC", "LUAD", "PRAD", "SKCM"]:
    # New cell lines
    vals = cl_alpha_df.loc[cl_alpha_df["cancer_type"] == ct_code, "alpha_neg"].values
    # Old cell lines mapped to this type
    for old_name, old_ct in old_cancer_map.items():
        if old_ct == ct_code:
            old_vals = old_cl.loc[old_cl["Sample"] == old_name, "Alpha"].values
            vals = np.concatenate([vals, old_vals])
    # GBM cell lines
    if ct_code == "GBM":
        vals = np.concatenate([vals, gbm_cl_alpha["alpha_neg"].values])
    if len(vals) > 0:
        mapped_key = "GBMx" if ct_code == "GBM" else ct_code
        cl_type_alpha[mapped_key] = np.mean(vals)


# === Step 1: Extract TCGA P(s) from contact txt files ===
def compute_ps_from_contacts(sample_prefix, tcga_dir, resolution):
    """Read all chr contact files for a sample, compute genome-wide P(s)."""
    all_contacts = {}
    for chrom_num in range(1, 23):
        pattern = os.path.join(tcga_dir, f"{sample_prefix}*_chr{chrom_num}_contacts.txt")
        files = glob.glob(pattern)
        if not files:
            continue
        for f in files:
            try:
                data = np.loadtxt(f, dtype=int)
            except:
                continue
            if data.ndim == 1:
                data = data.reshape(1, -1)
            for row in data:
                if len(row) < 3:
                    continue
                bi, bj, count = int(row[0]), int(row[1]), int(row[2])
                if bi == bj or count <= 0:
                    continue
                dist = abs(bj - bi) * resolution
                if dist > 0:
                    if dist not in all_contacts:
                        all_contacts[dist] = []
                    all_contacts[dist].append(count)

    if len(all_contacts) < 20:
        return None

    distances = np.array(sorted(all_contacts.keys()))
    mean_counts = np.array([np.mean(all_contacts[d]) for d in distances])

    # Log-bin
    log_d = np.log10(distances)
    bins = np.linspace(log_d.min(), log_d.max(), 150)
    bc, bm = [], []
    for i in range(len(bins) - 1):
        mask = (log_d >= bins[i]) & (log_d < bins[i + 1])
        if mask.sum() > 0:
            bc.append(10 ** ((bins[i] + bins[i + 1]) / 2))
            bm.append(np.mean(mean_counts[mask]))
    df = pd.DataFrame({"genomic_distance": bc, "contact_freq": bm})
    return normalise_ps(df)


# Find all individual TCGA samples for the 8 matched types
MATCHED = ["BLCA", "BRCA", "COAD", "GBMx", "LIHC", "LUAD", "PRAD", "SKCM"]

print("Extracting TCGA P(s) curves from contact files...")
tcga_ps = {}
for ct in MATCHED:
    chr1_files = glob.glob(os.path.join(TCGA_DIR, f"{ct}-*_chr1_contacts.txt"))
    sample_prefixes = []
    for f in chr1_files:
        base = os.path.basename(f)
        prefix = base.rsplit("_chr", 1)[0]
        sample_prefixes.append(prefix)

    tcga_ps[ct] = []
    for sp in sample_prefixes:
        short_name = sp.split("-")[1][:8] if "-" in sp else sp
        print(f"  {ct} / {short_name}...", end=" ", flush=True)
        ps = compute_ps_from_contacts(sp, TCGA_DIR, RESOLUTION)
        if ps is not None:
            ps["sample"] = short_name
            tcga_ps[ct].append(ps)
            print("OK")
        else:
            print("SKIP")

    # Also try pooled file
    pooled_chr1 = os.path.join(TCGA_DIR, f"{ct}_H3K27ac.allValidPairs.hic_chr1_contacts.txt")
    if os.path.exists(pooled_chr1):
        prefix = f"{ct}_H3K27ac.allValidPairs.hic"
        print(f"  {ct} / pooled...", end=" ", flush=True)
        ps = compute_ps_from_contacts(prefix, TCGA_DIR, RESOLUTION)
        if ps is not None:
            ps["sample"] = f"{ct}_pooled"
            tcga_ps[ct].append(ps)
            print("OK")
        else:
            print("SKIP")

    print(f"  {ct}: {len(tcga_ps[ct])} P(s) curves")

# === Step 2: Load cell line P(s) and normalise ===
ps_cl = pd.read_csv(os.path.join(R, "cellline_ps_curves.csv"))
ps_gbm_cl = pd.read_csv(os.path.join(R, "gbm_cellline_ps_curves.csv"))
ps_cl_all = pd.concat([ps_cl, ps_gbm_cl], ignore_index=True)

# Normalise each cell line sample
cl_normed = []
for samp in ps_cl_all["sample"].unique():
    sub = ps_cl_all[ps_cl_all["sample"] == samp].copy()
    sub = normalise_ps(sub)
    cl_normed.append(sub)
ps_cl_all = pd.concat(cl_normed, ignore_index=True)

# Load GBM patient P(s) and normalise
ps_gbm_pat = pd.read_csv(os.path.join(R, "gbm_patient_ps_curves.csv"))
gbm_normed = []
for samp in ps_gbm_pat["sample"].unique():
    sub = ps_gbm_pat[ps_gbm_pat["sample"] == samp].copy()
    sub = normalise_ps(sub)
    gbm_normed.append(sub)
ps_gbm_pat = pd.concat(gbm_normed, ignore_index=True)

# === Step 3: Plot 8-panel figure ===
TYPE_COLS = {
    "BLCA": "#C62828", "BRCA": "#AD1457", "COAD": "#6A1B9A", "GBMx": "#1565C0",
    "LIHC": "#2E7D32", "LUAD": "#00838F", "PRAD": "#4E342E", "SKCM": "#E65100",
}
DISPLAY = {
    "BLCA": "Bladder (BLCA)", "BRCA": "Breast (BRCA)", "COAD": "Colon (COAD)",
    "GBMx": "GBM", "LIHC": "Liver (LIHC)", "LUAD": "Lung Adeno (LUAD)",
    "PRAD": "Prostate (PRAD)", "SKCM": "Melanoma (SKCM)",
}

plt.rcParams.update({"font.family": "sans-serif", "font.size": 9,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 300, "savefig.dpi": 300, "savefig.bbox": "tight"})

fig, axes = plt.subplots(2, 4, figsize=(18, 9), sharex=True, sharey=True)

for idx, ct in enumerate(MATCHED):
    ax = axes.flatten()[idx]
    col = TYPE_COLS.get(ct, "gray")
    ct_cl = "GBM" if ct == "GBMx" else ct

    # --- TCGA tumour curves ---
    tumour_curves = tcga_ps.get(ct, [])

    # For GBM, also add the 30 patient samples
    if ct == "GBMx":
        gbm_samps = ps_gbm_pat["sample"].unique()
        for samp in gbm_samps:
            sub = ps_gbm_pat[ps_gbm_pat["sample"] == samp].sort_values("genomic_distance")
            tumour_curves.append(sub[["genomic_distance", "contact_freq"]].assign(sample=samp))

    n_tumour = len(tumour_curves)
    for tc in tumour_curves:
        tc_sorted = tc.sort_values("genomic_distance")
        ax.loglog(tc_sorted["genomic_distance"] / 1e6, tc_sorted["contact_freq"],
                  color=col, alpha=0.25, lw=0.7)

    # Tumour mean
    if n_tumour > 0:
        all_t = pd.concat(tumour_curves)
        grid = np.logspace(np.log10(all_t["genomic_distance"].min()),
                           np.log10(all_t["genomic_distance"].max()), 100)
        interps = []
        for tc in tumour_curves:
            tc_s = tc.sort_values("genomic_distance")
            interps.append(np.interp(np.log10(grid),
                np.log10(tc_s["genomic_distance"].values),
                np.log10(tc_s["contact_freq"].values)))
        mean_t = np.mean(interps, axis=0)
        t_alpha = tcga_type_alpha.get(ct, None)
        t_label = f"Tumour (n={n_tumour})" + (f", α={t_alpha:.2f}" if t_alpha is not None else "")
        ax.loglog(grid / 1e6, 10**mean_t, color=col, lw=2.5, alpha=0.9, label=t_label)

    # --- Cell line curves ---
    cl_samps = ps_cl_all.loc[ps_cl_all["cancer_type"] == ct_cl, "sample"].unique()
    n_cl = len(cl_samps)
    for samp in cl_samps:
        sub = ps_cl_all[(ps_cl_all["cancer_type"] == ct_cl) & (ps_cl_all["sample"] == samp)].sort_values("genomic_distance")
        ax.loglog(sub["genomic_distance"] / 1e6, sub["contact_freq"],
                  color="gray", alpha=0.3, lw=0.7)

    # Cell line mean
    if n_cl > 0:
        type_ps = ps_cl_all[ps_cl_all["cancer_type"] == ct_cl]
        grid_c = np.logspace(np.log10(type_ps["genomic_distance"].min()),
                             np.log10(type_ps["genomic_distance"].max()), 100)
        interps_c = []
        for samp in cl_samps:
            sub = type_ps[type_ps["sample"] == samp].sort_values("genomic_distance")
            interps_c.append(np.interp(np.log10(grid_c),
                np.log10(sub["genomic_distance"].values),
                np.log10(sub["contact_freq"].values)))
        mean_c = np.mean(interps_c, axis=0)
        c_alpha = cl_type_alpha.get(ct, None)
        c_label = f"Cell line (n={n_cl})" + (f", α={c_alpha:.2f}" if c_alpha is not None else "")
        ax.loglog(grid_c / 1e6, 10**mean_c, color="black", lw=2, ls="--", alpha=0.8,
                  label=c_label)

    # Δα annotation box
    t_a = tcga_type_alpha.get(ct, None)
    c_a = cl_type_alpha.get(ct, None)
    if t_a is not None and c_a is not None:
        delta = c_a - t_a
        ax.text(0.97, 0.95, r"$\Delta\alpha$" + f" = {delta:+.2f}",
                transform=ax.transAxes, fontsize=9, fontweight="bold",
                ha="right", va="top", color=col,
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=col, alpha=0.8))

    # Reference
    xr = np.logspace(-1, 2, 50)
    ax.loglog(xr, 0.3 * xr**(-1.0), color="gray", ls=":", alpha=0.2, lw=1)
    ax.axvspan(0.1, 10, alpha=0.04, color="gray")

    ax.set_title(f"{DISPLAY[ct]}\nTumour n={n_tumour}, CL n={n_cl}", fontsize=10, fontweight="bold")
    ax.legend(fontsize=7, loc="lower left", framealpha=0.8)
    ax.set_xlim(0.01, 100)

for ax in axes[1, :]:
    ax.set_xlabel("Genomic distance (Mb)", fontsize=11)
for ax in axes[:, 0]:
    ax.set_ylabel(r"Contact probability P(s)", fontsize=11)

fig.suptitle("P(s) curves by cancer type — tumours vs cell lines",
             fontsize=14, fontweight="bold", y=1.02)
fig.tight_layout()
for e in ("png", "pdf"):
    fig.savefig(os.path.join(R, f"fig5_8panel_ps.{e}"))
    print(f"Saved fig5_8panel_ps.{e}")
plt.close()
print("Done")