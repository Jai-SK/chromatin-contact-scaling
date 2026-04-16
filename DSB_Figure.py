#!/usr/bin/env python3
"""
Figure — Cohesin depletion & DSB disrupts contact scaling
Single panel: contact frequency vs genomic distance (8 conditions)
Data: *_frequencies.csv from DIvA/OHT experiments
"""
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial','Helvetica','DejaVu Sans'],
    'font.size': 10, 'axes.linewidth': 0.8, 'axes.spines.top': False,
    'axes.spines.right': False, 'pdf.fonttype': 42, 'ps.fonttype': 42,
    'savefig.dpi': 300, 'savefig.bbox': 'tight',
})

DATA = Path('/mnt/user-data/uploads')
OUT = Path('/mnt/user-data/outputs/thesis_figures')
OUT.mkdir(parents=True, exist_ok=True)

CONDITIONS = [
    ('Control (DIvA)',     'Control_DIvA_frequencies.csv',     '#2CA02C', '-',  2.2),
    ('Control (OHT)',      'Control_OHT_frequencies.csv',      '#2CA02C', '--', 1.4),
    ('Cohesin KD (DIvA)',  'Cohesin_KD_DIvA_frequencies.csv',  '#D62728', '-',  2.2),
    ('Cohesin KD (OHT)',   'Cohesin_KD_OHT_frequencies.csv',   '#D62728', '--', 1.4),
    ('DSB (DIvA) Rep A',   'DSB_DIvA_RepA_frequencies.csv',    '#FF7F0E', '-',  1.7),
    ('DSB (DIvA) Rep B',   'DSB_DIvA_RepB_frequencies.csv',    '#FF7F0E', ':',  1.3),
    ('DSB (OHT) Rep A',    'DSB_OHT_RepA_frequencies.csv',     '#E6AB02', '-',  1.7),
    ('DSB (OHT) Rep B',    'DSB_OHT_RepB_frequencies.csv',     '#E6AB02', ':',  1.3),
]

perturb = {}
for name, fname, col, ls, lw in CONDITIONS:
    df = pd.read_csv(DATA / fname)
    df['dist_mb'] = df['distance_bp'] / 1e6
    mask = (df['dist_mb'] >= 0.1) & (df['dist_mb'] <= 10) & (df['mean_count'] > 0)
    if mask.sum() >= 3:
        x = np.log10(df.loc[mask, 'dist_mb'].values)
        y = np.log10(df.loc[mask, 'mean_count'].values)
        slope, intercept, r, p, se = stats.linregress(x, y)
        alpha, r2 = slope, r**2
    else:
        alpha, r2 = np.nan, np.nan
    perturb[name] = {'df': df, 'alpha': alpha, 'r2': r2, 'col': col, 'ls': ls, 'lw': lw}

# ── Figure ──
fig, ax = plt.subplots(figsize=(8, 6))

for name, fname, col, ls, lw in CONDITIONS:
    d = perturb[name]
    df = d['df']
    alpha = d['alpha']
    mask = df['mean_count'] > 0
    ax.loglog(df.loc[mask, 'dist_mb'], df.loc[mask, 'mean_count'],
              ls=ls, color=col, lw=lw, alpha=0.85,
              label=f'{name} (α={alpha:.2f})')

ax.axvspan(0.1, 10, alpha=0.05, color='#888', zorder=0)
ax.set_xlabel('Genomic distance (Mb)', fontsize=12)
ax.set_ylabel('Contact frequency', fontsize=12)
ax.set_title('Cohesin depletion disrupts contact scaling', fontsize=14,
             fontweight='bold', loc='left', pad=10)
ax.legend(fontsize=8, loc='upper right', frameon=True, fancybox=False, edgecolor='#CCC')

ctrl_a = np.mean([perturb['Control (DIvA)']['alpha'], perturb['Control (OHT)']['alpha']])
ckd_a = np.mean([perturb['Cohesin KD (DIvA)']['alpha'], perturb['Cohesin KD (OHT)']['alpha']])
dsb_a = np.mean([perturb[k]['alpha'] for k in perturb if 'DSB' in k])
ax.text(0.03, 0.55,
        f'Control: α = {ctrl_a:.2f}\nCohesin KD: α = {ckd_a:.2f}\nDSB: α = {dsb_a:.2f}\n'
        f'Δα (Control→KD) = {ckd_a-ctrl_a:.2f}\nΔα (Control→DSB) = {dsb_a-ctrl_a:.2f}',
        transform=ax.transAxes, ha='left', va='top', fontsize=9,
        bbox=dict(boxstyle='round,pad=0.4', fc='white', ec='#CCC', alpha=0.9))

plt.tight_layout()
for ext in ['png', 'pdf']:
    fig.savefig(OUT / f'Figure_perturbation.{ext}', dpi=300, facecolor='white')
plt.close()
print("Done: Figure_perturbation")