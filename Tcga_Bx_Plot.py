#!/usr/bin/env python3
"""
TCGA box plot — all tumour samples, GBM merged into single group
GBMx (TCGA n=3) + GBM patients (GSE229962 n=30) = GBM (n=33)

Reads: TCGA_powerlaw_results.csv, gbm_patient_alpha_results.csv,
       cellline_alpha_results.csv, DISSERTATION_master_alpha_table.csv,
       gbm_cellline_alpha_results.csv
Output: fig_tcga_boxplot.png/.pdf
"""
import os, platform, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

R = "/mnt/c/Users/jaipa/OneDrive/Documents/Dissertation" if platform.system()=="Linux" else r"C:\Users\jaipa\OneDrive\Documents\Dissertation"

plt.rcParams.update({'font.family':'sans-serif','font.size':10,'axes.linewidth':0.8,
    'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':300,'savefig.bbox':'tight'})

# Load TCGA
tcga = pd.read_csv(os.path.join(R,"TCGA_HiChIP_hic","TCGA_powerlaw_results.csv"))
tcga['alpha'] = -tcga['Alpha']
# Rename GBMx -> GBM
tcga['Cancer_Type'] = tcga['Cancer_Type'].replace({'GBMx': 'GBM'})

# Load GBM patients — label as GBM (same group)
gbm_pat = pd.read_csv(os.path.join(R,"gbm_patient_alpha_results.csv"))
gbm_df = pd.DataFrame({'Cancer_Type': ['GBM']*len(gbm_pat), 'alpha': -gbm_pat['alpha'].values})

# Combine
tcga_all = pd.concat([tcga[['Cancer_Type','alpha']], gbm_df], ignore_index=True)

# Cell line range
cl_new = pd.read_csv(os.path.join(R,"cellline_alpha_results.csv"))
master = pd.read_csv(os.path.join(R,"DISSERTATION_master_alpha_table.csv"))
gbm_cl = pd.read_csv(os.path.join(R,"gbm_cellline_alpha_results.csv"))
old_cl = master[master['Dataset']!='TCGA']
all_cl_a = np.concatenate([old_cl['Alpha'].values, -cl_new['alpha'].values, -gbm_cl['alpha'].values])

CT_COLOURS = {
    'BLCA':'#E41A1C','BRCA':'#FF7F00','COAD':'#FFD700','ESCA':'#4DAF4A',
    'GBM':'#377EB8','KIRC':'#984EA3','LIHC':'#A65628','LUAD':'#F781BF',
    'LUSC':'#00CED1','PRAD':'#66C2A5','SKCM':'#FC8D62','STAD':'#8DA0CB',
    'THCA':'#E78AC3','UCEC':'#A6D854',
}

ct_stats = tcga_all.groupby('Cancer_Type')['alpha'].agg(['median','count']).sort_values('median')
ct_order = ct_stats.index.tolist()

fig, ax = plt.subplots(figsize=(15, 5.5))

for i, ct in enumerate(ct_order):
    vals = tcga_all[tcga_all['Cancer_Type']==ct]['alpha'].values
    col = CT_COLOURS.get(ct,'#888888'); n=len(vals)
    if n>=2:
        bp = ax.boxplot([vals],positions=[i],widths=0.5,patch_artist=True,showfliers=False,zorder=2,
            boxprops=dict(facecolor=col,alpha=0.3,edgecolor=col,lw=1),
            medianprops=dict(color='black',lw=1.5),
            whiskerprops=dict(color=col,lw=0.8),capprops=dict(color=col,lw=0.8))
    rng = np.random.default_rng(42+i)
    jit = rng.uniform(-0.15,0.15,n)
    ax.scatter(np.full(n,i)+jit,vals,s=30,color=col,alpha=0.8,edgecolors='white',linewidth=0.4,zorder=3)

# Cell line range
ax.axhspan(all_cl_a.min(),all_cl_a.max(),alpha=0.08,color='#4393C3',zorder=0)
ax.text(0.5,all_cl_a.max()+0.008,'cell line α range',fontsize=8,color='#4393C3',ha='left',va='bottom',style='italic')

# Steeper decay arrow
yl=ax.get_ylim()
ax.annotate('',xy=(len(ct_order)-0.3,yl[0]+0.01),xytext=(len(ct_order)-0.3,yl[0]+0.12),
    arrowprops=dict(arrowstyle='->',color='#888888',lw=1))
ax.text(len(ct_order)-0.3,yl[0]+0.13,'steeper\ndecay',fontsize=7,color='#888888',ha='center',va='bottom')

ct_labels = [f"{ct}\n(n={int(ct_stats.loc[ct,'count'])})" for ct in ct_order]
ax.set_xticks(range(len(ct_order)))
ax.set_xticklabels(ct_labels,fontsize=9,fontweight='bold')
ax.set_ylabel(r'Power-law exponent  $\alpha$',fontsize=11)
ax.set_title('Contact decay scaling across tumour samples',fontsize=13,fontweight='bold',loc='left',pad=10)
ax.text(0.99,0.03,f'N = {len(tcga_all)} tumour samples\n{len(ct_order)} cancer types',
    transform=ax.transAxes,ha='right',va='bottom',fontsize=9,
    bbox=dict(boxstyle='round,pad=0.3',fc='white',ec='#CCCCCC',alpha=0.9))

plt.tight_layout()
for e in['png','pdf']:
    fig.savefig(os.path.join(R,f'fig_tcga_boxplot.{e}'),dpi=300,facecolor='white')
    print(f"Saved fig_tcga_boxplot.{e}")
plt.close()
print("Done")
