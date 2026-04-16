# Chromatin Contact Scaling Exponents Across Primary Tumours and Cancer Cell Lines

**Author:** Jaipal Singh Kooner  
**Institution:** School of Biosciences, University of Surrey  
**Year:** 2026  

---

## Overview

This repository contains all analysis code for the MSci dissertation:

> *Cancer-wide analysis reveals systematic divergence in chromatin contact scaling between primary tumours and cell line models*

We systematically compare the power-law contact scaling exponent α between 78 primary tumour samples (14 cancer types) and 44 cancer cell lines across eight tissue-matched cancer types. Cell lines show consistently steeper P(s) decay than matched tumours in seven of eight cancer types (median Δα = -0.18; Mann-Whitney U p = 1.81 × 10-5; Cohen's d = 0.87).

---

## Repository Structure

```
.
├── README.md
├── requirements.txt
│
├── pipeline/
│   ├── hic_to_cool.py            # Convert .hic / .mcool / allValidPairs to .cool
│   ├── compute_ps.py             # Compute P(s) using cooltools.expected_cis
│   ├── fit_alpha.py              # OLS power-law regression (0.1–10 Mb window)
│   └── batch_process.py          # Run full pipeline across all samples
│
├── analysis/
│   ├── compare_tumour_cellline.py  # Mann-Whitney U, Welch t, Cohen's d
│   ├── tissue_matched.py           # Per-cancer-type Δα with Bonferroni correction
│   └── cohesin_perturbation.py     # Reanalysis of OHT/DIvA cohesin depletion data
│
├── figures/
│   ├── fig1_tumour_boxplot.py      # Fig 1A: Tumour α by cancer type
│   ├── fig2_boxplot_by_type.py     # Fig 1B: Cell line α by cancer type
│   ├── fig3_ps_curves.py           # Fig 3: P(s) curves for cell lines
│   ├── fig4_tissue_matched_ps.py   # Fig 4: Tissue-matched P(s) comparison
│   └── fig5_cohesin_dsb.py         # Fig 5: Cohesin depletion / DSB scaling
│
└── data/
    └── README_data.md              # Data sources and download instructions
```

---

## Data Sources

| Dataset | Description | Access |
|---------|-------------|--------|
| TCGA Pan-Cancer HiChIP | 48 primary tumour samples, H3K27ac HiChIP | [GDC Portal](https://portal.gdc.cancer.gov) |
| GSE229962 | 30 primary GBM Hi-C samples | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE229962) |
| GSE118629 | PRAD urological cell lines | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118629) |
| GSE92804 | Respiratory cell lines | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92804) |
| GSE280831 | GBM glioma stem cell lines | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE280831) |
| GSE253407 | Melanoma cell lines | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253407) |
| GSE185069 | Pancreatic cell lines | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185069) |
| GSE133928 | COAD allValidPairs | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE133928) |
| ENCODE | Reference cell lines (GM12878, K562 etc.) | [ENCODE Portal](https://www.encodeproject.org) |

---

## Installation

### Requirements

```bash
pip install -r requirements.txt
```

### `requirements.txt`

```
cooler>=0.9.3
cooltools>=0.7.1
hic2cool>=0.8.3
hicstraw>=1.3.1
numpy>=1.24
pandas>=2.0
scipy>=1.11
matplotlib>=3.7
seaborn>=0.12
statsmodels>=0.14
```

### Conda environment (recommended)

```bash
conda create -n hic-scaling python=3.10
conda activate hic-scaling
pip install -r requirements.txt
```

---

## Usage

### 1. Process raw Hi-C data

```bash
# Convert .hic file to .cool
python pipeline/hic_to_cool.py \
    --input sample.hic \
    --output sample.cool \
    --resolution 50000

# Compute P(s) contact probability curve
python pipeline/compute_ps.py \
    --cool sample.cool \
    --output sample_ps.csv

# Fit power-law scaling exponent alpha
python pipeline/fit_alpha.py \
    --ps_file sample_ps.csv \
    --fit_min 0.1 \
    --fit_max 10 \
    --output sample_alpha.csv
```

### 2. Run full batch pipeline

```bash
python pipeline/batch_process.py \
    --input_dir /path/to/cool_files \
    --output_dir /path/to/results
```

### 3. Statistical analysis

```bash
# Tumour vs cell line comparison
python analysis/compare_tumour_cellline.py \
    --tumour_alpha DISSERTATION_master_alpha_table.csv \
    --cellline_alpha cellline_alpha_results.csv

# Tissue-matched Δα with Bonferroni correction
python analysis/tissue_matched.py
```

### 4. Generate figures

```bash
cd /path/to/Dissertation

python figures/fig1_tumour_boxplot.py
python figures/fig2_boxplot_by_type.py
python figures/fig3_ps_curves.py
python figures/fig4_tissue_matched_ps.py
python figures/fig5_cohesin_dsb.py
```

All figures are saved as both `.png` and `.pdf` in the working directory.

---

## Methods Summary

### P(s) Estimation

Genome-wide contact probability P(s) was estimated using `cooltools.expected_cis` (v0.7.1):

```
P(d) = (1/N(d)) * Σ C(i, i+d)
```

where N(d) is the count of valid bin pairs at diagonal offset d. Genomic bin size: 40 kb or 50 kb.

### Power-Law Fitting

The scaling exponent α was extracted by OLS regression of the linearised model:

```
log10 P(s) = α · log10(s) + log10(A)
```

Fitting range: **0.1–10 Mb** (mesoscale window isolating TAD and loop extrusion dynamics).

### Statistical Tests

| Test | Purpose |
|------|---------|
| Mann-Whitney U | Non-parametric comparison of tumour vs cell line α distributions |
| Welch's t-test | Parametric comparison (unequal variance) |
| Cohen's d | Effect size (scale-independent) |
| Bonferroni correction | Multiple testing across k=8 tissue-matched comparisons |

### Key Results

| Comparison | Result |
|-----------|--------|
| Overall p-value | 1.81 × 10⁻⁵ (Mann-Whitney U) |
| Effect size | Cohen's d = 0.87 |
| Median Δα | −0.18 |
| Cancer types with negative Δα | 7 of 8 |
| GBM exception | Δα = +0.26 |

---

## Biophysical Model

Based on the Fractal Globule polymer model, the spatial distance scaling follows:

```
d(s) ∝ s^(−α/3)
```

| Regime | α | d(s) scaling |
|--------|---|-------------|
| Fractal globule | -1.0 | s^(1/3) |
| Equilibrium globule | -1.5 | s^(1/2) |
| Loop-enriched / mitotic | -0.5 | s^(1/6) |

---

## Citation

If you use this code, please cite:

```
Kooner JS. Cancer-wide analysis reveals systematic divergence in chromatin 
contact scaling between primary tumours and cell line models. 
MSci Dissertation, University of Surrey, 2026.
```

---

## Acknowledgements

Supervised by Dr Bingxin Lu, School of Biosciences, University of Surrey.

---

## Licence

MIT License — free to use and adapt with attribution.
