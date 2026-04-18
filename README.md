# Chromatin Contact Scaling Exponents Across Primary Tumours and Cancer Cell Lines

**Author:** Jaipal Singh Kooner
**Institution:** School of Biosciences, University of Surrey
**Year:** 2026

---

## Overview

This repository contains all analysis code for the MSci dissertation:

> *Cancer-wide analysis reveals systematic divergence in chromatin contact scaling between primary tumours and cell line models*

We systematically compare the power-law contact scaling exponent α between 78 primary tumour samples (14 cancer types) and 44 cancer cell lines across eight tissue-matched cancer types. Cell lines show consistently steeper P(s) decay than matched tumours in seven of eight cancer types (median Δα = -0.18; Mann-Whitney U p = 1.81 × 10⁻⁵; Cohen's d = 0.87).

---

## Repository Contents

All scripts live in the repository root. They are grouped below by role in the analysis pipeline.

### Hi-C data processing

Convert raw Hi-C data (`.hic`, `.mcool`, `.allValidPairs`) into `.cool` format, compute P(s) via `cooltools.expected_cis`, and fit the power-law scaling exponent α via OLS regression across the 0.1-10 Mb window.

| Script | Purpose |
| --- | --- |
| `process_all_cellline_hic.py` | Pipeline for processing the full cell line cohort (ENCODE + tissue-matched lines) |
| `process_gbm_celllines.py` | Dedicated pipeline for GBM glioma stem cell lines (GSE280831) and GBM_P4CC3 |
| `Cell-line_data.py` | Cell line α collation and formatting |

### Data aggregation

| Script | Purpose |
| --- | --- |
| `build_master_alpha_table.py` | Consolidates all tumour and cell line α values into the master alpha table used for downstream analysis |

### Statistical analysis

| Script | Purpose |
| --- | --- |
| `8_sample_comparison.py` | Tissue-matched tumour vs cell line comparison across the eight cancer types (Mann-Whitney U, Welch's t, Cohen's d, Bonferroni correction for k = 8) |
| `Tcga_Bx_Plot.py` | TCGA Pan-Cancer Atlas α exploratory analysis and boxplot |

### Figure generation

| Script | Manuscript figure | Description |
| --- | --- | --- |
| `fig1_violin.py` | Fig 1 | Violin / box plot of α for tumours vs cell lines across all cancer types |
| `fig2_boxplot_by_type.py` | Fig 2 | Cell line α by cancer type |
| `fig4_ps_curves.py` | Fig 3 | P(s) curves for all cell line samples with mean per-type overlay |
| `fig3_tissue_matched_v2.py` | Fig 4 | Tissue-matched P(s) curve comparison, eight-panel log-log plot |
| `DSB_Figure.py` | Fig 5 | Cohesin depletion and DNA double-strand break scaling analysis |

> **Note on script names:** `fig3_tissue_matched_v2.py` and `fig4_ps_curves.py` are named for historical draft order, not the final manuscript figure numbers. The correct manuscript mapping is given above.

### Supplementary data

| File | Description |
| --- | --- |
| `S1_Table.xlsx` | Cell line provenance and per-sample scaling exponents referenced as S1 Table in the manuscript |

---

## Data Sources

| Dataset | Description | Access |
| --- | --- | --- |
| TCGA Pan-Cancer HiChIP | 48 primary tumour samples, H3K27ac HiChIP | [GDC Portal](https://portal.gdc.cancer.gov) |
| GSE229962 | 30 primary GBM Hi-C samples | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE229962) |
| GSE118629 | PRAD urological cell lines | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118629) |
| GSE92804 | Respiratory cell lines (A549) | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92804) |
| GSE280831 | GBM glioma stem cell lines | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE280831) |
| GSE253407 | Melanoma cell lines | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253407) |
| GSE185069 | Pancreatic cell lines | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185069) |
| GSE133928 | COAD allValidPairs (SW480) | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE133928) |
| ENCODE | Reference cell lines (GM12878, K562, IMR90, HMEC, NHEK, HUVEC, HepG2, HCT116) | [ENCODE Portal](https://www.encodeproject.org) |

---

## Installation

### Conda environment (recommended)

```bash
conda create -n hic-scaling python=3.10
conda activate hic-scaling
pip install cooler cooltools hic2cool hicstraw numpy pandas scipy matplotlib seaborn statsmodels openpyxl
```

### Minimum package versions

| Package | Version |
| --- | --- |
| cooler | ≥ 0.9.3 |
| cooltools | ≥ 0.7.1 |
| hic2cool | ≥ 0.8.3 |
| hicstraw | ≥ 1.3.1 |
| numpy | ≥ 1.24 |
| pandas | ≥ 2.0 |
| scipy | ≥ 1.11 |
| matplotlib | ≥ 3.7 |
| seaborn | ≥ 0.12 |
| statsmodels | ≥ 0.14 |
| openpyxl | ≥ 3.1 |

---

## Usage

### 1. Process raw Hi-C data

Run the cell line pipeline, which converts `.hic` / `.mcool` / `.allValidPairs` inputs to `.cool`, computes P(s), and fits α across 0.1-10 Mb:

```bash
python process_all_cellline_hic.py
python process_gbm_celllines.py
python Cell-line_data.py
```

### 2. Build master α table

```bash
python build_master_alpha_table.py
```

### 3. Statistical analysis

```bash
python 8_sample_comparison.py
python Tcga_Bx_Plot.py
```

### 4. Generate figures

```bash
python fig1_violin.py            # Fig 1: Tumour vs cell line α (overall)
python fig2_boxplot_by_type.py   # Fig 2: Cell line α by cancer type
python fig4_ps_curves.py         # Fig 3: P(s) curves for cell lines
python fig3_tissue_matched_v2.py # Fig 4: Tissue-matched P(s) comparison
python DSB_Figure.py             # Fig 5: Cohesin depletion / DSB
```

All figures are saved as both `.png` and `.pdf` in the working directory.

---

## Methods Summary

### P(s) estimation

Genome-wide contact probability P(s) was estimated using `cooltools.expected_cis` (v0.7.1):

```
P(d) = (1/N(d)) * Σ C(i, i+d)
```

where N(d) is the count of valid bin pairs at diagonal offset d. Genomic bin size: 40 kb or 50 kb.

### Power-law fitting

The scaling exponent α was extracted by OLS regression of the linearised model:

```
log10 P(s) = α · log10(s) + log10(A)
```

Fitting range: **0.1-10 Mb** (mesoscale window isolating TAD and loop extrusion dynamics).

### Statistical tests

| Test | Purpose |
| --- | --- |
| Mann-Whitney U | Non-parametric comparison of tumour vs cell line α distributions |
| Welch's t-test | Parametric comparison (unequal variance) |
| Cohen's d | Effect size (scale-independent) |
| Bonferroni correction | Multiple testing across k = 8 tissue-matched comparisons |

### Key results

| Comparison | Result |
| --- | --- |
| Overall p-value | 1.81 × 10⁻⁵ (Mann-Whitney U) |
| Effect size | Cohen's d = 0.87 |
| Median Δα | −0.18 |
| Cancer types with negative Δα | 7 of 8 |
| GBM exception | Δα = +0.26 |

---

## Biophysical Model

Based on the fractal globule polymer model, spatial distance scales with genomic separation as:

```
d(s) ∝ s^(−α/3)
```

| Regime | α | d(s) scaling |
| --- | --- | --- |
| Fractal globule | −1.0 | s^(1/3) |
| Equilibrium globule | −1.5 | s^(1/2) |
| Loop-enriched / mitotic | −0.5 | s^(1/6) |

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
