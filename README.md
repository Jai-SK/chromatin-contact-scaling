# Chromatin Contact Scaling Exponents Across Primary Tumours and Cancer Cell Lines

**Author:** Jaipal Singh Kooner
**Supervisor:** Dr Bingxin Lu
**Institution:** School of Biosciences, University of Surrey
**Year:** 2026

Analysis code for the MSci dissertation:

> *Cancer-wide analysis reveals systematic divergence in chromatin contact scaling between primary tumours and cell line models*

We compare the power-law contact scaling exponent α, where P(s) ∝ s^α and α is reported as negative (more negative means steeper decay), between 78 primary tumour samples (14 cancer types) and 44 cancer cell lines across eight tissue-matched cancer types. In seven of the eight types, cell lines show steeper P(s) decay than matched tumours (median Δα = −0.18; Mann-Whitney U p = 1.81 × 10⁻⁵; Cohen's d = 0.87). GBM is the exception (Δα = +0.26).

---

## Repository layout

```
chromatin-contact-scaling/
├── preprocess/      # numbered steps: download Hi-C, compute P(s), fit α, build master table
├── plotting/        # one script per figure (figN_*), plus supplementary plots (supp_*)
├── data/
│   └── S1_Table.xlsx   # supplementary table S1: cell line provenance and per-sample α
├── LICENSE
└── README.md
```

### `preprocess/`

| Script | What it does | Main outputs |
| --- | --- | --- |
| `01_download_reference_hic.py` | Downloads and validates reference `.hic` files (GM12878, K562, IMR90, HMEC, HUVEC, NHEK from GSE63525; HCT116 and HepG2 from 4DN). Supports `--dry-run`, `--output` and `--logs`. | `.hic` files, `logs/results_*.csv` |
| `02_process_cellline_hic.py` | Main cell line pipeline for BLCA, BRCA, LUAD, LIHC, PRAD, SKCM, GBM and COAD. Converts `.hic`, `.mcool`, `.cool`, dense matrices and `allValidPairs` inputs, computes P(s) with `cooltools.expected_cis`, log-bins the curve (200 bins) and fits α by OLS over 100 kb to 10 Mb. | `cellline_alpha_results.csv`, `cellline_ps_curves.csv` |
| `03_process_gbm_celllines.py` | Handles the GBM glioma stem cell lines from GSE280831 (GSC_320, GSC_728, GSC_810, GSC_387). It uses `hicstraw` because hic2cool cannot read these v9 `.hic` files, and applies the same P(s) and fitting method. | `gbm_cellline_alpha_results.csv`, `gbm_cellline_ps_curves.csv` |
| `04_build_master_alpha_table.py` | Combines TCGA tumour α, GSE185069 cell line α and published reference values (Rao et al. 2014) into one table. | `DISSERTATION_master_alpha_table.csv` |

### `plotting/`

| Script | Figure | Output files |
| --- | --- | --- |
Each script saves its figure as `<script name>.png` and `<script name>.pdf`.

| Script | Figure |
| --- | --- |
| `fig1_violin_tumour_vs_cellline.py` | Fig 1: tumour vs cell line α across all cancer types. **Also computes the main statistics** (Mann-Whitney U, Welch's t, Cohen's d). |
| `fig2_cellline_alpha_by_type.py` | Fig 2: cell line α by cancer type |
| `fig3_cellline_ps_curves.py` | Fig 3: P(s) curves for all cell lines, with a per-type mean overlay |
| `fig4_tissue_matched.py` | Fig 4: tissue-matched tumour vs cell line comparison |
| `fig5_dsb_cohesin.py` | Fig 5: effect of cohesin knockdown and DNA double-strand breaks (DIvA/OHT) on contact frequency vs distance |
| `supp_8panel_ps_comparison.py` | Eight-panel plot of tumour vs cell line P(s), one panel per matched cancer type, with α annotated on each panel |
| `supp_tcga_alpha_boxplot.py` | Exploratory boxplot of TCGA tumour α by cancer type, with cell lines alongside |

---

## Before you run anything: configure paths

The scripts do **not** take input paths on the command line. Each script sets a data root near the top:

```python
# Linux / WSL
ROOT = "/mnt/c/Users/jaipa/OneDrive/Documents/Dissertation"
# Windows
ROOT = r"C:\Users\jaipa\OneDrive\Documents\Dissertation"
```

The variable is called `ROOT`, `R` or `DATA_ROOT`, depending on the script. Before running a script, set it to your own data directory. All intermediate CSVs and figures are read from and written to this directory.

`fig5_dsb_cohesin.py` is the exception. It reads from `DATA = /mnt/user-data/uploads` and writes to `OUT = /mnt/user-data/outputs/thesis_figures`, so edit both.

Expected layout under `ROOT`, which `02_process_cellline_hic.py` relies on:

```
ROOT/
├── GSE148079/          GSE148079_{HT1376,RT4,SCABER,SW780}_HiC.mcool     (BLCA)
├── GSE167154/          {T47D,ZR7530,HCC1954,HCC70,BT549}.hic             (BRCA)
├── GSE92804/           GSE92804_ENCFF*_chromatin_interactions.hic        (LUAD, A549)
├── GSE226216/          GSE226216_{Hep3B,Huh1,Huh7,SNU449}_*.hic          (LIHC)
├── GSE253407_analysis/ {LC676,…,LC821,GBM_P4CC3}.cool                    (SKCM, GBM)
├── GSE133928/          SW480_allValidPairs.txt                           (COAD)
├── GSE280831_GBM_celllines/  GSM86062{79..82}_B6-S*.hic                  (GBM GSCs)
├── GSE185069_analysis/ cancer_powerlaw_results.csv
├── TCGA_HiChIP_hic/    TCGA_powerlaw_results.csv, *_chr*_contacts.txt
└── GSE118629_22Rv1_HiC_40k.normalized.matrix.txt                         (PRAD)
```

### Inputs not produced by this repository

Some downstream scripts need files that come from steps whose code is not included here. You must supply these yourself:

| File | Used by | Origin |
| --- | --- | --- |
| `TCGA_HiChIP_hic/TCGA_powerlaw_results.csv` | `04_build_master_alpha_table.py`, most figures | TCGA HiChIP processing (`process_tcga_hichip.py`, not included) |
| `TCGA_HiChIP_hic/*_chr*_contacts.txt` | `supp_8panel_ps_comparison.py` | Per-chromosome contacts dumped from the TCGA `.hic` files |
| `gbm_patient_alpha_results.csv`, `gbm_patient_ps_curves.csv` | `fig4_tissue_matched.py`, `supp_8panel_ps_comparison.py`, `supp_tcga_alpha_boxplot.py` | Processing of the GSE229962 primary GBM samples (not included) |
| `GSE185069_analysis/cancer_powerlaw_results.csv` | `04_build_master_alpha_table.py` | Processing of GSE185069 (not included) |
| `*_frequencies.csv` (8 DIvA/OHT conditions) | `fig5_dsb_cohesin.py` | Cohesin knockdown and DSB contact-frequency tables |

---

## Installation

```bash
conda create -n hic-scaling python=3.10
conda activate hic-scaling
pip install "cooler>=0.9.3" "cooltools>=0.7.1" "hic2cool>=0.8.3" "hicstraw>=1.3.1" \
            "numpy>=1.24" "pandas>=2.0" "scipy>=1.11" "matplotlib>=3.7" requests openpyxl
```

`hic2cool` is called as a command-line tool by `02_process_cellline_hic.py`, so it must be on your `PATH`. `openpyxl` is needed only to open `data/S1_Table.xlsx` with pandas.

---

## Running the pipeline

Run the scripts in this order from the repository root, after setting the paths as described above.

```bash
# 1. (Optional) download reference .hic files
python preprocess/01_download_reference_hic.py --dry-run          # check URLs only
python preprocess/01_download_reference_hic.py --output <ROOT>    # download

# 2. Compute P(s) and fit α for cell lines
python preprocess/02_process_cellline_hic.py
python preprocess/03_process_gbm_celllines.py

# 3. Build the master α table (needs TCGA and GSE185069 results, see above)
python preprocess/04_build_master_alpha_table.py

# 4. Figures and statistics
python plotting/fig1_violin_tumour_vs_cellline.py   # also prints MWU / Welch / Cohen's d
python plotting/fig2_cellline_alpha_by_type.py
python plotting/fig3_cellline_ps_curves.py
python plotting/fig4_tissue_matched.py
python plotting/fig5_dsb_cohesin.py
python plotting/supp_8panel_ps_comparison.py
python plotting/supp_tcga_alpha_boxplot.py
```

Each figure is saved to `ROOT` as a 300 dpi `.png` and a `.pdf`.

---

## Methods summary

**P(s) estimation.** Genome-wide cis contact probability is estimated from intra-chromosomal contacts averaged at each diagonal offset d:

```
P(d) = (1 / N(d)) · Σ_i C(i, i+d)
```

N(d) is the number of valid bin pairs at offset d. Most samples use `cooltools.expected_cis`. The GSE280831 GBM lines use an equivalent `hicstraw` implementation on KR-normalised contacts. Bin sizes are 40 or 50 kb for most datasets and 10 kb for GSE253407.

**Power-law fit.** The curve is log-binned into 200 bins and fitted by OLS in log-log space over **100 kb to 10 Mb**, the mesoscale range covering TADs and loop extrusion:

```
log10 P(s) = −α · log10(s) + log10(A)
```

Fits need at least 5 points. R² is reported for each sample.

> **Sign convention:** the per-sample results files (`cellline_alpha_results.csv`, `gbm_cellline_alpha_results.csv`) store α as a **positive** number (α = −slope). `DISSERTATION_master_alpha_table.csv` stores the **negative** slope, where more negative means steeper decay. The figure scripts negate the positive values when they load them, so all figures and reported statistics use the negative convention.

**Statistics** (`fig1_violin_tumour_vs_cellline.py`): two-sided Mann-Whitney U, Welch's t-test and Cohen's d (pooled SD) for tumours vs cell lines.

### Key results

| Comparison | Result |
| --- | --- |
| Mann-Whitney U, tumour vs cell line | p = 1.81 × 10⁻⁵ |
| Effect size | Cohen's d = 0.87 |
| Median Δα (cell line − tumour, negative-α convention) | −0.18 |
| Cancer types where cell lines decay more steeply | 7 of 8 |
| GBM exception | Δα = +0.26 |

### Polymer interpretation

| Regime | Slope of P(s) | Spatial distance d(s) |
| --- | --- | --- |
| Loop-enriched / mitotic-like | ≈ −0.5 | ∝ s^(1/6) |
| Fractal globule | ≈ −1.0 | ∝ s^(1/3) |
| Equilibrium globule | ≈ −1.5 | ∝ s^(1/2) |

---

## Data sources

| Dataset | Description | Access |
| --- | --- | --- |
| TCGA Pan-Cancer HiChIP | 48 primary tumour samples, H3K27ac HiChIP | [GDC Portal](https://portal.gdc.cancer.gov) |
| GSE229962 | 30 primary GBM Hi-C samples | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE229962) |
| GSE148079 | Bladder cancer cell lines (BLCA) | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148079) |
| GSE167154 | Breast cancer cell lines (BRCA) | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167154) |
| GSE92804 | A549 lung adenocarcinoma (LUAD) | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92804) |
| GSE226216 | Liver cancer cell lines (LIHC) | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE226216) |
| GSE118629 | 22Rv1 prostate cancer (PRAD) | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118629) |
| GSE253407 | Droplet-HiC: melanoma lines (SKCM) and GBM_P4CC3 | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253407) |
| GSE133928 | SW480 colorectal (COAD), allValidPairs | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE133928) |
| GSE280831 | GBM glioma stem cell lines | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE280831) |
| GSE185069 | Pancreatic cancer cell lines | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185069) |
| GSE63525 / 4DN | Reference lines (GM12878, K562, IMR90, HMEC, NHEK, HUVEC; HCT116, HepG2) | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525), [4DN](https://data.4dnucleome.org) |

---

## Citation

```
Kooner JS. Cancer-wide analysis reveals systematic divergence in chromatin
contact scaling between primary tumours and cell line models.
MSci Dissertation, University of Surrey, 2026.
```

## Licence

MIT. See [LICENSE](LICENSE).
