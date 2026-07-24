# Integrity Risk and Income Analysis

## Overview

This repository contains the data, extraction scripts, analysis code, and supplementary materials for a **cross-sectional study based on secondary data** examining global patterns in the [SCImago IRIS](https://www.scimagoiris.com/) (Integrity Risk Indicators by SCImago) framework.

The study addresses three objectives:

1. **Global distribution** — Characterize the worldwide distribution of IRIS Overall scores and structural risk categories across 5,475 higher education institutions in 151 countries.
2. **Exploratory institutional profiles** — Identify multivariate indicator configurations using PCA and k-means clustering on the nine IRIS indicators.
3. **National context** — Evaluate associations between IRIS outcomes and country-level research system capacity (R&D expenditure, researchers per million, and World Bank income group) through ecological correlations, weighted regression, and multilevel modeling.

---

## Repository Structure

```
Integrity-Risk-Income-Analysis/
│
├── data/
│   └── Scimago_IRIS_Index_Data.csv          # Extracted SCImago IRIS dataset (5,475 institutions)
│
├── extraction data/
│   ├── integrity_risk_extraction_data.py     # Python web-scraping script (Google Colab)
│   └── Integrity_Risk_Extraction_Data.ipynb  # Jupyter Notebook version of the extraction
│
├── src/
│   └── code.R                                # Complete reproducible R analysis pipeline
│
├── supplementary material/
│   ├── Supplementary Dataset S1.csv          # Institution-level SCImago IRIS dataset
│   ├── Supplementary Dataset S2.csv          # Profile assignments and transformed indicator values
│   ├── Supplementary Figure S1.pdf            # Global distribution and country representation
│   ├── Supplementary Table S1.csv            # Overall dataset characteristics and IRIS categories
│   ├── Supplementary Table S2.csv            # Country-level IRIS summaries
│   ├── Supplementary Table S3.csv            # PCA variance decomposition and loadings
│   ├── Supplementary Table S4.csv            # K-means diagnostics and profile characteristics
│   ├── Supplementary Table S5.csv            # World Bank matching and analytical coverage
│   ├── Supplementary Table S6.csv            # Correlations, regression, and influence analyses
│   ├── Supplementary Table S7.csv            # Multilevel models and sensitivity analyses
│   └── Supplementary Table S8.csv            # Variable dictionary and analytical definitions
│
├── .gitignore
├── LICENSE                                   # MIT License
└── README.md
```

---

## Data Sources

| Source | Description | Access |
|--------|-------------|--------|
| **SCImago IRIS** | Institutional integrity risk indicators for higher education institutions worldwide. Indicators are z-score standardized, with positive values indicating above-average risk exposure. | [scimagoiris.com](https://www.scimagoiris.com/) |
| **World Bank WDI** | Country-level R&D expenditure (% of GDP) and researchers in R&D per million people. Income group classifications (low, lower-middle, upper-middle, high income). | [databank.worldbank.org](https://databank.worldbank.org/source/world-development-indicators) |

### SCImago IRIS Variables

The nine individual IRIS risk indicators used in the analysis are:

| Indicator | Description |
|-----------|-------------|
| Multi-affiliation | Proportion of output with multi-institutional affiliations |
| Retracted output | Proportion of retracted publications |
| Self-citation | Proportion of self-citations |
| Discontinued journals output | Proportion of output in discontinued journals |
| Hyperauthored output | Proportion of output with an unusually large number of authors |
| Leadership impact gap | Gap between leadership and non-leadership citation impact |
| Hyperprolific authors | Proportion of output by hyperprolific authors |
| Institutional journal output | Proportion of output in journals managed by the same institution |
| Redundant output | Proportion of output flagged as potentially redundant |

---

## Analytical Pipeline

The R analysis pipeline (`src/code.R`) is organized into three blocks:

### Block 1 — Global Descriptive Overview

- Data cleaning, type validation, and missing value assessment
- Global summary statistics (Table 1)
- Distribution of IRIS structural risk categories (very low, low, medium, significant)
- Top 20 countries by number of institutions
- Country-level aggregated summaries
- Top 20 countries by proportion of institutions in the significant category (≥10 institutions filter)
- Choropleth maps of mean IRIS Overall score and percentage of institutions in the significant category

### Block 2 — Exploratory Institutional Profiles

- Winsorization (1st–99th percentile) and z-score standardization of the nine IRIS indicators
- Spearman correlation matrix among indicators
- Principal Component Analysis (PCA) on the standardized indicator matrix
- K-means clustering diagnostics (k = 2–8; elbow method + silhouette width)
- Final k = 5 clustering solution with profile assignment
- Profile summary table with dominant indicators
- Profile heatmap, PCA projection, and IRIS structural risk category distribution by profile

### Block 3 — National Context and Modeling

- World Bank indicator download via the `WDI` package (most recent non-missing values from 2015 onward)
- IRIS–World Bank matching using ISO3 country codes
- Descriptive analysis by income group
- Four prespecified Spearman correlations with 10,000 bootstrap resamples, asymptotic P values, and Holm adjustment
- Weighted country-level linear regression with conventional model-based inference
- Cook's distance, leverage, standardized residuals, and leave-one-country-out analyses
- Multilevel linear model with a country-level random intercept
- Sensitivity analyses excluding the verified maximum IRIS Overall score and restricting the sample to countries with ≥10 institutions
- Exploratory institutional profile distribution across World Bank income groups

---

## Reproduction Instructions

### Prerequisites

- **R** ≥ 4.4.2
- **Python** ≥ 3.8 (only for data extraction; the extracted dataset is already provided)

### R Packages

The following R packages are required:

```
tidyverse, janitor, broom, lme4, broom.mixed, performance,
WDI, sf, rnaturalearth, rnaturalearthdata, cluster,
patchwork, scales, viridis
```

### Steps to Reproduce

1. **Clone the repository**
   ```bash
   git clone https://github.com/AndyPaaai/Integrity-Risk-Income-Analysis.git
   cd Integrity-Risk-Income-Analysis
   ```

2. **Run the analysis**
   From the repository root, run:
   ```bash
   Rscript src/code.R
   ```
   Alternatively, open `src/code.R` in RStudio and run the entire script. The pipeline will:
   - Validate the required packages and input data
   - Create the `outputs_revised/` directory and its subdirectories
   - Generate manuscript tables, figures, supplementary tables, supplementary datasets, and session information

   The script also accepts explicit paths:
   ```bash
   Rscript src/code.R \
     data/Scimago_IRIS_Index_Data.csv \
     data/block3_wdi_raw_download.csv \
     outputs_revised
   ```

   The second argument is an optional fixed World Bank WDI snapshot. If it is omitted or unavailable, the script downloads the indicators and saves a new snapshot. Exact reproduction of the published country-level results requires the archived snapshot used in the study.

3. **Data extraction** *(optional)*
   If you want to re-extract the SCImago IRIS data from the website, run the Python extraction script:
   ```bash
   cd "extraction data"
   python integrity_risk_extraction_data.py
   ```
   > **Note:** The extraction depends on the current structure of the SCImago IRIS website and may break if the site layout changes.

### Output

By default, generated outputs are saved to the `outputs_revised/` directory (not tracked by Git):

- `outputs_revised/tables/` — Summary statistics, model results, diagnostics, and intermediate tables
- `outputs_revised/figures/` — Main and supplementary figures in PNG and PDF formats
- `outputs_revised/supplementary/` — Supplementary Tables S1–S8 as sectioned CSV files
- `outputs_revised/supplementary_data/` — Cleaned Supplementary Data S1 and S2
- `outputs_revised/analysis_session_info.txt` — R session and package information

---

## Supplementary Materials

The `supplementary material/` directory contains all supplementary datasets, tables, and the supplementary figure referenced in the manuscript:

| File | Content |
|------|---------|
| **Supplementary Dataset S1** | Cleaned institution-level SCImago IRIS dataset for 5,475 institutions |
| **Supplementary Dataset S2** | Institution-level profile assignments, original cluster labels, and winsorized indicator values for 5,470 complete cases |
| **Supplementary Figure S1** | Global distribution of IRIS structural risk categories and country representation in the dataset |
| **Supplementary Table S1** | Overall characteristics of the institutional dataset and distribution across IRIS structural risk categories |
| **Supplementary Table S2** | Country-level summary of the SCImago IRIS dataset |
| **Supplementary Table S3** | PCA variance decomposition, cumulative variance, eigenvalues, and indicator loadings |
| **Supplementary Table S4** | K-means diagnostics and characteristics of the five exploratory institutional profiles |
| **Supplementary Table S5** | IRIS–World Bank matching, contextual-data coverage, indicator years, and geographical map matching |
| **Supplementary Table S6** | Country-level Spearman associations, conventional weighted regression, leave-one-country-out estimates, and influence diagnostics |
| **Supplementary Table S7** | Multilevel model estimates, model fit, and influence sensitivity analyses |
| **Supplementary Table S8** | Variable dictionary and analytical definitions for raw and derived variables |

Supplementary tables containing multiple analytical panels use a `section` column to preserve the original workbook structure in CSV format.

---

## License

This project is licensed under the [MIT License](LICENSE).

© 2026 Andy A. Acosta-Monterrosa
