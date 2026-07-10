# Integrity Risk and Income Analysis

## Overview

This repository contains the data, extraction scripts, analysis code, and supplementary materials for a **cross-sectional meta-research study** examining global patterns of institutional research integrity risk using the [SCImago IRIS](https://www.scimagoiris.com/) (Integrity Risk Indicators by SCImago) platform.

The study addresses three objectives:

1. **Global distribution** — Characterize the worldwide distribution of institutional integrity risk across 5,475 higher education institutions in 151 countries.
2. **Institutional risk profiles** — Identify multivariate risk profiles using PCA and k-means clustering on the nine IRIS indicators.
3. **National context** — Evaluate associations between institutional integrity risk and country-level research system capacity (R&D expenditure, researchers per million, World Bank income group) through ecological correlations, weighted regression, and multilevel modeling.

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
│   └── code.R                               # Full R analysis pipeline (Blocks 1–3)
│
├── supplementary material/
│   ├── Supplementary Dataset S1.csv          # Institution-level SCImago IRIS dataset
│   ├── Supplementary Dataset S2.csv          # Institution-level profile assignments (k-means)
│   ├── Supplementary Table S1.csv            # Variable dictionary
│   ├── Supplementary Table S2.csv            # Country-level IRIS summaries
│   ├── Supplementary Table S3.xlsx           # PCA variance decomposition and loadings
│   ├── Supplementary Table S4.xlsx           # K-means diagnostics and profile centroids
│   ├── Supplementary Table S5.xlsx           # Country-level correlations and regression
│   ├── Supplementary Table S6.xlsx           # Multilevel model estimates
│   └── Supplementary Table S7.xlsx           # Matching and data coverage summaries
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
- Distribution of IRIS risk categories (very low, low, medium, significant)
- Top 20 countries by number of institutions
- Country-level aggregated summaries
- Top 20 countries by proportion of significant-risk institutions (≥10 institutions filter)
- Choropleth maps: mean Overall score and % significant risk by country

### Block 2 — Institutional Risk Profiles

- Winsorization (1st–99th percentile) and z-score standardization of the nine IRIS indicators
- Spearman correlation matrix among indicators
- Principal Component Analysis (PCA) on the standardized indicator matrix
- K-means clustering diagnostics (k = 2–8; elbow method + silhouette width)
- Final k = 5 clustering solution with profile assignment
- Profile summary table with dominant indicators
- Profile heatmap, PCA projection by profile, and IRIS category distribution by profile

### Block 3 — National Context and Modeling

- World Bank indicator download via the `WDI` package (most recent non-missing values from 2015 onward)
- IRIS–World Bank matching using ISO3 country codes
- Descriptive analysis by income group (boxplots, scatterplots with R&D and researchers)
- Spearman/Pearson correlations between national capacity and IRIS risk
- Weighted country-level linear regression (mean Overall ~ R&D + researchers + income group + institutional count)
- Multilevel linear model (institution-level Overall ~ output + SIR rank + R&D + researchers + income group, random intercept by country)
- Institutional risk profile distribution across income groups

---

## Reproduction Instructions

### Prerequisites

- **R** ≥ 4.4.0
- **Python** ≥ 3.8 (only for data extraction; the extracted dataset is already provided)

### R Packages

The following R packages are required and will be installed automatically if missing:

```
tidyverse, janitor, readr, stringr, forcats, scales, sf,
rnaturalearth, rnaturalearthdata, viridis, patchwork,
cluster, WDI, countrycode, lme4, broom.mixed, performance
```

### Steps to Reproduce

1. **Clone the repository**
   ```bash
   git clone https://github.com/AndyPaaai/Integrity-Risk-Income-Analysis.git
   cd Integrity-Risk-Income-Analysis
   ```

2. **Run the analysis**
   Open `src/code.R` in RStudio (or any R environment) and run the entire script. The pipeline will:
   - Install any missing R packages
   - Create the `outputs/` directory with `tables/` and `figures/` subdirectories
   - Generate all tables (CSV) and figures (PNG/PDF)

3. **Data extraction** *(optional)*
   If you want to re-extract the SCImago IRIS data from the website, run the Python extraction script:
   ```bash
   cd "extraction data"
   python integrity_risk_extraction_data.py
   ```
   > **Note:** The extraction depends on the current structure of the SCImago IRIS website and may break if the site layout changes.

### Output

All generated outputs are saved to the `outputs/` directory (not tracked by Git):

- `outputs/tables/` — CSV files with summary statistics, model results, and intermediate tables
- `outputs/figures/` — PNG and PDF files with all figures reported in the manuscript

---

## Supplementary Materials

The `supplementary material/` directory contains all supplementary datasets and tables referenced in the manuscript:

| File | Content |
|------|---------|
| **Dataset S1** | Full institution-level SCImago IRIS dataset (CSV) |
| **Dataset S2** | Institution-level k-means profile assignments (CSV) |
| **Table S1** | Variable dictionary with indicator descriptions |
| **Table S2** | Country-level IRIS summaries |
| **Table S3** | PCA variance decomposition and indicator loadings |
| **Table S4** | K-means clustering diagnostics and profile centroids |
| **Table S5** | Country-level correlation and weighted regression results |
| **Table S6** | Multilevel model fixed and random effects estimates |
| **Table S7** | Data matching and coverage summaries (IRIS–World Bank–Natural Earth) |

---

## License

This project is licensed under the [MIT License](LICENSE).

© 2026 Andy A. Acosta-Monterrosa
