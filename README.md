# Integrity Risk and Income Analysis

## Overview

This repository contains the archived data snapshots, extraction scripts, analysis code, and supplementary materials for a **cross-sectional study based on secondary data** examining global patterns in the [SCImago IRIS](https://www.scimagoiris.com/) (Integrity Risk Indicators by SCImago) framework.

The study addresses three objectives:

1. **Global distribution** — Characterize the worldwide distribution of IRIS Overall scores and structural risk categories across 5,475 higher education institutions in 151 countries.
2. **Exploratory institutional profiles** — Identify multivariate indicator configurations using PCA and k-means clustering on the nine IRIS indicators.
3. **National context** — Evaluate associations between IRIS outcomes and country-level research system capacity (R&D expenditure, researchers per million, and World Bank income group) through ecological correlations, weighted regression, and multilevel modeling.

---

## Repository Structure

```text
Integrity-Risk-Income-Analysis/
│
├── data/
│   ├── Scimago_IRIS_Index_Data.csv           # Archived SCImago IRIS snapshot used in the study (5,475 institutions)
│   └── block3_wdi_raw_download.csv            # Archived World Bank WDI snapshot used in the study
│
├── extraction data/
│   ├── integrity_risk_extraction_data.py      # Optional Python re-extraction script
│   └── Integrity_Risk_Extraction_Data.ipynb   # Notebook version of the extraction workflow
│
├── src/
│   └── code.R                                 # Complete reproducible R analysis pipeline
│
├── supplementary material/
│   ├── Supplementary Dataset S1.csv           # Institution-level SCImago IRIS dataset
│   ├── Supplementary Dataset S2.csv           # Profile assignments and transformed indicator values
│   ├── Supplementary Figure S1.pdf            # Global distribution and country representation
│   ├── Supplementary Table S1.csv             # Dataset characteristics, IRIS categories, and category-boundary check
│   ├── Supplementary Table S2.csv             # Country-level IRIS summaries
│   ├── Supplementary Table S3.csv             # PCA variance decomposition and loadings
│   ├── Supplementary Table S4.csv             # K-means profiles plus descriptive AHP decomposition of Profile 5
│   ├── Supplementary Table S5.csv             # World Bank matching and analytical coverage
│   ├── Supplementary Table S6.csv             # Correlations, predictor-correlation check, regression, and influence analyses
│   ├── Supplementary Table S7.csv             # Multilevel models, sensitivities, and extreme-observation AHP decomposition
│   └── Supplementary Table S8.csv             # Variable dictionary and analytical definitions
│
├── analysis_session_info.txt                  # Session/package versions from the archived-input reproduction run
├── .gitignore
├── LICENSE                                    # MIT License for repository code
└── README.md
```

---

## Archived Data Snapshots and Exact Reproduction

Exact reproduction of the study uses the two archived input files in `data/`:

- `data/Scimago_IRIS_Index_Data.csv`: SCImago IRIS institutional data extracted on **4 March 2026**. This is the 5,475-institution dataset used in the manuscript.
- `data/block3_wdi_raw_download.csv`: fixed World Bank WDI snapshot used for the country-level analyses. The World Bank resources were accessed on **8 July 2026**; the archived file contains World Bank source metadata with `lastupdated = 2026-07-01`. For each contextual indicator, the analysis selects the most recent non-missing value from 2015 onward. In the study data, selected years ranged from 2015 to 2024.

The archived snapshots, rather than a fresh web download, should be used for exact reproduction. SCImago IRIS and World Bank data may change after the study date.

Re-running the SCImago extraction script is therefore **optional** and should not be used to replace the archived study dataset when reproducing the published analyses. Likewise, if the fixed WDI snapshot is not supplied, `src/code.R` can request current WDI data, but the result may differ from the archived analysis.

---

## Data Sources

| Source | Description | Access |
|---|---|---|
| **SCImago IRIS** | Institutional research-integrity risk indicators. The study uses the archived institutional snapshot described above. | [scimagoiris.com](https://www.scimagoiris.com/) |
| **World Bank WDI** | Country-level R&D expenditure (% GDP), researchers in R&D per million people, and associated country metadata. | [World Development Indicators](https://databank.worldbank.org/source/world-development-indicators) |
| **World Bank Country and Lending Groups** | Income-group metadata used for contextual analyses. | [World Bank country classifications](https://datahelpdesk.worldbank.org/knowledgebase/articles/906519) |

### SCImago IRIS Variables

The nine individual IRIS indicators used in the analysis are:

| Indicator | Description |
|---|---|
| Multi-affiliation | Proportion of output in which at least one author lists multiple institutional affiliations |
| Retracted output | Proportion of institutional output that has been formally retracted |
| Self-citation | Proportion of citations to institutional output originating from the same institution |
| Discontinued journals output | Proportion of output published in journals discontinued from major bibliographic databases |
| Hyperauthored output | Proportion of output with exceptionally large author lists |
| Leadership impact gap | Difference between the institution's overall normalized impact and the impact of output in which it holds a leading-authorship position |
| Hyperprolific authors | Proportion of institutional contributions linked to authors exceeding the IRIS annual productivity threshold |
| Institutional journal output | Proportion of output published in journals owned by the same institution |
| Redundant output | Proportion of output flagged for unusually high bibliographic overlap with other works by the same authors within a short period |

`Overall` and the structural risk category (`very low`, `low`, `medium`, `significant`) are reported as supplied by SCImago IRIS. In particular, **“significant” is a SCImago-provided structural risk category; it is not defined in this repository by a Tukey outlier rule or by a hypothesis-test P value.**

### AHP weights used for descriptive decomposition

SCImago's IRIS methodology combines the nine original indicators using Analytic Hierarchy Process (AHP) weights. The repository uses these weights only to reproduce/decompose the supplied `Overall` score:

| Indicator | AHP weight |
|---|---:|
| Multi-affiliation | 0.10 |
| Retracted output | 0.30 |
| Self-citation | 0.05 |
| Discontinued journals output | 0.20 |
| Hyperauthored output | 0.05 |
| Leadership impact gap | 0.05 |
| Hyperprolific authors | 0.10 |
| Institutional journal output | 0.10 |
| Redundant output | 0.05 |

Methodology: https://www.scimagoiris.com/integrity-risk-methodology.php

These weights sum to 1.00. They are applied to the **original supplied indicator values** for post hoc decomposition and rounding-tolerant validation only. They are not applied to winsorized (`*_w`) or restandardized clustering variables and do not change PCA or k-means.

---

## Analytical Pipeline

The R analysis pipeline (`src/code.R`) is organized into three analytical blocks.

### Block 1 — Global Descriptive Overview

- Data cleaning, type validation, and missing-value assessment
- Global summary statistics
- Distribution of the four SCImago-supplied structural risk categories
- Independent check of the observed `medium`/`significant` boundary and comparison with a Tukey upper fence
- Country-level aggregated summaries and representation plots

### Block 2 — Exploratory Institutional Profiles

- Winsorization (1st–99th percentile) and z-score standardization of the nine IRIS indicators
- Spearman correlation matrix among indicators
- PCA on the winsorized and standardized indicators
- K-means diagnostics for k = 2–8
- Final k = 5 clustering solution
- Profile summaries, centroids, dominant indicators, and risk-category distributions
- **Post hoc descriptive AHP decomposition of Profile 5 using the original IRIS indicator values and the published indicator weights**

The AHP weights are used only for descriptive decomposition of the supplied IRIS `Overall` score. They are **not** applied to the winsorized or standardized clustering variables, and PCA/k-means remain unweighted.

### Block 3 — National Context and Modeling

- Fixed World Bank snapshot linkage by ISO3 code
- Descriptive analysis by income group
- Four prespecified Spearman correlations with 10,000 bootstrap resamples and Holm adjustment
- Pearson correlation matrix among the three continuous weighted-regression predictors
- Weighted country-level linear regression with conventional model-based inference
- Cook's distance, leverage, standardized residuals, and leave-one-country-out analyses
- Multilevel linear model with a country random intercept
- Sensitivities excluding the verified maximum IRIS Overall score and restricting to countries represented by at least 10 institutions
- Descriptive AHP decomposition of the verified extreme Overall observation

The code does **not** replace the SCImago-supplied `Overall` score with a reconstructed score. A weighted reconstruction is used only as a reproducibility check. Small differences (maximum approximately 0.0008 in the archived data) are expected because the publicly available indicators and `Overall` values are rounded; a tolerance of 0.001 is used.

---

## Reproduction Instructions

### Prerequisites

- **R 4.4.2** (the manuscript analysis version)
- **Python ≥ 3.8** only if the optional SCImago re-extraction workflow is used

### R Packages

```text
tidyverse, janitor, broom, lme4, broom.mixed, performance,
WDI, sf, rnaturalearth, rnaturalearthdata, cluster,
patchwork, scales, viridis
```

The analysis script writes `outputs_revised/analysis_session_info.txt` at the end of a successful run. This records the R session and package versions in the environment in which the script was actually executed. No pre-existing session-info file is required to run the repository. For the versioned repository/Zenodo release, copy that generated file to the repository root as `analysis_session_info.txt` so the exact software environment used for the released outputs is archived.

### Exact Reproduction

From the repository root:

```bash
Rscript src/code.R \
  data/Scimago_IRIS_Index_Data.csv \
  data/block3_wdi_raw_download.csv \
  outputs_revised
```

Opening `src/code.R` in RStudio and running the whole script with the same archived inputs is also valid.

The pipeline creates the output directories, runs internal validation checks, and generates manuscript tables, figures, supplementary tables, supplementary datasets, and `analysis_session_info.txt`.

### Optional SCImago Re-extraction

```bash
cd "extraction data"
python integrity_risk_extraction_data.py
```

The script writes its CSV output to the **current working directory** (`extraction data/` when run exactly as above), not automatically to `data/`. A newly extracted file may differ from the 4 March 2026 archived study snapshot because the live SCImago IRIS site can change. Do not overwrite `data/Scimago_IRIS_Index_Data.csv` if exact reproduction is the goal.

---

## Generated Output

By default, outputs are saved to `outputs_revised/` and are not intended to be tracked as an additional duplicate of every intermediate object:

- `outputs_revised/tables/` — summary statistics, model results, diagnostics, and intermediate tables
- `outputs_revised/figures/` — main and supplementary figures
- `outputs_revised/supplementary/` — Supplementary Tables S1–S8 as sectioned CSV files
- `outputs_revised/supplementary_data/` — cleaned Supplementary Datasets S1 and S2
- `outputs_revised/analysis_session_info.txt` — R session and package information generated by the run

For the versioned repository deposit, the final Supplementary Tables S1–S8 and Supplementary Datasets S1–S2 referenced by the manuscript are retained in `supplementary material/`. The `analysis_session_info.txt` produced by the same exact-reproduction run is also retained at the repository root.

---

## Supplementary Materials

| File | Content |
|---|---|
| **Supplementary Dataset S1** | Cleaned institution-level SCImago IRIS dataset for 5,475 institutions |
| **Supplementary Dataset S2** | Institution-level profile assignments, original cluster labels, and winsorized indicator values for 5,470 complete cases |
| **Supplementary Figure S1** | Global distribution of IRIS structural risk categories and country representation |
| **Supplementary Table S1** | Overall dataset characteristics, SCImago structural risk categories, and independent category-boundary/Tukey comparison |
| **Supplementary Table S2** | Country-level summary of the SCImago IRIS dataset |
| **Supplementary Table S3** | PCA variance decomposition, cumulative variance, eigenvalues, and indicator loadings |
| **Supplementary Table S4** | K-means diagnostics and five-profile characteristics, including the descriptive AHP decomposition of Profile 5 |
| **Supplementary Table S5** | IRIS–World Bank matching, contextual-data coverage, indicator years, and geographical map matching |
| **Supplementary Table S6** | Country-level Spearman associations, Pearson correlations among continuous regression predictors, conventional weighted regression, leave-one-country-out estimates, and influence diagnostics |
| **Supplementary Table S7** | Multilevel model estimates, model fit, sensitivity analyses, influential-observation details, and the AHP decomposition of the verified maximum Overall observation |
| **Supplementary Table S8** | Variable dictionary and analytical definitions, including `ahp_weight` and `weighted_indicator_contribution` |

Supplementary tables containing multiple analytical panels use a `section` column to preserve the panel structure in CSV format.

---

## Licensing and Third-Party Data Terms

### Repository code

Original code in this repository is released under the [MIT License](LICENSE).

### SCImago IRIS data

SCImago IRIS is a third-party data source and is **not relicensed under the repository's MIT License**. SCImago states on the IRIS website that the information shown there may be used for non-commercial purposes provided that SCImago IRIS is cited. Users should consult the current SCImago IRIS terms before reuse or redistribution:

https://www.scimagoiris.com/

### World Bank WDI data

The archived WDI snapshot is third-party World Bank data and is **not relicensed under MIT**. World Development Indicators is currently listed by the World Bank under **Creative Commons Attribution 4.0 (CC BY 4.0)**, subject to the World Bank's dataset terms and any indicator-specific restrictions:

https://datacatalog.worldbank.org/search/dataset/0037712/world-development-indicators

https://www.worldbank.org/ext/en/legal/terms-conditions/datasets

Users are responsible for complying with the applicable source-data licenses and attribution requirements.

---

© 2026 Andy A. Acosta-Monterrosa
