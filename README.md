# Global patterns of institutional research integrity risk

## Project Description
This repository contains the data and code for the meta-research analysis of the **SCImago IRIS Index** (Integrity Risk Index) evaluated across different income groups and countries. The objective of this study is to characterize global patterns of institutional research integrity risk and to examine whether institutional risk patterns are associated with country-level research system context.

Research integrity is increasingly understood as a system-level issue shaped by institutional incentives, publication cultures, research evaluation practices, and national scientific capacity. We conducted a cross-sectional meta-research analysis of 5,475 higher education institutions from 151 countries, linking institutional-level data from SCImago IRIS with country-level contextual data from the World Bank (e.g., R&D expenditure, researchers in R&D per million people, and income group).

The original integrity indicators come from the official platform: [SCImago IRIS](https://www.scimagoiris.com/).

## Repository Contents

The repository is organized into the following folders:

- `data/`: Contains the main raw dataset `Scimago_IRIS_Index_Data.csv` used for the analysis.
- `extraction data/`: Contains the Python scripts (`integrity_risk_extraction_data.py`) and Jupyter Notebooks (`Integrity_Risk_Extraction_Data.ipynb`) used for the initial web scraping and API extraction of the IRIS data.
- `src/`: Contains the main analytic pipeline (`code.R`).
- `supplementary material/`: Contains all supplementary datasets and tables (S1 to S7) referenced in the study.
- `outputs/` *(Generated automatically)*: All tables and high-resolution figures produced by the script will be saved here.

## How to Reproduce the Analysis

The complete analysis pipeline, from data cleaning to statistical modeling (PCA, K-means clustering, Multilevel modeling) and visualization, is contained in a single R script. 

To reproduce the analysis:
1. Clone this repository to your local machine.
2. Open your R console or RStudio and ensure your working directory is set to the root of this repository.
3. Open `src/code.R` and run it entirely (e.g., using `source("src/code.R")`).
4. The script will automatically:
   - Install any missing R packages.
   - Read the dataset from `data/Scimago_IRIS_Index_Data.csv`.
   - Create the `outputs/tables/` and `outputs/figures/` folders if they don't exist.
   - Run Block 1 (Global Patterns), Block 2 (Risk Profiles), and Block 3 (National Context).
   - Save all generated tables (`.csv`) and figures (`.png`/`.pdf`) directly into the `outputs/` folder.

## License
Refer to the `LICENSE` file for the terms of use of this codebase.
