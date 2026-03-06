# Integrity Risk Income Analysis

## Project Description
This repository hosts a descriptive study and cross-sectional analysis of the **SCImago IRIS Index** (Integrity Risk Index), evaluated across different income groups and countries. The main objective of the project is to quantitatively explore and understand how scientific integrity risk is distributed and varies globally, depending on the income level and socioeconomic characteristics of each nation.

The original data and metrics come from the official platform: [SCImago IRIS](https://www.scimagoiris.com/).

## Current Status and Methodology
So far, the project has completed its initial phase of data collection and structuring:
- **Data Extraction (Web Scraping / API):** Python code was developed (using the Google Colab environment) to query, extract, and clean the index data directly from the SCImago IRIS website.
- **Storage and Structuring:** The raw extracted data has been consolidated and saved locally in CSV format (`Scimago_IRIS_Index_Data.csv`). This structured database will serve as the main input for future statistical analysis, modeling, and cross-sectional visualizations.

## Repository Contents
- `integrity_risk_extraction_data.py` / `Integrity_Risk_Extraction_Data.ipynb`: Python scripts and Jupyter Notebooks used for the initial extraction and processing of information from the web.
- `Scimago_IRIS_Index_Data.csv`: The final processed dataset resulting from the extraction.
- `README.md`: This file, which documents the purpose, evolution, and general structure of the project.
