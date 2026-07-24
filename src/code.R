# ==============================================================================
# SCImago IRIS: complete reproducible analysis 
# ==============================================================================
#
# Target software: R 4.4.2
#
# Recommended command:
# Rscript SCImago_IRIS_complete_revised_analysis.R \
#   data/Scimago_IRIS_Index_Data.csv \
#   data/block3_wdi_raw_download.csv \
#   outputs
#
# Arguments:
#   1. Institution-level SCImago IRIS CSV.
#   2. Fixed World Bank WDI snapshot CSV. If omitted or unavailable, the script
#      downloads the indicators and writes a new snapshot; exact reproduction of
#      the published results requires the archived snapshot.
#   3. Output directory.
#
# Statistical decisions retained in the revised manuscript:
#   - Four prespecified Spearman correlations.
#   - Two-sided asymptotic P values.
#   - Percentile 95% CIs from 10,000 country-level bootstrap resamples.
#   - Holm adjustment across the four correlation tests.
#   - Conventional model-based inference for the weighted regression.
#   - Leave-one-country-out and influence diagnostics.
#   - Primary multilevel model retaining the verified Overall = 37.599 value.
#   - Sensitivities excluding that institution and restricting to countries
#     represented by at least 10 institutions.
#
# The script writes Supplementary Tables S1-S8 as sectioned rectangular CSV
# files. For tables with several analytical panels, the `section` column
# identifies the original panel or worksheet.

# ------------------------------------------------------------------------------
# 0. Setup
# ------------------------------------------------------------------------------

required_packages <- c(
  "tidyverse", "janitor", "broom", "lme4", "broom.mixed", "performance",
  "WDI", "sf", "rnaturalearth", "rnaturalearthdata", "cluster", "patchwork",
  "scales", "viridis"
)

missing_packages <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages) > 0) {
  stop(
    "Install these packages before running the analysis: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(broom)
  library(lme4)
  library(broom.mixed)
  library(performance)
  library(WDI)
  library(sf)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(cluster)
  library(patchwork)
  library(scales)
  library(viridis)
})

options(
  stringsAsFactors = FALSE,
  contrasts = c("contr.treatment", "contr.poly")
)
RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")

args <- commandArgs(trailingOnly = TRUE)

first_existing <- function(paths) {
  hits <- paths[file.exists(paths)]
  if (length(hits) == 0) NA_character_ else hits[[1]]
}

iris_file <- if (length(args) >= 1) {
  args[[1]]
} else {
  first_existing(c(
    "data/Scimago_IRIS_Index_Data.csv",
    "Scimago_IRIS_Index_Data.csv",
    "upload/Scimago_IRIS_Index_Data.csv"
  ))
}

wdi_snapshot_file <- if (length(args) >= 2) {
  args[[2]]
} else {
  first_existing(c(
    "data/block3_wdi_raw_download.csv",
    "block3_wdi_raw_download.csv",
    "upload/block3_wdi_raw_download.csv"
  ))
}

output_root <- if (length(args) >= 3) args[[3]] else "outputs_revised"
table_dir <- file.path(output_root, "tables")
figure_dir <- file.path(output_root, "figures")
supp_dir <- file.path(output_root, "supplementary")
data_dir <- file.path(output_root, "supplementary_data")

walk(
  c(output_root, table_dir, figure_dir, supp_dir, data_dir),
  ~ dir.create(.x, recursive = TRUE, showWarnings = FALSE)
)

if (is.na(iris_file) || !file.exists(iris_file)) {
  stop(
    "The SCImago IRIS input CSV was not found. Supply it as the first argument.",
    call. = FALSE
  )
}

message("SCImago IRIS input: ", normalizePath(iris_file))
if (!is.na(wdi_snapshot_file) && file.exists(wdi_snapshot_file)) {
  message("Fixed World Bank snapshot: ", normalizePath(wdi_snapshot_file))
} else {
  warning(
    "No fixed World Bank snapshot was found. A new WDI download will be used; ",
    "this may not exactly reproduce the archived results."
  )
}

# ------------------------------------------------------------------------------
# 1. Helper functions and labels
# ------------------------------------------------------------------------------

risk_indicator_vars <- c(
  "multi_affiliation",
  "retracted_output",
  "self_citation",
  "discontinued_journals_output",
  "hyperauthored_output",
  "leadership_impact_gap",
  "hyperprolific_authors",
  "institutional_journal_output",
  "redundant_output"
)

numeric_vars <- c("sir_rank", "overall", "output", risk_indicator_vars)

indicator_labels <- c(
  multi_affiliation = "Multi-affiliation",
  retracted_output = "Retracted output",
  self_citation = "Self-citation",
  discontinued_journals_output = "Discontinued journals output",
  hyperauthored_output = "Hyperauthored output",
  leadership_impact_gap = "Leadership impact gap",
  hyperprolific_authors = "Hyperprolific authors",
  institutional_journal_output = "Institutional journal output",
  redundant_output = "Redundant output"
)

iris_colors <- c(
  "very low" = "#9EC753",
  "low" = "#E3CE47",
  "medium" = "#F49D4A",
  "significant" = "#D95C4F"
)

profile_colors <- c(
  "Profile 1" = "#0072B2",
  "Profile 2" = "#D55E00",
  "Profile 3" = "#009E73",
  "Profile 4" = "#CC79A7",
  "Profile 5" = "#E69F00"
)

valid_income_groups <- c(
  "Low income",
  "Lower middle income",
  "Upper middle income",
  "High income"
)

theme_paper <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(family = "sans", color = "black"),
      plot.title = element_text(
        face = "bold", size = base_size + 2, margin = margin(b = 8)
      ),
      plot.subtitle = element_text(
        size = base_size, color = "grey30", margin = margin(b = 8)
      ),
      plot.caption = element_text(
        size = base_size - 2, color = "grey45", hjust = 0,
        margin = margin(t = 10)
      ),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      panel.grid = element_blank()
    )
}

z_score <- function(x) as.numeric(scale(x))

winsorize_vec <- function(x, probs = c(0.01, 0.99)) {
  limits <- quantile(x, probs = probs, na.rm = TRUE, type = 7)
  pmin(pmax(x, limits[[1]]), limits[[2]])
}

bootstrap_spearman <- function(x, y, seed, replicates = 10000) {
  complete <- complete.cases(x, y)
  x <- x[complete]
  y <- y[complete]
  n <- length(x)

  set.seed(seed)
  estimates <- replicate(
    replicates,
    {
      index <- sample.int(n, size = n, replace = TRUE)
      suppressWarnings(cor(x[index], y[index], method = "spearman"))
    }
  )

  tibble(
    ci_low = unname(quantile(estimates, 0.025, na.rm = TRUE, type = 7)),
    ci_high = unname(quantile(estimates, 0.975, na.rm = TRUE, type = 7))
  )
}

write_sectioned_csv <- function(
  table_number,
  title,
  description,
  notes,
  sections,
  file
) {
  combined <- imap_dfr(
    sections,
    function(section_data, section_name) {
      section_data <- as_tibble(section_data, .name_repair = "unique")
      section_data |>
        mutate(
          supplementary_table = paste0("S", table_number),
          section = section_name,
          section_row = row_number(),
          title = if_else(row_number() == 1L, title, NA_character_),
          description = if_else(row_number() == 1L, description, NA_character_),
          notes = if_else(row_number() == 1L, paste(notes, collapse = " | "), NA_character_),
          .before = 1
        )
    }
  ) |>
    relocate(title, description, notes, .after = last_col())

  write_csv(combined, file, na = "")
  invisible(combined)
}

fixed_effect_table <- function(model, model_label = NULL, reference_label) {
  output <- broom.mixed::tidy(model, effects = "fixed", conf.int = FALSE) |>
    transmute(
      term,
      estimate,
      se = std.error,
      z_statistic = estimate / se,
      p_value = 2 * pnorm(abs(z_statistic), lower.tail = FALSE),
      ci_low = estimate - qnorm(0.975) * se,
      ci_high = estimate + qnorm(0.975) * se
    )

  if (!is.null(model_label)) {
    output <- output |> mutate(model = model_label, .before = 1)
  }
  output |> mutate(income_group_reference = reference_label, .before = 1)
}

random_effect_table <- function(model, model_label) {
  as.data.frame(VarCorr(model)) |>
    as_tibble() |>
    transmute(
      model = model_label,
      random_component = if_else(
        grp == "Residual", "Residual", "Country random intercept"
      ),
      variance = vcov,
      sd = sdcor
    )
}

model_fit_row <- function(model, model_label, data) {
  r2_values <- performance::r2_nakagawa(model)
  icc_value <- performance::icc(model)

  tibble(
    model = model_label,
    institutions_n = nrow(data),
    countries_n = n_distinct(data$country),
    marginal_r2 = unname(r2_values$R2_marginal),
    conditional_r2 = unname(r2_values$R2_conditional),
    adjusted_icc = unname(icc_value$ICC_adjusted)
  )
}

rename_fixed_terms <- function(data) {
  term_labels <- c(
    "(Intercept)" = "(Intercept)",
    "output_log_z" = "Log scientific output, standardized",
    "sir_rank_log_z" = "Log SIR Rank, standardized",
    "rd_gdp_z" = "R&D expenditure, standardized",
    "researchers_pm_log_z" = "Log researcher density, standardized",
    "income_group_wbLower middle income" = "Lower-middle-income group",
    "income_group_wbUpper middle income" = "Upper-middle-income group",
    "income_group_wbHigh income" = "High-income group"
  )
  data |> mutate(term = recode(term, !!!term_labels))
}

# ------------------------------------------------------------------------------
# 2. Import, clean, and validate the institution-level dataset
# ------------------------------------------------------------------------------

raw_data <- read_csv(
  iris_file,
  na = c("", "NA", "background: #ffffff; color: #000000;"),
  show_col_types = FALSE
) |>
  clean_names()

missing_columns <- setdiff(
  c("id", "institution", "country", "risk", numeric_vars),
  names(raw_data)
)
if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ", paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

data_clean <- raw_data |>
  mutate(
    across(all_of(numeric_vars), ~ parse_number(as.character(.x))),
    institution = str_squish(as.character(institution)),
    country = str_squish(as.character(country)),
    risk = str_to_lower(str_squish(as.character(risk))),
    risk = factor(
      risk,
      levels = c("very low", "low", "medium", "significant"),
      ordered = TRUE
    )
  )

stopifnot(
  nrow(data_clean) == 5475,
  n_distinct(data_clean$institution) == 5475,
  n_distinct(data_clean$country) == 151,
  abs(max(data_clean$overall, na.rm = TRUE) - 37.599) < 1e-9
)

write_csv(
  data_clean,
  file.path(data_dir, "Supplementary_Data_S1.csv"),
  na = ""
)

missing_summary <- data_clean |>
  summarise(across(everything(), ~ sum(is.na(.x)))) |>
  pivot_longer(everything(), names_to = "variable", values_to = "missing_n") |>
  arrange(desc(missing_n), variable)

write_csv(missing_summary, file.path(table_dir, "missing_values_summary.csv"))

# ------------------------------------------------------------------------------
# 3. Descriptive analyses and country summaries
# ------------------------------------------------------------------------------

global_summary <- data_clean |>
  summarise(
    n_institutions = n(),
    n_unique_institutions = n_distinct(institution),
    n_countries = n_distinct(country),
    overall_mean = mean(overall, na.rm = TRUE),
    overall_sd = sd(overall, na.rm = TRUE),
    overall_median = median(overall, na.rm = TRUE),
    overall_iqr = IQR(overall, na.rm = TRUE),
    overall_min = min(overall, na.rm = TRUE),
    overall_max = max(overall, na.rm = TRUE),
    output_mean = mean(output, na.rm = TRUE),
    output_sd = sd(output, na.rm = TRUE),
    output_median = median(output, na.rm = TRUE),
    output_iqr = IQR(output, na.rm = TRUE),
    output_min = min(output, na.rm = TRUE),
    output_max = max(output, na.rm = TRUE)
  )

s1_panel_a <- tribble(
  ~Characteristic, ~Value,
  "Institutions, n", global_summary$n_institutions,
  "Unique institutions, n", global_summary$n_unique_institutions,
  "Countries or territories, n", global_summary$n_countries,
  "Mean IRIS Overall score", global_summary$overall_mean,
  "SD IRIS Overall score", global_summary$overall_sd,
  "Median IRIS Overall score", global_summary$overall_median,
  "IQR IRIS Overall score", global_summary$overall_iqr,
  "Minimum IRIS Overall score", global_summary$overall_min,
  "Maximum IRIS Overall score", global_summary$overall_max,
  "Mean scientific output", global_summary$output_mean,
  "SD scientific output", global_summary$output_sd,
  "Median scientific output", global_summary$output_median,
  "IQR scientific output", global_summary$output_iqr,
  "Minimum scientific output", global_summary$output_min,
  "Maximum scientific output", global_summary$output_max
)

risk_distribution <- data_clean |>
  group_by(risk) |>
  summarise(
    institutions_n = n(),
    institutions_percent = n() / nrow(data_clean) * 100,
    mean_overall = mean(overall, na.rm = TRUE),
    median_overall = median(overall, na.rm = TRUE),
    median_output = median(output, na.rm = TRUE),
    .groups = "drop"
  )

s1_panel_b <- risk_distribution |>
  transmute(
    `IRIS structural risk category` = str_to_title(as.character(risk)),
    `Institutions, n` = institutions_n,
    `Institutions (%)` = institutions_percent,
    `Mean IRIS Overall score` = mean_overall,
    `Median IRIS Overall score` = median_overall,
    `Median scientific output` = median_output
  )

country_summary <- data_clean |>
  group_by(country) |>
  summarise(
    n_institutions = n(),
    mean_overall = mean(overall, na.rm = TRUE),
    median_overall = median(overall, na.rm = TRUE),
    sd_overall = sd(overall, na.rm = TRUE),
    pct_very_low = mean(risk == "very low", na.rm = TRUE) * 100,
    pct_low = mean(risk == "low", na.rm = TRUE) * 100,
    pct_medium = mean(risk == "medium", na.rm = TRUE) * 100,
    pct_significant = mean(risk == "significant", na.rm = TRUE) * 100,
    pct_medium_or_significant =
      mean(risk %in% c("medium", "significant"), na.rm = TRUE) * 100,
    mean_output = mean(output, na.rm = TRUE),
    median_output = median(output, na.rm = TRUE),
    across(
      all_of(risk_indicator_vars),
      ~ mean(.x, na.rm = TRUE),
      .names = "mean_{.col}"
    ),
    .groups = "drop"
  ) |>
  arrange(desc(n_institutions), country)

s2_country_summary <- country_summary |>
  transmute(
    `Country or territory code` = country,
    `Institutions, n` = n_institutions,
    `Mean IRIS Overall score` = mean_overall,
    `Median IRIS Overall score` = median_overall,
    `SD IRIS Overall score` = sd_overall,
    `Very low (%)` = pct_very_low,
    `Low (%)` = pct_low,
    `Medium (%)` = pct_medium,
    `Significant (%)` = pct_significant,
    `Medium or significant (%)` = pct_medium_or_significant,
    `Mean scientific output` = mean_output,
    `Median scientific output` = median_output,
    `Mean multi-affiliation` = mean_multi_affiliation,
    `Mean retracted output` = mean_retracted_output,
    `Mean self-citation` = mean_self_citation,
    `Mean discontinued journals output` = mean_discontinued_journals_output,
    `Mean hyperauthored output` = mean_hyperauthored_output,
    `Mean leadership impact gap` = mean_leadership_impact_gap,
    `Mean hyperprolific authors` = mean_hyperprolific_authors,
    `Mean institutional journal output` = mean_institutional_journal_output,
    `Mean redundant output` = mean_redundant_output
  )

top_countries <- country_summary |>
  arrange(desc(n_institutions), country) |>
  slice_head(n = 20)

top_significant <- country_summary |>
  filter(n_institutions >= 10) |>
  arrange(desc(pct_significant), desc(n_institutions), country) |>
  slice_head(n = 20)

write_csv(global_summary, file.path(table_dir, "global_summary.csv"))
write_csv(country_summary, file.path(table_dir, "country_summary.csv"))

# ------------------------------------------------------------------------------
# 4. PCA, clustering, and exploratory institutional profiles
# ------------------------------------------------------------------------------

block2_data <- data_clean |>
  select(
    id, institution, country, sir_rank, overall, risk, output,
    all_of(risk_indicator_vars)
  )

block2_complete <- block2_data |>
  filter(if_all(all_of(risk_indicator_vars), ~ !is.na(.x)))

stopifnot(nrow(block2_complete) == 5470)

inclusion_summary <- tibble(
  item = c(
    "Institutions in cleaned dataset",
    "Institutions with complete values for all nine indicators",
    "Institutions excluded because at least one indicator was missing"
  ),
  institutions_n = c(
    nrow(data_clean),
    nrow(block2_complete),
    nrow(data_clean) - nrow(block2_complete)
  )
)

winsor_thresholds <- map_dfr(
  risk_indicator_vars,
  function(variable) {
    tibble(
      indicator_variable = variable,
      indicator = indicator_labels[[variable]],
      p01 = quantile(block2_complete[[variable]], 0.01, na.rm = TRUE, type = 7),
      p99 = quantile(block2_complete[[variable]], 0.99, na.rm = TRUE, type = 7)
    )
  }
)

block2_winsor <- block2_complete |>
  mutate(
    across(
      all_of(risk_indicator_vars),
      winsorize_vec,
      .names = "{.col}_w"
    )
  )

winsor_vars <- paste0(risk_indicator_vars, "_w")
indicator_matrix <- block2_winsor |>
  select(all_of(winsor_vars)) |>
  as.matrix()
colnames(indicator_matrix) <- risk_indicator_vars
indicator_matrix_scaled <- scale(indicator_matrix)

scaling_parameters <- tibble(
  indicator_variable = risk_indicator_vars,
  indicator = unname(indicator_labels[risk_indicator_vars]),
  scaling_mean = as.numeric(attr(indicator_matrix_scaled, "scaled:center")),
  scaling_sd = as.numeric(attr(indicator_matrix_scaled, "scaled:scale"))
)

indicator_correlations <- cor(
  indicator_matrix,
  method = "spearman",
  use = "pairwise.complete.obs"
) |>
  as.data.frame() |>
  rownames_to_column("indicator_1") |>
  pivot_longer(-indicator_1, names_to = "indicator_2", values_to = "spearman_rho")

pca_fit <- prcomp(indicator_matrix_scaled, center = FALSE, scale. = FALSE)

pca_variance <- tibble(
  principal_component = paste0("PC", seq_along(pca_fit$sdev)),
  eigenvalue = pca_fit$sdev^2,
  variance_explained = eigenvalue / sum(eigenvalue),
  cumulative_variance = cumsum(variance_explained)
)

pca_loadings <- as_tibble(pca_fit$rotation, rownames = "indicator_variable") |>
  mutate(
    indicator = recode(indicator_variable, !!!indicator_labels),
    .after = indicator_variable
  )

pca_scores <- as_tibble(pca_fit$x) |>
  mutate(
    id = block2_winsor$id,
    institution = block2_winsor$institution,
    country = block2_winsor$country,
    overall = block2_winsor$overall,
    risk = block2_winsor$risk,
    output = block2_winsor$output,
    .before = 1
  )

set.seed(2026)
k_grid <- 2:8
silhouette_sample_size <- min(3000, nrow(indicator_matrix_scaled))
silhouette_sample <- sample(
  seq_len(nrow(indicator_matrix_scaled)),
  silhouette_sample_size
)
silhouette_dist <- dist(indicator_matrix_scaled[silhouette_sample, , drop = FALSE])

k_diagnostics <- map_dfr(
  k_grid,
  function(k) {
    set.seed(2026 + k)
    candidate <- kmeans(
      indicator_matrix_scaled,
      centers = k,
      nstart = 50,
      iter.max = 200
    )
    sil <- silhouette(candidate$cluster[silhouette_sample], silhouette_dist)
    tibble(
      k = k,
      total_withinss = candidate$tot.withinss,
      average_silhouette = mean(sil[, "sil_width"])
    )
  }
)

set.seed(2026)
final_kmeans <- kmeans(
  indicator_matrix_scaled,
  centers = 5,
  nstart = 100,
  iter.max = 500
)

clustered_raw <- block2_winsor |>
  mutate(cluster_raw = final_kmeans$cluster)

cluster_order <- clustered_raw |>
  group_by(cluster_raw) |>
  summarise(
    profile_mean_overall = mean(overall, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(profile_mean_overall) |>
  mutate(profile = paste0("Profile ", row_number()))

block2_clustered <- clustered_raw |>
  left_join(cluster_order, by = "cluster_raw") |>
  mutate(profile = factor(profile, levels = paste0("Profile ", 1:5)))

cluster_summary <- block2_clustered |>
  group_by(profile) |>
  summarise(
    n = n(),
    percent = n() / nrow(block2_clustered) * 100,
    countries = n_distinct(country),
    mean_overall = mean(overall, na.rm = TRUE),
    median_overall = median(overall, na.rm = TRUE),
    mean_output = mean(output, na.rm = TRUE),
    median_output = median(output, na.rm = TRUE),
    pct_very_low = mean(risk == "very low", na.rm = TRUE) * 100,
    pct_low = mean(risk == "low", na.rm = TRUE) * 100,
    pct_medium = mean(risk == "medium", na.rm = TRUE) * 100,
    pct_significant = mean(risk == "significant", na.rm = TRUE) * 100,
    .groups = "drop"
  )

scaled_with_profile <- as_tibble(indicator_matrix_scaled) |>
  mutate(id = block2_winsor$id, .before = 1) |>
  left_join(block2_clustered |> select(id, profile), by = "id")

profile_centroids <- scaled_with_profile |>
  group_by(profile) |>
  summarise(
    across(all_of(risk_indicator_vars), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

profile_centroids_long <- profile_centroids |>
  pivot_longer(
    all_of(risk_indicator_vars),
    names_to = "indicator_variable",
    values_to = "mean_z"
  ) |>
  mutate(
    indicator = recode(indicator_variable, !!!indicator_labels),
    indicator = factor(
      indicator,
      levels = unname(indicator_labels[risk_indicator_vars])
    )
  )

dominant_indicators <- profile_centroids_long |>
  group_by(profile) |>
  arrange(desc(mean_z), .by_group = TRUE) |>
  slice_head(n = 3) |>
  ungroup()

risk_by_profile <- block2_clustered |>
  count(profile, risk, name = "n") |>
  group_by(profile) |>
  mutate(within_profile_percent = n / sum(n) * 100) |>
  ungroup()

supplementary_data_s2 <- block2_clustered |>
  select(
    id, institution, country, sir_rank, overall, risk, output,
    all_of(risk_indicator_vars), all_of(winsor_vars),
    cluster_raw, profile, profile_mean_overall
  )

write_csv(
  supplementary_data_s2,
  file.path(data_dir, "Supplementary_Data_S2.csv"),
  na = ""
)

write_csv(pca_scores, file.path(table_dir, "pca_scores.csv"))
write_csv(indicator_correlations, file.path(table_dir, "indicator_correlations.csv"))
write_csv(scaling_parameters, file.path(table_dir, "scaling_parameters.csv"))

# ------------------------------------------------------------------------------
# 5. World Bank linkage and analytical coverage
# ------------------------------------------------------------------------------

start_year <- 2015
end_year <- 2024
wdi_indicators <- c(
  rd_gdp = "GB.XPD.RSDV.GD.ZS",
  researchers_pm = "SP.POP.SCIE.RD.P6"
)

if (!is.na(wdi_snapshot_file) && file.exists(wdi_snapshot_file)) {
  wdi_raw <- read_csv(wdi_snapshot_file, show_col_types = FALSE) |>
    clean_names()
} else {
  wdi_raw <- WDI::WDI(
    country = "all",
    indicator = wdi_indicators,
    start = start_year,
    end = end_year,
    extra = TRUE,
    cache = NULL
  ) |>
    as_tibble() |>
    clean_names()
  wdi_snapshot_file <- file.path(table_dir, "block3_wdi_raw_download.csv")
  write_csv(wdi_raw, wdi_snapshot_file)
}

wdi_countries <- wdi_raw |>
  filter(region != "Aggregates", !is.na(iso3c))

country_metadata <- wdi_countries |>
  arrange(iso3c, desc(year)) |>
  group_by(iso3c) |>
  summarise(
    country_name_wb = first(country),
    region_wb = first(region),
    income_group_wb = first(income),
    lending_group_wb = first(lending),
    .groups = "drop"
  )

wdi_latest_long <- wdi_countries |>
  select(iso3c, year, rd_gdp, researchers_pm) |>
  pivot_longer(
    c(rd_gdp, researchers_pm),
    names_to = "indicator",
    values_to = "value"
  ) |>
  filter(!is.na(value)) |>
  arrange(iso3c, indicator, desc(year)) |>
  group_by(iso3c, indicator) |>
  slice_head(n = 1) |>
  ungroup()

wdi_latest_values <- wdi_latest_long |>
  select(iso3c, indicator, value) |>
  pivot_wider(names_from = indicator, values_from = value)

wdi_latest_years <- wdi_latest_long |>
  select(iso3c, indicator, year) |>
  mutate(indicator = paste0(indicator, "_year")) |>
  pivot_wider(names_from = indicator, values_from = year)

country_context <- country_metadata |>
  left_join(wdi_latest_values, by = "iso3c") |>
  left_join(wdi_latest_years, by = "iso3c")

country_context_iris <- country_summary |>
  left_join(country_context, by = c("country" = "iso3c"))

match_summary <- tibble(
  item = c(
    "IRIS countries",
    "IRIS countries matched with World Bank metadata",
    "IRIS countries with R&D expenditure data",
    "IRIS countries with researchers per million data",
    "IRIS countries with both R&D and researchers data"
  ),
  countries_n = c(
    n_distinct(country_summary$country),
    sum(!is.na(country_context_iris$country_name_wb)),
    sum(!is.na(country_context_iris$rd_gdp)),
    sum(!is.na(country_context_iris$researchers_pm)),
    sum(
      !is.na(country_context_iris$rd_gdp) &
        !is.na(country_context_iris$researchers_pm)
    )
  )
)

country_analysis <- country_context_iris |>
  filter(n_institutions >= 10) |>
  mutate(
    income_group_wb = if_else(
      income_group_wb %in% valid_income_groups,
      income_group_wb,
      NA_character_
    ),
    income_group_wb = factor(income_group_wb, levels = valid_income_groups),
    researchers_pm_log = log1p(researchers_pm)
  )

country_analysis_coverage <- tibble(
  item = c(
    "Countries in the ≥10-institution analytical subset",
    "Institutions represented in that subset",
    "Countries with World Bank income group",
    "Countries with R&D expenditure",
    "Countries with researchers per million",
    "Countries with both R&D expenditure and researcher density"
  ),
  n = c(
    nrow(country_analysis),
    sum(country_analysis$n_institutions),
    sum(!is.na(country_analysis$income_group_wb)),
    sum(!is.na(country_analysis$rd_gdp)),
    sum(!is.na(country_analysis$researchers_pm)),
    sum(!is.na(country_analysis$rd_gdp) & !is.na(country_analysis$researchers_pm))
  )
)

unmatched_iris_wb <- country_context_iris |>
  filter(is.na(country_name_wb)) |>
  transmute(
    `IRIS code` = country,
    `Institutions, n` = n_institutions,
    `Mean IRIS Overall score` = mean_overall,
    `Significant (%)` = pct_significant
  ) |>
  arrange(desc(`Institutions, n`))

# Map matching is calculated before Supplementary Table S5 is assembled.
world_map <- rnaturalearth::ne_countries(
  scale = "medium",
  returnclass = "sf"
) |>
  filter(name != "Antarctica") |>
  mutate(
    join_iso3 = case_when(
      name == "Kosovo" ~ "XKX",
      !is.na(iso_a3) & iso_a3 != "-99" ~ iso_a3,
      !is.na(adm0_a3) & adm0_a3 != "-99" ~ adm0_a3,
      TRUE ~ NA_character_
    )
  ) |>
  select(name, iso_a3, adm0_a3, join_iso3, geometry)

map_data <- world_map |>
  left_join(country_summary, by = c("join_iso3" = "country"))

unmatched_map_codes <- country_summary |>
  anti_join(world_map, by = c("country" = "join_iso3")) |>
  arrange(desc(n_institutions))

map_match_summary <- tibble(
  item = c(
    "Countries in SCImago IRIS country_summary",
    "Countries matched with map",
    "Countries unmatched after correction"
  ),
  countries_n = c(
    n_distinct(country_summary$country),
    n_distinct(map_data$join_iso3[!is.na(map_data$mean_overall)]),
    nrow(unmatched_map_codes)
  )
)

# ------------------------------------------------------------------------------
# 6. Country correlations, weighted regression, and influence analyses
# ------------------------------------------------------------------------------

correlation_specs <- tribble(
  ~predictor, ~outcome, ~x_variable, ~y_variable, ~bootstrap_seed,
  "R&D expenditure", "Mean IRIS Overall score",
  "rd_gdp", "mean_overall", 2026,
  "Researcher density", "Mean IRIS Overall score",
  "researchers_pm", "mean_overall", 2027,
  "R&D expenditure", "Institutions in the significant category (%)",
  "rd_gdp", "pct_significant", 2028,
  "Researcher density", "Institutions in the significant category (%)",
  "researchers_pm", "pct_significant", 2029
)

country_correlations <- pmap_dfr(
  correlation_specs,
  function(predictor, outcome, x_variable, y_variable, bootstrap_seed) {
    pair <- country_analysis |>
      select(all_of(c(x_variable, y_variable))) |>
      drop_na()

    test <- suppressWarnings(
      cor.test(
        pair[[x_variable]],
        pair[[y_variable]],
        method = "spearman",
        exact = FALSE,
        alternative = "two.sided"
      )
    )
    ci <- bootstrap_spearman(
      pair[[x_variable]],
      pair[[y_variable]],
      seed = bootstrap_seed,
      replicates = 10000
    )

    tibble(
      predictor = predictor,
      outcome = outcome,
      n = nrow(pair),
      spearman_rho = unname(test$estimate),
      ci_low = ci$ci_low,
      ci_high = ci$ci_high,
      p_unadjusted = test$p.value
    )
  }
) |>
  mutate(p_holm = p.adjust(p_unadjusted, method = "holm"))

country_model_data <- country_analysis |>
  filter(
    !is.na(mean_overall),
    !is.na(rd_gdp),
    !is.na(researchers_pm),
    !is.na(income_group_wb)
  ) |>
  mutate(
    rd_gdp_z = z_score(rd_gdp),
    researchers_pm_log_z = z_score(log1p(researchers_pm)),
    n_institutions_log_z = z_score(log1p(n_institutions)),
    income_group_wb = droplevels(income_group_wb),
    income_group_wb = relevel(income_group_wb, ref = "Lower middle income")
  )

stopifnot(nrow(country_model_data) == 59)

country_model_overall <- lm(
  mean_overall ~
    rd_gdp_z +
    researchers_pm_log_z +
    income_group_wb +
    n_institutions_log_z,
  data = country_model_data,
  weights = n_institutions
)

regression_term_labels <- c(
  "(Intercept)" = "(Intercept)",
  "rd_gdp_z" = "R&D expenditure, standardized",
  "researchers_pm_log_z" = "Log researcher density, standardized",
  "income_group_wbUpper middle income" = "Upper-middle-income group",
  "income_group_wbHigh income" = "High-income group",
  "n_institutions_log_z" = "Log number of institutions, standardized"
)

weighted_regression <- broom::tidy(
  country_model_overall,
  conf.int = TRUE,
  conf.level = 0.95
) |>
  transmute(
    term = recode(term, !!!regression_term_labels),
    estimate,
    model_based_se = std.error,
    t_statistic = statistic,
    p_value = p.value,
    ci_low = conf.low,
    ci_high = conf.high
  )

regression_glance <- broom::glance(country_model_overall)
income_counts <- country_model_data |>
  count(income_group_wb, name = "countries_n")

regression_summary <- tibble(
  model_statistic = c(
    "Countries, n",
    "Residual degrees of freedom",
    "R²",
    "Adjusted R²",
    "F statistic",
    "Model P value",
    "Lower-middle-income countries, n",
    "Upper-middle-income countries, n",
    "High-income countries, n"
  ),
  value = c(
    nrow(country_model_data),
    df.residual(country_model_overall),
    regression_glance$r.squared,
    regression_glance$adj.r.squared,
    regression_glance$statistic,
    regression_glance$p.value,
    income_counts$countries_n[
      income_counts$income_group_wb == "Lower middle income"
    ],
    income_counts$countries_n[
      income_counts$income_group_wb == "Upper middle income"
    ],
    income_counts$countries_n[
      income_counts$income_group_wb == "High income"
    ]
  )
)

influence_diagnostics <- country_model_data |>
  transmute(
    country_code = country,
    institutions_n = n_institutions,
    cooks_distance = cooks.distance(country_model_overall),
    leverage = hatvalues(country_model_overall),
    standardized_weighted_residual = rstandard(country_model_overall)
  )

loo_estimates <- map_dfr(
  seq_len(nrow(country_model_data)),
  function(index) {
    refit <- update(
      country_model_overall,
      data = country_model_data[-index, , drop = FALSE]
    )
    tibble(
      excluded_country = country_model_data$country[[index]],
      rd_coefficient = unname(coef(refit)[["rd_gdp_z"]]),
      researcher_density_coefficient =
        unname(coef(refit)[["researchers_pm_log_z"]])
    )
  }
)

primary_coefficients <- coef(country_model_overall)
loo_summary <- tibble(
  term = c(
    "R&D expenditure, standardized",
    "Log researcher density, standardized"
  ),
  primary_estimate = c(
    primary_coefficients[["rd_gdp_z"]],
    primary_coefficients[["researchers_pm_log_z"]]
  ),
  minimum_loo_estimate = c(
    min(loo_estimates$rd_coefficient),
    min(loo_estimates$researcher_density_coefficient)
  ),
  maximum_loo_estimate = c(
    max(loo_estimates$rd_coefficient),
    max(loo_estimates$researcher_density_coefficient)
  ),
  country_largest_change = c(
    loo_estimates$excluded_country[which.max(
      abs(loo_estimates$rd_coefficient - primary_coefficients[["rd_gdp_z"]])
    )],
    loo_estimates$excluded_country[which.max(
      abs(
        loo_estimates$researcher_density_coefficient -
          primary_coefficients[["researchers_pm_log_z"]]
      )
    )]
  )
)

country_analysis_coverage <- bind_rows(
  country_analysis_coverage,
  tibble(
    item = "Countries in the complete-case weighted regression",
    n = nrow(country_model_data)
  )
)

write_csv(country_correlations, file.path(table_dir, "country_correlations.csv"))
write_csv(weighted_regression, file.path(table_dir, "weighted_regression.csv"))
write_csv(loo_estimates, file.path(table_dir, "leave_one_country_out.csv"))
write_csv(
  influence_diagnostics,
  file.path(table_dir, "weighted_regression_influence.csv")
)

# ------------------------------------------------------------------------------
# 7. Institution-level multilevel model and sensitivity analyses
# ------------------------------------------------------------------------------

institution_context <- data_clean |>
  left_join(country_context, by = c("country" = "iso3c")) |>
  mutate(
    income_group_wb = if_else(
      income_group_wb %in% valid_income_groups,
      income_group_wb,
      NA_character_
    ),
    income_group_wb = factor(income_group_wb, levels = valid_income_groups),
    output_log_z = z_score(log1p(output)),
    sir_rank_log_z = z_score(log1p(sir_rank)),
    rd_gdp_z = z_score(rd_gdp),
    researchers_pm_log_z = z_score(log1p(researchers_pm))
  )

institution_model_data <- institution_context |>
  filter(
    !is.na(overall),
    !is.na(country),
    !is.na(output_log_z),
    !is.na(sir_rank_log_z),
    !is.na(rd_gdp_z),
    !is.na(researchers_pm_log_z),
    !is.na(income_group_wb)
  ) |>
  mutate(income_group_wb = droplevels(income_group_wb))

stopifnot(
  nrow(institution_model_data) == 4891,
  n_distinct(institution_model_data$country) == 110
)

mixed_formula <- overall ~
  output_log_z +
  sir_rank_log_z +
  rd_gdp_z +
  researchers_pm_log_z +
  income_group_wb +
  (1 | country)

lmer_control <- lmerControl(
  optimizer = "bobyqa",
  optCtrl = list(maxfun = 200000),
  check.conv.singular = "warning"
)

primary_mixed_model <- lmer(
  mixed_formula,
  data = institution_model_data,
  REML = FALSE,
  control = lmer_control
)

max_overall_index <- which.max(institution_model_data$overall)
verified_maximum <- institution_model_data[max_overall_index, ]

multilevel_influence_diagnostics <- institution_model_data |>
  transmute(
    id,
    institution,
    country,
    overall,
    conditional_residual = resid(primary_mixed_model),
    standardized_conditional_residual =
      conditional_residual / sigma(primary_mixed_model)
  ) |>
  arrange(desc(abs(standardized_conditional_residual)))

max_residual_id <- multilevel_influence_diagnostics$id[[1]]
stopifnot(max_residual_id == verified_maximum$id)

write_csv(
  multilevel_influence_diagnostics,
  file.path(table_dir, "multilevel_influence_diagnostics.csv")
)

without_maximum_data <- institution_model_data[-max_overall_index, , drop = FALSE]
without_maximum_model <- lmer(
  mixed_formula,
  data = without_maximum_data,
  REML = FALSE,
  control = lmer_control
)

eligible_min10 <- country_summary |>
  filter(n_institutions >= 10) |>
  pull(country)

min10_model_data <- institution_model_data |>
  filter(country %in% eligible_min10) |>
  mutate(
    income_group_wb = droplevels(income_group_wb),
    income_group_wb = relevel(income_group_wb, ref = "Lower middle income")
  )

min10_mixed_model <- lmer(
  mixed_formula,
  data = min10_model_data,
  REML = FALSE,
  control = lmer_control
)

model_fit_summary <- bind_rows(
  model_fit_row(
    primary_mixed_model,
    "Primary complete-case model",
    institution_model_data
  ),
  model_fit_row(
    without_maximum_model,
    "Excluding verified maximum Overall score",
    without_maximum_data
  ),
  model_fit_row(
    min10_mixed_model,
    "Countries represented by at least 10 institutions",
    min10_model_data
  )
)

primary_fixed_effects <- fixed_effect_table(
  primary_mixed_model,
  reference_label = "Low income"
) |>
  rename_fixed_terms()

sensitivity_fixed_effects <- bind_rows(
  fixed_effect_table(
    without_maximum_model,
    model_label = "Excluding verified maximum Overall score",
    reference_label = "Low income"
  ),
  fixed_effect_table(
    min10_mixed_model,
    model_label = "Countries represented by at least 10 institutions",
    reference_label = "Lower middle income"
  )
) |>
  rename_fixed_terms()

random_effects <- bind_rows(
  random_effect_table(
    primary_mixed_model,
    "Primary complete-case model"
  ),
  random_effect_table(
    without_maximum_model,
    "Excluding verified maximum Overall score"
  ),
  random_effect_table(
    min10_mixed_model,
    "Countries represented by at least 10 institutions"
  )
)

influential_observation <- institution_model_data |>
  arrange(desc(overall)) |>
  slice_head(n = 2) |>
  transmute(
    status_in_primary_model = c(
      "Verified maximum; retained in the primary model",
      "Second-highest IRIS Overall score"
    ),
    institution_id = id,
    institution,
    country_code = country,
    iris_overall_score = overall,
    retracted_output
  )

write_csv(
  primary_fixed_effects,
  file.path(table_dir, "multilevel_primary_fixed_effects.csv")
)
write_csv(
  sensitivity_fixed_effects,
  file.path(table_dir, "multilevel_sensitivity_fixed_effects.csv")
)
write_csv(
  model_fit_summary,
  file.path(table_dir, "multilevel_model_fit_summary.csv")
)
write_csv(
  random_effects,
  file.path(table_dir, "multilevel_random_effects.csv")
)

# ------------------------------------------------------------------------------
# 8. Profile distribution by World Bank income group
# ------------------------------------------------------------------------------

profiles_context <- block2_clustered |>
  left_join(country_context, by = c("country" = "iso3c")) |>
  mutate(
    income_group_wb = if_else(
      income_group_wb %in% valid_income_groups,
      income_group_wb,
      NA_character_
    ),
    income_group_wb = factor(income_group_wb, levels = valid_income_groups)
  )

profile_income_distribution <- profiles_context |>
  filter(!is.na(income_group_wb), !is.na(profile)) |>
  count(income_group_wb, profile, name = "n") |>
  group_by(income_group_wb) |>
  mutate(percent = n / sum(n) * 100) |>
  ungroup()

stopifnot(sum(profile_income_distribution$n) == 5360)
write_csv(
  profile_income_distribution,
  file.path(table_dir, "profiles_by_income_group.csv")
)

# ------------------------------------------------------------------------------
# 9. Supplementary Tables S1-S8
# ------------------------------------------------------------------------------

s3_variance <- pca_variance |>
  transmute(
    `Principal component` = principal_component,
    `Eigenvalue` = eigenvalue,
    `Variance explained` = variance_explained,
    `Cumulative variance` = cumulative_variance
  )

s3_loadings <- pca_loadings |>
  transmute(
    `Indicator variable` = indicator_variable,
    `Indicator` = indicator,
    across(starts_with("PC"))
  )

s4_inclusion <- inclusion_summary |>
  transmute(`Item` = item, `Institutions, n` = institutions_n)

s4_parameters <- tribble(
  ~Parameter, ~Value, ~Definition,
  "Analytical sample", "5470",
  "Institutions with complete values for all nine indicators",
  "Winsorization", "1st and 99th percentiles",
  "Applied separately to each indicator",
  "Standardization", "z score", "Applied after winsorization",
  "Random seed", "2026", "Used for k-means and silhouette sampling",
  "Candidate k values", "2–8",
  "Evaluated using total within-cluster sum of squares and silhouette width",
  "Silhouette sample", "3000",
  "Institutions sampled for diagnostic silhouette calculations",
  "Diagnostic nstart", "50", "Random initializations for k = 2–8",
  "Diagnostic maximum iterations", "200",
  "Maximum iterations for candidate k = 2–8 models",
  "Final k", "5", "Retained for substantive interpretability",
  "Final nstart", "100",
  "Random initializations for the final k = 5 solution",
  "Final maximum iterations", "500",
  "Maximum iterations for the final k = 5 model"
)

s4_k_diagnostics <- k_diagnostics |>
  transmute(
    `Number of profiles (k)` = k,
    `Total within-cluster sum of squares` = total_withinss,
    `Average silhouette width` = average_silhouette
  )

s4_profile_summary <- cluster_summary |>
  transmute(
    `Exploratory profile` = as.character(profile),
    `Institutions, n` = n,
    `Institutions (%)` = percent,
    `Countries or territories, n` = countries,
    `Mean IRIS Overall score` = mean_overall,
    `Median IRIS Overall score` = median_overall,
    `Mean scientific output` = mean_output,
    `Median scientific output` = median_output,
    `Very low (%)` = pct_very_low,
    `Low (%)` = pct_low,
    `Medium (%)` = pct_medium,
    `Significant (%)` = pct_significant
  )

s4_centroids <- profile_centroids |>
  rename_with(
    ~ unname(indicator_labels[.x]),
    all_of(risk_indicator_vars)
  ) |>
  rename(`Exploratory profile` = profile) |>
  mutate(`Exploratory profile` = as.character(`Exploratory profile`))

s4_dominant <- dominant_indicators |>
  transmute(
    `Exploratory profile` = as.character(profile),
    `Indicator variable` = indicator_variable,
    `Mean standardized value` = mean_z,
    `Indicator` = as.character(indicator)
  )

s4_risk_by_profile <- risk_by_profile |>
  transmute(
    `Exploratory profile` = as.character(profile),
    `IRIS structural risk category` = str_to_title(as.character(risk)),
    `Institutions, n` = n,
    `Within-profile percentage` = within_profile_percent
  )

s4_winsor <- winsor_thresholds |>
  transmute(
    `Indicator variable` = indicator_variable,
    `Indicator` = indicator,
    `1st percentile` = p01,
    `99th percentile` = p99
  )

regression_country_codes <- country_model_data$country
s5_context_by_country <- country_context_iris |>
  transmute(
    `IRIS code` = country,
    `World Bank country or economy` = country_name_wb,
    `Institutions, n` = n_institutions,
    `Matched to World Bank` = if_else(!is.na(country_name_wb), "Yes", "No"),
    `World Bank region` = region_wb,
    `World Bank income group` = income_group_wb,
    `R&D expenditure (% GDP)` = rd_gdp,
    `R&D year` = rd_gdp_year,
    `Researchers per million` = researchers_pm,
    `Researchers year` = researchers_pm_year,
    `≥10-institution subset` = if_else(n_institutions >= 10, "Yes", "No"),
    `Both indicators available` = if_else(
      !is.na(rd_gdp) & !is.na(researchers_pm), "Yes", "No"
    ),
    `Weighted regression sample` = if_else(
      country %in% regression_country_codes, "Yes", "No"
    )
  )

s5_context_years <- full_join(
  country_context_iris |>
    filter(!is.na(rd_gdp)) |>
    count(rd_gdp_year, name = "rd_n") |>
    rename(year = rd_gdp_year),
  country_context_iris |>
    filter(!is.na(researchers_pm)) |>
    count(researchers_pm_year, name = "researcher_n") |>
    rename(year = researchers_pm_year),
  by = "year"
) |>
  arrange(year) |>
  transmute(
    `Indicator year` = year,
    `IRIS-matched countries with R&D value, n` = replace_na(rd_n, 0L),
    `IRIS-matched countries with researcher value, n` =
      replace_na(researcher_n, 0L)
  )

s5_unmatched_map <- unmatched_map_codes |>
  transmute(
    `Country or territory code` = country,
    `Institutions, n` = n_institutions,
    `Mean IRIS Overall score` = mean_overall,
    `Median IRIS Overall score` = median_overall,
    `SD IRIS Overall score` = sd_overall,
    `Very low (%)` = pct_very_low,
    `Low (%)` = pct_low,
    `Medium (%)` = pct_medium,
    `Significant (%)` = pct_significant,
    `Medium or significant (%)` = pct_medium_or_significant,
    `Mean scientific output` = mean_output,
    `Median scientific output` = median_output,
    `Mean multi-affiliation` = mean_multi_affiliation,
    `Mean retracted output` = mean_retracted_output,
    `Mean self-citation` = mean_self_citation,
    `Mean discontinued journals output` = mean_discontinued_journals_output,
    `Mean hyperauthored output` = mean_hyperauthored_output,
    `Mean leadership impact gap` = mean_leadership_impact_gap,
    `Mean hyperprolific authors` = mean_hyperprolific_authors,
    `Mean institutional journal output` = mean_institutional_journal_output,
    `Mean redundant output` = mean_redundant_output
  )

s6_correlations <- country_correlations |>
  transmute(
    `Predictor` = predictor,
    `Outcome` = outcome,
    `n` = n,
    `Spearman rho` = spearman_rho,
    `95% CI low` = ci_low,
    `95% CI high` = ci_high,
    `Unadjusted P value` = p_unadjusted,
    `Holm-adjusted P value` = p_holm
  )

s6_regression <- weighted_regression |>
  transmute(
    `Term` = term,
    `Estimate` = estimate,
    `Model-based SE` = model_based_se,
    `t statistic` = t_statistic,
    `P value` = p_value,
    `95% CI low` = ci_low,
    `95% CI high` = ci_high
  )

s6_regression_summary <- regression_summary |>
  transmute(`Model statistic` = model_statistic, `Value` = value)

s6_loo_summary <- loo_summary |>
  transmute(
    `Term` = term,
    `Primary estimate` = primary_estimate,
    `Minimum leave-one-country-out estimate` = minimum_loo_estimate,
    `Maximum leave-one-country-out estimate` = maximum_loo_estimate,
    `Country producing largest absolute change` = country_largest_change
  )

s6_loo_estimates <- loo_estimates |>
  transmute(
    `Excluded country` = excluded_country,
    `R&D coefficient` = rd_coefficient,
    `Researcher-density coefficient` = researcher_density_coefficient
  )

s6_influence <- influence_diagnostics |>
  transmute(
    `Country code` = country_code,
    `Institutions, n` = institutions_n,
    `Cook's distance` = cooks_distance,
    `Leverage` = leverage,
    `Standardized weighted residual` = standardized_weighted_residual
  )

s7_fit <- model_fit_summary |>
  transmute(
    `Model` = model,
    `Institutions, n` = institutions_n,
    `Countries, n` = countries_n,
    `Marginal R²` = marginal_r2,
    `Conditional R²` = conditional_r2,
    `Adjusted ICC` = adjusted_icc
  )

s7_primary_fixed <- primary_fixed_effects |>
  transmute(
    `Income-group reference` = income_group_reference,
    `Term` = term,
    `Estimate` = estimate,
    `SE` = se,
    `z statistic` = z_statistic,
    `P value` = p_value,
    `95% CI low` = ci_low,
    `95% CI high` = ci_high
  )

s7_sensitivity_fixed <- sensitivity_fixed_effects |>
  transmute(
    `Model` = model,
    `Income-group reference` = income_group_reference,
    `Term` = term,
    `Estimate` = estimate,
    `SE` = se,
    `z statistic` = z_statistic,
    `P value` = p_value,
    `95% CI low` = ci_low,
    `95% CI high` = ci_high
  )

s7_random <- random_effects |>
  transmute(
    `Model` = model,
    `Random component` = random_component,
    `Variance` = variance,
    `SD` = sd
  )

s7_influential <- influential_observation |>
  transmute(
    `Status in primary model` = status_in_primary_model,
    `Institution ID` = institution_id,
    `Institution` = institution,
    `Country code` = country_code,
    `IRIS Overall score` = iris_overall_score,
    `Retracted output` = retracted_output
  )

# The variable dictionary is explicit so the analytical meaning of each raw and
# derived field remains versioned with the code.
s8_dictionary <- tribble(
  ~Variable, ~Definition, ~Source, ~`Analytical level`, ~`Used in`, ~`Transformation or coding`,
  "id", "Institution identifier extracted from SCImago IRIS", "SCImago IRIS",
  "Institution", "Data cleaning; institutional analyses; profile assignment", "Raw",
  "institution", "Name of the higher education institution", "SCImago IRIS",
  "Institution", "Descriptive analyses; institutional profile assignment",
  "Whitespace cleaned; raw text",
  "country", "Country or territory code used by SCImago IRIS, harmonized where needed for World Bank and map matching",
  "SCImago IRIS", "Institution/Country",
  "All analyses; country aggregation; World Bank linkage",
  "Raw; used as ISO3-style join key",
  "sir_rank", "SCImago Institutions Ranking position", "SCImago IRIS",
  "Institution", "Descriptive analyses; multilevel model",
  "Raw; log1p-transformed and standardized in multilevel model",
  "overall", "Institution-level IRIS Overall score", "SCImago IRIS",
  "Institution",
  "Descriptive analyses; country summaries; weighted regression; multilevel model; profile characterization",
  "Raw",
  "risk", "IRIS structural risk category", "SCImago IRIS", "Institution",
  "Descriptive analyses; profile characterization",
  "Ordered categorical: very low, low, medium, significant",
  "output", "Institutional scientific output", "SCImago IRIS", "Institution",
  "Descriptive analyses; multilevel model",
  "Raw; log1p-transformed and standardized in multilevel model",
  "multi_affiliation", "IRIS indicator capturing multi-affiliation patterns",
  "SCImago IRIS", "Institution", "PCA; k-means; profile characterization",
  "Winsorized at p1/p99 and standardized for PCA/clustering",
  "retracted_output", "IRIS indicator capturing output associated with retracted publications",
  "SCImago IRIS", "Institution", "PCA; k-means; profile characterization",
  "Winsorized at p1/p99 and standardized for PCA/clustering",
  "self_citation", "IRIS indicator capturing institutional self-citation patterns",
  "SCImago IRIS", "Institution", "PCA; k-means; profile characterization",
  "Winsorized at p1/p99 and standardized for PCA/clustering",
  "discontinued_journals_output", "IRIS indicator capturing output in discontinued journals",
  "SCImago IRIS", "Institution", "PCA; k-means; profile characterization",
  "Winsorized at p1/p99 and standardized for PCA/clustering",
  "hyperauthored_output", "IRIS indicator capturing output in hyperauthored publications",
  "SCImago IRIS", "Institution", "PCA; k-means; profile characterization",
  "Winsorized at p1/p99 and standardized for PCA/clustering",
  "leadership_impact_gap", "IRIS indicator capturing gaps between leadership and impact patterns",
  "SCImago IRIS", "Institution", "PCA; k-means; profile characterization",
  "Winsorized at p1/p99 and standardized for PCA/clustering",
  "hyperprolific_authors", "IRIS indicator capturing hyperprolific author patterns",
  "SCImago IRIS", "Institution", "PCA; k-means; profile characterization",
  "Winsorized at p1/p99 and standardized for PCA/clustering",
  "institutional_journal_output", "IRIS indicator capturing output in institutional journals",
  "SCImago IRIS", "Institution", "PCA; k-means; profile characterization",
  "Winsorized at p1/p99 and standardized for PCA/clustering",
  "redundant_output", "IRIS indicator capturing redundant or potentially overlapping output patterns",
  "SCImago IRIS", "Institution", "PCA; k-means; profile characterization",
  "Winsorized at p1/p99 and standardized for PCA/clustering",
  "profile", "Exploratory institutional profile assigned by k-means clustering",
  "Derived from analysis", "Institution",
  "Profile characterization; distribution by income group",
  "Profile 1 to Profile 5; ordered by profile mean IRIS Overall score",
  "cluster_raw", "Original numerical k-means cluster label before profile ordering",
  "Derived from analysis", "Institution", "Internal clustering output", "Raw k-means cluster ID",
  "profile_mean_overall", "Mean IRIS Overall score of the exploratory profile to which the institution was assigned",
  "Derived from analysis", "Profile/Institution",
  "Ordering and labeling of profiles; Supplementary Data S2",
  "Profile-level mean repeated for assigned institutions",
  "n_institutions", "Number of institutions represented within each country",
  "Derived from SCImago IRIS", "Country",
  "Country summaries; map inclusion; weighted regression",
  "Raw count; log1p-transformed and standardized in weighted regression",
  "mean_overall", "Mean institution-level IRIS Overall score within each country",
  "Derived from SCImago IRIS", "Country",
  "Country analysis; maps; correlations; weighted regression", "Raw country mean",
  "median_overall", "Median institution-level IRIS Overall score within each country",
  "Derived from SCImago IRIS", "Country", "Country descriptive summaries", "Raw country median",
  "pct_significant", "Percentage of institutions in the IRIS significant structural risk category within each country",
  "Derived from SCImago IRIS", "Country",
  "Country analysis; maps; correlations", "Raw percentage",
  "income_group_wb", "World Bank income group classification", "World Bank",
  "Country", "Contextual analyses; profile distribution; regression models",
  "Low, lower middle, upper middle, or high income",
  "region_wb", "World Bank region classification", "World Bank", "Country",
  "Coverage and contextual description", "Categorical",
  "rd_gdp", "R&D expenditure as percentage of GDP",
  "World Bank World Development Indicators", "Country",
  "Correlations; weighted regression; multilevel model",
  "Most recent non-missing value from 2015 onward; standardized in models",
  "rd_gdp_year", "Year corresponding to the selected R&D expenditure value",
  "World Bank World Development Indicators", "Country", "Coverage and reproducibility", "Raw year",
  "researchers_pm", "Researchers in R&D per million people",
  "World Bank World Development Indicators", "Country",
  "Correlations; weighted regression; multilevel model",
  "Most recent non-missing value from 2015 onward; log1p-transformed and standardized",
  "researchers_pm_year", "Year corresponding to the selected researcher-density value",
  "World Bank World Development Indicators", "Country", "Coverage and reproducibility", "Raw year"
)

winsor_dictionary <- tibble(
  Variable = winsor_vars,
  Definition = paste0(
    "Winsorized ",
    str_to_lower(unname(indicator_labels[risk_indicator_vars])),
    " value"
  ),
  Source = "Derived from SCImago IRIS",
  `Analytical level` = "Institution",
  `Used in` = "PCA; k-means clustering; profile characterization",
  `Transformation or coding` =
    "Capped at the indicator-specific 1st and 99th percentiles"
)

model_dictionary <- tribble(
  ~Variable, ~Definition, ~Source, ~`Analytical level`, ~`Used in`, ~`Transformation or coding`,
  "rd_gdp_z", "Standardized R&D expenditure as percentage of GDP",
  "Derived from World Bank WDI", "Country",
  "Weighted regression; multilevel model", "z score",
  "researchers_pm_log_z", "Standardized log-transformed researcher density",
  "Derived from World Bank WDI", "Country",
  "Weighted regression; multilevel model", "log(1 + x), then z score",
  "n_institutions_log_z", "Standardized log-transformed number of represented institutions",
  "Derived from SCImago IRIS", "Country", "Weighted regression",
  "log(1 + x), then z score",
  "output_log_z", "Standardized log-transformed institutional scientific output",
  "Derived from SCImago IRIS", "Institution", "Multilevel model",
  "log(1 + x), then z score",
  "sir_rank_log_z", "Standardized log-transformed SIR Rank",
  "Derived from SCImago IRIS", "Institution", "Multilevel model",
  "log(1 + x), then z score"
)

s8_dictionary <- bind_rows(s8_dictionary, winsor_dictionary, model_dictionary)

write_sectioned_csv(
  1,
  "Supplementary Table S1. Overall characteristics of the institutional dataset and distribution across IRIS structural risk categories.",
  "Panel A summarizes the full institution-level dataset. Panel B presents counts, percentages, IRIS Overall scores, and scientific output across the four IRIS structural risk categories.",
  c(
    "IRIS Overall score and structural risk category are reported as supplied by SCImago IRIS.",
    "Percentages use all 5,475 institutions as the denominator.",
    "Scientific output is the institutional output variable supplied by SCImago IRIS."
  ),
  list(
    `Panel A - Global` = s1_panel_a,
    `Panel B - IRIS categories` = s1_panel_b
  ),
  file.path(supp_dir, "Supplementary_Table_S1.csv")
)

write_sectioned_csv(
  2,
  "Supplementary Table S2. Country-level summary of the SCImago IRIS dataset.",
  "Country- and territory-level counts, IRIS Overall statistics, structural risk categories, scientific output, and means of the nine individual IRIS indicators.",
  c(
    "The table includes all 151 country or territory codes in the institution-level dataset.",
    "Country summaries describe the institutions represented in IRIS and are not national rankings."
  ),
  list(`Country summary` = s2_country_summary),
  file.path(supp_dir, "Supplementary_Table_S2.csv")
)

write_sectioned_csv(
  3,
  "Supplementary Table S3. Principal component analysis of the nine SCImago IRIS indicators.",
  "Variance explained, cumulative variance, eigenvalues, and loadings from PCA of winsorized and standardized indicators.",
  c(
    "PCA used the 5,470 institutions with complete values for all nine indicators.",
    "The sign of a principal-component loading is arbitrary; relative magnitudes and directions within a component are interpreted."
  ),
  list(
    `Variance explained` = s3_variance,
    `PCA loadings` = s3_loadings
  ),
  file.path(supp_dir, "Supplementary_Table_S3.csv")
)

write_sectioned_csv(
  4,
  "Supplementary Table S4. K-means diagnostics and characteristics of the five exploratory institutional profiles.",
  "Diagnostics for k = 2–8, profile summaries and centroids, dominant indicators, structural risk category distributions, and winsorization thresholds.",
  c(
    "Profiles are exploratory indicator configurations, not definitive institutional classes.",
    "IRIS Overall score and structural risk category were not clustering inputs.",
    "Diagnostic models used iter.max = 200; the final k = 5 model used iter.max = 500."
  ),
  list(
    `Inclusion summary` = s4_inclusion,
    `Clustering parameters` = s4_parameters,
    `K diagnostics` = s4_k_diagnostics,
    `Profile summary` = s4_profile_summary,
    `Profile centroids z` = s4_centroids,
    `Dominant indicators` = s4_dominant,
    `IRIS categories by profile` = s4_risk_by_profile,
    `Winsor thresholds` = s4_winsor
  ),
  file.path(supp_dir, "Supplementary_Table_S4.csv")
)

write_sectioned_csv(
  5,
  "Supplementary Table S5. Matching of SCImago IRIS countries with World Bank metadata and contextual-data coverage.",
  "IRIS–World Bank matching, analytical coverage, indicator years, unmatched codes, and geographical map matching.",
  c(
    "Contextual indicators use the most recent non-missing value from 2015 onward.",
    "The weighted regression used 59 complete-case countries.",
    "Unmatched geographical codes are absent from the mapped layer but retained in non-spatial analyses."
  ),
  list(
    `IRIS-WB matching` = match_summary |>
      transmute(`Item` = item, `Countries or territories, n` = countries_n),
    `Country analysis coverage` = country_analysis_coverage |>
      transmute(`Item` = item, `n` = n),
    `Contextual data by country` = s5_context_by_country,
    `Context years` = s5_context_years,
    `Unmatched IRIS-WB` = unmatched_iris_wb,
    `Map matching summary` = map_match_summary |>
      transmute(`Item` = item, `Countries or territories, n` = countries_n),
    `Unmatched map codes` = s5_unmatched_map
  ),
  file.path(supp_dir, "Supplementary_Table_S5.csv")
)

write_sectioned_csv(
  6,
  "Supplementary Table S6. Country-level associations and conventional weighted regression estimates.",
  "Prespecified Spearman correlations with bootstrap CIs and unadjusted and Holm-adjusted P values, plus conventional weighted regression, leave-one-country-out estimates, and influence diagnostics.",
  c(
    "Correlation CIs use 10,000 country-level bootstrap resamples.",
    "The four P values were adjusted as one family using the Holm procedure.",
    "Only conventional model-based regression standard errors and confidence intervals are reported.",
    "Countries were weighted by their number of represented IRIS institutions."
  ),
  list(
    `Spearman correlations` = s6_correlations,
    `Weighted regression` = s6_regression,
    `Regression summary` = s6_regression_summary,
    `LOO summary` = s6_loo_summary,
    `LOO estimates` = s6_loo_estimates,
    `Influence diagnostics` = s6_influence
  ),
  file.path(supp_dir, "Supplementary_Table_S6.csv")
)

write_sectioned_csv(
  7,
  "Supplementary Table S7. Multilevel model estimates, model fit, and influence sensitivity analyses.",
  "Primary random-intercept model and sensitivity analyses excluding the verified maximum IRIS Overall score and restricting to countries represented by at least 10 institutions.",
  c(
    "The verified Overall = 37.599 observation was retained in the primary model.",
    "Models were fitted by maximum likelihood.",
    "The sensitivity analysis was performed after identifying the influential observation and is not described as prespecified."
  ),
  list(
    `Model fit summary` = s7_fit,
    `Primary fixed effects` = s7_primary_fixed,
    `Sensitivity fixed effects` = s7_sensitivity_fixed,
    `Random effects` = s7_random,
    `Influential observation` = s7_influential
  ),
  file.path(supp_dir, "Supplementary_Table_S7.csv")
)

write_sectioned_csv(
  8,
  "Supplementary Table S8. Variable dictionary and analytical definitions.",
  "Variable names, definitions, sources, analytical levels, uses, and transformations for raw and derived variables.",
  c(
    "IRIS Overall score, IRIS structural risk category, and exploratory institutional profile are distinct constructs.",
    "The nine individual indicators were clustering inputs; Overall score and structural risk category were not.",
    "profile_mean_overall is a profile-level value repeated for assigned institutions and is distinct from the country-level mean_overall variable."
  ),
  list(`Variable dictionary` = s8_dictionary),
  file.path(supp_dir, "Supplementary_Table_S8.csv")
)

# ------------------------------------------------------------------------------
# 10. Main manuscript tables
# ------------------------------------------------------------------------------

main_table_1 <- s4_profile_summary

main_table_2 <- country_analysis |>
  filter(!is.na(income_group_wb)) |>
  group_by(income_group_wb) |>
  summarise(
    countries_n = n(),
    mean_overall = mean(mean_overall, na.rm = TRUE),
    mean_pct_significant = mean(pct_significant, na.rm = TRUE),
    mean_rd_gdp = mean(rd_gdp, na.rm = TRUE),
    rd_countries_n = sum(!is.na(rd_gdp)),
    mean_researchers_pm = mean(researchers_pm, na.rm = TRUE),
    researcher_countries_n = sum(!is.na(researchers_pm)),
    .groups = "drop"
  )

write_csv(main_table_1, file.path(table_dir, "Main_Table_1.csv"))
write_csv(main_table_2, file.path(table_dir, "Main_Table_2.csv"))

# ------------------------------------------------------------------------------
# 11. Figures with revised terminology
# ------------------------------------------------------------------------------

# Supplementary Figure S1
fig_s1a <- risk_distribution |>
  ggplot(aes(x = risk, y = institutions_percent, fill = risk)) +
  geom_col(width = 0.65, color = "black", linewidth = 0.25) +
  geom_text(
    aes(label = paste0(institutions_n, "\n", round(institutions_percent, 1), "%")),
    vjust = -0.25,
    size = 3.1,
    lineheight = 0.9
  ) +
  scale_fill_manual(values = iris_colors, guide = "none") +
  scale_x_discrete(labels = c(
    "very low" = "Very low", "low" = "Low",
    "medium" = "Medium", "significant" = "Significant"
  )) +
  scale_y_continuous(
    limits = c(0, max(risk_distribution$institutions_percent) * 1.22),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = "IRIS structural risk category",
    y = "Institutions (%)",
    title = "a  IRIS structural risk categories"
  ) +
  theme_paper()

fig_s1b <- top_countries |>
  mutate(country = fct_reorder(country, n_institutions)) |>
  ggplot(aes(country, n_institutions)) +
  geom_col(fill = "#3C5488", width = 0.65) +
  coord_flip() +
  labs(
    x = NULL, y = "Institutions, n",
    title = "b  Countries or territories with the most institutions"
  ) +
  theme_paper()

fig_s1c <- top_significant |>
  mutate(country = fct_reorder(country, pct_significant)) |>
  ggplot(aes(country, pct_significant)) +
  geom_col(fill = "#D95C4F", width = 0.65) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Institutions in the significant category (%)",
    title = "c  Highest percentages in the significant category"
  ) +
  theme_paper()

supplementary_figure_s1 <- fig_s1a / fig_s1b / fig_s1c +
  plot_annotation(
    title = "Global distribution and country representation in the SCImago IRIS dataset"
  )

ggsave(
  file.path(figure_dir, "Supplementary_Figure_S1.png"),
  supplementary_figure_s1,
  width = 9,
  height = 15,
  dpi = 400
)
ggsave(
  file.path(figure_dir, "Supplementary_Figure_S1.pdf"),
  supplementary_figure_s1,
  width = 9,
  height = 15
)

# Figure 1: geographical distribution
map_data_min10 <- map_data |>
  mutate(
    mean_overall_plot = if_else(n_institutions >= 10, mean_overall, NA_real_),
    pct_significant_plot = if_else(
      n_institutions >= 10, pct_significant, NA_real_
    )
  )

figure_1a <- ggplot(map_data_min10) +
  geom_sf(aes(fill = mean_overall_plot), color = "white", linewidth = 0.15) +
  scale_fill_viridis_c(
    option = "inferno",
    na.value = "grey88",
    name = "Mean IRIS Overall score"
  ) +
  labs(title = "a  Mean IRIS Overall score by country") +
  theme_void(base_size = 11) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

figure_1b <- ggplot(map_data_min10) +
  geom_sf(aes(fill = pct_significant_plot), color = "white", linewidth = 0.15) +
  scale_fill_viridis_c(
    option = "inferno",
    na.value = "grey88",
    labels = label_percent(scale = 1),
    name = "Institutions in significant category (%)"
  ) +
  labs(
    title =
      "b  Institutions in the IRIS significant structural risk category (%)"
  ) +
  theme_void(base_size = 11) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

figure_1 <- figure_1a / figure_1b
ggsave(
  file.path(figure_dir, "Figure_1.png"),
  figure_1,
  width = 10,
  height = 10,
  dpi = 400
)
ggsave(
  file.path(figure_dir, "Figure_1.pdf"),
  figure_1,
  width = 10,
  height = 10
)

# Figure 2: profile heatmap
figure_2 <- profile_centroids_long |>
  ggplot(aes(indicator, profile, fill = mean_z)) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(
    aes(label = round(mean_z, 2), color = abs(mean_z) > 1.2),
    size = 3.2
  ) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "white"), guide = "none") +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-2, 2),
    oob = squish,
    name = "Mean z score"
  ) +
  labs(
    title = "Indicator configurations across five exploratory institutional profiles",
    x = NULL,
    y = "Exploratory institutional profile"
  ) +
  theme_paper() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_blank(),
    axis.ticks = element_blank()
  )

ggsave(
  file.path(figure_dir, "Figure_2.png"),
  figure_2,
  width = 10,
  height = 6,
  dpi = 400
)
ggsave(
  file.path(figure_dir, "Figure_2.pdf"),
  figure_2,
  width = 10,
  height = 6
)

# Figure 3: structural categories and PCA distribution
figure_3a <- risk_by_profile |>
  ggplot(aes(profile, within_profile_percent, fill = risk)) +
  geom_col(position = "fill", color = "black", linewidth = 0.25) +
  scale_fill_manual(
    values = iris_colors,
    name = "IRIS structural risk category"
  ) +
  scale_y_continuous(labels = label_percent()) +
  labs(
    title = "a  IRIS structural risk categories",
    x = "Exploratory institutional profile",
    y = "Proportion of institutions"
  ) +
  theme_paper()

figure_3b_data <- pca_scores |>
  left_join(block2_clustered |> select(id, profile), by = "id")

figure_3b <- figure_3b_data |>
  ggplot(aes(PC1, PC2, color = profile)) +
  geom_point(alpha = 0.75, size = 1.4) +
  scale_color_manual(
    values = profile_colors,
    name = "Exploratory institutional profile"
  ) +
  labs(
    title = "b  Principal-component distribution",
    x = paste0(
      "PC1 (", percent(pca_variance$variance_explained[[1]], accuracy = 0.1), ")"
    ),
    y = paste0(
      "PC2 (", percent(pca_variance$variance_explained[[2]], accuracy = 0.1), ")"
    )
  ) +
  theme_paper()

figure_3 <- figure_3a / figure_3b
ggsave(
  file.path(figure_dir, "Figure_3.png"),
  figure_3,
  width = 9,
  height = 12,
  dpi = 400
)
ggsave(
  file.path(figure_dir, "Figure_3.pdf"),
  figure_3,
  width = 9,
  height = 12
)

# Figure 4: national context
figure_4a <- country_analysis |>
  filter(!is.na(income_group_wb)) |>
  ggplot(aes(income_group_wb, mean_overall)) +
  geom_boxplot(outlier.shape = NA, fill = "grey95", width = 0.5) +
  geom_point(
    aes(size = n_institutions),
    color = "#3C5488",
    alpha = 0.65,
    position = position_jitter(width = 0.14, height = 0)
  ) +
  labs(
    title = "a  World Bank income group",
    x = NULL,
    y = "Mean IRIS Overall score",
    size = "Institutions"
  ) +
  theme_paper() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

figure_4b <- country_analysis |>
  filter(!is.na(rd_gdp), !is.na(mean_overall)) |>
  ggplot(aes(rd_gdp, mean_overall)) +
  geom_point(aes(size = n_institutions), color = "#3C5488", alpha = 0.65) +
  geom_smooth(method = "lm", se = TRUE, color = "#D95C4F", fill = "grey80") +
  labs(
    title = "b  R&D expenditure",
    x = "R&D expenditure (% GDP)",
    y = "Mean IRIS Overall score",
    size = "Institutions"
  ) +
  theme_paper()

figure_4c <- country_analysis |>
  filter(!is.na(researchers_pm), !is.na(mean_overall)) |>
  ggplot(aes(researchers_pm, mean_overall)) +
  geom_point(aes(size = n_institutions), color = "#3C5488", alpha = 0.65) +
  geom_smooth(method = "lm", se = TRUE, color = "#D95C4F", fill = "grey80") +
  scale_x_log10(labels = label_comma()) +
  labs(
    title = "c  Researcher density",
    x = "Researchers in R&D per million people (log scale)",
    y = "Mean IRIS Overall score",
    size = "Institutions"
  ) +
  theme_paper()

figure_4 <- figure_4a / figure_4b / figure_4c
ggsave(
  file.path(figure_dir, "Figure_4.png"),
  figure_4,
  width = 9,
  height = 15,
  dpi = 400
)
ggsave(
  file.path(figure_dir, "Figure_4.pdf"),
  figure_4,
  width = 9,
  height = 15
)

# Figure 5: profile composition by income group
figure_5 <- profile_income_distribution |>
  ggplot(aes(income_group_wb, percent, fill = profile)) +
  geom_col(position = "fill", color = "black", linewidth = 0.3, width = 0.65) +
  scale_fill_manual(
    values = profile_colors,
    name = "Institutional profile"
  ) +
  scale_y_continuous(labels = label_percent()) +
  labs(
    title = "Exploratory institutional profiles by World Bank income group",
    x = NULL,
    y = "Proportion of institutions"
  ) +
  theme_paper() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

ggsave(
  file.path(figure_dir, "Figure_5.png"),
  figure_5,
  width = 8,
  height = 6,
  dpi = 400
)
ggsave(
  file.path(figure_dir, "Figure_5.pdf"),
  figure_5,
  width = 8,
  height = 6
)

# ------------------------------------------------------------------------------
# 12. Final checks and session information
# ------------------------------------------------------------------------------

expected_fit <- model_fit_summary |>
  select(model, adjusted_icc, marginal_r2, conditional_r2)

print(global_summary)
print(country_correlations)
print(weighted_regression)
print(expected_fit)

stopifnot(
  sum(as.character(data_clean$risk) == "significant", na.rm = TRUE) == 365,
  n_distinct(block2_clustered$profile) == 5,
  nrow(country_analysis) == 69,
  sum(country_analysis$n_institutions) == 5213,
  nrow(min10_model_data) == 4718,
  n_distinct(min10_model_data$country) == 59
)

capture.output(
  sessionInfo(),
  file = file.path(output_root, "analysis_session_info.txt")
)

message("Analysis completed. Outputs written to: ", normalizePath(output_root))
