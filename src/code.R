# ============================================================
# INITIAL SETUP AND DATA IMPORT
# ============================================================

# 1. Load required packages ----------------------------------
# (Install if not present)
packages <- c("tidyverse", "janitor", "readr", "stringr", "forcats", 
              "scales", "sf", "rnaturalearth", "rnaturalearthdata", 
              "viridis", "patchwork", "cluster", "WDI", "countrycode", 
              "lme4", "broom.mixed", "performance")

installed <- packages %in% rownames(installed.packages())
if (any(!installed)) {
  install.packages(packages[!installed])
}

invisible(lapply(packages, library, character.only = TRUE))

# 2. Output directories setup --------------------------------
base_dir <- "outputs"
dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(base_dir, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(base_dir, "figures"), showWarnings = FALSE, recursive = TRUE)

# 3. Graphical Theme (Q1: Science/Nature style) --------------
theme_paper <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(family = "sans", color = "black"),
      plot.title = element_text(face = "bold", size = base_size + 2, margin = margin(b = 10)),
      plot.subtitle = element_text(size = base_size, color = "grey30", margin = margin(b = 10)),
      plot.caption = element_text(size = base_size - 3, color = "grey50", hjust = 0, margin = margin(t = 15)),
      axis.title = element_text(face = "bold", color = "black"),
      axis.text = element_text(color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = base_size),
      legend.background = element_blank(),
      plot.margin = margin(15, 15, 15, 15)
    )
}

# 4. Data Import ---------------------------------------------
data <- read_csv(
  "data/Scimago_IRIS_Index_Data.csv",
  na = c("", "NA", "background: #ffffff; color: #000000;"),
  show_col_types = FALSE
) |> 
  clean_names()

# ============================================================
# BLOCK 1
# Global overview of institutional research integrity risk
# Dataset: SCImago IRIS
# Objeto de entrada esperado: data
# ============================================================

# Paquetes ----------------------------------------------------


# Carpetas de salida -----------------------------------------

# Define base path

# Creamos la carpeta principal 'outputs' y sus subcarpetas necesarias

# Graphical Theme (Q1: Science/Nature style) ------


# ============================================================
# 1. Basic cleaning and type checking
# ============================================================

# Revisar nombres
names(data)
glimpse(data)

# Expected numeric variables
numeric_vars <- c(
  "sir_rank", "overall", "output", "multi_affiliation", "retracted_output",
  "self_citation", "discontinued_journals_output", "hyperauthored_output",
  "leadership_impact_gap", "hyperprolific_authors", "institutional_journal_output",
  "redundant_output"
)

# Check for missing columns
missing_numeric_vars <- setdiff(numeric_vars, names(data))

if (length(missing_numeric_vars) > 0) {
  stop(
    paste(
      "These expected columns are not in the dataset:",
      paste(missing_numeric_vars, collapse = ", ")
    )
  )
}

# Cleaning
data_clean <- data |>
  mutate(
    across(
      all_of(numeric_vars),
      ~ parse_number(as.character(.x))
    ),
    institution = str_squish(as.character(institution)),
    country = str_squish(as.character(country)),
    risk = str_squish(str_to_lower(as.character(risk))),
    risk = factor(
      risk,
      levels = c("very low", "low", "medium", "significant"),
      ordered = TRUE
    )
  )

# Save clean dataset
write_csv(
  data_clean,
  file.path(base_dir, "tables/block1_data_clean.csv")
)

# Missing values summary
missing_summary <- data_clean |>
  summarise(across(everything(), ~ sum(is.na(.x)))) |>
  pivot_longer(
    cols = everything(),
    names_to = "variable",
    values_to = "missing_n"
  ) |>
  arrange(desc(missing_n))


write_csv(
  missing_summary,
  file.path(base_dir, "tables/block1_missing_summary.csv")
)

# ============================================================
# 2. Global dataset summary
# ============================================================

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


write_csv(
  global_summary,
  file.path(base_dir, "tables/block1_global_summary.csv")
)

# Summary Table 1 for manuscript
table_1_global <- tibble(
  indicator = c(
    "Institutions", "Unique institutions", "Countries",
    "Mean Overall", "SD Overall", "Median Overall", "IQR Overall", "Minimum Overall", "Maximum Overall",
    "Mean Output", "SD Output", "Median Output", "IQR Output", "Minimum Output", "Maximum Output"
  ),
  value = c(
    global_summary$n_institutions, global_summary$n_unique_institutions, global_summary$n_countries,
    round(global_summary$overall_mean, 3), round(global_summary$overall_sd, 3), round(global_summary$overall_median, 3),
    round(global_summary$overall_iqr, 3), round(global_summary$overall_min, 3), round(global_summary$overall_max, 3),
    round(global_summary$output_mean, 1), round(global_summary$output_sd, 1), round(global_summary$output_median, 1),
    round(global_summary$output_iqr, 1), round(global_summary$output_min, 1), round(global_summary$output_max, 1)
  )
)


write_csv(
  table_1_global,
  file.path(base_dir, "tables/table_1_global_summary.csv")
)

# ============================================================
# 3. Global distribution by risk category
# ============================================================

risk_distribution <- data_clean |>
  count(risk, name = "n") |>
  mutate(
    percent = n / sum(n) * 100,
    percent_label = paste0(round(percent, 1), "%"),
    risk_label = str_to_title(as.character(risk))
  )


write_csv(
  risk_distribution,
  file.path(base_dir, "tables/block1_risk_distribution.csv")
)

# Figure: risk distribution
fig_risk_distribution <- risk_distribution |>
  ggplot(aes(x = risk, y = percent, fill = risk)) +
  geom_col(width = 0.6, alpha = 0.9, color = "black", linewidth = 0.3) +
  geom_text(
    aes(label = paste0(n, "\n", percent_label)),
    vjust = -0.3,
    size = 3.5,
    lineheight = 0.9
  ) +
  # Official SCImago IRIS color palette from website
  scale_fill_manual(
    values = c(
      "very low" = "#9EC753",     # Verde
      "low" = "#E3CE47",          # Amarillo
      "medium" = "#F49D4A",       # Naranja
      "significant" = "#D95C4F"   # Rojo
    ),
    guide = "none"
  ) +
  scale_x_discrete(
    labels = c(
      "very low" = "Very low",
      "low" = "Low",
      "medium" = "Medium",
      "significant" = "Significant"
    )
  ) +
  scale_y_continuous(
    labels = label_percent(scale = 1),
    limits = c(0, max(risk_distribution$percent, na.rm = TRUE) * 1.25),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = "Integrity Risk Category",
    y = "Proportion of Institutions (%)",
    title = "Distribution of institutional integrity risk categories",
    caption = "Source: SCImago IRIS dataset."
  ) +
  theme_paper()


ggsave(
  filename = "outputs/figures/block1_risk_distribution.png",
  plot = fig_risk_distribution,
  width = 7, height = 5, dpi = 300
)

# ============================================================
# 4. Countries with the most institutions
# ============================================================

top_countries_by_institutions <- data_clean |>
  count(country, name = "n_institutions") |>
  arrange(desc(n_institutions)) |>
  slice_head(n = 20)

write_csv(
  top_countries_by_institutions,
  file.path(base_dir, "tables/block1_top20_countries_by_institutions.csv")
)

fig_top_countries <- top_countries_by_institutions |>
  mutate(country = fct_reorder(country, n_institutions)) |>
  # Solid scientific blue. Q1 journals avoid gradients unless mapping useful continuous variables.
  ggplot(aes(x = country, y = n_institutions)) +
  geom_col(width = 0.65, fill = "#3C5488FF", alpha = 0.9) +
  geom_text(
    aes(label = n_institutions),
    hjust = -0.2,
    size = 3.5
  ) +
  coord_flip() +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.15))
  ) +
  labs(
    x = NULL,
    y = "Number of institutions",
    title = "Top 20 countries by number of institutions",
    caption = "Source: SCImago IRIS dataset."
  ) +
  theme_paper() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank()) # Cleaner look for horizontal bars


ggsave(
  filename = "outputs/figures/block1_top20_countries_by_institutions.png",
  plot = fig_top_countries,
  width = 7, height = 6, dpi = 300
)

# ============================================================
# 5. Aggregated summary by country
# ============================================================

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
    pct_medium_or_significant = mean(risk %in% c("medium", "significant"), na.rm = TRUE) * 100,
    mean_output = mean(output, na.rm = TRUE),
    median_output = median(output, na.rm = TRUE),
    mean_multi_affiliation = mean(multi_affiliation, na.rm = TRUE),
    mean_retracted_output = mean(retracted_output, na.rm = TRUE),
    mean_self_citation = mean(self_citation, na.rm = TRUE),
    mean_discontinued_journals_output = mean(discontinued_journals_output, na.rm = TRUE),
    mean_hyperauthored_output = mean(hyperauthored_output, na.rm = TRUE),
    mean_leadership_impact_gap = mean(leadership_impact_gap, na.rm = TRUE),
    mean_hyperprolific_authors = mean(hyperprolific_authors, na.rm = TRUE),
    mean_institutional_journal_output = mean(institutional_journal_output, na.rm = TRUE),
    mean_redundant_output = mean(redundant_output, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(desc(n_institutions))

write_csv(
  country_summary,
  file.path(base_dir, "tables/block1_country_summary.csv")
)

# ============================================================
# 6. Countries with highest proportion of significant risk
# ============================================================

countries_significant_risk_min10 <- country_summary |>
  filter(n_institutions >= 10) |>
  arrange(desc(pct_significant)) |>
  slice_head(n = 20)

write_csv(
  countries_significant_risk_min10,
  file.path(base_dir, "tables/block1_top20_countries_significant_risk_min10.csv")
)

fig_significant_risk_countries <- countries_significant_risk_min10 |>
  mutate(country = fct_reorder(country, pct_significant)) |>
  # Solid red to highlight "significant risk"
  ggplot(aes(x = country, y = pct_significant)) +
  geom_col(width = 0.65, fill = "#E64B35FF", alpha = 0.9) +
  geom_text(
    aes(label = paste0(round(pct_significant, 1), "%")),
    hjust = -0.2,
    size = 3.5
  ) +
  coord_flip() +
  scale_y_continuous(
    labels = label_percent(scale = 1),
    expand = expansion(mult = c(0, 0.15))
  ) +
  labs(
    x = NULL,
    y = "Institutions with significant risk (%)",
    title = "Top countries by proportion of institutions with significant risk",
    subtitle = "Included countries with ≥ 10 institutions",
    caption = "Source: SCImago IRIS dataset."
  ) +
  theme_paper() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())


ggsave(
  filename = "outputs/figures/block1_top20_countries_significant_risk_min10.png",
  plot = fig_significant_risk_countries,
  width = 8, height = 6, dpi = 300
)

# ============================================================
# 7. World map with ISO3 corrections
# ============================================================

world_raw <- ne_countries(
  scale = "medium",
  returnclass = "sf"
) |>
  filter(name != "Antarctica")

world_map <- world_raw |>
  mutate(
    join_iso3 = case_when(
      !is.na(iso_a3) & iso_a3 != "-99" ~ iso_a3,
      !is.na(adm0_a3) & adm0_a3 != "-99" ~ adm0_a3,
      TRUE ~ NA_character_
    )
  ) |>
  mutate(
    join_iso3 = case_when(
      name == "Kosovo" ~ "XKX",
      TRUE ~ join_iso3
    )
  ) |>
  select(name, iso_a3, adm0_a3, join_iso3, geometry)

map_data <- world_map |>
  left_join(
    country_summary,
    by = c("join_iso3" = "country")
  )

unmatched_countries <- country_summary |>
  anti_join(
    world_map,
    by = c("country" = "join_iso3")
  ) |>
  arrange(desc(n_institutions))

write_csv(unmatched_countries, file.path(base_dir, "tables/block1_unmatched_countries_map_after_fix.csv"))

saveRDS(map_data, file.path(base_dir, "tables/block1_map_data_fixed.rds"))

map_check_summary <- tibble(
  item = c(
    "Countries in SCImago IRIS country_summary",
    "Countries matched with map",
    "Countries unmatched after correction"
  ),
  n = c(
    n_distinct(country_summary$country),
    n_distinct(map_data$join_iso3[!is.na(map_data$mean_overall)]),
    nrow(unmatched_countries)
  )
)

write_csv(map_check_summary, file.path(base_dir, "tables/block1_map_check_summary.csv"))

# ============================================================
# 8. Recommended map for manuscript
# ============================================================

map_data_min10 <- map_data |>
  mutate(
    mean_overall_plot = if_else(
      n_institutions >= 10,
      mean_overall,
      NA_real_
    ),
    pct_significant_plot = if_else(
      n_institutions >= 10,
      pct_significant,
      NA_real_
    )
  )

saveRDS(map_data_min10, file.path(base_dir, "tables/block1_map_data_min10_fixed.rds"))

fig_world_mean_overall_min10_publication <- map_data_min10 |>
  ggplot() +
  geom_sf(
    aes(fill = mean_overall_plot),
    color = "white", # White borders for a polished look
    linewidth = 0.2
  ) +
  # Changed to inferno as requested
  scale_fill_viridis_c(
    option = "inferno",
    na.value = "grey90",
    name = "Mean Overall"
  ) +
  labs(
    title = "Mean institutional integrity risk score by country",
    subtitle = "Included countries with ≥ 10 institutions",
    caption = "Source: SCImago IRIS dataset. Countries with fewer than 10 institutions or without data are shown in grey."
  ) +
  theme_void(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(2, "cm"),
    plot.title = element_text(face = "bold", size = 14, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 11, color = "grey30", margin = margin(b = 15)),
    plot.caption = element_text(size = 9, color = "grey50", hjust = 0, margin = margin(t = 15)),
    legend.title = element_text(face = "bold"),
    plot.margin = margin(15, 15, 15, 15)
  )


ggsave(
  filename = "outputs/figures/figure_world_mean_overall_min10_publication.png",
  plot = fig_world_mean_overall_min10_publication,
  width = 10, height = 6, dpi = 300
)
ggsave(
  filename = "outputs/figures/figure_world_mean_overall_min10_publication.pdf",
  plot = fig_world_mean_overall_min10_publication,
  width = 10, height = 6
)

# ============================================================
# 9. Recommended map for supplement or secondary figure
# ============================================================

fig_world_pct_significant_min10_supplement <- map_data_min10 |>
  ggplot() +
  geom_sf(
    aes(fill = pct_significant_plot),
    color = "white",
    linewidth = 0.2
  ) +
  scale_fill_viridis_c(
    option = "inferno",
    na.value = "grey90",
    labels = label_percent(scale = 1),
    name = "Significant risk"
  ) +
  labs(
    title = "Institutions classified as significant integrity risk by country",
    subtitle = "Included countries with ≥ 10 institutions",
    caption = "Source: SCImago IRIS dataset. Countries with fewer than 10 institutions or without data are shown in grey."
  ) +
  theme_void(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(2, "cm"),
    plot.title = element_text(face = "bold", size = 14, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 11, color = "grey30", margin = margin(b = 15)),
    plot.caption = element_text(size = 9, color = "grey50", hjust = 0, margin = margin(t = 15)),
    legend.title = element_text(face = "bold"),
    plot.margin = margin(15, 15, 15, 15)
  )


ggsave(
  filename = "outputs/figures/figure_s1_world_pct_significant_min10.png",
  plot = fig_world_pct_significant_min10_supplement,
  width = 10, height = 6, dpi = 300
)
ggsave(
  filename = "outputs/figures/figure_s1_world_pct_significant_min10.pdf",
  plot = fig_world_pct_significant_min10_supplement,
  width = 10, height = 6
)

# ============================================================
# 11. Export list of generated main files
# ============================================================

generated_files <- tibble(
  type = c("clean_data", "table", "table", "table", "table", "table", "table"),
  file = c(
    file.path(base_dir, "tables/block1_data_clean.csv"),
    file.path(base_dir, "tables/block1_missing_summary.csv"),
    file.path(base_dir, "tables/block1_global_summary.csv"),
    file.path(base_dir, "tables/table_1_global_summary.csv"),
    file.path(base_dir, "tables/block1_risk_distribution.csv"),
    file.path(base_dir, "tables/block1_country_summary.csv"),
    file.path(base_dir, "tables/block1_unmatched_countries_map_after_fix.csv")
  )
)
write_csv(generated_files, file.path(base_dir, "tables/block1_generated_files.csv"))

# End of Block 1

# ============================================================
# BLOCK 2: Institutional risk profiles
# PCA + clustering + profile heatmap
# Base dataset: data_clean
# ============================================================


# Instalar si no lo tienes
}


# ------------------------------------------------------------
# Output directories and Q1 Graphical Theme
# ------------------------------------------------------------


# Create output folders

# Q1 Graphical Theme

# Official SCImago IRIS risk palette
iris_colors <- c(
  "very low" = "#9EC753",
  "low" = "#E3CE47",
  "medium" = "#F49D4A",
  "significant" = "#D95C4F"
)

# ------------------------------------------------------------
# 1. Define IRIS indicators for profile analysis
# ------------------------------------------------------------

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

# ------------------------------------------------------------
# 2. Prepare analytical dataset
# ------------------------------------------------------------

block2_data <- data_clean |>
  mutate(
    risk = factor(
      as.character(risk),
      levels = c("very low", "low", "medium", "significant"),
      ordered = TRUE
    )
  ) |>
  select(
    id, institution, country, sir_rank, overall, risk, output, all_of(risk_indicator_vars)
  )

block2_missing <- block2_data |>
  summarise(across(all_of(risk_indicator_vars), ~ sum(is.na(.x)))) |>
  pivot_longer(
    cols = everything(),
    names_to = "indicator",
    values_to = "missing_n"
  ) |>
  mutate(
    indicator_label = recode(indicator, !!!indicator_labels)
  ) |>
  arrange(desc(missing_n))

write_csv(
  block2_missing,
  file.path(base_dir, "tables/block2_missing_indicators.csv")
)

block2_complete <- block2_data |>
  filter(if_all(all_of(risk_indicator_vars), ~ !is.na(.x)))

block2_inclusion_summary <- tibble(
  item = c(
    "Institutions in data_clean",
    "Institutions with complete IRIS indicators",
    "Institutions excluded due to missing IRIS indicators"
  ),
  n = c(
    nrow(data_clean),
    nrow(block2_complete),
    nrow(data_clean) - nrow(block2_complete)
  )
)

write_csv(
  block2_inclusion_summary,
  file.path(base_dir, "tables/block2_inclusion_summary.csv")
)

# ------------------------------------------------------------
# 3. Winsorization
# ------------------------------------------------------------

winsorize_vec <- function(x, probs = c(0.01, 0.99)) {
  qs <- quantile(x, probs = probs, na.rm = TRUE, type = 7)
  x <- pmax(x, qs[[1]], na.rm = TRUE)
  x <- pmin(x, qs[[2]], na.rm = TRUE)
  return(x)
}

winsor_thresholds <- map_dfr(
  risk_indicator_vars,
  function(v) {
    tibble(
      indicator = v,
      indicator_label = indicator_labels[[v]],
      p01 = quantile(block2_complete[[v]], 0.01, na.rm = TRUE),
      p99 = quantile(block2_complete[[v]], 0.99, na.rm = TRUE)
    )
  }
)

write_csv(
  winsor_thresholds,
  file.path(base_dir, "tables/block2_winsor_thresholds.csv")
)

block2_winsor <- block2_complete |>
  mutate(
    across(
      all_of(risk_indicator_vars),
      ~ winsorize_vec(.x),
      .names = "{.col}_w"
    )
  )

winsor_vars <- paste0(risk_indicator_vars, "_w")

indicator_matrix <- block2_winsor |>
  select(all_of(winsor_vars)) |>
  as.matrix()

colnames(indicator_matrix) <- risk_indicator_vars
indicator_matrix_scaled <- scale(indicator_matrix)

scaling_params <- tibble(
  indicator = colnames(indicator_matrix_scaled),
  indicator_label = indicator_labels[indicator],
  mean_used_for_scaling = attr(indicator_matrix_scaled, "scaled:center"),
  sd_used_for_scaling = attr(indicator_matrix_scaled, "scaled:scale")
)

write_csv(
  scaling_params,
  file.path(base_dir, "tables/block2_scaling_parameters.csv")
)

indicator_scaled_df <- as_tibble(indicator_matrix_scaled) |>
  mutate(
    id = block2_winsor$id,
    institution = block2_winsor$institution,
    country = block2_winsor$country,
    overall = block2_winsor$overall,
    risk = block2_winsor$risk,
    output = block2_winsor$output
  ) |>
  relocate(id, institution, country, overall, risk, output)

write_csv(
  indicator_scaled_df,
  file.path(base_dir, "tables/block2_scaled_indicators.csv")
)

# ------------------------------------------------------------
# 4. Correlation matrix between indicators
# ------------------------------------------------------------

cor_mat <- cor(
  indicator_matrix,
  method = "spearman",
  use = "pairwise.complete.obs"
)

cor_df <- as.data.frame(cor_mat) |>
  rownames_to_column("indicator_1") |>
  pivot_longer(
    cols = -indicator_1,
    names_to = "indicator_2",
    values_to = "spearman_r"
  ) |>
  mutate(
    indicator_1_label = recode(indicator_1, !!!indicator_labels),
    indicator_2_label = recode(indicator_2, !!!indicator_labels)
  )

write_csv(
  cor_df,
  file.path(base_dir, "tables/block2_spearman_correlations.csv")
)

fig_cor_heatmap <- cor_df |>
  ggplot(aes(x = indicator_1_label, y = indicator_2_label, fill = spearman_r)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient2(
    low = "#3C5488FF", # Scientific blue
    mid = "white",
    high = "#E64B35FF", # Rojo científico
    midpoint = 0,
    limits = c(-1, 1),
    name = "Spearman r"
  ) +
  labs(
    title = "Correlation structure among SCImago IRIS risk indicators",
    x = NULL, y = NULL,
    caption = "Source: SCImago IRIS dataset."
  ) +
  theme_paper() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_blank(),   # Eliminar líneas de eje para heatmaps
    axis.ticks = element_blank(),
    legend.key.width = unit(1.5, "cm")
  )


ggsave(
  filename = file.path(base_dir, "figures/block2_correlation_heatmap.png"),
  plot = fig_cor_heatmap,
  width = 8, height = 7, dpi = 300
)

# ------------------------------------------------------------
# 5. PCA on standardized IRIS indicators
# ------------------------------------------------------------

pca_fit <- prcomp(
  indicator_matrix_scaled,
  center = FALSE,
  scale. = FALSE
)

pca_variance <- tibble(
  pc = paste0("PC", seq_along(pca_fit$sdev)),
  eigenvalue = pca_fit$sdev^2,
  variance_explained = eigenvalue / sum(eigenvalue),
  cumulative_variance = cumsum(variance_explained)
)

write_csv(
  pca_variance,
  file.path(base_dir, "tables/block2_pca_variance.csv")
)

pca_loadings <- as_tibble(
  pca_fit$rotation,
  rownames = "indicator"
) |>
  mutate(
    indicator_label = recode(indicator, !!!indicator_labels)
  ) |>
  relocate(indicator, indicator_label)

write_csv(
  pca_loadings,
  file.path(base_dir, "tables/block2_pca_loadings.csv")
)

pca_scores <- as_tibble(pca_fit$x) |>
  mutate(
    id = block2_winsor$id,
    institution = block2_winsor$institution,
    country = block2_winsor$country,
    overall = block2_winsor$overall,
    risk = block2_winsor$risk,
    output = block2_winsor$output
  ) |>
  relocate(id, institution, country, overall, risk, output)

write_csv(
  pca_scores,
  file.path(base_dir, "tables/block2_pca_scores.csv")
)

# Scree plot
fig_pca_scree <- pca_variance |>
  mutate(pc = factor(pc, levels = pc)) |>
  ggplot(aes(x = pc, y = variance_explained)) +
  geom_col(fill = "#3C5488FF", width = 0.6) +
  geom_line(aes(group = 1), color = "black", linewidth = 0.8) +
  geom_point(size = 3, color = "black", fill = "white", shape = 21, stroke = 1) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.1))
  ) +
  labs(
    title = "Variance explained by principal components",
    x = "Principal component",
    y = "Variance explained",
    caption = "PCA based on winsorized and standardized SCImago IRIS risk indicators."
  ) +
  theme_paper()


ggsave(
  filename = file.path(base_dir, "figures/block2_pca_scree_plot.png"),
  plot = fig_pca_scree,
  width = 7, height = 5, dpi = 300
)

# PCA scatter plot (Official IRIS Colors)
fig_pca_risk <- pca_scores |>
  ggplot(aes(x = PC1, y = PC2, color = risk)) +
  geom_point(alpha = 0.7, size = 1.8, stroke = 0) +
  scale_color_manual(values = iris_colors) +
  labs(
    title = "PCA of institutional risk indicators by IRIS risk category",
    x = paste0("PC1 (", percent(pca_variance$variance_explained[1], accuracy = 0.1), ")"),
    y = paste0("PC2 (", percent(pca_variance$variance_explained[2], accuracy = 0.1), ")"),
    color = "Risk category",
    caption = "PCA based on winsorized and standardized SCImago IRIS risk indicators."
  ) +
  theme_paper() +
  theme(legend.position = "bottom")


ggsave(
  filename = file.path(base_dir, "figures/block2_pca_by_risk_category.png"),
  plot = fig_pca_risk,
  width = 8, height = 6, dpi = 300
)

# ------------------------------------------------------------
# 6. Diagnostics to choose number of clusters
# ------------------------------------------------------------

set.seed(2026)
k_grid <- 2:8
silhouette_sample_size <- min(3000, nrow(indicator_matrix_scaled))
silhouette_sample <- sample(seq_len(nrow(indicator_matrix_scaled)), silhouette_sample_size)
dist_silhouette <- dist(indicator_matrix_scaled[silhouette_sample, ])

k_diagnostics <- map_dfr(
  k_grid,
  function(k) {
    set.seed(2026 + k)
    km <- kmeans(indicator_matrix_scaled, centers = k, nstart = 50, iter.max = 200)
    sil <- silhouette(km$cluster[silhouette_sample], dist_silhouette)
    tibble(k = k, total_withinss = km$tot.withinss, avg_silhouette = mean(sil[, 3]))
  }
)

write_csv(
  k_diagnostics,
  file.path(base_dir, "tables/block2_kmeans_k_diagnostics.csv")
)

fig_k_diagnostics_wss <- k_diagnostics |>
  ggplot(aes(x = k, y = total_withinss)) +
  geom_line(color = "#3C5488FF", linewidth = 0.8) +
  geom_point(size = 3, color = "#3C5488FF") +
  scale_x_continuous(breaks = k_grid) +
  labs(
    title = "Elbow method for k-means clustering",
    x = "Number of clusters",
    y = "Total within-cluster sum of squares"
  ) +
  theme_paper()

fig_k_diagnostics_silhouette <- k_diagnostics |>
  ggplot(aes(x = k, y = avg_silhouette)) +
  geom_line(color = "#E64B35FF", linewidth = 0.8) +
  geom_point(size = 3, color = "#E64B35FF") +
  scale_x_continuous(breaks = k_grid) +
  labs(
    title = "Average silhouette by number of clusters",
    x = "Number of clusters",
    y = "Average silhouette"
  ) +
  theme_paper()

fig_k_diagnostics <- fig_k_diagnostics_wss / fig_k_diagnostics_silhouette

ggsave(
  filename = file.path(base_dir, "figures/block2_k_diagnostics.png"),
  plot = fig_k_diagnostics,
  width = 8, height = 8, dpi = 300
)

# ------------------------------------------------------------
# 7. Final clustering
# ------------------------------------------------------------

k_final <- 5
set.seed(2026)

km_final <- kmeans(
  indicator_matrix_scaled, centers = k_final, nstart = 100, iter.max = 500
)

block2_clustered_raw <- block2_winsor |> mutate(cluster_raw = km_final$cluster)

cluster_order <- block2_clustered_raw |>
  group_by(cluster_raw) |>
  summarise(mean_overall = mean(overall, na.rm = TRUE), .groups = "drop") |>
  arrange(mean_overall) |>
  mutate(profile = paste0("Profile ", row_number()))

block2_clustered <- block2_clustered_raw |>
  left_join(cluster_order, by = "cluster_raw") |>
  mutate(profile = factor(profile, levels = paste0("Profile ", seq_len(k_final))))

write_csv(
  block2_clustered,
  file.path(base_dir, "tables/block2_cluster_assignments.csv")
)

# ------------------------------------------------------------
# 8-10. Summary of profiles and dominant indicators
# ------------------------------------------------------------

cluster_summary <- block2_clustered |>
  group_by(profile) |>
  summarise(
    n = n(), percent = n / nrow(block2_clustered) * 100,
    countries = n_distinct(country), mean_overall = mean(overall, na.rm = TRUE),
    median_overall = median(overall, na.rm = TRUE), mean_output = mean(output, na.rm = TRUE),
    median_output = median(output, na.rm = TRUE),
    pct_very_low = mean(risk == "very low", na.rm = TRUE) * 100,
    pct_low = mean(risk == "low", na.rm = TRUE) * 100,
    pct_medium = mean(risk == "medium", na.rm = TRUE) * 100,
    pct_significant = mean(risk == "significant", na.rm = TRUE) * 100,
    .groups = "drop"
  )

write_csv(cluster_summary, file.path(base_dir, "tables/block2_cluster_summary.csv"))

scaled_indicators_with_profile <- as_tibble(indicator_matrix_scaled) |>
  mutate(id = block2_winsor$id) |>
  left_join(block2_clustered |> select(id, profile), by = "id")

cluster_indicator_profiles_z <- scaled_indicators_with_profile |>
  group_by(profile) |>
  summarise(across(all_of(risk_indicator_vars), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

write_csv(cluster_indicator_profiles_z, file.path(base_dir, "tables/block2_cluster_indicator_profiles_z.csv"))

cluster_indicator_profiles_long <- cluster_indicator_profiles_z |>
  pivot_longer(cols = all_of(risk_indicator_vars), names_to = "indicator", values_to = "mean_z") |>
  mutate(
    indicator_label = recode(indicator, !!!indicator_labels),
    indicator_label = factor(indicator_label, levels = indicator_labels[risk_indicator_vars])
  )

write_csv(cluster_indicator_profiles_long, file.path(base_dir, "tables/block2_cluster_indicator_profiles_long.csv"))

cluster_dominant_indicators <- cluster_indicator_profiles_long |>
  group_by(profile) |> arrange(desc(mean_z), .by_group = TRUE) |> slice_head(n = 3) |> ungroup()

write_csv(cluster_dominant_indicators, file.path(base_dir, "tables/block2_cluster_dominant_indicators.csv"))

cluster_dominant_summary <- cluster_dominant_indicators |>
  group_by(profile) |>
  summarise(
    dominant_indicators = paste0(indicator_label, " (z = ", round(mean_z, 2), ")", collapse = "; "),
    .groups = "drop"
  ) |>
  left_join(cluster_summary, by = "profile") |>
  select(profile, n, percent, mean_overall, pct_significant, dominant_indicators)

write_csv(cluster_dominant_summary, file.path(base_dir, "tables/block2_cluster_dominant_summary.csv"))

# ------------------------------------------------------------
# 11. Institutional profiles heatmap
# ------------------------------------------------------------

fig_cluster_heatmap <- cluster_indicator_profiles_long |>
  ggplot(aes(x = indicator_label, y = profile, fill = mean_z)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(
    # Evaluate contrast directly inside aes() for safety
    aes(
      label = round(mean_z, 2),
      color = abs(mean_z) > 1.2 
    ),
    size = 3.5
  ) +
  # Force text colors based on the previous evaluation
  scale_color_manual(
    values = c("FALSE" = "black", "TRUE" = "white"),
    guide = "none" # Hide this sub-legend
  ) +
  scale_fill_gradient2(
    low = "#2166AC",       
    mid = "white", 
    high = "#B2182B",      
    midpoint = 0, 
    limits = c(-2, 2),
    # Force specific breaks in the legend
    breaks = c(-2, -1, 0, 1, 2),
    # Modify visual labels adding "less than or equal" and "greater than or equal"
    labels = c("\u2264 -2", "-1", "0", "1", "\u2265 2"), 
    oob = scales::squish,  
    name = "Mean z-score"
  ) +
  labs(
    title = "Institutional integrity risk profiles based on SCImago IRIS indicators",
    subtitle = paste0("K-means clustering with k = ", k_final),
    x = NULL, y = NULL,
    caption = "Indicators were winsorized and standardized. Color scale is capped at z = ±2 for visual contrast."
  ) +
  theme_paper() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.key.width = unit(2, "cm")
  )


ggsave(
  filename = file.path(base_dir, "figures/block2_cluster_profile_heatmap.png"),
  plot = fig_cluster_heatmap,
  width = 10, height = 6, dpi = 300
)

# ------------------------------------------------------------
# 12. PCA colored by profiles
# ------------------------------------------------------------

pca_scores_clustered <- pca_scores |>
  left_join(block2_clustered |> select(id, profile), by = "id")

fig_pca_clusters <- pca_scores_clustered |>
  ggplot(aes(x = PC1, y = PC2, color = profile)) +
  # Slightly increased alpha to 0.8 for more solid colors
  geom_point(alpha = 0.8, size = 1.8, stroke = 0) +
  # Okabe-Ito palette: High contrast, recommended by Nature Methods
  scale_color_manual(values = c(
    "Profile 1" = "#0072B2", # Vibrant dark blue
    "Profile 2" = "#D55E00", # Tile red / Vermilion
    "Profile 3" = "#009E73", # Intense emerald green
    "Profile 4" = "#CC79A7", # Soft purple / Magenta
    "Profile 5" = "#E69F00"  # Golden orange
  )) +
  labs(
    title = "PCA projection of institutional risk profiles",
    x = paste0("PC1 (", percent(pca_variance$variance_explained[1], accuracy = 0.1), ")"),
    y = paste0("PC2 (", percent(pca_variance$variance_explained[2], accuracy = 0.1), ")"),
    color = "Risk profile",
    caption = "PCA based on winsorized and standardized SCImago IRIS risk indicators."
  ) +
  theme_paper() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )


  ggsave(
  filename = file.path(base_dir, "figures/block2_pca_clusters.png"),
  plot = fig_pca_clusters,
  width = 8, height = 6, dpi = 300
)

# ------------------------------------------------------------
# 13. Distribution of IRIS categories by profile
# ------------------------------------------------------------

cluster_risk_distribution <- block2_clustered |>
  count(profile, risk, name = "n") |>
  group_by(profile) |>
  mutate(percent = n / sum(n) * 100) |>
  ungroup()

write_csv(
  cluster_risk_distribution,
  file.path(base_dir, "tables/block2_cluster_risk_distribution.csv")
)

fig_cluster_risk_distribution <- cluster_risk_distribution |>
  ggplot(aes(x = profile, y = percent, fill = risk)) +
  geom_col(position = "fill", width = 0.6, color = "black", linewidth = 0.3) +
  scale_y_continuous(labels = percent_format(), expand = expansion(mult = c(0, 0))) +
  scale_fill_manual(values = iris_colors) + # Official SCImago palette
  labs(
    title = "IRIS risk categories across institutional profiles",
    x = NULL, y = "Proportion of institutions",
    fill = "Risk category",
    caption = "Profiles were estimated using k-means clustering on IRIS risk indicators."
  ) +
  theme_paper() +
  theme(legend.position = "bottom")


ggsave(
  filename = file.path(base_dir, "figures/block2_risk_categories_by_profile.png"),
  plot = fig_cluster_risk_distribution,
  width = 8, height = 6, dpi = 300
)

# ============================================================
# BLOCK 3: National context and scientific capacity
# ============================================================

)


}


# ------------------------------------------------------------
# Output directories and Q1 Graphical Theme
# ------------------------------------------------------------




# ------------------------------------------------------------
# 1. Download World Bank national indicators
# ------------------------------------------------------------

start_year <- 2015
end_year <- as.integer(format(Sys.Date(), "%Y")) - 1

wdi_indicators <- c(
  rd_gdp = "GB.XPD.RSDV.GD.ZS",
  researchers_pm = "SP.POP.SCIE.RD.P6"
)

wdi_raw <- WDI::WDI(
  country = "all",
  indicator = wdi_indicators,
  start = start_year,
  end = end_year,
  extra = TRUE,
  cache = NULL
)

write_csv(
  wdi_raw,
  file.path(base_dir, "tables/block3_wdi_raw_download.csv")
)

# ------------------------------------------------------------
# 2. Latest available data by country and indicator
# ------------------------------------------------------------

wdi_countries <- wdi_raw |>
  filter(
    region != "Aggregates",
    !is.na(iso3c)
  )

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
    cols = c(rd_gdp, researchers_pm),
    names_to = "indicator",
    values_to = "value"
  ) |>
  filter(!is.na(value)) |>
  arrange(iso3c, indicator, desc(year)) |>
  group_by(iso3c, indicator) |>
  slice(1) |>
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
  left_join(wdi_latest_years, by = "iso3c") |>
  mutate(
    researchers_pm_log = log1p(researchers_pm),
    rd_gdp_z = as.numeric(scale(rd_gdp)),
    researchers_pm_log_z = as.numeric(scale(researchers_pm_log))
  )

write_csv(
  country_context,
  file.path(base_dir, "tables/block3_country_context_worldbank_latest.csv")
)

# ------------------------------------------------------------
# 3. IRIS + World Bank matching coverage
# ------------------------------------------------------------

country_context_iris <- country_summary |>
  left_join(country_context, by = c("country" = "iso3c"))

block3_match_summary <- tibble(
  item = c(
    "IRIS countries",
    "IRIS countries matched with World Bank metadata",
    "IRIS countries with R&D expenditure data",
    "IRIS countries with researchers per million data",
    "IRIS countries with both R&D and researchers data"
  ),
  n = c(
    n_distinct(country_summary$country),
    sum(!is.na(country_context_iris$country_name_wb)),
    sum(!is.na(country_context_iris$rd_gdp)),
    sum(!is.na(country_context_iris$researchers_pm)),
    sum(!is.na(country_context_iris$rd_gdp) & !is.na(country_context_iris$researchers_pm))
  )
)

write_csv(
  block3_match_summary,
  file.path(base_dir, "tables/block3_match_summary.csv")
)

block3_unmatched_countries <- country_context_iris |>
  filter(is.na(country_name_wb)) |>
  select(country, n_institutions, mean_overall, pct_significant) |>
  arrange(desc(n_institutions))

write_csv(
  block3_unmatched_countries,
  file.path(base_dir, "tables/block3_unmatched_iris_worldbank.csv")
)

# ------------------------------------------------------------
# 4. Country dataset for contextual analysis
# ------------------------------------------------------------

country_analysis <- country_context_iris |>
  filter(n_institutions >= 10) |>
  mutate(
    income_group_wb = factor(
      income_group_wb,
      levels = c("Low income", "Lower middle income", "Upper middle income", "High income")
    ),
    researchers_pm_log = log1p(researchers_pm)
  )

country_analysis_summary <- country_analysis |>
  summarise(
    countries_n = n(),
    institutions_n = sum(n_institutions, na.rm = TRUE),
    countries_with_income_group = sum(!is.na(income_group_wb)),
    countries_with_rd_gdp = sum(!is.na(rd_gdp)),
    countries_with_researchers = sum(!is.na(researchers_pm)),
    countries_with_both_rd_and_researchers = sum(!is.na(rd_gdp) & !is.na(researchers_pm))
  )

write_csv(
  country_analysis_summary,
  file.path(base_dir, "tables/block3_country_analysis_summary.csv")
)

fig_block3_income_overall <- country_analysis |>
  filter(!is.na(income_group_wb)) |>
  ggplot(aes(x = income_group_wb, y = mean_overall)) +
  geom_boxplot(outlier.alpha = 0, fill = "grey95", color = "black", width = 0.5) +
  geom_point(
    aes(size = n_institutions),
    alpha = 0.6,
    color = "#3C5488FF", # Scientific blue
    position = position_jitter(width = 0.15, height = 0)
  ) +
  scale_size_continuous(name = "Institutions") +
  labs(
    title = "Country-level mean IRIS Overall score by World Bank income group",
    x = NULL, y = "Mean Overall score"
  ) +
  theme_paper() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))


fig_block3_rd_overall <- country_analysis |>
  filter(!is.na(rd_gdp), !is.na(mean_overall)) |>
  ggplot(aes(x = rd_gdp, y = mean_overall)) +
  geom_point(aes(size = n_institutions), alpha = 0.65, color = "#3C5488FF") +
  geom_smooth(method = "lm", se = TRUE, color = "#E64B35FF", fill = "grey80") +
  scale_size_continuous(name = "Institutions") +
  labs(
    title = "R&D expenditure and country-level mean IRIS Overall score",
    x = "R&D expenditure (% of GDP)", y = "Mean Overall score"
  ) +
  theme_paper()


fig_block3_researchers_overall <- country_analysis |>
  filter(!is.na(researchers_pm), !is.na(mean_overall)) |>
  ggplot(aes(x = researchers_pm, y = mean_overall)) +
  geom_point(aes(size = n_institutions), alpha = 0.65, color = "#3C5488FF") +
  geom_smooth(method = "lm", se = TRUE, color = "#E64B35FF", fill = "grey80") +
  scale_x_log10(labels = comma_format()) +
  scale_size_continuous(name = "Institutions") +
  labs(
    title = "Researchers in R&D and country-level mean IRIS Overall score",
    x = "Researchers in R&D per million people, log scale", y = "Mean Overall score"
  ) +
  theme_paper()


# Complete panel
fig_block3_context <- fig_block3_income_overall / fig_block3_rd_overall / fig_block3_researchers_overall +
  plot_annotation(
    title = "Institutional integrity risk and national research system context",
    caption = "Source: SCImago IRIS dataset and World Bank World Development Indicators. Countries with ≥ 10 institutions.",
    theme = theme(
      plot.title = element_text(family = "sans", face = "bold", size = 15, margin = margin(b = 10)),
      plot.caption = element_text(family = "sans", size = 9, color = "grey50", hjust = 0, margin = margin(t = 15))
    )
  )


ggsave(
  filename = file.path(base_dir, "figures/figure_block3_national_context.png"),
  plot = fig_block3_context,
  width = 9, height = 14, dpi = 300
)

# ------------------------------------------------------------
# 8. Correlations between aggregate risk and national variables
# ------------------------------------------------------------

country_correlations <- country_analysis |>
  summarise(
    n_rd = sum(!is.na(rd_gdp) & !is.na(mean_overall)),
    spearman_rd_mean_overall = cor(rd_gdp, mean_overall, method = "spearman", use = "complete.obs"),
    pearson_rd_mean_overall = cor(rd_gdp, mean_overall, method = "pearson", use = "complete.obs"),
    
    n_researchers = sum(!is.na(researchers_pm) & !is.na(mean_overall)),
    spearman_researchers_mean_overall = cor(researchers_pm, mean_overall, method = "spearman", use = "complete.obs"),
    pearson_researchers_mean_overall = cor(researchers_pm, mean_overall, method = "pearson", use = "complete.obs"),
    
    n_rd_pctsig = sum(!is.na(rd_gdp) & !is.na(pct_significant)),
    spearman_rd_pct_significant = cor(rd_gdp, pct_significant, method = "spearman", use = "complete.obs"),
    
    n_researchers_pctsig = sum(!is.na(researchers_pm) & !is.na(pct_significant)),
    spearman_researchers_pct_significant = cor(researchers_pm, pct_significant, method = "spearman", use = "complete.obs")
  )

write_csv(
  country_correlations,
  file.path(base_dir, "tables/block3_country_correlations.csv")
)

# ------------------------------------------------------------
# 9. Country models weighted by number of institutions
# ------------------------------------------------------------

country_model_data <- country_analysis |>
  filter(
    !is.na(mean_overall), !is.na(rd_gdp), !is.na(researchers_pm), !is.na(income_group_wb)
  ) |>
  mutate(
    rd_gdp_z = as.numeric(scale(rd_gdp)),
    researchers_pm_log_z = as.numeric(scale(log1p(researchers_pm))),
    n_institutions_log_z = as.numeric(scale(log1p(n_institutions)))
  )

country_model_overall <- lm(
  mean_overall ~ rd_gdp_z + researchers_pm_log_z + income_group_wb + n_institutions_log_z,
  data = country_model_data,
  weights = n_institutions
)

country_model_overall_table <- broom::tidy(
  country_model_overall,
  conf.int = TRUE
)

write_csv(
  country_model_overall_table,
  file.path(base_dir, "tables/block3_country_weighted_lm_mean_overall.csv")
)

# ------------------------------------------------------------
# 10. Institutional dataset with national variables
# ------------------------------------------------------------

institution_context <- data_clean |>
  left_join(country_context, by = c("country" = "iso3c")) |>
  mutate(
    output_log = log1p(output),
    sir_rank_log = log1p(sir_rank),
    output_log_z = as.numeric(scale(output_log)),
    sir_rank_log_z = as.numeric(scale(sir_rank_log)),
    rd_gdp_z = as.numeric(scale(rd_gdp)),
    researchers_pm_log = log1p(researchers_pm),
    researchers_pm_log_z = as.numeric(scale(researchers_pm_log)),
    income_group_wb = factor(
      income_group_wb,
      levels = c("Low income", "Lower middle income", "Upper middle income", "High income")
    )
  )

institution_model_data <- institution_context |>
  filter(
    !is.na(overall), !is.na(country), !is.na(output_log_z), !is.na(sir_rank_log_z),
    !is.na(rd_gdp_z), !is.na(researchers_pm_log_z), !is.na(income_group_wb)
  )

model_data_summary <- tibble(
  institutions = nrow(institution_model_data),
  countries = n_distinct(institution_model_data$country)
)

write_csv(
  model_data_summary,
  file.path(base_dir, "tables/block3_multilevel_model_data_summary.csv")
)

# ------------------------------------------------------------
# 11. Multilevel model with Overall as outcome
# ------------------------------------------------------------

model_overall_mixed <- lmer(
  overall ~ output_log_z + sir_rank_log_z +
    rd_gdp_z + researchers_pm_log_z +
    income_group_wb +
    (1 | country),
  data = institution_model_data,
  REML = FALSE
)

model_overall_mixed_table <- broom.mixed::tidy(
  model_overall_mixed, effects = "fixed", conf.int = TRUE
)

write_csv(
  model_overall_mixed_table,
  file.path(base_dir, "tables/block3_multilevel_model_overall_fixed_effects.csv")
)

model_overall_mixed_random <- broom.mixed::tidy(
  model_overall_mixed, effects = "ran_pars", conf.int = TRUE
)

write_csv(
  model_overall_mixed_random,
  file.path(base_dir, "tables/block3_multilevel_model_overall_random_effects.csv")
)

model_r2 <- performance::r2(model_overall_mixed)
write_csv(
  as_tibble(model_r2),
  file.path(base_dir, "tables/block3_multilevel_model_overall_r2.csv")
)

model_icc <- performance::icc(model_overall_mixed)
write_csv(
  as_tibble(model_icc),
  file.path(base_dir, "tables/block3_multilevel_model_overall_icc.csv")
)

# ------------------------------------------------------------
# 12. Institutional profiles by income group
# ------------------------------------------------------------

profiles_context <- block2_clustered |>
  left_join(country_context, by = c("country" = "iso3c")) |>
  mutate(
    income_group_wb = factor(
      income_group_wb,
      levels = c("Low income", "Lower middle income", "Upper middle income", "High income")
    )
  )

profile_income_distribution <- profiles_context |>
  filter(!is.na(income_group_wb), !is.na(profile)) |>
  count(income_group_wb, profile, name = "n") |>
  group_by(income_group_wb) |>
  mutate(percent = n / sum(n) * 100) |>
  ungroup()

write_csv(
  profile_income_distribution,
  file.path(base_dir, "tables/block3_profile_distribution_by_income_group.csv")
)

fig_block3_profiles_income <- profile_income_distribution |>
  ggplot(aes(x = income_group_wb, y = percent, fill = profile)) +
  # Return to solid black contour (color = "black") with fine line (linewidth = 0.4)
  geom_col(position = "fill", color = "black", linewidth = 0.4, width = 0.6) +
  scale_y_continuous(labels = percent_format(), expand = expansion(mult = c(0, 0))) +
  # Super vivid, contrasting, and professional color palette (No Pastels)
  scale_fill_manual(values = c(
    "Profile 1" = "#1F77B4", # Vibrant strong blue
    "Profile 2" = "#D62728", # Solid red
    "Profile 3" = "#2CA02C", # Bright green
    "Profile 4" = "#9467BD", # Intense purple
    "Profile 5" = "#FF7F0E"  # Bright orange
  )) +
  labs(
    title = "Institutional risk profiles by World Bank income group",
    x = NULL, y = "Proportion of institutions",
    fill = "Risk profile",
    caption = "Profiles were estimated using k-means clustering on SCImago IRIS risk indicators."
  ) +
  theme_paper() +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "bottom"
  )


# Graphs are still set to be saved manually.
ggsave(
  filename = file.path(base_dir, "figures/block3_profiles_by_income_group.png"),
  plot = fig_block3_profiles_income,
  width = 8, height = 6, dpi = 300
)