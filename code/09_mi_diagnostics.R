################################################################################
# 09_mi_diagnostics.R
#
# Purpose:
#   Diagnostic checks for understanding why LP estimates using multiple-imputed
#   BAC may differ from estimates using raw measured BAC.
#
# Main question:
#   Are missing-BAC cases that get filled in by MI systematically less
#   alcohol-involved than observed raw-BAC cases?
#
# Diagnostics:
#   1. Driver/person-level drunk-share comparisons:
#        observed raw BAC vs missing raw BAC filled by MI
#      for thresholds BAC >= 0.05, 0.08, and 0.10.
#
#   2. Same comparison by:
#        one-car vs two-car crashes,
#        day vs night,
#        weekday vs weekend,
#        rural/urban,
#        state,
#        police-reported drinking, if available,
#        injury severity, if available.
#
#   3. LP crash-type count and mechanism-ratio comparisons:
#        DD, DS, SS, D, S
#        DD / DS
#        DS / SS
#        DD / SS
#        D / S
#        D / (D + S)
#        (DD + 0.5 DS) / (DD + DS + SS)
#        [D / (D + S)] / [(DD + 0.5 DS) / (DD + DS + SS)]
#
# Inputs:
#   data/processed/driver_miper_long_2015_2024.csv
#
#   data/processed/script07_mi_lp_mle/{state,rur_urb,location_weekend}/
#     raw_crash_type_counts.csv
#     mi_crash_type_counts_by_imputation.csv
#
# Outputs:
#   data/processed/script09_mi_diagnostics/...
#   figures/script09_mi_diagnostics/...
################################################################################

rm(list = ls())

library(tidyverse)
library(scales)

# ------------------------------------------------------------------------------
# 0. User options
# ------------------------------------------------------------------------------

mi_driver_path <- "data/processed/driver_miper_long_2015_2024.csv"

script07_root <- "data/processed/script07_mi_lp_mle"

out_dir <- "data/processed/script09_mi_diagnostics"
fig_dir <- "figures/script09_mi_diagnostics"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Script 07 location-spec folders to use.
location_spec_folders <- c("location_weekend", "rur_urb", "state")

# Main comparison:
#   B = raw BAC, 95%-complete state-years
#   D = MI BAC, all state-years
raw_sample_main <- "raw_bac_complete_95_state_years"
mi_sample_main  <- "mi_all_state_years"

# Also include C vs B in the LP crash-type diagnostic because C uses MI but keeps
# the same 95%-complete state-years as B.
mi_sample_same_95 <- "mi_same_95_state_years"

# Preferred sparse-cell filter from script 07.
main_min_n2 <- 5

# Preferred compact location specification for the most useful display tables
# and figures.
preferred_location_spec <- "LOCATION_WEEKEND"
preferred_cell_spec <- "rur_urb_state_weekend"

# ------------------------------------------------------------------------------
# 1. Helper functions
# ------------------------------------------------------------------------------

safe_divide <- function(num, den) {
  ifelse(is.na(den) | den == 0, NA_real_, num / den)
}

standardize_location_spec <- function(x) {
  case_when(
    str_to_lower(x) == "location_weekend" ~ "LOCATION_WEEKEND",
    str_to_lower(x) == "rur_urb" ~ "RUR_URB",
    str_to_lower(x) == "state" ~ "STATE",
    TRUE ~ as.character(x)
  )
}

first_existing_col <- function(df, candidates) {
  out <- candidates[candidates %in% names(df)]
  if (length(out) == 0) return(NA_character_)
  out[1]
}

make_optional_character_col <- function(df, source_col, default_value = "not available") {
  if (!is.na(source_col) && source_col %in% names(df)) {
    as.character(df[[source_col]])
  } else {
    rep(default_value, nrow(df))
  }
}

# ------------------------------------------------------------------------------
# 2. Load driver-level MI file
# ------------------------------------------------------------------------------

if (!file.exists(mi_driver_path)) {
  stop(
    "Missing input file: ", mi_driver_path,
    "\nPlease run 06_prep_mi_data.R first."
  )
}

driver_mi <- read_csv(mi_driver_path, show_col_types = FALSE) %>%
  rename_with(toupper)

cat("Loaded driver-level MI file:\n")
cat("  Rows:", comma(nrow(driver_mi)), "\n")
cat("  Distinct imputations:", paste(sort(unique(driver_mi$IMPUTATION)), collapse = ", "), "\n\n")

required_cols <- c(
  "YEAR_FILE", "STATE", "ST_CASE", "VEH_NO", "PER_NO",
  "VE_FORMS", "HOUR", "IMPUTATION", "RUR_URB",
  "BAC_GE_05", "BAC_GE_08", "BAC_GE_10",
  "MI_BAC_GE_05", "MI_BAC_GE_08", "MI_BAC_GE_10"
)

missing_required <- setdiff(required_cols, names(driver_mi))

if (length(missing_required) > 0) {
  stop(
    "driver_mi is missing required columns: ",
    paste(missing_required, collapse = ", ")
  )
}

threshold_specs <- tibble(
  threshold = c(0.05, 0.08, 0.10),
  threshold_label = c("BAC >= 0.05", "BAC >= 0.08", "BAC >= 0.10"),
  raw_definition_var = c("BAC_GE_05", "BAC_GE_08", "BAC_GE_10"),
  mi_definition_var  = c("MI_BAC_GE_05", "MI_BAC_GE_08", "MI_BAC_GE_10")
)

# Optional columns. These may or may not exist depending on how script 06 created
# the driver-level file.
police_drinking_col <- first_existing_col(
  driver_mi,
  c(
    "DRINKING",
    "DRINKING_DRIVER",
    "DRINKING_REPORTED",
    "POLICE_REPORTED_DRINKING",
    "DRINKING_INV"
  )
)

injury_col <- first_existing_col(
  driver_mi,
  c(
    "INJ_SEV",
    "INJURY_SEVERITY",
    "MAX_SEV",
    "PER_TYP_INJ_SEV"
  )
)

cat("Optional diagnostic columns found:\n")
cat("  police_drinking_col:", ifelse(is.na(police_drinking_col), "none", police_drinking_col), "\n")
cat("  injury_col:", ifelse(is.na(injury_col), "none", injury_col), "\n\n")

# ------------------------------------------------------------------------------
# 3. Create driver/person-level diagnostic variables
# ------------------------------------------------------------------------------

driver_mi <- driver_mi %>%
  mutate(
    PERSON_ID_09 = paste(YEAR_FILE, STATE, ST_CASE, VEH_NO, PER_NO, sep = "_"),
    CRASH_ID_09 = paste(YEAR_FILE, STATE, ST_CASE, sep = "_"),
    
    CRASH_SIZE_09 = case_when(
      VE_FORMS == 1 ~ "one-car",
      VE_FORMS == 2 ~ "two-car",
      TRUE ~ "other"
    ),
    
    DAY_NIGHT_09 = case_when(
      is.na(HOUR) ~ "missing hour",
      HOUR >= 18 | HOUR <= 5 ~ "night: 6pm-5:59am",
      HOUR >= 6 & HOUR <= 17 ~ "day: 6am-5:59pm",
      TRUE ~ "missing/invalid hour"
    ),
    
    WEEKEND_09 = case_when(
      "WEEKEND_LP" %in% names(.) & as.logical(WEEKEND_LP) ~ "weekend",
      "WEEKEND_LP" %in% names(.) & !as.logical(WEEKEND_LP) ~ "weekday",
      "DAY_WEEK" %in% names(.) & DAY_WEEK %in% c(1, 7) ~ "weekend",
      "DAY_WEEK" %in% names(.) & !is.na(DAY_WEEK) ~ "weekday",
      TRUE ~ "not available"
    ),
    
    RUR_URB_GROUP_09 = case_when(
      RUR_URB == 1 ~ "rural",
      RUR_URB == 2 ~ "urban",
      RUR_URB %in% c(6, 8, 9) ~ "other/unknown",
      TRUE ~ "missing"
    ),
    
    STATE_09 = as.character(STATE),
    ALL_09 = "all drivers"
  ) %>%
  filter(VE_FORMS %in% c(1, 2))

# Add optional columns after mutate so this is robust if those source columns
# are not present.
driver_mi$POLICE_REPORTED_DRINKING_09 <- make_optional_character_col(
  driver_mi,
  police_drinking_col
)

driver_mi$INJURY_SEVERITY_09 <- make_optional_character_col(
  driver_mi,
  injury_col
)

# Raw measured-BAC variables are duplicated across imputations. Keep one row per
# person/crash for observed raw BAC summaries.
driver_raw_one_row <- driver_mi %>%
  filter(IMPUTATION == min(IMPUTATION, na.rm = TRUE))

cat("Rows used for driver-level diagnostics:\n")
cat("  Raw one-row copy:", comma(nrow(driver_raw_one_row)), "\n")
cat("  MI long rows:", comma(nrow(driver_mi)), "\n\n")

# ------------------------------------------------------------------------------
# 4. Driver-level observed raw BAC vs missing-BAC / MI-filled diagnostics
# ------------------------------------------------------------------------------

diagnostic_group_vars <- tribble(
  ~group_var,                         ~group_label,
  "ALL_09",                           "all drivers",
  "CRASH_SIZE_09",                    "one-car vs two-car",
  "DAY_NIGHT_09",                     "day vs night",
  "WEEKEND_09",                       "weekday vs weekend",
  "RUR_URB_GROUP_09",                 "rural/urban",
  "STATE_09",                         "state",
  "POLICE_REPORTED_DRINKING_09",      "police-reported drinking",
  "INJURY_SEVERITY_09",               "injury severity"
)

summarize_observed_vs_imputed_missing <- function(raw_one_row,
                                                  mi_long,
                                                  raw_var,
                                                  mi_var,
                                                  threshold,
                                                  threshold_label,
                                                  group_var,
                                                  group_label) {
  
  observed_raw <- raw_one_row %>%
    filter(!is.na(.data[[raw_var]])) %>%
    mutate(
      dimension = group_label,
      dimension_level = as.character(.data[[group_var]]),
      bac_source_group = "Observed raw BAC"
    ) %>%
    group_by(dimension, dimension_level, bac_source_group) %>%
    summarize(
      threshold = threshold,
      threshold_label = threshold_label,
      n_driver_rows = n(),
      n_distinct_drivers = n_distinct(PERSON_ID_09),
      n_distinct_crashes = n_distinct(CRASH_ID_09),
      drunk_share = mean(.data[[raw_var]] == 1, na.rm = TRUE),
      drunk_share_min_across_imputations = NA_real_,
      drunk_share_max_across_imputations = NA_real_,
      .groups = "drop"
    )
  
  imputed_missing_by_imp <- mi_long %>%
    filter(is.na(.data[[raw_var]]), !is.na(.data[[mi_var]])) %>%
    mutate(
      dimension = group_label,
      dimension_level = as.character(.data[[group_var]]),
      bac_source_group = "Missing raw BAC, MI-filled"
    ) %>%
    group_by(imputation = IMPUTATION, dimension, dimension_level, bac_source_group) %>%
    summarize(
      n_driver_rows_imp = n(),
      n_distinct_drivers_imp = n_distinct(PERSON_ID_09),
      n_distinct_crashes_imp = n_distinct(CRASH_ID_09),
      drunk_share_imp = mean(.data[[mi_var]] == 1, na.rm = TRUE),
      .groups = "drop"
    )
  
  imputed_missing <- imputed_missing_by_imp %>%
    group_by(dimension, dimension_level, bac_source_group) %>%
    summarize(
      threshold = threshold,
      threshold_label = threshold_label,
      n_driver_rows = mean(n_driver_rows_imp, na.rm = TRUE),
      n_distinct_drivers = mean(n_distinct_drivers_imp, na.rm = TRUE),
      n_distinct_crashes = mean(n_distinct_crashes_imp, na.rm = TRUE),
      drunk_share = mean(drunk_share_imp, na.rm = TRUE),
      drunk_share_min_across_imputations = min(drunk_share_imp, na.rm = TRUE),
      drunk_share_max_across_imputations = max(drunk_share_imp, na.rm = TRUE),
      .groups = "drop"
    )
  
  bind_rows(observed_raw, imputed_missing) %>%
    mutate(
      raw_definition_var = raw_var,
      mi_definition_var = mi_var,
      .before = 1
    )
}

driver_bac_group_diagnostics_long <- map_dfr(seq_len(nrow(threshold_specs)), function(t) {
  map_dfr(seq_len(nrow(diagnostic_group_vars)), function(g) {
    summarize_observed_vs_imputed_missing(
      raw_one_row = driver_raw_one_row,
      mi_long = driver_mi,
      raw_var = threshold_specs$raw_definition_var[t],
      mi_var = threshold_specs$mi_definition_var[t],
      threshold = threshold_specs$threshold[t],
      threshold_label = threshold_specs$threshold_label[t],
      group_var = diagnostic_group_vars$group_var[g],
      group_label = diagnostic_group_vars$group_label[g]
    )
  })
}) %>%
  mutate(
    drunk_share_pct = 100 * drunk_share,
    drunk_share_min_pct = 100 * drunk_share_min_across_imputations,
    drunk_share_max_pct = 100 * drunk_share_max_across_imputations
  ) %>%
  arrange(dimension, dimension_level, threshold, bac_source_group)

driver_bac_group_comparison <- driver_bac_group_diagnostics_long %>%
  select(
    dimension, dimension_level, threshold, threshold_label, bac_source_group,
    n_driver_rows, n_distinct_drivers, n_distinct_crashes, drunk_share
  ) %>%
  pivot_wider(
    names_from = bac_source_group,
    values_from = c(n_driver_rows, n_distinct_drivers, n_distinct_crashes, drunk_share),
    names_glue = "{.value}_{bac_source_group}"
  ) %>%
  rename_with(~ str_replace_all(.x, "[ ,\\-]+", "_")) %>%
  mutate(
    drunk_share_gap_imputed_missing_minus_observed =
      drunk_share_Missing_raw_BAC_MI_filled - drunk_share_Observed_raw_BAC,
    drunk_share_ratio_imputed_missing_to_observed =
      safe_divide(
        drunk_share_Missing_raw_BAC_MI_filled,
        drunk_share_Observed_raw_BAC
      )
  ) %>%
  arrange(dimension, dimension_level, threshold)

driver_bac_group_summary <- driver_bac_group_comparison %>%
  filter(dimension == "all drivers")

driver_bac_group_summary_display <- driver_bac_group_summary %>%
  transmute(
    threshold_label,
    observed_raw_n = round(n_driver_rows_Observed_raw_BAC, 1),
    imputed_missing_n_avg_across_imputations =
      round(n_driver_rows_Missing_raw_BAC_MI_filled, 1),
    observed_raw_drunk_share =
      percent(drunk_share_Observed_raw_BAC, accuracy = 0.1),
    imputed_missing_drunk_share =
      percent(drunk_share_Missing_raw_BAC_MI_filled, accuracy = 0.1),
    imputed_minus_observed_percentage_points =
      round(100 * drunk_share_gap_imputed_missing_minus_observed, 2),
    imputed_to_observed_ratio =
      round(drunk_share_ratio_imputed_missing_to_observed, 3)
  )

write_csv(
  driver_bac_group_diagnostics_long,
  file.path(out_dir, "diagnostic_driver_bac_group_by_covariates_long.csv")
)

write_csv(
  driver_bac_group_comparison,
  file.path(out_dir, "diagnostic_driver_bac_group_by_covariates.csv")
)

write_csv(
  driver_bac_group_summary,
  file.path(out_dir, "diagnostic_driver_bac_group_summary.csv")
)

write_csv(
  driver_bac_group_summary_display,
  file.path(out_dir, "diagnostic_driver_bac_group_summary_display.csv")
)

# ------------------------------------------------------------------------------
# 5. Driver-level figures
# ------------------------------------------------------------------------------

plot_driver_summary <- driver_bac_group_diagnostics_long %>%
  filter(dimension == "all drivers") %>%
  mutate(
    bac_source_group = factor(
      bac_source_group,
      levels = c("Observed raw BAC", "Missing raw BAC, MI-filled")
    ),
    threshold_label = factor(
      threshold_label,
      levels = threshold_specs$threshold_label
    )
  )

p_observed_vs_imputed <- ggplot(
  plot_driver_summary,
  aes(
    x = threshold_label,
    y = drunk_share,
    group = bac_source_group,
    shape = bac_source_group,
    linetype = bac_source_group
  )
) +
  geom_line() +
  geom_point(size = 2.5) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "Observed raw BAC vs imputed-missing BAC drunk share",
    subtitle = "Observed group uses non-missing raw BAC; imputed-missing group averages MI indicators across 10 imputations",
    x = "BAC threshold",
    y = "Share above threshold",
    shape = "BAC group",
    linetype = "BAC group",
    caption = "If the imputed-missing group has a lower drunk share, MI can mechanically add more sober-classified cases."
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(fig_dir, "observed_vs_imputed_missing_drunk_share_by_threshold.png"),
  p_observed_vs_imputed,
  width = 9,
  height = 5,
  dpi = 300
)

plot_driver_crash_size <- driver_bac_group_diagnostics_long %>%
  filter(dimension == "one-car vs two-car") %>%
  mutate(
    bac_source_group = factor(
      bac_source_group,
      levels = c("Observed raw BAC", "Missing raw BAC, MI-filled")
    ),
    threshold_label = factor(
      threshold_label,
      levels = threshold_specs$threshold_label
    )
  )

p_observed_vs_imputed_crash_size <- ggplot(
  plot_driver_crash_size,
  aes(
    x = threshold_label,
    y = drunk_share,
    group = bac_source_group,
    shape = bac_source_group,
    linetype = bac_source_group
  )
) +
  geom_line() +
  geom_point(size = 2.5) +
  facet_wrap(~ dimension_level) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "Observed raw BAC vs imputed-missing BAC drunk share, by crash size",
    subtitle = "This checks whether MI adds relatively more sober one-car or two-car cases",
    x = "BAC threshold",
    y = "Share above threshold",
    shape = "BAC group",
    linetype = "BAC group"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(fig_dir, "observed_vs_imputed_missing_drunk_share_by_crash_size.png"),
  p_observed_vs_imputed_crash_size,
  width = 10,
  height = 5.5,
  dpi = 300
)

# ------------------------------------------------------------------------------
# 6. LP crash-type count and mechanism-ratio diagnostics
# ------------------------------------------------------------------------------

read_script07_crash_counts_one_folder <- function(folder_slug) {
  folder <- file.path(script07_root, folder_slug)
  
  raw_path <- file.path(folder, "raw_crash_type_counts.csv")
  mi_path  <- file.path(folder, "mi_crash_type_counts_by_imputation.csv")
  
  if (!file.exists(raw_path) || !file.exists(mi_path)) {
    warning("Skipping crash-count diagnostics for ", folder_slug, " because expected files are missing.")
    return(list(raw = tibble(), mi = tibble()))
  }
  
  raw_counts <- read_csv(raw_path, show_col_types = FALSE) %>%
    mutate(source_folder = folder_slug)
  
  mi_counts <- read_csv(mi_path, show_col_types = FALSE) %>%
    mutate(source_folder = folder_slug)
  
  list(raw = raw_counts, mi = mi_counts)
}

add_lp_mechanism_ratios <- function(df) {
  df %>%
    mutate(
      two_car_complete = DD + DS + SS,
      one_car_complete = D + S,
      
      ratio_DD_to_DS = safe_divide(DD, DS),
      ratio_DS_to_SS = safe_divide(DS, SS),
      ratio_DD_to_SS = safe_divide(DD, SS),
      ratio_D_to_S = safe_divide(D, S),
      
      one_car_drunk_share = safe_divide(D, D + S),
      two_car_drunk_involvement_share = safe_divide(
        DD + 0.5 * DS,
        DD + DS + SS
      ),
      one_car_to_two_car_drunk_involvement_ratio = safe_divide(
        one_car_drunk_share,
        two_car_drunk_involvement_share
      )
    )
}

all_crash_count_reads <- map(location_spec_folders, read_script07_crash_counts_one_folder)

raw_crash_counts <- map_dfr(all_crash_count_reads, "raw") %>%
  mutate(location_spec = standardize_location_spec(location_spec))

mi_crash_counts <- map_dfr(all_crash_count_reads, "mi") %>%
  mutate(location_spec = standardize_location_spec(location_spec))

if (nrow(raw_crash_counts) == 0 || nrow(mi_crash_counts) == 0) {
  warning("No script 07 crash-count outputs found for LP crash-type diagnostics.")
} else {
  
  crash_count_types <- c("DD", "DS", "SS", "D", "S")
  
  ratio_types <- c(
    "ratio_DD_to_DS",
    "ratio_DS_to_SS",
    "ratio_DD_to_SS",
    "ratio_D_to_S",
    "one_car_drunk_share",
    "two_car_drunk_involvement_share",
    "one_car_to_two_car_drunk_involvement_ratio"
  )
  
  ratio_labels <- tibble(
    metric = ratio_types,
    metric_label = c(
      "DD / DS",
      "DS / SS",
      "DD / SS",
      "D / S",
      "D / (D + S)",
      "(DD + 0.5 DS) / (DD + DS + SS)",
      "[D / (D + S)] / [(DD + 0.5 DS) / (DD + DS + SS)]"
    ),
    metric_family = c(
      "two-car ratio",
      "two-car ratio",
      "two-car ratio",
      "one-car ratio",
      "one-car drunk share",
      "two-car drunk involvement share",
      "one-car vs two-car comparison"
    )
  )
  
  raw_crash_metrics <- raw_crash_counts %>%
    filter(
      analysis_sample == raw_sample_main,
      min_n2 == main_min_n2
    ) %>%
    add_lp_mechanism_ratios() %>%
    select(
      analysis_sample, analysis_label, threshold, threshold_label,
      location_spec, cell_spec, cell_spec_label, min_n2, min_n2_label,
      all_of(crash_count_types), all_of(ratio_types)
    )
  
  mi_crash_metrics_by_imp <- mi_crash_counts %>%
    filter(
      analysis_sample %in% c(mi_sample_same_95, mi_sample_main),
      min_n2 == main_min_n2
    ) %>%
    add_lp_mechanism_ratios() %>%
    select(
      analysis_sample, analysis_label, imputation, threshold, threshold_label,
      location_spec, cell_spec, cell_spec_label, min_n2, min_n2_label,
      all_of(crash_count_types), all_of(ratio_types)
    )
  
  # --------------------------------------------------------------------------
  # 6a. Count diagnostics: DD / DS / SS / D / S
  # --------------------------------------------------------------------------
  
  raw_crash_counts_long <- raw_crash_metrics %>%
    select(
      analysis_sample, analysis_label, threshold, threshold_label,
      location_spec, cell_spec, cell_spec_label, min_n2, min_n2_label,
      all_of(crash_count_types)
    ) %>%
    pivot_longer(
      cols = all_of(crash_count_types),
      names_to = "metric",
      values_to = "raw_value"
    ) %>%
    mutate(
      metric_label = metric,
      metric_family = "LP crash-type count"
    )
  
  mi_crash_counts_avg_long <- mi_crash_metrics_by_imp %>%
    group_by(
      analysis_sample, analysis_label, threshold, threshold_label,
      location_spec, cell_spec, cell_spec_label, min_n2, min_n2_label
    ) %>%
    summarize(
      across(all_of(crash_count_types), ~ mean(.x, na.rm = TRUE)),
      n_imputations = n_distinct(imputation),
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols = all_of(crash_count_types),
      names_to = "metric",
      values_to = "mi_value_avg"
    ) %>%
    mutate(
      metric_label = metric,
      metric_family = "LP crash-type count"
    ) %>%
    rename(
      mi_analysis_sample = analysis_sample,
      mi_analysis_label = analysis_label
    )
  
  lp_crash_count_diagnostics <- mi_crash_counts_avg_long %>%
    inner_join(
      raw_crash_counts_long,
      by = c(
        "threshold", "threshold_label", "location_spec", "cell_spec",
        "cell_spec_label", "min_n2", "min_n2_label",
        "metric", "metric_label", "metric_family"
      )
    ) %>%
    mutate(
      comparison_label = paste0(mi_analysis_label, " minus ", analysis_label),
      difference_mi_minus_raw = mi_value_avg - raw_value,
      pct_difference_mi_minus_raw = ifelse(
        raw_value > 0,
        difference_mi_minus_raw / raw_value,
        NA_real_
      )
    )
  
  # --------------------------------------------------------------------------
  # 6b. Ratio diagnostics
  # --------------------------------------------------------------------------
  
  raw_ratio_long <- raw_crash_metrics %>%
    select(
      analysis_sample, analysis_label, threshold, threshold_label,
      location_spec, cell_spec, cell_spec_label, min_n2, min_n2_label,
      all_of(ratio_types)
    ) %>%
    pivot_longer(
      cols = all_of(ratio_types),
      names_to = "metric",
      values_to = "raw_value"
    ) %>%
    left_join(ratio_labels, by = "metric")
  
  mi_ratio_avg_long <- mi_crash_metrics_by_imp %>%
    group_by(
      analysis_sample, analysis_label, threshold, threshold_label,
      location_spec, cell_spec, cell_spec_label, min_n2, min_n2_label
    ) %>%
    summarize(
      across(all_of(ratio_types), ~ mean(.x, na.rm = TRUE)),
      n_imputations = n_distinct(imputation),
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols = all_of(ratio_types),
      names_to = "metric",
      values_to = "mi_value_avg"
    ) %>%
    left_join(ratio_labels, by = "metric") %>%
    rename(
      mi_analysis_sample = analysis_sample,
      mi_analysis_label = analysis_label
    )
  
  lp_ratio_diagnostics <- mi_ratio_avg_long %>%
    inner_join(
      raw_ratio_long,
      by = c(
        "threshold", "threshold_label", "location_spec", "cell_spec",
        "cell_spec_label", "min_n2", "min_n2_label",
        "metric", "metric_label", "metric_family"
      )
    ) %>%
    mutate(
      comparison_label = paste0(mi_analysis_label, " minus ", analysis_label),
      difference_mi_minus_raw = mi_value_avg - raw_value,
      pct_difference_mi_minus_raw = ifelse(
        is.finite(raw_value) & raw_value != 0,
        difference_mi_minus_raw / raw_value,
        NA_real_
      )
    )
  
  # --------------------------------------------------------------------------
  # 6c. Combined mechanism table
  # --------------------------------------------------------------------------
  
  lp_mechanism_diagnostics <- bind_rows(
    lp_crash_count_diagnostics,
    lp_ratio_diagnostics
  ) %>%
    mutate(
      direction_interpretation = case_when(
        metric %in% c("SS", "S") & difference_mi_minus_raw > 0 ~
          "MI adds more sober-classified crashes",
        
        metric %in% c("DD", "D") & difference_mi_minus_raw > 0 ~
          "MI adds more drunk-classified crashes",
        
        metric == "DS" & difference_mi_minus_raw > 0 ~
          "MI adds more mixed drunk/sober two-car crashes",
        
        metric %in% c("ratio_D_to_S", "one_car_drunk_share") &
          difference_mi_minus_raw < 0 ~
          "MI lowers one-car drunk involvement",
        
        metric %in% c("ratio_DS_to_SS", "ratio_DD_to_SS") &
          difference_mi_minus_raw < 0 ~
          "MI lowers two-car drunk involvement relative to sober two-car crashes",
        
        metric == "one_car_to_two_car_drunk_involvement_ratio" &
          difference_mi_minus_raw < 0 ~
          "MI lowers one-car drunk involvement relative to two-car drunk involvement",
        
        TRUE ~ "check context"
      )
    ) %>%
    select(
      comparison_label,
      raw_analysis_sample = analysis_sample,
      raw_analysis_label = analysis_label,
      mi_analysis_sample,
      mi_analysis_label,
      threshold, threshold_label,
      location_spec, cell_spec, cell_spec_label,
      min_n2, min_n2_label,
      metric, metric_label, metric_family,
      raw_value, mi_value_avg,
      difference_mi_minus_raw,
      pct_difference_mi_minus_raw,
      n_imputations,
      direction_interpretation
    ) %>%
    arrange(
      location_spec,
      cell_spec,
      threshold,
      mi_analysis_sample,
      factor(
        metric,
        levels = c(
          "DD", "DS", "SS", "D", "S",
          "ratio_DD_to_DS",
          "ratio_DS_to_SS",
          "ratio_DD_to_SS",
          "ratio_D_to_S",
          "one_car_drunk_share",
          "two_car_drunk_involvement_share",
          "one_car_to_two_car_drunk_involvement_ratio"
        )
      )
    )
  
  lp_mechanism_diagnostics_display <- lp_mechanism_diagnostics %>%
    mutate(
      raw_value = round(raw_value, 4),
      mi_value_avg = round(mi_value_avg, 4),
      difference_mi_minus_raw = round(difference_mi_minus_raw, 4),
      pct_difference_mi_minus_raw = percent(pct_difference_mi_minus_raw, accuracy = 0.1)
    )
  
  write_csv(
    lp_mechanism_diagnostics,
    file.path(out_dir, "diagnostic_lp_mechanism_counts_and_ratios_raw_vs_mi.csv")
  )
  
  write_csv(
    lp_mechanism_diagnostics_display,
    file.path(out_dir, "diagnostic_lp_mechanism_counts_and_ratios_raw_vs_mi_display.csv")
  )
  
  # Backward-compatible count-only files.
  write_csv(
    lp_crash_count_diagnostics,
    file.path(out_dir, "diagnostic_lp_crash_counts_raw_vs_mi.csv")
  )
  
  write_csv(
    lp_crash_count_diagnostics %>%
      mutate(
        raw_value = round(raw_value, 1),
        mi_value_avg = round(mi_value_avg, 1),
        difference_mi_minus_raw = round(difference_mi_minus_raw, 1),
        pct_difference_mi_minus_raw = percent(pct_difference_mi_minus_raw, accuracy = 0.1)
      ),
    file.path(out_dir, "diagnostic_lp_crash_counts_raw_vs_mi_display.csv")
  )
  
  # Compact preferred-spec tables.
  preferred_lp_mechanism_diagnostics <- lp_mechanism_diagnostics %>%
    filter(
      location_spec == preferred_location_spec,
      cell_spec == preferred_cell_spec,
      mi_analysis_sample == mi_sample_main
    )
  
  preferred_lp_mechanism_diagnostics_display <- preferred_lp_mechanism_diagnostics %>%
    mutate(
      raw_value = round(raw_value, 4),
      mi_value_avg = round(mi_value_avg, 4),
      difference_mi_minus_raw = round(difference_mi_minus_raw, 4),
      pct_difference_mi_minus_raw = percent(pct_difference_mi_minus_raw, accuracy = 0.1)
    )
  
  write_csv(
    preferred_lp_mechanism_diagnostics,
    file.path(out_dir, "diagnostic_preferred_lp_mechanism_counts_and_ratios_raw_vs_mi.csv")
  )
  
  write_csv(
    preferred_lp_mechanism_diagnostics_display,
    file.path(out_dir, "diagnostic_preferred_lp_mechanism_counts_and_ratios_raw_vs_mi_display.csv")
  )
  
  # Wide preferred-spec table that is easier to read quickly.
  preferred_lp_mechanism_wide_display <- preferred_lp_mechanism_diagnostics %>%
    select(
      threshold_label,
      metric_label,
      metric_family,
      raw_value,
      mi_value_avg,
      difference_mi_minus_raw,
      pct_difference_mi_minus_raw,
      direction_interpretation
    ) %>%
    mutate(
      raw_value = round(raw_value, 4),
      mi_value_avg = round(mi_value_avg, 4),
      difference_mi_minus_raw = round(difference_mi_minus_raw, 4),
      pct_difference_mi_minus_raw = percent(pct_difference_mi_minus_raw, accuracy = 0.1)
    ) %>%
    arrange(
      threshold_label,
      factor(
        metric_label,
        levels = c(
          "DD", "DS", "SS", "D", "S",
          "DD / DS",
          "DS / SS",
          "DD / SS",
          "D / S",
          "D / (D + S)",
          "(DD + 0.5 DS) / (DD + DS + SS)",
          "[D / (D + S)] / [(DD + 0.5 DS) / (DD + DS + SS)]"
        )
      )
    )
  
  write_csv(
    preferred_lp_mechanism_wide_display,
    file.path(out_dir, "diagnostic_preferred_lp_mechanism_wide_display.csv")
  )
  
  # --------------------------------------------------------------------------
  # 6d. LP crash-type count figure
  # --------------------------------------------------------------------------
  
  lp_crash_plot_df <- preferred_lp_mechanism_diagnostics %>%
    filter(metric %in% crash_count_types) %>%
    mutate(
      metric = factor(metric, levels = crash_count_types),
      threshold_label = factor(threshold_label, levels = threshold_specs$threshold_label)
    )
  
  if (nrow(lp_crash_plot_df) > 0) {
    p_lp_crash_counts <- ggplot(
      lp_crash_plot_df,
      aes(
        x = metric,
        y = difference_mi_minus_raw
      )
    ) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_col() +
      facet_wrap(~ threshold_label, scales = "free_y") +
      labs(
        title = "LP crash-type count shift: MI all-state-years minus raw 95%-complete",
        subtitle = "Preferred LOCATION_WEEKEND specification: rural/urban x state x weekend, n2 >= 5",
        x = "LP crash type",
        y = "Average MI count minus raw count",
        caption = "Positive S or SS shifts mean MI adds relatively more sober-classified crashes; positive D or DD shifts mean more drunk-classified crashes."
      ) +
      theme_minimal(base_size = 12) +
      theme(panel.grid.minor = element_blank())
    
    ggsave(
      file.path(fig_dir, "lp_crash_type_counts_raw_vs_mi.png"),
      p_lp_crash_counts,
      width = 10,
      height = 5.5,
      dpi = 300
    )
  }
  
  # --------------------------------------------------------------------------
  # 6e. LP mechanism-ratio figure
  # --------------------------------------------------------------------------
  
  lp_ratio_plot_df <- preferred_lp_mechanism_diagnostics %>%
    filter(metric %in% ratio_types) %>%
    mutate(
      metric_label = factor(
        metric_label,
        levels = ratio_labels$metric_label
      ),
      threshold_label = factor(threshold_label, levels = threshold_specs$threshold_label)
    )
  
  if (nrow(lp_ratio_plot_df) > 0) {
    p_lp_ratio_shifts <- ggplot(
      lp_ratio_plot_df,
      aes(
        x = difference_mi_minus_raw,
        y = metric_label
      )
    ) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      geom_point(size = 2.5) +
      facet_wrap(~ threshold_label, scales = "free_x") +
      labs(
        title = "LP mechanism-ratio shift: MI all-state-years minus raw 95%-complete",
        subtitle = "Preferred LOCATION_WEEKEND specification: rural/urban x state x weekend, n2 >= 5",
        x = "Average MI ratio minus raw ratio",
        y = NULL,
        caption = "Negative D/S or D/(D+S) shifts suggest MI lowers one-car drunk involvement; negative DS/SS or DD/SS shifts suggest MI adds relatively more sober two-car crashes."
      ) +
      theme_minimal(base_size = 12) +
      theme(panel.grid.minor = element_blank())
    
    ggsave(
      file.path(fig_dir, "lp_mechanism_ratio_shifts_raw_vs_mi.png"),
      p_lp_ratio_shifts,
      width = 11,
      height = 6,
      dpi = 300
    )
  }
}
