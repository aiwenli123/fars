################################################################################
# 08_regression_analysis.R
#
# Purpose:
#   Exploratory regression analysis of the cell-level LP implied exposure ratio
#   produced by script 07.
#
#   Main outcomes:
#     1. log(N_hat_c), where N_hat_c is the LP-implied drinking-to-sober
#        driver exposure ratio in cell c.
#     2. Delta_c = log(N_hat_{c, MI}) - log(N_hat_{c, raw}), comparing
#        Sample D to Sample B for matched cells.
#
#   Main samples:
#     B. Raw BAC, 95% complete state-years
#     D. MI BAC, all state-years
#
#   Main restriction:
#     n2 >= 5 only, matching the preferred sparse-cell filter from script 07.
#
# Notes:
#   - This script runs the regressions separately within each location-spec family
#     rather than pooling STATE, RUR_URB, and LOCATION_WEEKEND cells together.
#   - The Delta_c comparison here is B vs. D, so it reflects both the MI BAC
#     replacement and the change from restricted state-years to all state-years.
#     If you later want a cleaner MI-only comparison within the same state-years,
#     change mi_sample_for_delta from "mi_all_state_years" to
#     "mi_same_95_state_years".
#
# Inputs:
#   data/processed/script07_mi_lp_mle/{state,rur_urb,location_weekend}/
#     raw_cell_N_hat_values.csv
#     mi_cell_N_hat_values_by_imputation.csv
#
# Outputs:
#   data/processed/script08_regression_analysis/...
#   figures/script08_regression_analysis/...
################################################################################

rm(list = ls())

library(tidyverse)
library(scales)

# ------------------------------------------------------------------------------
# 0. User options
# ------------------------------------------------------------------------------

script07_root <- "data/processed/script07_mi_lp_mle"
out_dir <- "data/processed/script08_regression_analysis"
fig_dir <- "figures/script08_regression_analysis"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Script 07 location-spec folders to analyze. These are deliberately analyzed
# separately because their cell definitions are not comparable observations.
location_spec_folders <- c("location_weekend", "rur_urb", "state")

# Main samples requested.
raw_sample_main <- "raw_bac_complete_95_state_years"  # Sample B
mi_sample_main  <- "mi_all_state_years"               # Sample D

# Main sparse-cell filter requested.
main_min_n2 <- 5

# Which MI sample to use for the Delta_c comparison.
# Main requested version is D - B.
mi_sample_for_delta <- mi_sample_main
raw_sample_for_delta <- raw_sample_main

# Distribution check heuristic. These are only diagnostic flags, not formal tests.
skew_flag_cutoff <- 1
mean_median_ratio_cutoff <- 2

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

skewness_simple <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3 || sd(x) == 0) return(NA_real_)
  mean(((x - mean(x)) / sd(x))^3)
}

coef_table_from_lm <- function(fit) {
  coef_mat <- summary(fit)$coefficients
  tibble(
    term = rownames(coef_mat),
    estimate = coef_mat[, "Estimate"],
    std_error = coef_mat[, "Std. Error"],
    statistic = coef_mat[, "t value"],
    p_value = coef_mat[, "Pr(>|t|)"],
    n_obs = nobs(fit),
    r_squared = summary(fit)$r.squared,
    adj_r_squared = summary(fit)$adj.r.squared
  )
}

pool_lm_terms_rubin <- function(coef_list) {
  # coef_list should contain one coefficient table per imputation.
  # This pools terms that are present in at least one imputation.
  coef_df <- bind_rows(coef_list)
  if (nrow(coef_df) == 0) return(tibble())
  
  coef_df %>%
    group_by(term) %>%
    summarize(
      m_total = n_distinct(imputation),
      m_est = sum(is.finite(estimate)),
      qbar = mean(estimate, na.rm = TRUE),
      W = mean(std_error^2, na.rm = TRUE),
      B = var(estimate, na.rm = TRUE),
      avg_n_obs = mean(n_obs, na.rm = TRUE),
      avg_r_squared = mean(r_squared, na.rm = TRUE),
      avg_adj_r_squared = mean(adj_r_squared, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      B = ifelse(is.na(B), 0, B),
      total_var = W + (1 + 1 / m_total) * B,
      std_error = sqrt(total_var),
      statistic = qbar / std_error,
      # Large-sample normal approximation. This is enough for exploratory tables.
      p_value = 2 * pnorm(abs(statistic), lower.tail = FALSE),
      estimate = qbar
    ) %>%
    select(
      term, estimate, std_error, statistic, p_value,
      m_total, m_est, avg_n_obs, avg_r_squared, avg_adj_r_squared,
      W, B, total_var
    )
}

# Parse script 07 cell_id values into usable regression covariates.
# This intentionally uses factor variables because the goal is diagnostic /
# explanatory variation across cells, not causal interpretation.
add_cell_covariates <- function(df) {
  df %>%
    mutate(
      location_spec = standardize_location_spec(location_spec),
      cell_id = as.character(cell_id),
      HOUR = NA_character_,
      PERIOD = NA_character_,
      WEEKEND = NA_character_,
      RUR_URB_GROUP = NA_character_,
      STATE = NA_character_
    ) %>%
    mutate(
      # ------------------------------
      # LOCATION_WEEKEND specs
      # ------------------------------
      RUR_URB_GROUP = case_when(
        location_spec == "LOCATION_WEEKEND" & cell_spec == "rur_urb" ~ cell_id,
        location_spec == "LOCATION_WEEKEND" & cell_spec == "rur_urb_state" ~ str_match(cell_id, "^(rural|urban)_(.+)$")[, 2],
        location_spec == "LOCATION_WEEKEND" & cell_spec == "rur_urb_weekend" ~ str_match(cell_id, "^(rural|urban)_(TRUE|FALSE)$")[, 2],
        location_spec == "LOCATION_WEEKEND" & cell_spec == "rur_urb_state_weekend" ~ str_match(cell_id, "^(rural|urban)_(.+)_(TRUE|FALSE)$")[, 2],
        TRUE ~ RUR_URB_GROUP
      ),
      STATE = case_when(
        location_spec == "LOCATION_WEEKEND" & cell_spec == "state" ~ cell_id,
        location_spec == "LOCATION_WEEKEND" & cell_spec == "rur_urb_state" ~ str_match(cell_id, "^(rural|urban)_(.+)$")[, 3],
        location_spec == "LOCATION_WEEKEND" & cell_spec == "state_weekend" ~ str_match(cell_id, "^(.+)_(TRUE|FALSE)$")[, 2],
        location_spec == "LOCATION_WEEKEND" & cell_spec == "rur_urb_state_weekend" ~ str_match(cell_id, "^(rural|urban)_(.+)_(TRUE|FALSE)$")[, 3],
        TRUE ~ STATE
      ),
      WEEKEND = case_when(
        location_spec == "LOCATION_WEEKEND" & cell_spec == "rur_urb_weekend" ~ str_match(cell_id, "^(rural|urban)_(TRUE|FALSE)$")[, 3],
        location_spec == "LOCATION_WEEKEND" & cell_spec == "state_weekend" ~ str_match(cell_id, "^(.+)_(TRUE|FALSE)$")[, 3],
        location_spec == "LOCATION_WEEKEND" & cell_spec == "rur_urb_state_weekend" ~ str_match(cell_id, "^(rural|urban)_(.+)_(TRUE|FALSE)$")[, 4],
        TRUE ~ WEEKEND
      ),
      
      # ------------------------------
      # RUR_URB specs
      # ------------------------------
      HOUR = case_when(
        location_spec == "RUR_URB" & cell_spec %in% c(
          "hour", "hour_weekend", "hour_rur_urb", "hour_rur_urb_weekend",
          "hour_rur_urb_period", "hour_rur_urb_period_weekend"
        ) ~ str_match(cell_id, "^([^_]+)")[, 2],
        TRUE ~ HOUR
      ),
      RUR_URB_GROUP = case_when(
        location_spec == "RUR_URB" & cell_spec == "hour_rur_urb" ~ str_match(cell_id, "^[^_]+_(rural|urban)$")[, 2],
        location_spec == "RUR_URB" & cell_spec == "hour_rur_urb_weekend" ~ str_match(cell_id, "^[^_]+_(rural|urban)_(TRUE|FALSE)$")[, 2],
        location_spec == "RUR_URB" & cell_spec == "hour_rur_urb_period" ~ str_match(cell_id, "^[^_]+_(rural|urban)_(pre_covid|covid|post_covid)$")[, 2],
        location_spec == "RUR_URB" & cell_spec == "hour_rur_urb_period_weekend" ~ str_match(cell_id, "^[^_]+_(rural|urban)_(pre_covid|covid|post_covid)_(TRUE|FALSE)$")[, 2],
        TRUE ~ RUR_URB_GROUP
      ),
      PERIOD = case_when(
        location_spec == "RUR_URB" & cell_spec == "hour_rur_urb_period" ~ str_match(cell_id, "^[^_]+_(rural|urban)_(pre_covid|covid|post_covid)$")[, 3],
        location_spec == "RUR_URB" & cell_spec == "hour_rur_urb_period_weekend" ~ str_match(cell_id, "^[^_]+_(rural|urban)_(pre_covid|covid|post_covid)_(TRUE|FALSE)$")[, 3],
        TRUE ~ PERIOD
      ),
      WEEKEND = case_when(
        location_spec == "RUR_URB" & cell_spec == "hour_weekend" ~ str_match(cell_id, "^[^_]+_(TRUE|FALSE)$")[, 2],
        location_spec == "RUR_URB" & cell_spec == "hour_rur_urb_weekend" ~ str_match(cell_id, "^[^_]+_(rural|urban)_(TRUE|FALSE)$")[, 3],
        location_spec == "RUR_URB" & cell_spec == "hour_rur_urb_period_weekend" ~ str_match(cell_id, "^[^_]+_(rural|urban)_(pre_covid|covid|post_covid)_(TRUE|FALSE)$")[, 4],
        TRUE ~ WEEKEND
      ),
      
      # ------------------------------
      # STATE specs
      # ------------------------------
      HOUR = case_when(
        location_spec == "STATE" & cell_spec %in% c(
          "hour", "hour_period", "hour_period_weekend",
          "hour_state_period", "hour_state_period_weekend"
        ) ~ str_match(cell_id, "^([^_]+)")[, 2],
        TRUE ~ HOUR
      ),
      PERIOD = case_when(
        location_spec == "STATE" & cell_spec == "hour_period" ~ str_match(cell_id, "^[^_]+_(pre_covid|covid|post_covid)$")[, 2],
        location_spec == "STATE" & cell_spec == "hour_period_weekend" ~ str_match(cell_id, "^[^_]+_(pre_covid|covid|post_covid)_(TRUE|FALSE)$")[, 2],
        location_spec == "STATE" & cell_spec == "hour_state_period" ~ str_match(cell_id, "^[^_]+_([^_]+)_(pre_covid|covid|post_covid)$")[, 3],
        location_spec == "STATE" & cell_spec == "hour_state_period_weekend" ~ str_match(cell_id, "^[^_]+_([^_]+)_(pre_covid|covid|post_covid)_(TRUE|FALSE)$")[, 3],
        TRUE ~ PERIOD
      ),
      STATE = case_when(
        location_spec == "STATE" & cell_spec == "hour_state_period" ~ str_match(cell_id, "^[^_]+_([^_]+)_(pre_covid|covid|post_covid)$")[, 2],
        location_spec == "STATE" & cell_spec == "hour_state_period_weekend" ~ str_match(cell_id, "^[^_]+_([^_]+)_(pre_covid|covid|post_covid)_(TRUE|FALSE)$")[, 2],
        TRUE ~ STATE
      ),
      WEEKEND = case_when(
        location_spec == "STATE" & cell_spec == "hour_period_weekend" ~ str_match(cell_id, "^[^_]+_(pre_covid|covid|post_covid)_(TRUE|FALSE)$")[, 3],
        location_spec == "STATE" & cell_spec == "hour_state_period_weekend" ~ str_match(cell_id, "^[^_]+_([^_]+)_(pre_covid|covid|post_covid)_(TRUE|FALSE)$")[, 4],
        TRUE ~ WEEKEND
      )
    ) %>%
    mutate(
      HOUR = factor(HOUR),
      PERIOD = factor(PERIOD, levels = c("pre_covid", "covid", "post_covid")),
      WEEKEND = factor(WEEKEND, levels = c("FALSE", "TRUE")),
      RUR_URB_GROUP = factor(RUR_URB_GROUP, levels = c("urban", "rural")),
      STATE = factor(STATE)
    )
}

make_formula_for_available_covariates <- function(df, outcome_var) {
  candidate_vars <- c("RUR_URB_GROUP", "WEEKEND", "PERIOD", "HOUR", "STATE")
  usable_vars <- candidate_vars[sapply(candidate_vars, function(v) {
    v %in% names(df) && n_distinct(df[[v]][!is.na(df[[v]])]) >= 2
  })]
  
  if (length(usable_vars) == 0) return(NULL)
  as.formula(paste(outcome_var, "~", paste(usable_vars, collapse = " + ")))
}

fit_one_lm <- function(df, outcome_var, weighted = FALSE) {
  formula_obj <- make_formula_for_available_covariates(df, outcome_var)
  if (is.null(formula_obj)) return(NULL)
  
  model_df <- df %>%
    filter(is.finite(.data[[outcome_var]]))
  
  # Drop rows with missing model variables.
  needed_vars <- all.vars(formula_obj)
  model_df <- model_df %>% drop_na(all_of(needed_vars))
  
  # Need enough observations after factor expansion.
  if (nrow(model_df) < 5) return(NULL)
  
  if (weighted && "n2" %in% names(model_df)) {
    model_df <- model_df %>% filter(is.finite(n2), n2 > 0)
    if (nrow(model_df) < 5) return(NULL)
    fit <- lm(formula_obj, data = model_df, weights = n2)
  } else {
    fit <- lm(formula_obj, data = model_df)
  }
  
  list(fit = fit, formula = deparse(formula_obj))
}

run_raw_regressions <- function(df, outcome_var, outcome_label) {
  split_df <- df %>%
    group_by(location_spec, cell_spec, cell_spec_label, threshold, threshold_label, min_n2, min_n2_label) %>%
    group_split()
  
  map_dfr(split_df, function(d) {
    meta <- d %>% slice(1)
    map_dfr(c(FALSE, TRUE), function(weighted) {
      fit_obj <- fit_one_lm(d, outcome_var = outcome_var, weighted = weighted)
      if (is.null(fit_obj)) return(tibble())
      coef_table_from_lm(fit_obj$fit) %>%
        mutate(
          analysis_sample = meta$analysis_sample,
          analysis_label = meta$analysis_label,
          outcome = outcome_label,
          weighted_by_n2 = weighted,
          formula = fit_obj$formula,
          location_spec = meta$location_spec,
          cell_spec = meta$cell_spec,
          cell_spec_label = meta$cell_spec_label,
          threshold = meta$threshold,
          threshold_label = meta$threshold_label,
          min_n2 = meta$min_n2,
          min_n2_label = meta$min_n2_label,
          .before = 1
        )
    })
  })
}

run_mi_regressions_pooled <- function(df, outcome_var, outcome_label) {
  split_df <- df %>%
    group_by(location_spec, cell_spec, cell_spec_label, threshold, threshold_label, min_n2, min_n2_label) %>%
    group_split()
  
  map_dfr(split_df, function(d_all) {
    meta <- d_all %>% slice(1)
    map_dfr(c(FALSE, TRUE), function(weighted) {
      coef_list <- d_all %>%
        group_by(imputation) %>%
        group_split() %>%
        map(function(d_imp) {
          imp <- unique(d_imp$imputation)
          fit_obj <- fit_one_lm(d_imp, outcome_var = outcome_var, weighted = weighted)
          if (is.null(fit_obj)) return(tibble())
          coef_table_from_lm(fit_obj$fit) %>%
            mutate(imputation = imp, formula = fit_obj$formula)
        })
      
      pooled <- pool_lm_terms_rubin(coef_list)
      if (nrow(pooled) == 0) return(tibble())
      
      pooled %>%
        mutate(
          analysis_sample = meta$analysis_sample,
          analysis_label = meta$analysis_label,
          outcome = outcome_label,
          weighted_by_n2 = weighted,
          formula = paste(unique(bind_rows(coef_list)$formula), collapse = " | "),
          location_spec = meta$location_spec,
          cell_spec = meta$cell_spec,
          cell_spec_label = meta$cell_spec_label,
          threshold = meta$threshold,
          threshold_label = meta$threshold_label,
          min_n2 = meta$min_n2,
          min_n2_label = meta$min_n2_label,
          .before = 1
        )
    })
  })
}

# ------------------------------------------------------------------------------
# 2. Read script 07 cell-level N-hat outputs
# ------------------------------------------------------------------------------

read_one_location_folder <- function(folder_slug) {
  folder <- file.path(script07_root, folder_slug)
  raw_path <- file.path(folder, "raw_cell_N_hat_values.csv")
  mi_path <- file.path(folder, "mi_cell_N_hat_values_by_imputation.csv")
  
  if (!file.exists(raw_path) || !file.exists(mi_path)) {
    warning("Skipping ", folder_slug, " because expected script 07 outputs were not found.")
    return(list(raw = tibble(), mi = tibble()))
  }
  
  raw <- read_csv(raw_path, show_col_types = FALSE) %>%
    mutate(source_folder = folder_slug)
  mi <- read_csv(mi_path, show_col_types = FALSE) %>%
    mutate(source_folder = folder_slug)
  
  list(raw = raw, mi = mi)
}

all_read <- map(location_spec_folders, read_one_location_folder)

raw_nhat <- map_dfr(all_read, "raw")
mi_nhat <- map_dfr(all_read, "mi")

if (nrow(raw_nhat) == 0 || nrow(mi_nhat) == 0) {
  stop("Could not find the needed script 07 N-hat outputs. Run script 07 for the desired location-spec folders first.")
}

required_cols <- c(
  "analysis_sample", "analysis_label", "threshold", "threshold_label",
  "location_spec", "cell_spec", "cell_spec_label", "min_n2", "min_n2_label",
  "cell_id", "n2", "N_hat_drinking_to_sober", "frac_drinking_hat"
)

missing_raw <- setdiff(required_cols, names(raw_nhat))
missing_mi <- setdiff(c(required_cols, "imputation"), names(mi_nhat))

if (length(missing_raw) > 0) stop("raw_cell_N_hat_values.csv missing columns: ", paste(missing_raw, collapse = ", "))
if (length(missing_mi) > 0) stop("mi_cell_N_hat_values_by_imputation.csv missing columns: ", paste(missing_mi, collapse = ", "))

raw_main <- raw_nhat %>%
  filter(analysis_sample == raw_sample_main, min_n2 == main_min_n2) %>%
  mutate(
    log_N_hat = log(N_hat_drinking_to_sober),
    logit_drunk_share_hat = qlogis(pmin(pmax(frac_drinking_hat, 1e-6), 1 - 1e-6)),
    data_source_for_regression = "raw_B"
  ) %>%
  add_cell_covariates()

mi_main <- mi_nhat %>%
  filter(analysis_sample == mi_sample_main, min_n2 == main_min_n2) %>%
  mutate(
    log_N_hat = log(N_hat_drinking_to_sober),
    logit_drunk_share_hat = qlogis(pmin(pmax(frac_drinking_hat, 1e-6), 1 - 1e-6)),
    data_source_for_regression = "mi_D"
  ) %>%
  add_cell_covariates()

# ------------------------------------------------------------------------------
# 3. Distribution checks for N_hat before deciding whether log transform helps
# ------------------------------------------------------------------------------

# For MI, use both imputation-level summaries and an across-imputation cell average
# for easier display. The regression itself still pools imputation-level models.
mi_cell_avg_for_distribution <- mi_main %>%
  group_by(
    analysis_sample, analysis_label, threshold, threshold_label,
    location_spec, cell_spec, cell_spec_label, min_n2, min_n2_label, cell_id
  ) %>%
  summarize(
    N_hat_drinking_to_sober = mean(N_hat_drinking_to_sober, na.rm = TRUE),
    log_N_hat = mean(log_N_hat, na.rm = TRUE),
    n2 = mean(n2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(data_source_for_regression = "mi_D_cell_avg")

raw_for_distribution <- raw_main %>%
  select(
    analysis_sample, analysis_label, threshold, threshold_label,
    location_spec, cell_spec, cell_spec_label, min_n2, min_n2_label,
    cell_id, N_hat_drinking_to_sober, log_N_hat, n2, data_source_for_regression
  )

distribution_input <- bind_rows(raw_for_distribution, mi_cell_avg_for_distribution)

distribution_summary <- distribution_input %>%
  group_by(
    data_source_for_regression, analysis_sample, analysis_label,
    location_spec, cell_spec, cell_spec_label, threshold, threshold_label, min_n2, min_n2_label
  ) %>%
  summarize(
    n_cells = n(),
    n2_total = sum(n2, na.rm = TRUE),
    N_min = min(N_hat_drinking_to_sober, na.rm = TRUE),
    N_p25 = quantile(N_hat_drinking_to_sober, 0.25, na.rm = TRUE),
    N_median = median(N_hat_drinking_to_sober, na.rm = TRUE),
    N_mean = mean(N_hat_drinking_to_sober, na.rm = TRUE),
    N_p75 = quantile(N_hat_drinking_to_sober, 0.75, na.rm = TRUE),
    N_p90 = quantile(N_hat_drinking_to_sober, 0.90, na.rm = TRUE),
    N_p95 = quantile(N_hat_drinking_to_sober, 0.95, na.rm = TRUE),
    N_max = max(N_hat_drinking_to_sober, na.rm = TRUE),
    N_skewness = skewness_simple(N_hat_drinking_to_sober),
    N_mean_median_ratio = safe_divide(mean(N_hat_drinking_to_sober, na.rm = TRUE), median(N_hat_drinking_to_sober, na.rm = TRUE)),
    logN_min = min(log_N_hat, na.rm = TRUE),
    logN_p25 = quantile(log_N_hat, 0.25, na.rm = TRUE),
    logN_median = median(log_N_hat, na.rm = TRUE),
    logN_mean = mean(log_N_hat, na.rm = TRUE),
    logN_p75 = quantile(log_N_hat, 0.75, na.rm = TRUE),
    logN_p90 = quantile(log_N_hat, 0.90, na.rm = TRUE),
    logN_p95 = quantile(log_N_hat, 0.95, na.rm = TRUE),
    logN_max = max(log_N_hat, na.rm = TRUE),
    logN_skewness = skewness_simple(log_N_hat),
    .groups = "drop"
  ) %>%
  mutate(
    log_transform_recommended = abs(N_skewness) > skew_flag_cutoff | N_mean_median_ratio > mean_median_ratio_cutoff,
    log_transform_reason = case_when(
      abs(N_skewness) > skew_flag_cutoff & N_mean_median_ratio > mean_median_ratio_cutoff ~ "N is skewed and mean is much larger than median",
      abs(N_skewness) > skew_flag_cutoff ~ "N is skewed",
      N_mean_median_ratio > mean_median_ratio_cutoff ~ "mean of N is much larger than median",
      TRUE ~ "N distribution not strongly skewed by simple heuristic"
    )
  ) %>%
  arrange(location_spec, cell_spec, threshold, data_source_for_regression)

write_csv(distribution_summary, file.path(out_dir, "N_hat_distribution_summary.csv"))

# Plot distributions (for diagnostics)
for (loc in unique(distribution_input$location_spec)) {
  plot_df <- distribution_input %>% filter(location_spec == loc)
  loc_slug <- str_to_lower(loc)
  
  p_N <- ggplot(plot_df, aes(x = N_hat_drinking_to_sober)) +
    geom_histogram(bins = 30) +
    facet_grid(threshold_label ~ cell_spec_label + data_source_for_regression, scales = "free_y") +
    labs(
      title = paste0("Distribution of N-hat: ", loc),
      x = "N-hat: drinking-to-sober exposure ratio",
      y = "Number of cells"
    ) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(
    filename = file.path(fig_dir, paste0("N_hat_distribution_", loc_slug, ".png")),
    plot = p_N, width = 16, height = 8, dpi = 300
  )
  
  p_logN <- ggplot(plot_df, aes(x = log_N_hat)) +
    geom_histogram(bins = 30) +
    facet_grid(threshold_label ~ cell_spec_label + data_source_for_regression, scales = "free_y") +
    labs(
      title = paste0("Distribution of log(N-hat): ", loc),
      x = "log(N-hat)",
      y = "Number of cells"
    ) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(
    filename = file.path(fig_dir, paste0("log_N_hat_distribution_", loc_slug, ".png")),
    plot = p_logN, width = 16, height = 8, dpi = 300
  )
}

# ------------------------------------------------------------------------------
# 4. Regress log(N_hat_c) on available cell characteristics
# ------------------------------------------------------------------------------

raw_logN_regressions <- run_raw_regressions(
  df = raw_main,
  outcome_var = "log_N_hat",
  outcome_label = "log_N_hat_raw_B"
)

mi_logN_regressions_pooled <- run_mi_regressions_pooled(
  df = mi_main,
  outcome_var = "log_N_hat",
  outcome_label = "log_N_hat_mi_D_pooled"
)

write_csv(raw_logN_regressions, file.path(out_dir, "regression_log_N_hat_raw_B.csv"))
write_csv(mi_logN_regressions_pooled, file.path(out_dir, "regression_log_N_hat_mi_D_pooled.csv"))

# ------------------------------------------------------------------------------
# 5. Delta_c = log(N_hat_MI) - log(N_hat_raw) for matched cells
# ------------------------------------------------------------------------------

raw_for_delta <- raw_nhat %>%
  filter(analysis_sample == raw_sample_for_delta, min_n2 == main_min_n2) %>%
  transmute(
    threshold, threshold_label, location_spec = standardize_location_spec(location_spec),
    cell_spec, cell_spec_label, min_n2, min_n2_label, cell_id,
    log_N_raw = log(N_hat_drinking_to_sober),
    N_hat_raw = N_hat_drinking_to_sober,
    n2_raw = n2
  )

mi_for_delta <- mi_nhat %>%
  filter(analysis_sample == mi_sample_for_delta, min_n2 == main_min_n2) %>%
  transmute(
    imputation, threshold, threshold_label, location_spec = standardize_location_spec(location_spec),
    cell_spec, cell_spec_label, min_n2, min_n2_label, cell_id,
    log_N_mi = log(N_hat_drinking_to_sober),
    N_hat_mi = N_hat_drinking_to_sober,
    n2_mi = n2
  )

delta_cell_level <- mi_for_delta %>%
  inner_join(
    raw_for_delta,
    by = c(
      "threshold", "threshold_label", "location_spec", "cell_spec", "cell_spec_label",
      "min_n2", "min_n2_label", "cell_id"
    )
  ) %>%
  mutate(
    delta_log_N_hat = log_N_mi - log_N_raw,
    ratio_N_hat_mi_to_raw = N_hat_mi / N_hat_raw,
    n2 = n2_mi,
    analysis_sample = paste0(mi_sample_for_delta, " minus ", raw_sample_for_delta),
    analysis_label = "D. MI BAC, all state-years minus B. Raw BAC, 95% complete"
  ) %>%
  add_cell_covariates()

write_csv(delta_cell_level, file.path(out_dir, "delta_log_N_hat_cell_level_B_to_D.csv"))

delta_summary <- delta_cell_level %>%
  group_by(location_spec, cell_spec, cell_spec_label, threshold, threshold_label, min_n2, min_n2_label) %>%
  summarize(
    n_matched_cell_imputations = n(),
    n_imputations = n_distinct(imputation),
    avg_matched_cells_per_imputation = n_matched_cell_imputations / n_imputations,
    mean_delta = mean(delta_log_N_hat, na.rm = TRUE),
    median_delta = median(delta_log_N_hat, na.rm = TRUE),
    p25_delta = quantile(delta_log_N_hat, 0.25, na.rm = TRUE),
    p75_delta = quantile(delta_log_N_hat, 0.75, na.rm = TRUE),
    mean_ratio_N_mi_to_raw = mean(ratio_N_hat_mi_to_raw, na.rm = TRUE),
    median_ratio_N_mi_to_raw = median(ratio_N_hat_mi_to_raw, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(location_spec, cell_spec, threshold)

write_csv(delta_summary, file.path(out_dir, "delta_log_N_hat_summary_B_to_D.csv"))

delta_regressions_pooled <- run_mi_regressions_pooled(
  df = delta_cell_level,
  outcome_var = "delta_log_N_hat",
  outcome_label = "delta_log_N_hat_D_minus_B_pooled"
)

write_csv(delta_regressions_pooled, file.path(out_dir, "regression_delta_log_N_hat_B_to_D_pooled.csv"))

# ------------------------------------------------------------------------------
# 6. Combined cleaned tables for easier review
# ------------------------------------------------------------------------------

combined_regression_results <- bind_rows(
  raw_logN_regressions %>% mutate(model_family = "raw_log_N_hat"),
  mi_logN_regressions_pooled %>% mutate(model_family = "mi_log_N_hat_pooled"),
  delta_regressions_pooled %>% mutate(model_family = "delta_log_N_hat_pooled")
) %>%
  mutate(
    stars = case_when(
      is.na(p_value) ~ "",
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      p_value < 0.10 ~ ".",
      TRUE ~ ""
    )
  ) %>%
  relocate(model_family, .before = 1) %>%
  arrange(location_spec, cell_spec, threshold, model_family, weighted_by_n2, term)

combined_regression_results_display <- combined_regression_results %>%
  mutate(
    estimate = round(estimate, 4),
    std_error = round(std_error, 4),
    statistic = round(statistic, 3),
    p_value = round(p_value, 4),
    r_squared = ifelse("r_squared" %in% names(.), round(r_squared, 3), NA_real_),
    adj_r_squared = ifelse("adj_r_squared" %in% names(.), round(adj_r_squared, 3), NA_real_),
    avg_r_squared = ifelse("avg_r_squared" %in% names(.), round(avg_r_squared, 3), NA_real_),
    avg_adj_r_squared = ifelse("avg_adj_r_squared" %in% names(.), round(avg_adj_r_squared, 3), NA_real_)
  )

write_csv(combined_regression_results, file.path(out_dir, "combined_regression_results.csv"))
write_csv(combined_regression_results_display, file.path(out_dir, "combined_regression_results_display.csv"))

# A compact model key to quickly see which regressions actually ran.
model_key <- combined_regression_results %>%
  distinct(
    model_family, outcome, analysis_label, location_spec, cell_spec, cell_spec_label,
    threshold, threshold_label, min_n2, min_n2_label, weighted_by_n2, formula
  ) %>%
  arrange(location_spec, cell_spec, threshold, model_family, weighted_by_n2)

write_csv(model_key, file.path(out_dir, "model_key.csv"))

# ------------------------------------------------------------------------------
# 7. figures for main regression results
# ------------------------------------------------------------------------------

# These figures are meant to be summaries of the regression
# tables. They focus on weighted models and selected interpretable coefficients.
# Coefficients are plotted as percent changes in N_hat:
#   100 * (exp(beta) - 1)
# because the regression outcome is log(N_hat).

analysis_fig_dir <- file.path(fig_dir, "analysis_figures")
dir.create(analysis_fig_dir, recursive = TRUE, showWarnings = FALSE)

# Helper: clean selected coefficient terms
clean_coef_for_plot <- function(df) {
  df %>%
    filter(
      weighted_by_n2 == TRUE,
      term %in% c(
        "RUR_URB_GROUPrural",
        "WEEKENDTRUE",
        "PERIODcovid",
        "PERIODpost_covid"
      )
    ) %>%
    mutate(
      term_label = case_when(
        term == "RUR_URB_GROUPrural" ~ "Rural",
        term == "WEEKENDTRUE" ~ "Weekend",
        term == "PERIODcovid" ~ "COVID period",
        term == "PERIODpost_covid" ~ "Post-COVID period",
        TRUE ~ term
      ),
      sample_label_short = case_when(
        model_family == "raw_log_N_hat" ~ "B. Raw BAC, 95% complete",
        model_family == "mi_log_N_hat_pooled" ~ "D. MI BAC, all state-years",
        model_family == "delta_log_N_hat_pooled" ~ "D - B",
        TRUE ~ model_family
      ),
      estimate_pct = 100 * (exp(estimate) - 1),
      ci_low_pct = 100 * (exp(estimate - 1.96 * std_error) - 1),
      ci_high_pct = 100 * (exp(estimate + 1.96 * std_error) - 1)
    )
}

# ------------------------------------------------------------------------------
# Figure 1: Preferred LOCATION_WEEKEND model
# ------------------------------------------------------------------------------

# This is probably the cleanest single figure
# It uses the preferred lower-dimensional location specification:
# rural/urban x state x weekend.
# It shows how rural and weekend cells differ in implied drinking exposure.

preferred_location_weekend_plot_df <- combined_regression_results %>%
  filter(
    model_family %in% c("raw_log_N_hat", "mi_log_N_hat_pooled"),
    location_spec == "LOCATION_WEEKEND",
    cell_spec == "rur_urb_state_weekend"
  ) %>%
  clean_coef_for_plot() %>%
  filter(term_label %in% c("Rural", "Weekend")) %>%
  mutate(
    term_label = factor(term_label, levels = c("Weekend", "Rural")),
    sample_label_short = factor(
      sample_label_short,
      levels = c("B. Raw BAC, 95% complete", "D. MI BAC, all state-years")
    )
  )

if (nrow(preferred_location_weekend_plot_df) > 0) {
  p_location_weekend <- ggplot(
    preferred_location_weekend_plot_df,
    aes(
      x = estimate_pct,
      y = term_label,
      xmin = ci_low_pct,
      xmax = ci_high_pct,
      shape = sample_label_short
    )
  ) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_pointrange(
      position = position_dodge(width = 0.55),
      size = 0.45
    ) +
    facet_wrap(~ threshold_label, nrow = 1) +
    labs(
      title = "Main regression results: LOCATION_WEEKEND specification",
      subtitle = "Weighted regressions of log(N-hat) using rural/urban x state x weekend cells",
      x = "Estimated percent difference in N-hat",
      y = NULL,
      shape = "Analysis sample",
      caption = "Points show coefficient estimates transformed as 100*(exp(beta)-1). Bars show approximate 95% CIs."
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
  
  ggsave(
    filename = file.path(analysis_fig_dir, "location_weekend_coefficients.png"),
    plot = p_location_weekend,
    width = 11,
    height = 5,
    dpi = 300
  )
}

# ------------------------------------------------------------------------------
# Figure 2: Robustness across cell-specification families
# ------------------------------------------------------------------------------

# This figure compares the same main coefficient types across representative
# models from each cell-specification family. It is useful for showing that
# the weekend result is stable, while rural/urban appears where available.

representative_models <- tribble(
  ~location_spec,        ~cell_spec,                         ~representative_model_label,
  "LOCATION_WEEKEND",    "rur_urb_state_weekend",             "Location/weekend: rural/urban x state x weekend",
  "RUR_URB",             "hour_rur_urb_period_weekend",       "Rural/urban: hour x rural/urban x period x weekend",
  "STATE",               "hour_period_weekend",              "State family: hour x period x weekend"
)

robustness_plot_df <- combined_regression_results %>%
  filter(model_family %in% c("raw_log_N_hat", "mi_log_N_hat_pooled")) %>%
  inner_join(
    representative_models,
    by = c("location_spec", "cell_spec")
  ) %>%
  clean_coef_for_plot() %>%
  filter(term_label %in% c("Weekend", "Rural", "COVID period", "Post-COVID period")) %>%
  mutate(
    term_label = factor(
      term_label,
      levels = c("Weekend", "Rural", "COVID period", "Post-COVID period")
    ),
    sample_label_short = factor(
      sample_label_short,
      levels = c("B. Raw BAC, 95% complete", "D. MI BAC, all state-years")
    ),
    representative_model_label = factor(
      representative_model_label,
      levels = representative_models$representative_model_label
    )
  )

if (nrow(robustness_plot_df) > 0) {
  p_robustness <- ggplot(
    robustness_plot_df,
    aes(
      x = estimate_pct,
      y = term_label,
      xmin = ci_low_pct,
      xmax = ci_high_pct,
      shape = sample_label_short
    )
  ) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_pointrange(
      position = position_dodge(width = 0.55),
      size = 0.35
    ) +
    facet_grid(representative_model_label ~ threshold_label, scales = "free_y", space = "free_y") +
    labs(
      title = "Regression coefficient summary across representative specifications",
      subtitle = "Weighted regressions of log(N-hat); coefficients shown as percent differences in N-hat",
      x = "Estimated percent difference in N-hat",
      y = NULL,
      shape = "Analysis sample",
      caption = "STATE panel does not include rural/urban because that representative model has no rural/urban variable."
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      strip.text.y = element_text(angle = 0)
    )
  
  ggsave(
    filename = file.path(analysis_fig_dir, "regression_coefficients_by_spec_family.png"),
    plot = p_robustness,
    width = 13,
    height = 8,
    dpi = 300
  )
}

# ------------------------------------------------------------------------------
# Figure 3: delta plot, D minus B
# ------------------------------------------------------------------------------

# This figure summarizes whether MI all-state-years changes log(N-hat) relative
# to raw BAC 95%-complete state-years in matched cells.
# Negative coefficient: D - B difference is smaller for that cell type.
# Positive coefficient: D - B difference is larger for that cell type.

delta_plot_df <- combined_regression_results %>%
  filter(
    model_family == "delta_log_N_hat_pooled",
    location_spec == "LOCATION_WEEKEND",
    cell_spec == "rur_urb_state_weekend"
  ) %>%
  clean_coef_for_plot() %>%
  filter(term_label %in% c("Rural", "Weekend")) %>%
  mutate(
    term_label = factor(term_label, levels = c("Weekend", "Rural"))
  )

if (nrow(delta_plot_df) > 0) {
  p_delta <- ggplot(
    delta_plot_df,
    aes(
      x = estimate,
      y = term_label,
      xmin = estimate - 1.96 * std_error,
      xmax = estimate + 1.96 * std_error
    )
  ) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_pointrange(size = 0.45) +
    facet_wrap(~ threshold_label, nrow = 1) +
    labs(
      title = "Matched-cell difference in log(N-hat): MI all-state-years minus raw 95%-complete",
      subtitle = "LOCATION_WEEKEND: rural/urban x state x weekend cells",
      x = "Coefficient in D - B regression",
      y = NULL,
      caption = "Outcome is log(N-hat_D) - log(N-hat_B). Bars show approximate 95% CIs."
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank()
    )
  
  ggsave(
    filename = file.path(analysis_fig_dir, "delta_log_N_hat_location_weekend.png"),
    plot = p_delta,
    width = 11,
    height = 5,
    dpi = 300
  )
}

# Save the plotting data just in case
write_csv(
  preferred_location_weekend_plot_df,
  file.path(out_dir, "figure_location_weekend_coefficients_data.csv")
)

write_csv(
  robustness_plot_df,
  file.path(out_dir, "figure_regression_coefficients_by_spec_family_data.csv")
)

write_csv(
  delta_plot_df,
  file.path(out_dir, "figure_delta_log_N_hat_location_weekend_data.csv")
)

cat("Analysis figures saved to:\n")
cat(analysis_fig_dir, "\n")


