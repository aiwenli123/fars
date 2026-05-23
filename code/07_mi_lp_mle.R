################################################################################
# 07_mi_lp_mle.R
#
# Purpose:
#   MI version of the LP table-2-style MLE.
#
#   Main comparison:
#     (1) Raw FARS measured BAC, restricted to BAC_COMPLETE_95 state-years
#     (2) MIPER multiple-imputation BAC, using all state-years / all rows
#
#   For thresholds BAC >= 0.05, 0.08, and 0.10:
#     - Recreate LP one- and two-vehicle crash types
#     - Run the two-car theta MLE and plug-in one-car lambda MLE
#     - For MI, estimate separately in each imputed dataset m = 1,...,10
#     - Combine MI estimates using Rubin's rules as outlined in Subramanian (2002)
#
# Main input:
#   data/processed/driver_miper_long_2015_2024.csv
#
# Main outputs:
#   data/processed/script07_mi_lp_mle/mi_lp_mle_imputation_estimates.csv
#   data/processed/script07_mi_lp_mle/mi_lp_mle_rubin_combined.csv
#   data/processed/script07_mi_lp_mle/raw_bac95_lp_mle_estimates.csv
#   data/processed/script07_mi_lp_mle/lp_mle_raw_vs_mi_comparison.csv
#   data/processed/script07_mi_lp_mle/lp_mle_raw_vs_mi_wide_theta.csv
#   data/processed/script07_mi_lp_mle/lp_mle_raw_vs_mi_wide_lambda.csv
#   data/processed/script07_mi_lp_mle/crash_type_counts_raw_vs_mi.csv
################################################################################

rm(list = ls())

library(tidyverse)
library(scales)

# ------------------------------------------------------------------------------
# 0. User options
# ------------------------------------------------------------------------------

mi_driver_path <- "data/processed/driver_miper_long_2015_2024.csv"

out_dir <- "data/processed/script07_mi_lp_mle"
fig_dir <- "figures/script07_mi_lp_mle"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Bounds for theta and lambda, matching script 04.
lower_risk <- 1
upper_risk <- 200

# Finite-difference step for curvature-based SEs on the log scale.
fd_step <- 1e-4

threshold_specs <- tibble(
  threshold = c(0.05, 0.08, 0.10),
  threshold_label = c("BAC >= 0.05", "BAC >= 0.08", "BAC >= 0.10"),
  raw_definition_var = c("BAC_GE_05", "BAC_GE_08", "BAC_GE_10"),
  mi_definition_var  = c("MI_BAC_GE_05", "MI_BAC_GE_08", "MI_BAC_GE_10")
)

cell_specs <- tibble(
  cell_spec = c(
    "all_data",
    "hour",
    "hour_period",
    "hour_period_weekend",
    "hour_state_period",
    "hour_state_period_weekend"
  ),
  cell_var = c(
    "CELL_ALL",
    "CELL_HOUR_07",
    "CELL_HOUR_PERIOD_07",
    "CELL_HOUR_PERIOD_WEEKEND_07",
    "CELL_HOUR_STATE_PERIOD_07",
    "CELL_HOUR_STATE_PERIOD_WEEKEND_07"
  ),
  cell_spec_label = c(
    "all data",
    "hour",
    "hour x period",
    "hour x period x weekend",
    "hour x state x period",
    "hour x state x period x weekend"
  )
)

min_n2_specs <- tibble(
  min_n2 = c(1, 5),
  min_n2_label = c("No n2 filter", "n2 >= 5")
)

# ------------------------------------------------------------------------------
# 1. Helper functions
# ------------------------------------------------------------------------------

safe_divide <- function(num, den) {
  ifelse(is.na(den) | den == 0, NA_real_, num / den)
}

format_est_se_flag <- function(est, se, status, digits = 3) {
  case_when(
    is.na(est) ~ NA_character_,
    !is.na(se) ~ paste0(round(est, digits), " (", round(se, digits), ")"),
    TRUE ~ paste0(round(est, digits), " [", status, "]")
  )
}

curvature_se_1d <- function(objective, x_hat, est_hat, h = 1e-4) {
  f0 <- objective(x_hat)
  fp <- objective(x_hat + h)
  fm <- objective(x_hat - h)

  second_deriv <- (fp - 2 * f0 + fm) / (h^2)

  if (!is.finite(second_deriv) || second_deriv <= 0) {
    return(list(
      se = NA_real_,
      second_deriv = second_deriv,
      status = "bad_curvature"
    ))
  }

  se_log <- sqrt(1 / second_deriv)
  se_original <- est_hat * se_log

  list(
    se = se_original,
    second_deriv = second_deriv,
    status = "ok"
  )
}

closed_form_theta_cell <- function(DD, DS, SS) {
  R <- safe_divide(DS^2, DD * SS)
  theta <- rep(NA_real_, length(R))

  valid_boundary <- !is.na(R) & DD > 0 & SS > 0 & R < 4
  valid_closed <- !is.na(R) & DD > 0 & SS > 0 & R >= 4

  theta[valid_boundary] <- 1

  disc <- R[valid_closed]^2 - 4 * R[valid_closed]
  theta[valid_closed] <- ((R[valid_closed] - 2) + sqrt(disc)) / 2

  theta
}

combine_mi <- function(est, se) {
  ok_est <- is.finite(est)
  ok_se <- is.finite(se)
  ok_both <- ok_est & ok_se

  m_total <- length(est)
  m_est <- sum(ok_est)
  m_used_for_W <- sum(ok_both)

  qbar <- mean(est[ok_est], na.rm = TRUE)
  W <- mean(se[ok_both]^2, na.rm = TRUE)
  B <- var(est[ok_est], na.rm = TRUE)

  # If only one finite estimate exists, var() is NA
  Tvar <- W + (1 + 1 / m_total) * B

  tibble(
    m_total = m_total,
    m_est = m_est,
    m_used_for_W = m_used_for_W,
    estimate = qbar,
    within_var = W,
    between_var = B,
    total_var = Tvar,
    se_total = sqrt(Tvar),
    ci_low = qbar - 1.96 * sqrt(Tvar),
    ci_high = qbar + 1.96 * sqrt(Tvar)
  )
}

# ------------------------------------------------------------------------------
# 2. Load MI long file and create analysis variables
# ------------------------------------------------------------------------------

if (!file.exists(mi_driver_path)) {
  stop("Missing input file: ", mi_driver_path, "\nPlease run 06_prep_mi_data.R first.")
}

driver_mi_long <- read_csv(mi_driver_path, show_col_types = FALSE) %>%
  rename_with(toupper)

cat("Loaded MI driver-level long file:\n")
cat("  Rows:", comma(nrow(driver_mi_long)), "\n")
cat("  Distinct imputations:", paste(sort(unique(driver_mi_long$IMPUTATION)), collapse = ", "), "\n")
cat("  Distinct driver rows per imputation should be about:", comma(nrow(driver_mi_long) / n_distinct(driver_mi_long$IMPUTATION)), "\n\n")

required_cols <- c(
  "YEAR_FILE", "YEAR", "STATE", "ST_CASE", "VEH_NO", "PER_NO",
  "VE_FORMS", "HOUR", "IMPUTATION",
  "BAC_COMPLETE_95", "BAC_GE_05", "BAC_GE_08", "BAC_GE_10",
  "MI_BAC_GE_05", "MI_BAC_GE_08", "MI_BAC_GE_10"
)

missing_required <- setdiff(required_cols, names(driver_mi_long))
if (length(missing_required) > 0) {
  stop("Missing required columns in driver_mi_long: ", paste(missing_required, collapse = ", "))
}

# Create period/weekend/cell variables for this script.
driver_mi_long <- driver_mi_long %>%
  mutate(
    PERIOD_07 = case_when(
      YEAR %in% 2015:2019 ~ "pre_covid",
      YEAR %in% 2020:2021 ~ "covid",
      YEAR %in% 2022:2024 ~ "post_covid",
      TRUE ~ NA_character_
    ),
    PERIOD_07 = factor(PERIOD_07, levels = c("pre_covid", "covid", "post_covid")),

    WEEKEND_07 = case_when(
      "WEEKEND_LP" %in% names(.) ~ as.logical(WEEKEND_LP),
      DAY_WEEK %in% c(1, 7) ~ TRUE,
      TRUE ~ FALSE
    ),

    CELL_ALL = "all_data",
    CELL_HOUR_07 = as.character(HOUR),
    CELL_HOUR_PERIOD_07 = paste(HOUR, PERIOD_07, sep = "_"),
    CELL_HOUR_PERIOD_WEEKEND_07 = paste(HOUR, PERIOD_07, WEEKEND_07, sep = "_"),
    CELL_HOUR_STATE_PERIOD_07 = paste(HOUR, STATE, PERIOD_07, sep = "_"),
    CELL_HOUR_STATE_PERIOD_WEEKEND_07 = paste(HOUR, STATE, PERIOD_07, WEEKEND_07, sep = "_")
  ) %>%
  filter(!is.na(PERIOD_07), VE_FORMS %in% c(1, 2))

# Raw measured-BAC comparison should not use the ten duplicated imputation rows.
# The raw BAC variables are identical across imputations, so keep imputation 1 only.
driver_raw_one_row <- driver_mi_long %>%
  filter(IMPUTATION == min(IMPUTATION, na.rm = TRUE))

cat("Rows for raw measured-BAC comparison after keeping one imputation copy:", comma(nrow(driver_raw_one_row)), "\n")
cat("Rows for MI analysis:", comma(nrow(driver_mi_long)), "\n\n")

# ------------------------------------------------------------------------------
# 3. Crash classification and cell-count construction
# ------------------------------------------------------------------------------

make_crash_types <- function(driver_data,
                             drinking_var,
                             threshold_label,
                             analysis_sample,
                             imputation_value = NA_integer_) {

  if (!(drinking_var %in% names(driver_data))) {
    stop("Could not find drinking/BAC indicator column: ", drinking_var)
  }

  driver_data %>%
    mutate(
      DRINKING_STATUS_07 = .data[[drinking_var]],
      CRASH_ID_07 = paste(YEAR_FILE, ST_CASE, sep = "_")
    ) %>%
    group_by(
      YEAR_FILE, YEAR, PERIOD_07, STATE, ST_CASE, VE_FORMS,
      HOUR, WEEKEND_07,
      CELL_ALL, CELL_HOUR_07, CELL_HOUR_PERIOD_07,
      CELL_HOUR_PERIOD_WEEKEND_07, CELL_HOUR_STATE_PERIOD_07,
      CELL_HOUR_STATE_PERIOD_WEEKEND_07
    ) %>%
    summarize(
      n_drivers = n(),
      n_drinking = sum(DRINKING_STATUS_07 == 1, na.rm = TRUE),
      n_sober = sum(DRINKING_STATUS_07 == 0, na.rm = TRUE),
      n_missing = sum(is.na(DRINKING_STATUS_07)),
      .groups = "drop"
    ) %>%
    mutate(
      threshold_label = threshold_label,
      analysis_sample = analysis_sample,
      imputation = imputation_value,

      crash_type = case_when(
        VE_FORMS == 2 & n_drivers == 2 & n_missing == 0 & n_drinking == 2 & n_sober == 0 ~ "DD",
        VE_FORMS == 2 & n_drivers == 2 & n_missing == 0 & n_drinking == 1 & n_sober == 1 ~ "DS",
        VE_FORMS == 2 & n_drivers == 2 & n_missing == 0 & n_drinking == 0 & n_sober == 2 ~ "SS",

        VE_FORMS == 1 & n_drivers == 1 & n_missing == 0 & n_drinking == 1 & n_sober == 0 ~ "D",
        VE_FORMS == 1 & n_drivers == 1 & n_missing == 0 & n_drinking == 0 & n_sober == 1 ~ "S",

        TRUE ~ "missing_or_ambiguous"
      )
    )
}

make_cell_counts <- function(crash_types_df, cell_var) {
  crash_types_df %>%
    mutate(cell_id = .data[[cell_var]]) %>%
    group_by(cell_id) %>%
    summarize(
      DD = sum(crash_type == "DD"),
      DS = sum(crash_type == "DS"),
      SS = sum(crash_type == "SS"),
      D  = sum(crash_type == "D"),
      S  = sum(crash_type == "S"),

      missing_two_vehicle = sum(VE_FORMS == 2 & crash_type == "missing_or_ambiguous"),
      missing_one_vehicle = sum(VE_FORMS == 1 & crash_type == "missing_or_ambiguous"),

      total_two_vehicle = sum(VE_FORMS == 2),
      total_one_vehicle = sum(VE_FORMS == 1),
      .groups = "drop"
    ) %>%
    mutate(
      n2 = DD + DS + SS,
      n1 = D + S,
      total_complete_case = n1 + n2,

      R_cell = safe_divide(DS^2, DD * SS),
      theta_cell_closed = closed_form_theta_cell(DD, DS, SS)
    ) %>%
    filter(total_complete_case > 0)
}

summarize_crash_types <- function(crash_types_df) {
  crash_types_df %>%
    summarize(
      DD = sum(crash_type == "DD"),
      DS = sum(crash_type == "DS"),
      SS = sum(crash_type == "SS"),
      D  = sum(crash_type == "D"),
      S  = sum(crash_type == "S"),
      missing_or_ambiguous = sum(crash_type == "missing_or_ambiguous"),
      total_crashes = n(),
      complete_two_vehicle = DD + DS + SS,
      complete_one_vehicle = D + S,
      pct_missing_or_ambiguous = missing_or_ambiguous / total_crashes
    )
}

summarize_cell_diagnostics <- function(counts) {
  counts_two_available <- counts %>% filter(n2 > 0)

  if (nrow(counts_two_available) == 0) {
    return(tibble(
      n_cells_total = nrow(counts),
      n_theta_cells_available = 0,
      n_two_vehicle_complete_available = 0,
      n_two_vehicle_all = sum(counts$total_two_vehicle),
      pct_two_vehicle_missing_or_ambiguous =
        safe_divide(sum(counts$missing_two_vehicle), sum(counts$total_two_vehicle)),
      median_theta_cell_closed = NA_real_,
      p90_theta_cell_closed = NA_real_,
      p95_theta_cell_closed = NA_real_
    ))
  }

  tibble(
    n_cells_total = nrow(counts),
    n_theta_cells_available = nrow(counts_two_available),
    n_two_vehicle_complete_available = sum(counts_two_available$n2),
    n_two_vehicle_all = sum(counts$total_two_vehicle),
    pct_two_vehicle_missing_or_ambiguous =
      safe_divide(sum(counts$missing_two_vehicle), sum(counts$total_two_vehicle)),
    n_theta_cells_no_DD_available = sum(counts_two_available$DD == 0),
    n_theta_cells_no_DS_available = sum(counts_two_available$DS == 0),
    n_theta_cells_no_SS_available = sum(counts_two_available$SS == 0),
    pct_theta_cells_n2_le_1_available = mean(counts_two_available$n2 <= 1),
    pct_theta_cells_n2_le_2_available = mean(counts_two_available$n2 <= 2),
    pct_theta_cells_n2_le_5_available = mean(counts_two_available$n2 <= 5),
    median_theta_cell_closed = median(counts_two_available$theta_cell_closed, na.rm = TRUE),
    p90_theta_cell_closed = quantile(counts_two_available$theta_cell_closed, 0.90, na.rm = TRUE),
    p95_theta_cell_closed = quantile(counts_two_available$theta_cell_closed, 0.95, na.rm = TRUE)
  )
}

# ------------------------------------------------------------------------------
# 4. LP MLE functions from script 04
# ------------------------------------------------------------------------------

profile_eta_two_car <- function(theta, counts_two, max_iter = 100, tol = 1e-10) {
  DD <- counts_two$DD
  DS <- counts_two$DS
  SS <- counts_two$SS
  n2 <- counts_two$n2

  obs_D_drivers <- 2 * DD + DS
  obs_S_drivers <- 2 * SS + DS

  eta <- log((obs_D_drivers + 0.5) / (obs_S_drivers + 0.5))
  eta <- pmin(pmax(eta, -20), 20)

  for (iter in seq_len(max_iter)) {
    N <- exp(eta)

    Z <- theta * N^2 + (theta + 1) * N + 1
    A <- 2 * theta * N^2 + (theta + 1) * N
    Aprime <- 4 * theta * N^2 + (theta + 1) * N

    grad <- (2 * DD + DS) - n2 * (A / Z)
    hess <- -n2 * ((Aprime * Z - A^2) / Z^2)

    bad <- is.na(hess) | !is.finite(hess) | abs(hess) < 1e-10

    step <- grad / hess
    step[bad] <- 0
    step <- pmin(pmax(step, -5), 5)

    eta_new <- eta - step
    eta_new <- pmin(pmax(eta_new, -20), 20)

    if (max(abs(eta_new - eta), na.rm = TRUE) < tol) {
      eta <- eta_new
      break
    }

    eta <- eta_new
  }

  eta
}

two_car_loglik <- function(log_theta, counts_two, return_eta = FALSE) {
  theta <- exp(log_theta)

  if (!is.finite(theta) || theta <= 0) {
    if (return_eta) return(list(loglik = -Inf, theta = theta, eta = NA, N = NA))
    return(-Inf)
  }

  DD <- counts_two$DD
  DS <- counts_two$DS
  n2 <- counts_two$n2

  eta <- profile_eta_two_car(theta, counts_two)
  N <- exp(eta)

  Z <- theta * N^2 + (theta + 1) * N + 1

  ll <- DD * (log(theta) + 2 * eta) +
    DS * (log(theta + 1) + eta) -
    n2 * log(Z)

  total_ll <- sum(ll, na.rm = TRUE)

  if (return_eta) {
    return(list(
      loglik = total_ll,
      theta = theta,
      eta = eta,
      N = N
    ))
  }

  total_ll
}

fit_theta_simple <- function(counts,
                             min_n2 = 1,
                             lower_risk = 1,
                             upper_risk = 100,
                             fd_step = 1e-4) {

  counts_two_available <- counts %>% filter(n2 > 0)
  counts_two <- counts_two_available %>% filter(n2 >= min_n2)

  if (nrow(counts_two) == 0) {
    return(list(
      estimates = tibble(
        theta = NA_real_,
        theta_se = NA_real_,
        theta_logLik = NA_real_,
        theta_status = "no_two_car_data_after_min_n2_filter",
        theta_curvature = NA_real_,
        min_n2 = min_n2,
        n_theta_cells_available = nrow(counts_two_available),
        n_theta_cells_used = 0,
        n_two_vehicle_complete_available = sum(counts_two_available$n2),
        n_two_vehicle_complete_used = 0,
        pct_two_vehicle_complete_used = 0
      ),
      cell_prevalence = tibble()
    ))
  }

  lower_log <- log(lower_risk + 1e-8)
  upper_log <- log(upper_risk)

  objective <- function(log_theta) {
    val <- -two_car_loglik(log_theta, counts_two)
    if (!is.finite(val)) return(1e100)
    val
  }

  opt <- optimize(
    f = objective,
    interval = c(lower_log, upper_log),
    tol = 1e-10
  )

  log_theta_hat <- opt$minimum
  prof <- two_car_loglik(log_theta_hat, counts_two, return_eta = TRUE)
  theta_hat <- prof$theta

  near_lower <- abs(log_theta_hat - lower_log) < 1e-4
  near_upper <- abs(log_theta_hat - upper_log) < 1e-4

  if (near_lower) {
    theta_status <- "near_lower_bound"
    theta_se <- NA_real_
    theta_curvature <- NA_real_
  } else if (near_upper) {
    theta_status <- "near_upper_bound"
    theta_se <- NA_real_
    theta_curvature <- NA_real_
  } else {
    curv <- curvature_se_1d(
      objective = objective,
      x_hat = log_theta_hat,
      est_hat = theta_hat,
      h = fd_step
    )

    theta_status <- curv$status
    theta_se <- curv$se
    theta_curvature <- curv$second_deriv
  }

  cell_prevalence <- counts_two %>%
    mutate(
      min_n2 = min_n2,
      theta = theta_hat,
      eta_hat = prof$eta,
      N_hat_drinking_to_sober = prof$N,
      frac_drinking_hat = N_hat_drinking_to_sober / (1 + N_hat_drinking_to_sober)
    )

  n_two_available <- sum(counts_two_available$n2)
  n_two_used <- sum(counts_two$n2)

  estimates <- tibble(
    theta = theta_hat,
    theta_se = theta_se,
    theta_logLik = prof$loglik,
    theta_status = theta_status,
    theta_curvature = theta_curvature,
    min_n2 = min_n2,

    n_theta_cells_available = nrow(counts_two_available),
    n_theta_cells_used = nrow(counts_two),
    n_theta_cells_excluded_by_min_n2 = nrow(counts_two_available) - nrow(counts_two),

    n_two_vehicle_complete_available = n_two_available,
    n_two_vehicle_complete_used = n_two_used,
    pct_two_vehicle_complete_used = safe_divide(n_two_used, n_two_available),

    n_two_vehicle_all = sum(counts$total_two_vehicle),
    pct_two_vehicle_missing_or_ambiguous =
      safe_divide(sum(counts$missing_two_vehicle), sum(counts$total_two_vehicle)),

    implied_frac_drinking_weighted_two_car = weighted.mean(
      cell_prevalence$frac_drinking_hat,
      w = cell_prevalence$n2,
      na.rm = TRUE
    ),

    n_theta_cells_no_DD_used = sum(counts_two$DD == 0),
    n_theta_cells_no_DS_used = sum(counts_two$DS == 0),
    n_theta_cells_no_SS_used = sum(counts_two$SS == 0),

    median_theta_cell_closed_used = median(counts_two$theta_cell_closed, na.rm = TRUE),
    p90_theta_cell_closed_used = quantile(counts_two$theta_cell_closed, 0.90, na.rm = TRUE),
    p95_theta_cell_closed_used = quantile(counts_two$theta_cell_closed, 0.95, na.rm = TRUE)
  )

  list(
    estimates = estimates,
    cell_prevalence = cell_prevalence
  )
}

lambda_plugin_loglik <- function(log_lambda, plugin_counts) {
  lambda <- exp(log_lambda)
  if (!is.finite(lambda) || lambda <= 0) return(-Inf)

  D <- plugin_counts$D
  n1 <- plugin_counts$n1
  Nhat <- plugin_counts$N_hat_drinking_to_sober

  ll <- D * (log(lambda) + log(Nhat)) -
    n1 * log(1 + lambda * Nhat)

  sum(ll, na.rm = TRUE)
}

fit_lambda_simple <- function(counts,
                              cell_prevalence,
                              min_n2 = 1,
                              lower_risk = 1,
                              upper_risk = 100,
                              fd_step = 1e-4) {

  if (nrow(cell_prevalence) == 0) {
    return(tibble(
      lambda_plugin = NA_real_,
      lambda_plugin_se = NA_real_,
      lambda_plugin_logLik = NA_real_,
      lambda_plugin_status = "no_theta_prevalence",
      lambda_plugin_curvature = NA_real_,
      min_n2 = min_n2,
      n_lambda_cells = 0,
      n_one_vehicle_complete_used_for_lambda = 0,
      pct_one_vehicle_complete_used_for_lambda = NA_real_
    ))
  }

  plugin_counts <- counts %>%
    select(cell_id, D, S, n1) %>%
    inner_join(
      cell_prevalence %>% select(cell_id, N_hat_drinking_to_sober),
      by = "cell_id"
    ) %>%
    filter(
      n1 > 0,
      is.finite(N_hat_drinking_to_sober),
      N_hat_drinking_to_sober > 0
    )

  n_one_total <- sum(counts$n1)
  n_one_used <- sum(plugin_counts$n1)

  if (nrow(plugin_counts) == 0) {
    return(tibble(
      lambda_plugin = NA_real_,
      lambda_plugin_se = NA_real_,
      lambda_plugin_logLik = NA_real_,
      lambda_plugin_status = "no_one_car_data_with_Nhat",
      lambda_plugin_curvature = NA_real_,
      min_n2 = min_n2,
      n_lambda_cells = 0,
      n_one_vehicle_complete_total = n_one_total,
      n_one_vehicle_complete_used_for_lambda = 0,
      n_one_vehicle_complete_excluded_no_Nhat = n_one_total,
      pct_one_vehicle_complete_used_for_lambda = 0
    ))
  }

  lower_log <- log(lower_risk + 1e-8)
  upper_log <- log(upper_risk)

  objective <- function(log_lambda) {
    val <- -lambda_plugin_loglik(log_lambda, plugin_counts)
    if (!is.finite(val)) return(1e100)
    val
  }

  opt <- optimize(
    f = objective,
    interval = c(lower_log, upper_log),
    tol = 1e-10
  )

  log_lambda_hat <- opt$minimum
  lambda_hat <- exp(log_lambda_hat)
  ll_hat <- lambda_plugin_loglik(log_lambda_hat, plugin_counts)

  near_lower <- abs(log_lambda_hat - lower_log) < 1e-4
  near_upper <- abs(log_lambda_hat - upper_log) < 1e-4

  if (near_lower) {
    lambda_status <- "near_lower_bound"
    lambda_se <- NA_real_
    lambda_curvature <- NA_real_
  } else if (near_upper) {
    lambda_status <- "near_upper_bound"
    lambda_se <- NA_real_
    lambda_curvature <- NA_real_
  } else {
    curv <- curvature_se_1d(
      objective = objective,
      x_hat = log_lambda_hat,
      est_hat = lambda_hat,
      h = fd_step
    )

    lambda_status <- curv$status
    lambda_se <- curv$se
    lambda_curvature <- curv$second_deriv
  }

  tibble(
    lambda_plugin = lambda_hat,
    lambda_plugin_se = lambda_se,
    lambda_plugin_logLik = ll_hat,
    lambda_plugin_status = lambda_status,
    lambda_plugin_curvature = lambda_curvature,
    min_n2 = min_n2,

    n_lambda_cells = nrow(plugin_counts),
    n_one_vehicle_complete_total = n_one_total,
    n_one_vehicle_complete_used_for_lambda = n_one_used,
    n_one_vehicle_complete_excluded_no_Nhat = n_one_total - n_one_used,
    pct_one_vehicle_complete_used_for_lambda = safe_divide(n_one_used, n_one_total),

    n_lambda_cells_n1_le_1 = sum(plugin_counts$n1 <= 1),
    n_lambda_cells_n1_le_2 = sum(plugin_counts$n1 <= 2),
    n_lambda_cells_n1_le_5 = sum(plugin_counts$n1 <= 5),
    pct_lambda_cells_n1_le_1 = mean(plugin_counts$n1 <= 1),
    pct_lambda_cells_n1_le_2 = mean(plugin_counts$n1 <= 2),
    pct_lambda_cells_n1_le_5 = mean(plugin_counts$n1 <= 5)
  )
}

run_one_lp_spec <- function(driver_data,
                            drinking_var,
                            threshold,
                            threshold_label,
                            analysis_sample,
                            imputation_value,
                            cell_spec,
                            cell_var,
                            cell_spec_label,
                            min_n2,
                            min_n2_label) {

  cat(
    "Running:", analysis_sample,
    "| threshold:", threshold_label,
    "| imp:", imputation_value,
    "| cell:", cell_spec_label,
    "|", min_n2_label, "\n"
  )

  crash_types <- make_crash_types(
    driver_data = driver_data,
    drinking_var = drinking_var,
    threshold_label = threshold_label,
    analysis_sample = analysis_sample,
    imputation_value = imputation_value
  )

  counts <- make_cell_counts(crash_types, cell_var = cell_var)

  crash_summary <- summarize_crash_types(crash_types) %>%
    mutate(
      analysis_sample = analysis_sample,
      imputation = imputation_value,
      threshold = threshold,
      threshold_label = threshold_label,
      cell_spec = cell_spec,
      cell_spec_label = cell_spec_label,
      min_n2 = min_n2,
      min_n2_label = min_n2_label,
      .before = 1
    )

  cell_diag <- summarize_cell_diagnostics(counts) %>%
    mutate(
      analysis_sample = analysis_sample,
      imputation = imputation_value,
      threshold = threshold,
      threshold_label = threshold_label,
      cell_spec = cell_spec,
      cell_spec_label = cell_spec_label,
      min_n2 = min_n2,
      min_n2_label = min_n2_label,
      .before = 1
    )

  theta_fit <- fit_theta_simple(
    counts = counts,
    min_n2 = min_n2,
    lower_risk = lower_risk,
    upper_risk = upper_risk,
    fd_step = fd_step
  )

  lambda_fit <- fit_lambda_simple(
    counts = counts,
    cell_prevalence = theta_fit$cell_prevalence,
    min_n2 = min_n2,
    lower_risk = lower_risk,
    upper_risk = upper_risk,
    fd_step = fd_step
  )

  estimates <- bind_cols(theta_fit$estimates, lambda_fit %>% select(-min_n2)) %>%
    mutate(
      analysis_sample = analysis_sample,
      imputation = imputation_value,
      threshold = threshold,
      threshold_label = threshold_label,
      drinking_var = drinking_var,
      cell_spec = cell_spec,
      cell_spec_label = cell_spec_label,
      min_n2 = min_n2,
      min_n2_label = min_n2_label,
      .before = 1
    )

  list(
    estimates = estimates,
    crash_summary = crash_summary,
    cell_diagnostics = cell_diag
  )
}

# ------------------------------------------------------------------------------
# 5. Raw measured-BAC, BAC_COMPLETE_95 state-years only
# ------------------------------------------------------------------------------

raw_results <- list()
raw_crash_summaries <- list()
raw_cell_diagnostics <- list()
raw_counter <- 1

raw_bac95_data <- driver_raw_one_row %>%
  filter(BAC_COMPLETE_95 == TRUE)

cat("\nRaw measured-BAC 95%-complete sample:\n")
cat("  Driver rows:", comma(nrow(raw_bac95_data)), "\n")
cat("  Distinct crashes:", comma(n_distinct(paste(raw_bac95_data$YEAR_FILE, raw_bac95_data$ST_CASE))), "\n\n")

for (t in seq_len(nrow(threshold_specs))) {
  for (c in seq_len(nrow(cell_specs))) {
    for (m in seq_len(nrow(min_n2_specs))) {

      res <- run_one_lp_spec(
        driver_data = raw_bac95_data,
        drinking_var = threshold_specs$raw_definition_var[t],
        threshold = threshold_specs$threshold[t],
        threshold_label = threshold_specs$threshold_label[t],
        analysis_sample = "raw_bac_complete_95_state_years",
        imputation_value = NA_integer_,
        cell_spec = cell_specs$cell_spec[c],
        cell_var = cell_specs$cell_var[c],
        cell_spec_label = cell_specs$cell_spec_label[c],
        min_n2 = min_n2_specs$min_n2[m],
        min_n2_label = min_n2_specs$min_n2_label[m]
      )

      raw_results[[raw_counter]] <- res$estimates
      raw_crash_summaries[[raw_counter]] <- res$crash_summary
      raw_cell_diagnostics[[raw_counter]] <- res$cell_diagnostics
      raw_counter <- raw_counter + 1
    }
  }
}

raw_estimates <- bind_rows(raw_results)
raw_crash_counts <- bind_rows(raw_crash_summaries)
raw_cell_diagnostics <- bind_rows(raw_cell_diagnostics)

# ------------------------------------------------------------------------------
# 6. MI all state-years: run separately for imputations 1,...,10
# ------------------------------------------------------------------------------

mi_results <- list()
mi_crash_summaries <- list()
mi_cell_diagnostics <- list()
mi_counter <- 1

for (imp in sort(unique(driver_mi_long$IMPUTATION))) {
  mi_data_imp <- driver_mi_long %>% filter(IMPUTATION == imp)

  for (t in seq_len(nrow(threshold_specs))) {
    for (c in seq_len(nrow(cell_specs))) {
      for (m in seq_len(nrow(min_n2_specs))) {

        res <- run_one_lp_spec(
          driver_data = mi_data_imp,
          drinking_var = threshold_specs$mi_definition_var[t],
          threshold = threshold_specs$threshold[t],
          threshold_label = threshold_specs$threshold_label[t],
          analysis_sample = "mi_all_state_years",
          imputation_value = imp,
          cell_spec = cell_specs$cell_spec[c],
          cell_var = cell_specs$cell_var[c],
          cell_spec_label = cell_specs$cell_spec_label[c],
          min_n2 = min_n2_specs$min_n2[m],
          min_n2_label = min_n2_specs$min_n2_label[m]
        )

        mi_results[[mi_counter]] <- res$estimates
        mi_crash_summaries[[mi_counter]] <- res$crash_summary
        mi_cell_diagnostics[[mi_counter]] <- res$cell_diagnostics
        mi_counter <- mi_counter + 1
      }
    }
  }
}

mi_imputation_estimates <- bind_rows(mi_results)
mi_crash_counts_by_imp <- bind_rows(mi_crash_summaries)
mi_cell_diagnostics_by_imp <- bind_rows(mi_cell_diagnostics)

# ------------------------------------------------------------------------------
# 7. Combine MI estimates using Rubin's rules
# ------------------------------------------------------------------------------

mi_theta_combined <- mi_imputation_estimates %>%
  group_by(
    analysis_sample, threshold, threshold_label,
    cell_spec, cell_spec_label, min_n2, min_n2_label
  ) %>%
  group_modify(~ combine_mi(est = .x$theta, se = .x$theta_se)) %>%
  ungroup() %>%
  mutate(parameter = "theta")

mi_lambda_combined <- mi_imputation_estimates %>%
  group_by(
    analysis_sample, threshold, threshold_label,
    cell_spec, cell_spec_label, min_n2, min_n2_label
  ) %>%
  group_modify(~ combine_mi(est = .x$lambda_plugin, se = .x$lambda_plugin_se)) %>%
  ungroup() %>%
  mutate(parameter = "lambda")

mi_combined <- bind_rows(mi_theta_combined, mi_lambda_combined) %>%
  arrange(parameter, threshold, match(cell_spec, cell_specs$cell_spec), min_n2)

# Add average MI diagnostics/counts for context.
mi_crash_counts_avg <- mi_crash_counts_by_imp %>%
  group_by(
    analysis_sample, threshold, threshold_label,
    cell_spec, cell_spec_label, min_n2, min_n2_label
  ) %>%
  summarize(
    DD = mean(DD, na.rm = TRUE),
    DS = mean(DS, na.rm = TRUE),
    SS = mean(SS, na.rm = TRUE),
    D = mean(D, na.rm = TRUE),
    S = mean(S, na.rm = TRUE),
    missing_or_ambiguous = mean(missing_or_ambiguous, na.rm = TRUE),
    total_crashes = mean(total_crashes, na.rm = TRUE),
    complete_two_vehicle = mean(complete_two_vehicle, na.rm = TRUE),
    complete_one_vehicle = mean(complete_one_vehicle, na.rm = TRUE),
    pct_missing_or_ambiguous = mean(pct_missing_or_ambiguous, na.rm = TRUE),
    .groups = "drop"
  )

# ------------------------------------------------------------------------------
# 8. Build raw-vs-MI comparison tables
# ------------------------------------------------------------------------------

raw_theta_for_compare <- raw_estimates %>%
  transmute(
    analysis_sample,
    threshold,
    threshold_label,
    cell_spec,
    cell_spec_label,
    min_n2,
    min_n2_label,
    parameter = "theta",
    estimate = theta,
    se_total = theta_se,
    ci_low = theta - 1.96 * theta_se,
    ci_high = theta + 1.96 * theta_se,
    status = theta_status,
    within_var = NA_real_,
    between_var = NA_real_,
    total_var = theta_se^2,
    m_total = NA_integer_,
    m_est = NA_integer_,
    m_used_for_W = NA_integer_
  )

raw_lambda_for_compare <- raw_estimates %>%
  transmute(
    analysis_sample,
    threshold,
    threshold_label,
    cell_spec,
    cell_spec_label,
    min_n2,
    min_n2_label,
    parameter = "lambda",
    estimate = lambda_plugin,
    se_total = lambda_plugin_se,
    ci_low = lambda_plugin - 1.96 * lambda_plugin_se,
    ci_high = lambda_plugin + 1.96 * lambda_plugin_se,
    status = lambda_plugin_status,
    within_var = NA_real_,
    between_var = NA_real_,
    total_var = lambda_plugin_se^2,
    m_total = NA_integer_,
    m_est = NA_integer_,
    m_used_for_W = NA_integer_
  )

mi_for_compare <- mi_combined %>%
  mutate(
    status = case_when(
      m_est == m_total & m_used_for_W == m_total ~ "rubin_ok_all_10",
      m_est > 0 ~ "rubin_partial",
      TRUE ~ "rubin_no_finite_estimates"
    )
  ) %>%
  select(
    analysis_sample,
    threshold,
    threshold_label,
    cell_spec,
    cell_spec_label,
    min_n2,
    min_n2_label,
    parameter,
    estimate,
    se_total,
    ci_low,
    ci_high,
    status,
    within_var,
    between_var,
    total_var,
    m_total,
    m_est,
    m_used_for_W
  )

comparison_long <- bind_rows(
  raw_theta_for_compare,
  raw_lambda_for_compare,
  mi_for_compare
) %>%
  arrange(parameter, threshold, match(cell_spec, cell_specs$cell_spec), min_n2, analysis_sample)

comparison_differences <- comparison_long %>%
  select(
    threshold, threshold_label, cell_spec, cell_spec_label,
    min_n2, min_n2_label, parameter, analysis_sample, estimate, se_total
  ) %>%
  pivot_wider(
    names_from = analysis_sample,
    values_from = c(estimate, se_total)
  ) %>%
  mutate(
    mi_minus_raw = estimate_mi_all_state_years - estimate_raw_bac_complete_95_state_years,
    mi_div_raw = estimate_mi_all_state_years / estimate_raw_bac_complete_95_state_years
  )

comparison_long <- comparison_long %>%
  left_join(
    comparison_differences %>%
      select(
        threshold, cell_spec, min_n2, parameter,
        mi_minus_raw, mi_div_raw
      ),
    by = c("threshold", "cell_spec", "min_n2", "parameter")
  )

# Wider display tables.
comparison_wide_theta <- comparison_long %>%
  filter(parameter == "theta") %>%
  mutate(
    estimate_with_se = format_est_se_flag(estimate, se_total, status),
    row_label = paste(threshold_label, min_n2_label, analysis_sample, sep = " | ")
  ) %>%
  select(row_label, threshold_label, min_n2_label, analysis_sample, cell_spec_label, estimate_with_se) %>%
  pivot_wider(names_from = cell_spec_label, values_from = estimate_with_se)

comparison_wide_lambda <- comparison_long %>%
  filter(parameter == "lambda") %>%
  mutate(
    estimate_with_se = format_est_se_flag(estimate, se_total, status),
    row_label = paste(threshold_label, min_n2_label, analysis_sample, sep = " | ")
  ) %>%
  select(row_label, threshold_label, min_n2_label, analysis_sample, cell_spec_label, estimate_with_se) %>%
  pivot_wider(names_from = cell_spec_label, values_from = estimate_with_se)

crash_type_counts_raw_vs_mi <- bind_rows(
  raw_crash_counts %>% mutate(source_type = "raw"),
  mi_crash_counts_avg %>% mutate(source_type = "mi_avg_across_imputations")
) %>%
  arrange(threshold, match(cell_spec, cell_specs$cell_spec), min_n2, source_type)

# ------------------------------------------------------------------------------
# 9. Save outputs
# ------------------------------------------------------------------------------

write_csv(raw_estimates, file.path(out_dir, "raw_bac95_lp_mle_estimates.csv"))
write_csv(raw_crash_counts, file.path(out_dir, "raw_bac95_crash_type_counts.csv"))
write_csv(raw_cell_diagnostics, file.path(out_dir, "raw_bac95_cell_diagnostics.csv"))

write_csv(mi_imputation_estimates, file.path(out_dir, "mi_lp_mle_imputation_estimates.csv"))
write_csv(mi_combined, file.path(out_dir, "mi_lp_mle_rubin_combined.csv"))
write_csv(mi_crash_counts_by_imp, file.path(out_dir, "mi_crash_type_counts_by_imputation.csv"))
write_csv(mi_crash_counts_avg, file.path(out_dir, "mi_crash_type_counts_avg_across_imputations.csv"))
write_csv(mi_cell_diagnostics_by_imp, file.path(out_dir, "mi_cell_diagnostics_by_imputation.csv"))

write_csv(comparison_long, file.path(out_dir, "lp_mle_raw_vs_mi_comparison.csv"))
write_csv(comparison_differences, file.path(out_dir, "lp_mle_raw_vs_mi_differences.csv"))
write_csv(comparison_wide_theta, file.path(out_dir, "lp_mle_raw_vs_mi_wide_theta.csv"))
write_csv(comparison_wide_lambda, file.path(out_dir, "lp_mle_raw_vs_mi_wide_lambda.csv"))
write_csv(crash_type_counts_raw_vs_mi, file.path(out_dir, "crash_type_counts_raw_vs_mi.csv"))

# ------------------------------------------------------------------------------
# 10. Quick plots
# ------------------------------------------------------------------------------

plot_theta <- comparison_long %>%
  filter(parameter == "theta") %>%
  mutate(
    analysis_sample = recode(
      analysis_sample,
      raw_bac_complete_95_state_years = "Raw measured BAC, 95% state-years",
      mi_all_state_years = "MI BAC, all state-years"
    ),
    cell_spec_label = factor(cell_spec_label, levels = cell_specs$cell_spec_label),
    min_n2_label = factor(min_n2_label, levels = min_n2_specs$min_n2_label)
  ) %>%
  ggplot(aes(x = cell_spec_label, y = estimate, group = analysis_sample, linetype = analysis_sample)) +
  geom_line() +
  geom_point(size = 2) +
  facet_grid(threshold_label ~ min_n2_label, scales = "free_y") +
  coord_flip() +
  labs(
    title = "Theta: Raw BAC 95%-Complete State-Years vs MI All State-Years",
    x = "Cell specification",
    y = expression(theta),
    linetype = "Analysis sample"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  file.path(fig_dir, "theta_raw_bac95_vs_mi.png"),
  plot_theta,
  width = 14,
  height = 9,
  dpi = 300
)

plot_lambda <- comparison_long %>%
  filter(parameter == "lambda") %>%
  mutate(
    analysis_sample = recode(
      analysis_sample,
      raw_bac_complete_95_state_years = "Raw measured BAC, 95% state-years",
      mi_all_state_years = "MI BAC, all state-years"
    ),
    cell_spec_label = factor(cell_spec_label, levels = cell_specs$cell_spec_label),
    min_n2_label = factor(min_n2_label, levels = min_n2_specs$min_n2_label)
  ) %>%
  ggplot(aes(x = cell_spec_label, y = estimate, group = analysis_sample, linetype = analysis_sample)) +
  geom_line() +
  geom_point(size = 2) +
  facet_grid(threshold_label ~ min_n2_label, scales = "free_y") +
  coord_flip() +
  labs(
    title = "Lambda: Raw BAC 95%-Complete State-Years vs MI All State-Years",
    x = "Cell specification",
    y = expression(lambda),
    linetype = "Analysis sample"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  file.path(fig_dir, "lambda_raw_bac95_vs_mi.png"),
  plot_lambda,
  width = 14,
  height = 9,
  dpi = 300
)

