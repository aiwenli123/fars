################################################################################
# 04_lp_mle_period.R
#
# Purpose:
#   Simpler Levitt-Porter-style MLE with COVID-period cells:
#
#   1. Estimate theta using only two-car crashes: DD, DS, SS.
#   2. Estimate lambda in a second-stage plug-in model using one-car crashes
#      and N_hat from the two-car model.
#   3. Replace YEAR-based cells with PERIOD-based cells:
#        pre_covid  = 2015-2019
#        covid      = 2020-2021
#        post_covid = 2022-2024
#   4. Run each specification twice:
#        min_n2 = 1  : no real two-car cell-size filter except n2 > 0
#        min_n2 = 5  : keep only cells with at least 5 complete two-car crashes
#
# Main input:
#   data/processed/driver_analysis_lp_2015_2024.csv
#
# Main outputs:
#   data/processed/script04_simple_mle_period/...
################################################################################

rm(list = ls())

library(tidyverse)
library(scales)

# ------------------------------------------------------------------------------
# 0. User options
# ------------------------------------------------------------------------------

driver_path <- "data/processed/driver_analysis_lp_2015_2024.csv"

out_dir <- "data/processed/script04_simple_mle_period"
fig_dir <- "figures/script04_simple_mle_period"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Bounds for theta and lambda.
# If many estimates hit upper_risk, that is a diagnostic, not just a coding issue.
lower_risk <- 1
upper_risk <- 200

# Finite difference step size for curvature SEs on log scale.
fd_step <- 1e-4

# Cell-size variants for two-car theta estimation.
# min_n2 = 1 means all cells with at least one complete two-car crash are used.
# min_n2 = 5 means only cells with at least five complete two-car crashes are used.
min_n2_specs <- tibble(
  min_n2 = c(1, 5),
  min_n2_label = c("No n2 filter", "n2 >= 5")
)

# ------------------------------------------------------------------------------
# 0.5 Helper functions
# ------------------------------------------------------------------------------

safe_divide <- function(num, den) {
  ifelse(is.na(den) | den == 0, NA_real_, num / den)
}

format_est_se_flag <- function(est, se, status, digits = 3) {
  case_when(
    is.na(est) ~ NA_character_,
    status == "ok" & !is.na(se) ~ paste0(round(est, digits), " (", round(se, digits), ")"),
    TRUE ~ paste0(round(est, digits), " [", status, "]")
  )
}

curvature_se_1d <- function(objective, x_hat, est_hat, h = 1e-4) {
  # objective is negative log likelihood in terms of log parameter.
  # second derivative of negative log likelihood should be positive at minimum.
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

# ------------------------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------------------------

if (!file.exists(driver_path)) {
  stop("Missing input file: ", driver_path, "\nPlease run 01_clean_data.R first.")
}

driver_lp <- read_csv(driver_path, show_col_types = FALSE) %>%
  rename_with(toupper)

cat("Loaded driver-level LP sample:\n")
cat("  Driver rows:", comma(nrow(driver_lp)), "\n")
cat("  Distinct crashes:", comma(n_distinct(paste(driver_lp$YEAR_FILE, driver_lp$ST_CASE))), "\n\n")

# ------------------------------------------------------------------------------
# 2. Recreate key variables
# ------------------------------------------------------------------------------

driver_lp <- driver_lp %>%
  mutate(
    period = case_when(
      YEAR %in% 2015:2019 ~ "pre_covid",
      YEAR %in% 2020:2021 ~ "covid",
      YEAR %in% 2022:2024 ~ "post_covid",
      TRUE ~ NA_character_
    ),
    period = factor(period, levels = c("pre_covid", "covid", "post_covid")),
    
    police_drinking_complete = case_when(
      DRINKING == 1 ~ 1,
      DRINKING == 0 ~ 0,
      TRUE ~ NA_real_
    ),
    
    police_drinking_unknown_as_no = case_when(
      DRINKING == 1 ~ 1,
      TRUE ~ 0
    ),
    
    vehicle_dr_drink_complete = case_when(
      DR_DRINK == 1 ~ 1,
      DR_DRINK == 0 ~ 0,
      TRUE ~ NA_real_
    ),
    
    bac_valid = !is.na(ALC_RES) & ALC_RES >= 0 & ALC_RES < 995,
    bac_value = if_else(bac_valid, ALC_RES / 1000, NA_real_),
    
    bac_ge_08 = case_when(
      bac_valid & ALC_RES >= 80 ~ 1,
      bac_valid & ALC_RES < 80 ~ 0,
      TRUE ~ NA_real_
    ),
    
    bac_ge_10 = case_when(
      bac_valid & ALC_RES >= 100 ~ 1,
      bac_valid & ALC_RES < 100 ~ 0,
      TRUE ~ NA_real_
    ),
    
    weekend_lp = as.logical(WEEKEND_LP),
    
    # Period-based cell definitions.
    cell_all = "all_data",
    cell_hour = paste(HOUR, sep = "_"),
    cell_hour_period = paste(HOUR, period, sep = "_"),
    cell_hour_period_weekend = paste(HOUR, period, weekend_lp, sep = "_"),
    cell_hour_state_period = paste(HOUR, STATE, period, sep = "_"),
    cell_hour_state_period_weekend = paste(HOUR, STATE, period, weekend_lp, sep = "_")
  ) %>%
  filter(!is.na(period))

cat("Period counts, driver rows:\n")
print(driver_lp %>% count(period))
cat("\n")

# ------------------------------------------------------------------------------
# 2.5 Recreate BAC_COMPLETE_95 if needed
# ------------------------------------------------------------------------------

if (!"BAC_COMPLETE_95" %in% names(driver_lp)) {
  cat("BAC_COMPLETE_95 not found; recreating from driver-level file.\n")
  
  bac_completeness_state_year <- driver_lp %>%
    group_by(YEAR, STATE) %>%
    summarize(
      n_police_drinking_yes = sum(police_drinking_complete == 1, na.rm = TRUE),
      n_police_drinking_yes_bac_valid =
        sum(police_drinking_complete == 1 & bac_valid, na.rm = TRUE),
      bac_complete_among_police_drinkers = case_when(
        n_police_drinking_yes > 0 ~
          n_police_drinking_yes_bac_valid / n_police_drinking_yes,
        TRUE ~ NA_real_
      ),
      BAC_COMPLETE_95 = bac_complete_among_police_drinkers >= 0.95,
      .groups = "drop"
    )
  
  driver_lp <- driver_lp %>%
    left_join(
      bac_completeness_state_year %>%
        select(YEAR, STATE, BAC_COMPLETE_95),
      by = c("YEAR", "STATE")
    )
}

driver_lp <- driver_lp %>%
  mutate(BAC_COMPLETE_95 = if_else(is.na(BAC_COMPLETE_95), FALSE, BAC_COMPLETE_95))

# ------------------------------------------------------------------------------
# 2.6 Audit how many crashes are removed by BAC_COMPLETE_95 state-year filter
# ------------------------------------------------------------------------------

bac_complete_95_audit <- driver_lp %>%
  mutate(
    crash_id = paste(YEAR_FILE, ST_CASE, sep = "_")
  ) %>%
  group_by(BAC_COMPLETE_95) %>%
  summarize(
    n_driver_rows = n(),
    n_crashes = n_distinct(crash_id),
    n_one_vehicle_crashes = n_distinct(crash_id[VE_FORMS == 1]),
    n_two_vehicle_crashes = n_distinct(crash_id[VE_FORMS == 2]),
    .groups = "drop"
  ) %>%
  mutate(
    bac_sample_status = if_else(
      BAC_COMPLETE_95,
      "kept: BAC_COMPLETE_95 state-years",
      "removed: below 95% BAC-complete state-years"
    )
  ) %>%
  select(bac_sample_status, BAC_COMPLETE_95, everything())

bac_complete_95_audit_overall <- bac_complete_95_audit %>%
  summarize(
    n_driver_rows_total = sum(n_driver_rows),
    n_driver_rows_kept = sum(n_driver_rows[BAC_COMPLETE_95]),
    n_driver_rows_removed = sum(n_driver_rows[!BAC_COMPLETE_95]),
    pct_driver_rows_removed = n_driver_rows_removed / n_driver_rows_total,
    
    n_crashes_total = sum(n_crashes),
    n_crashes_kept = sum(n_crashes[BAC_COMPLETE_95]),
    n_crashes_removed = sum(n_crashes[!BAC_COMPLETE_95]),
    pct_crashes_removed = n_crashes_removed / n_crashes_total,
    
    n_one_vehicle_crashes_total = sum(n_one_vehicle_crashes),
    n_one_vehicle_crashes_kept = sum(n_one_vehicle_crashes[BAC_COMPLETE_95]),
    n_one_vehicle_crashes_removed = sum(n_one_vehicle_crashes[!BAC_COMPLETE_95]),
    pct_one_vehicle_crashes_removed =
      n_one_vehicle_crashes_removed / n_one_vehicle_crashes_total,
    
    n_two_vehicle_crashes_total = sum(n_two_vehicle_crashes),
    n_two_vehicle_crashes_kept = sum(n_two_vehicle_crashes[BAC_COMPLETE_95]),
    n_two_vehicle_crashes_removed = sum(n_two_vehicle_crashes[!BAC_COMPLETE_95]),
    pct_two_vehicle_crashes_removed =
      n_two_vehicle_crashes_removed / n_two_vehicle_crashes_total
  )

bac_complete_95_audit_by_year_state <- driver_lp %>%
  mutate(
    crash_id = paste(YEAR_FILE, ST_CASE, sep = "_")
  ) %>%
  group_by(YEAR, STATE, STATENAME, BAC_COMPLETE_95) %>%
  summarize(
    n_driver_rows = n(),
    n_crashes = n_distinct(crash_id),
    n_one_vehicle_crashes = n_distinct(crash_id[VE_FORMS == 1]),
    n_two_vehicle_crashes = n_distinct(crash_id[VE_FORMS == 2]),
    .groups = "drop"
  ) %>%
  arrange(YEAR, STATE)

cat("\nBAC_COMPLETE_95 audit:\n")
print(bac_complete_95_audit, width = Inf)

cat("\nOverall crashes removed by BAC_COMPLETE_95 filter:\n")
print(bac_complete_95_audit_overall)

write_csv(
  bac_complete_95_audit,
  file.path(out_dir, "bac_complete_95_audit.csv")
)

write_csv(
  bac_complete_95_audit_overall,
  file.path(out_dir, "bac_complete_95_audit_overall.csv")
)

write_csv(
  bac_complete_95_audit_by_year_state,
  file.path(out_dir, "bac_complete_95_audit_by_year_state.csv")
)

# ------------------------------------------------------------------------------
# 3. Specifications
# ------------------------------------------------------------------------------

alcohol_definitions <- tibble(
  alcohol_definition = c(
    "police_drinking_complete",
    "police_drinking_unknown_as_no",
    "vehicle_dr_drink_complete",
    "bac_ge_08",
    "bac_ge_10"
  ),
  alcohol_definition_label = c(
    "Police drinking, complete case",
    "Police drinking, unknown as no",
    "Vehicle DR_DRINK, complete case",
    "Observed BAC >= 0.08",
    "Observed BAC >= 0.10"
  )
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
    "cell_all",
    "cell_hour",
    "cell_hour_period",
    "cell_hour_period_weekend",
    "cell_hour_state_period",
    "cell_hour_state_period_weekend"
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

sample_variants <- tibble(
  sample_variant = c(
    "all_state_years",
    "bac_complete_95_state_years"
  ),
  sample_variant_label = c(
    "All state-years",
    "BAC-complete 95% state-years"
  )
)

# ------------------------------------------------------------------------------
# 4. Crash classification
# ------------------------------------------------------------------------------

make_crash_types_for_definition <- function(driver_data,
                                            definition_var,
                                            definition_label,
                                            sample_variant) {
  
  df <- driver_data
  
  if (sample_variant == "bac_complete_95_state_years") {
    df <- df %>% filter(BAC_COMPLETE_95)
  }
  
  df %>%
    mutate(drinking_status = .data[[definition_var]]) %>%
    group_by(
      YEAR_FILE,
      YEAR,
      period,
      STATE,
      STATENAME,
      ST_CASE,
      VE_FORMS,
      HOUR,
      weekend_lp,
      cell_all,
      cell_hour,
      cell_hour_period,
      cell_hour_period_weekend,
      cell_hour_state_period,
      cell_hour_state_period_weekend
    ) %>%
    summarize(
      n_drivers = n(),
      n_drinking = sum(drinking_status == 1, na.rm = TRUE),
      n_sober = sum(drinking_status == 0, na.rm = TRUE),
      n_missing = sum(is.na(drinking_status)),
      .groups = "drop"
    ) %>%
    mutate(
      alcohol_definition = definition_var,
      alcohol_definition_label = definition_label,
      sample_variant = sample_variant,
      
      crash_type = case_when(
        VE_FORMS == 2 & n_missing == 0 & n_drinking == 2 & n_sober == 0 ~ "DD",
        VE_FORMS == 2 & n_missing == 0 & n_drinking == 1 & n_sober == 1 ~ "DS",
        VE_FORMS == 2 & n_missing == 0 & n_drinking == 0 & n_sober == 2 ~ "SS",
        
        VE_FORMS == 1 & n_missing == 0 & n_drinking == 1 & n_sober == 0 ~ "D",
        VE_FORMS == 1 & n_missing == 0 & n_drinking == 0 & n_sober == 1 ~ "S",
        
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
      theta_cell_closed = closed_form_theta_cell(DD, DS, SS),
      
      sparse_bucket_n2 = case_when(
        n2 == 0 ~ "n2 = 0",
        n2 <= 1 ~ "n2 = 1",
        n2 <= 2 ~ "n2 = 2",
        n2 <= 5 ~ "3 <= n2 <= 5",
        n2 <= 10 ~ "6 <= n2 <= 10",
        TRUE ~ "n2 > 10"
      )
    ) %>%
    filter(total_complete_case > 0)
}

summarize_cell_diagnostics <- function(counts) {
  
  counts_two_available <- counts %>%
    filter(n2 > 0)
  
  if (nrow(counts_two_available) == 0) {
    return(tibble(
      n_cells_total = nrow(counts),
      n_theta_cells_available = 0,
      n_two_vehicle_complete_available = 0,
      n_two_vehicle_all = sum(counts$total_two_vehicle),
      pct_two_vehicle_missing_or_ambiguous =
        safe_divide(sum(counts$missing_two_vehicle), sum(counts$total_two_vehicle)),
      n_theta_cells_no_DD_available = NA_integer_,
      n_theta_cells_no_DS_available = NA_integer_,
      n_theta_cells_no_SS_available = NA_integer_,
      pct_theta_cells_n2_le_1_available = NA_real_,
      pct_theta_cells_n2_le_2_available = NA_real_,
      pct_theta_cells_n2_le_5_available = NA_real_,
      median_theta_cell_closed = NA_real_,
      p75_theta_cell_closed = NA_real_,
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
    p75_theta_cell_closed = quantile(counts_two_available$theta_cell_closed, 0.75, na.rm = TRUE),
    p90_theta_cell_closed = quantile(counts_two_available$theta_cell_closed, 0.90, na.rm = TRUE),
    p95_theta_cell_closed = quantile(counts_two_available$theta_cell_closed, 0.95, na.rm = TRUE)
  )
}

# ------------------------------------------------------------------------------
# 5. Two-car-only theta MLE
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
    if (return_eta) {
      return(list(loglik = -Inf, theta = theta, eta = NA, N = NA))
    }
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
  
  counts_two_available <- counts %>%
    filter(n2 > 0)
  
  counts_two <- counts_two_available %>%
    filter(n2 >= min_n2)
  
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
        n_theta_cells_excluded_by_min_n2 = nrow(counts_two_available),
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
    
    implied_frac_drinking_weighted_two_car =
      weighted.mean(
        cell_prevalence$frac_drinking_hat,
        w = cell_prevalence$n2,
        na.rm = TRUE
      ),
    
    n_theta_cells_n2_le_1_available = sum(counts_two_available$n2 <= 1),
    n_theta_cells_n2_le_2_available = sum(counts_two_available$n2 <= 2),
    n_theta_cells_n2_le_5_available = sum(counts_two_available$n2 <= 5),
    pct_theta_cells_n2_le_1_available = mean(counts_two_available$n2 <= 1),
    pct_theta_cells_n2_le_2_available = mean(counts_two_available$n2 <= 2),
    pct_theta_cells_n2_le_5_available = mean(counts_two_available$n2 <= 5),
    
    n_theta_cells_no_DD_available = sum(counts_two_available$DD == 0),
    n_theta_cells_no_DS_available = sum(counts_two_available$DS == 0),
    n_theta_cells_no_SS_available = sum(counts_two_available$SS == 0),
    
    n_theta_cells_no_DD_used = sum(counts_two$DD == 0),
    n_theta_cells_no_DS_used = sum(counts_two$DS == 0),
    n_theta_cells_no_SS_used = sum(counts_two$SS == 0),
    
    median_theta_cell_closed_available = median(counts_two_available$theta_cell_closed, na.rm = TRUE),
    p90_theta_cell_closed_available = quantile(counts_two_available$theta_cell_closed, 0.90, na.rm = TRUE),
    p95_theta_cell_closed_available = quantile(counts_two_available$theta_cell_closed, 0.95, na.rm = TRUE),
    
    median_theta_cell_closed_used = median(counts_two$theta_cell_closed, na.rm = TRUE),
    p90_theta_cell_closed_used = quantile(counts_two$theta_cell_closed, 0.90, na.rm = TRUE),
    p95_theta_cell_closed_used = quantile(counts_two$theta_cell_closed, 0.95, na.rm = TRUE)
  )
  
  list(
    estimates = estimates,
    cell_prevalence = cell_prevalence
  )
}

# ------------------------------------------------------------------------------
# 6. Plug-in lambda MLE
# ------------------------------------------------------------------------------

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
      cell_prevalence %>%
        select(cell_id, N_hat_drinking_to_sober),
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

# ------------------------------------------------------------------------------
# 7. Run one specification
# ------------------------------------------------------------------------------

run_one_spec <- function(driver_data,
                         definition_var,
                         definition_label,
                         sample_variant,
                         cell_spec,
                         cell_var,
                         cell_spec_label,
                         min_n2,
                         min_n2_label) {
  
  cat(
    "\nRunning:",
    definition_label,
    "| sample:", sample_variant,
    "| cell:", cell_spec_label,
    "|", min_n2_label,
    "\n"
  )
  
  crash_types <- make_crash_types_for_definition(
    driver_data = driver_data,
    definition_var = definition_var,
    definition_label = definition_label,
    sample_variant = sample_variant
  )
  
  counts <- make_cell_counts(crash_types, cell_var = cell_var)
  
  cell_diag <- summarize_cell_diagnostics(counts) %>%
    mutate(
      alcohol_definition = definition_var,
      alcohol_definition_label = definition_label,
      sample_variant = sample_variant,
      cell_spec = cell_spec,
      cell_spec_label = cell_spec_label,
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
  
  # Avoid duplicate min_n2 column from theta and lambda.
  lambda_fit_clean <- lambda_fit %>%
    select(-min_n2)
  
  estimates <- bind_cols(theta_fit$estimates, lambda_fit_clean) %>%
    mutate(
      alcohol_definition = definition_var,
      alcohol_definition_label = definition_label,
      sample_variant = sample_variant,
      cell_spec = cell_spec,
      cell_spec_label = cell_spec_label,
      min_n2 = min_n2,
      min_n2_label = min_n2_label,
      .before = 1
    )
  
  cell_prevalence <- theta_fit$cell_prevalence %>%
    mutate(
      alcohol_definition = definition_var,
      alcohol_definition_label = definition_label,
      sample_variant = sample_variant,
      cell_spec = cell_spec,
      cell_spec_label = cell_spec_label,
      min_n2 = min_n2,
      min_n2_label = min_n2_label,
      .before = 1
    )
  
  list(
    estimates = estimates,
    cell_prevalence = cell_prevalence,
    cell_diagnostics = cell_diag
  )
}

# ------------------------------------------------------------------------------
# 8. Run Table 2-style models
# ------------------------------------------------------------------------------

results <- list()
cell_prevalence_results <- list()
cell_diagnostic_results <- list()

counter <- 1

for (sv in sample_variants$sample_variant) {
  for (d in seq_len(nrow(alcohol_definitions))) {
    for (c in seq_len(nrow(cell_specs))) {
      for (m in seq_len(nrow(min_n2_specs))) {
        
        res <- run_one_spec(
          driver_data = driver_lp,
          definition_var = alcohol_definitions$alcohol_definition[d],
          definition_label = alcohol_definitions$alcohol_definition_label[d],
          sample_variant = sv,
          cell_spec = cell_specs$cell_spec[c],
          cell_var = cell_specs$cell_var[c],
          cell_spec_label = cell_specs$cell_spec_label[c],
          min_n2 = min_n2_specs$min_n2[m],
          min_n2_label = min_n2_specs$min_n2_label[m]
        )
        
        results[[counter]] <- res$estimates
        cell_prevalence_results[[counter]] <- res$cell_prevalence
        cell_diagnostic_results[[counter]] <- res$cell_diagnostics
        
        counter <- counter + 1
      }
    }
  }
}

table2_long <- bind_rows(results) %>%
  left_join(sample_variants, by = "sample_variant") %>%
  arrange(
    sample_variant,
    alcohol_definition,
    match(cell_spec, cell_specs$cell_spec),
    min_n2
  )

table2_cells <- bind_rows(cell_prevalence_results) %>%
  left_join(sample_variants, by = "sample_variant") %>%
  arrange(
    sample_variant,
    alcohol_definition,
    match(cell_spec, cell_specs$cell_spec),
    min_n2,
    cell_id
  )

cell_diagnostics <- bind_rows(cell_diagnostic_results) %>%
  distinct() %>%
  left_join(sample_variants, by = "sample_variant") %>%
  arrange(
    sample_variant,
    alcohol_definition,
    match(cell_spec, cell_specs$cell_spec)
  )

# ------------------------------------------------------------------------------
# 9. Save outputs
# ------------------------------------------------------------------------------

write_csv(
  table2_long,
  file.path(out_dir, "lp_period_table2_long.csv")
)

write_csv(
  table2_cells,
  file.path(out_dir, "lp_period_cell_prevalence_table2.csv")
)

write_csv(
  cell_diagnostics,
  file.path(out_dir, "lp_period_cell_diagnostics.csv")
)

table2_wide_theta <- table2_long %>%
  mutate(
    theta_with_se = format_est_se_flag(theta, theta_se, theta_status),
    row_label = paste(sample_variant_label, alcohol_definition_label, min_n2_label, sep = " | ")
  ) %>%
  select(
    row_label,
    sample_variant_label,
    alcohol_definition_label,
    min_n2_label,
    cell_spec_label,
    theta_with_se
  ) %>%
  pivot_wider(
    names_from = cell_spec_label,
    values_from = theta_with_se
  )

table2_wide_lambda <- table2_long %>%
  mutate(
    lambda_with_se = format_est_se_flag(lambda_plugin, lambda_plugin_se, lambda_plugin_status),
    row_label = paste(sample_variant_label, alcohol_definition_label, min_n2_label, sep = " | ")
  ) %>%
  select(
    row_label,
    sample_variant_label,
    alcohol_definition_label,
    min_n2_label,
    cell_spec_label,
    lambda_with_se
  ) %>%
  pivot_wider(
    names_from = cell_spec_label,
    values_from = lambda_with_se
  )

write_csv(
  table2_wide_theta,
  file.path(out_dir, "lp_period_table2_wide_theta.csv")
)

write_csv(
  table2_wide_lambda,
  file.path(out_dir, "lp_period_table2_wide_lambda.csv")
)

# ------------------------------------------------------------------------------
# 10. Plots
# ------------------------------------------------------------------------------

plot_theta <- table2_long %>%
  mutate(
    cell_spec_label = factor(cell_spec_label, levels = cell_specs$cell_spec_label),
    min_n2_label = factor(min_n2_label, levels = min_n2_specs$min_n2_label)
  ) %>%
  ggplot(aes(
    x = cell_spec_label,
    y = theta,
    group = min_n2_label,
    linetype = min_n2_label
  )) +
  geom_line() +
  geom_point(size = 2) +
  facet_wrap(~ sample_variant_label + alcohol_definition_label, scales = "free_y") +
  coord_flip() +
  labs(
    title = "Period-Based Two-Car LP MLE: Theta",
    subtitle = "Theta estimated from DD/DS/SS counts; period replaces year; compared with and without n2 >= 5 filter",
    x = "Cell specification",
    y = expression(theta),
    linetype = "Cell-size filter"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  file.path(fig_dir, "lp_period_theta_table2.png"),
  plot_theta,
  width = 30,
  height = 10,
  dpi = 300
)

plot_lambda <- table2_long %>%
  mutate(
    cell_spec_label = factor(cell_spec_label, levels = cell_specs$cell_spec_label),
    min_n2_label = factor(min_n2_label, levels = min_n2_specs$min_n2_label)
  ) %>%
  ggplot(aes(
    x = cell_spec_label,
    y = lambda_plugin,
    group = min_n2_label,
    linetype = min_n2_label
  )) +
  geom_line() +
  geom_point(size = 2) +
  facet_wrap(~ sample_variant_label + alcohol_definition_label, scales = "free_y") +
  coord_flip() +
  labs(
    title = "Period-Based Plug-in Lambda",
    subtitle = "Lambda estimated from one-car D/S counts using N_hat from two-car model",
    x = "Cell specification",
    y = expression(lambda),
    linetype = "Cell-size filter"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  file.path(fig_dir, "lp_period_lambda_table2.png"),
  plot_lambda,
  width = 30,
  height = 10,
  dpi = 300
)

plot_diagnostics <- table2_long %>%
  mutate(
    cell_spec_label = factor(cell_spec_label, levels = cell_specs$cell_spec_label),
    min_n2_label = factor(min_n2_label, levels = min_n2_specs$min_n2_label)
  ) %>%
  ggplot(aes(
    x = cell_spec_label,
    y = pct_two_vehicle_complete_used,
    group = min_n2_label,
    linetype = min_n2_label
  )) +
  geom_line() +
  geom_point(size = 2) +
  facet_wrap(~ sample_variant_label + alcohol_definition_label, scales = "free_y") +
  coord_flip() +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "Share of Complete Two-Car Crashes Used",
    subtitle = "Shows how much data is retained under the min_n2 filter",
    x = "Cell specification",
    y = "Percent of complete two-car crashes used",
    linetype = "Cell-size filter"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  file.path(fig_dir, "lp_period_pct_two_vehicle_used.png"),
  plot_diagnostics,
  width = 30,
  height = 10,
  dpi = 300
)


#######################################
### look at results
#######################################

rm(list=ls())
library(tidyverse)

driver_path <- "data/processed/driver_analysis_lp_2015_2024.csv"

out_dir <- "data/processed/script04_simple_mle_period"
fig_dir <- "figures/script04_simple_mle_period"

diagnostics <- read_csv(file.path(out_dir, "lp_period_cell_diagnostics.csv"))
table2_long <- read.csv(file.path(out_dir, "lp_period_table2_long.csv"))
table2_wide_theta <- read.csv(file.path(out_dir, "lp_period_table2_wide_theta.csv"))
table2_wide_lambda <- read.csv(file.path(out_dir, "lp_period_table2_wide_lambda.csv"))
cell_prevalence <- read.csv(file.path(out_dir, "lp_period_cell_prevalence_table2.csv")) # size of each cell

driver_lp <- read_csv(driver_path, show_col_types = FALSE) %>% rename_with(toupper)


# original crash type counts
crash_type_counts <- read.csv("data/processed/lp_crash_type_counts_2015_2024.csv")
onecar_crash_type_counts <- read.csv("data/processed/lp_one_vehicle_counts_2015_2024.csv")

# audit of BAC_COMPLETE_95 filter
bac_95_audit <- read.csv(file.path(out_dir, "bac_complete_95_audit.csv"))
bac_95_audit_overall <- read.csv(file.path(out_dir, "bac_complete_95_audit_overall.csv"))
bac_95_audit_by_year_state <- read.csv(file.path(out_dir, "bac_complete_95_audit_by_year_state.csv"))

var_names <- read.csv("data/processed/fars_variable_reference_for_replication.csv")
