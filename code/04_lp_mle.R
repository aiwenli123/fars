################################################################################
# 04_lp_mle.R
#
# Purpose:
#   Estimate full Levitt-Porter-style MLE models for drinking-driver relative risk.
#
#   This script estimates:
#     theta  = relative two-car fatal crash risk for drinking vs sober drivers
#     lambda = relative one-car fatal crash risk for drinking vs sober drivers
#
#   while allowing the drinking/sober driver prevalence ratio N to vary across
#   increasingly fine equal-mixing cells.
#
# Robustness upgrades in this version:
#   - Uses nlminb() as the main optimizer, with optim(Nelder-Mead) fallback.
#   - Only reports Hessian/delta-method SEs when the optimizer converges cleanly.
#   - Checks Hessian positive-definiteness and condition number.
#   - Flags suspiciously tiny SEs instead of reporting misleading precision.
#   - Adds sparse-cell diagnostics to help identify unstable specifications.
#   - Optionally computes profile-likelihood CIs for Table 2.
#
# Main input:
#   data/processed/driver_analysis_lp_2015_2024.csv
#
# Main outputs:
#   data/processed/script04_mle/lp_mle_table2_long.csv
#   data/processed/script04_mle/lp_mle_table2_wide_theta.csv
#   data/processed/script04_mle/lp_mle_table2_wide_lambda.csv
#   data/processed/script04_mle/lp_mle_table3_by_year.csv
#   data/processed/script04_mle/lp_mle_table4_by_hour.csv
#   data/processed/script04_mle/lp_mle_cell_prevalence_table2.csv
#
# Notes:
#   - Alcohol definitions:
#       1. Police drinking, complete case
#       2. Police drinking, unknown as no
#       3. Vehicle DR_DRINK
#       4. Observed BAC >= 0.08
#       5. Observed BAC >= 0.10
#
#   - Sample variants:
#       1. all_state_years
#       2. bac_complete_95_state_years
#
#   The second sample variant restricts to state-year pairs where BAC reporting
#   completeness among police-reported drinking drivers is at least 95% (same as what LP did). 
################################################################################

rm(list = ls())

library(tidyverse)
library(scales)

# ------------------------------------------------------------------------------
# 0. User options
# ------------------------------------------------------------------------------

run_table2 <- TRUE
run_table3_by_year <- FALSE
run_table4_by_hour <- FALSE

# Profile-likelihood CIs are more reliable than Hessian SEs, but slower.
compute_profile_ci_table2 <- TRUE
compute_profile_ci_table3 <- FALSE
compute_profile_ci_table4 <- FALSE

driver_path <- "data/processed/driver_analysis_lp_2015_2024.csv"

out_dir <- "data/processed/script04_mle"
fig_dir <- "figures/script04_mle"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 0.5 Helper functions
# ------------------------------------------------------------------------------

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

safe_divide <- function(num, den) {
  ifelse(is.na(den) | den == 0, NA_real_, num / den)
}

format_est_se <- function(est, se, se_status, digits = 3) {
  case_when(
    is.na(est) ~ NA_character_,
    se_status == "ok" & !is.na(se) ~ paste0(round(est, digits), " (", round(se, digits), ")"),
    TRUE ~ paste0(round(est, digits), " (SE unreliable: ", se_status, ")")
  )
}

format_est_ci <- function(est, lo, hi, digits = 3) {
  case_when(
    is.na(est) ~ NA_character_,
    !is.na(lo) & !is.na(hi) ~ paste0(
      round(est, digits),
      " [",
      round(lo, digits),
      ", ",
      round(hi, digits),
      "]"
    ),
    TRUE ~ paste0(round(est, digits), " [CI unavailable]")
  )
}

# ------------------------------------------------------------------------------
# 1. Load driver-level LP sample
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
# 2. Recreate / refresh key variables
# ------------------------------------------------------------------------------

driver_lp <- driver_lp %>%
  mutate(
    # Police drinking definition: complete case
    police_drinking_complete = case_when(
      DRINKING == 1 ~ 1,
      DRINKING == 0 ~ 0,
      TRUE ~ NA_real_
    ),
    
    # Police drinking definition: unknown / not reported as no
    police_drinking_unknown_as_no = case_when(
      DRINKING == 1 ~ 1,
      TRUE ~ 0
    ),
    
    # Vehicle-level driver drinking variable
    vehicle_dr_drink_complete = case_when(
      DR_DRINK == 1 ~ 1,
      DR_DRINK == 0 ~ 0,
      TRUE ~ NA_real_
    ),
    
    # BAC definitions
    # FARS ALC_RES appears to be BAC * 1000.
    # Example: ALC_RES = 80 means BAC = 0.080.
    # Special values 995-999 are treated as invalid / missing-style values.
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
    
    covid_period = case_when(
      YEAR %in% 2015:2019 ~ "2015-2019 pre-COVID",
      YEAR %in% 2020:2021 ~ "2020-2021 COVID",
      YEAR %in% 2022:2024 ~ "2022-2024 post-COVID",
      TRUE ~ NA_character_
    ),
    
    weekend_lp = as.logical(WEEKEND_LP),
    
    # Cell labels used for LP-style specifications
    cell_all = "all_data",
    cell_hour = paste(HOUR, sep = "_"),
    cell_hour_year = paste(HOUR, YEAR, sep = "_"),
    cell_hour_year_weekend = paste(HOUR, YEAR, weekend_lp, sep = "_"),
    cell_hour_state_year = paste(HOUR, STATE, YEAR, sep = "_"),
    cell_hour_state_year_weekend = paste(HOUR, STATE, YEAR, weekend_lp, sep = "_")
  )

# If BAC completeness flags are missing for any reason, recreate a simple version.
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
  mutate(
    BAC_COMPLETE_95 = if_else(is.na(BAC_COMPLETE_95), FALSE, BAC_COMPLETE_95)
  )

# ------------------------------------------------------------------------------
# 3. Define alcohol definitions, cell specifications, and sample variants
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
    "hour_year",
    "hour_year_weekend",
    "hour_state_year",
    "hour_state_year_weekend"
  ),
  cell_var = c(
    "cell_all",
    "cell_hour",
    "cell_hour_year",
    "cell_hour_year_weekend",
    "cell_hour_state_year",
    "cell_hour_state_year_weekend"
  ),
  cell_spec_label = c(
    "all data",
    "hour",
    "hour x year",
    "hour x year x weekend",
    "hour x state x year",
    "hour x state x year x weekend"
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
# 4. Collapse driver rows to crash type for one alcohol definition
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
    mutate(
      drinking_status = .data[[definition_var]]
    ) %>%
    group_by(
      YEAR_FILE,
      YEAR,
      STATE,
      STATENAME,
      ST_CASE,
      VE_FORMS,
      HOUR,
      weekend_lp,
      covid_period,
      cell_all,
      cell_hour,
      cell_hour_year,
      cell_hour_year_weekend,
      cell_hour_state_year,
      cell_hour_state_year_weekend
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

# ------------------------------------------------------------------------------
# 5. Make cell-level counts for MLE
# ------------------------------------------------------------------------------

make_cell_counts <- function(crash_types_df, cell_var) {
  
  crash_types_df %>%
    mutate(cell_id = .data[[cell_var]]) %>%
    filter(!is.na(cell_id)) %>%
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
      total_all = total_two_vehicle + total_one_vehicle,
      
      has_two_vehicle_info = n2 > 0,
      has_one_vehicle_info = n1 > 0,
      has_drinking_info = (DD + DS + D) > 0,
      has_sober_info = (SS + DS + S) > 0
    ) %>%
    filter(total_complete_case > 0)
}

# ------------------------------------------------------------------------------
# 6. MLE helper functions
# ------------------------------------------------------------------------------

# Model:
#
# Two-car crash probabilities within a cell:
#   P_DD = theta * N^2 / Z2
#   P_DS = (theta + 1) * N / Z2
#   P_SS = 1 / Z2
#   where Z2 = theta * N^2 + (theta + 1) * N + 1
#
# One-car crash probabilities within a cell:
#   Q_D = lambda * N / Z1
#   Q_S = 1 / Z1
#   where Z1 = lambda * N + 1
#
# theta and lambda are common across cells.
# N is cell-specific.
#
# This implementation profiles out cell-specific eta = log(N) using Newton updates.
# Then it optimizes theta and lambda over log scale.

profile_eta_newton <- function(theta, lambda, counts, max_iter = 100, tol = 1e-8) {
  
  DD <- counts$DD
  DS <- counts$DS
  SS <- counts$SS
  D  <- counts$D
  S  <- counts$S
  n2 <- counts$n2
  n1 <- counts$n1
  
  # Initial eta based on observed drinking/sober driver balance in classifiable crashes.
  obs_D_drivers <- 2 * DD + DS + D
  obs_S_drivers <- 2 * SS + DS + S
  
  eta <- log((obs_D_drivers + 0.5) / (obs_S_drivers + 0.5))
  eta <- pmin(pmax(eta, -20), 20)
  
  for (iter in seq_len(max_iter)) {
    N <- exp(eta)
    
    Z2 <- theta * N^2 + (theta + 1) * N + 1
    A2 <- 2 * theta * N^2 + (theta + 1) * N
    A2_prime <- 4 * theta * N^2 + (theta + 1) * N
    
    Z1 <- lambda * N + 1
    A1 <- lambda * N
    
    grad <- (2 * DD + DS + D) -
      n2 * (A2 / Z2) -
      n1 * (A1 / Z1)
    
    hess <- -
      n2 * ((A2_prime * Z2 - A2^2) / (Z2^2)) -
      n1 * (A1 / (Z1^2))
    
    bad_hess <- is.na(hess) | !is.finite(hess) | abs(hess) < 1e-10
    
    step <- grad / hess
    step[bad_hess] <- 0
    
    # Limit Newton steps to avoid huge jumps in sparse cells.
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

profile_loglik <- function(log_params, counts, return_eta = FALSE) {
  
  theta <- exp(log_params[1])
  lambda <- exp(log_params[2])
  
  if (!is.finite(theta) || !is.finite(lambda) || theta <= 0 || lambda <= 0) {
    if (return_eta) {
      return(list(
        loglik = -Inf,
        eta = rep(NA_real_, nrow(counts)),
        N = rep(NA_real_, nrow(counts)),
        theta = theta,
        lambda = lambda
      ))
    }
    return(-Inf)
  }
  
  DD <- counts$DD
  DS <- counts$DS
  SS <- counts$SS
  D  <- counts$D
  S  <- counts$S
  n2 <- counts$n2
  n1 <- counts$n1
  
  eta <- profile_eta_newton(theta, lambda, counts)
  N <- exp(eta)
  
  Z2 <- theta * N^2 + (theta + 1) * N + 1
  Z1 <- lambda * N + 1
  
  if (any(!is.finite(Z2)) || any(!is.finite(Z1)) || any(Z2 <= 0) || any(Z1 <= 0)) {
    if (return_eta) {
      return(list(
        loglik = -Inf,
        eta = eta,
        N = N,
        theta = theta,
        lambda = lambda
      ))
    }
    return(-Inf)
  }
  
  # Log likelihood excluding multinomial constants.
  ll_two <- DD * (log(theta) + 2 * eta) +
    DS * (log(theta + 1) + eta) -
    n2 * log(Z2)
  
  ll_one <- D * (log(lambda) + eta) -
    n1 * log(Z1)
  
  ll <- sum(ll_two + ll_one, na.rm = TRUE)
  
  if (return_eta) {
    return(list(
      loglik = ll,
      eta = eta,
      N = N,
      theta = theta,
      lambda = lambda
    ))
  }
  
  ll
}

# ------------------------------------------------------------------------------
# 6.5 Profile-likelihood confidence intervals
# ------------------------------------------------------------------------------

profile_ci_one_param <- function(par_hat,
                                 counts,
                                 which_param = c("theta", "lambda"),
                                 level = 0.95,
                                 lower_risk = 1,
                                 upper_risk = 100) {
  
  which_param <- match.arg(which_param)
  cutoff <- qchisq(level, df = 1)
  
  ll_hat <- profile_loglik(par_hat, counts)
  
  if (!is.finite(ll_hat)) {
    return(tibble(
      estimate = exp(ifelse(which_param == "theta", par_hat[1], par_hat[2])),
      ci_low = NA_real_,
      ci_high = NA_real_,
      ci_status = "ll_hat_not_finite"
    ))
  }
  
  lower_log_bound <- log(lower_risk + 1e-6)
  upper_log_bound <- log(upper_risk)
  mle_log <- if (which_param == "theta") par_hat[1] else par_hat[2]
  
  target_fun <- function(log_value) {
    
    if (which_param == "theta") {
      
      opt_other <- tryCatch(
        optimize(
          f = function(log_lambda) {
            par_tmp <- c(log_value, log_lambda)
            val <- -profile_loglik(par_tmp, counts)
            if (!is.finite(val)) return(1e100)
            val
          },
          interval = c(lower_log_bound, upper_log_bound)
        ),
        error = function(e) NULL
      )
      
      if (is.null(opt_other)) return(NA_real_)
      
      par_fixed <- c(log_value, opt_other$minimum)
      
    } else {
      
      opt_other <- tryCatch(
        optimize(
          f = function(log_theta) {
            par_tmp <- c(log_theta, log_value)
            val <- -profile_loglik(par_tmp, counts)
            if (!is.finite(val)) return(1e100)
            val
          },
          interval = c(lower_log_bound, upper_log_bound)
        ),
        error = function(e) NULL
      )
      
      if (is.null(opt_other)) return(NA_real_)
      
      par_fixed <- c(opt_other$minimum, log_value)
    }
    
    ll_fixed <- profile_loglik(par_fixed, counts)
    
    if (!is.finite(ll_fixed)) return(NA_real_)
    
    2 * (ll_hat - ll_fixed) - cutoff
  }
  
  get_root_safe <- function(lower, upper) {
    
    f_lower <- tryCatch(target_fun(lower), error = function(e) NA_real_)
    f_upper <- tryCatch(target_fun(upper), error = function(e) NA_real_)
    
    if (is.na(f_lower) || is.na(f_upper)) return(NA_real_)
    
    # Need a sign change for uniroot.
    if (f_lower * f_upper > 0) return(NA_real_)
    
    tryCatch(
      uniroot(
        f = target_fun,
        lower = lower,
        upper = upper,
        tol = 1e-5
      )$root,
      error = function(e) NA_real_
    )
  }
  
  lower_ci_log <- get_root_safe(lower_log_bound, mle_log)
  upper_ci_log <- get_root_safe(mle_log, upper_log_bound)
  
  ci_status <- case_when(
    !is.na(lower_ci_log) & !is.na(upper_ci_log) ~ "ok",
    is.na(lower_ci_log) & is.na(upper_ci_log) ~ "both_bounds_failed",
    is.na(lower_ci_log) ~ "lower_bound_failed",
    is.na(upper_ci_log) ~ "upper_bound_failed",
    TRUE ~ "unknown"
  )
  
  tibble(
    estimate = exp(mle_log),
    ci_low = exp(lower_ci_log),
    ci_high = exp(upper_ci_log),
    ci_status = ci_status
  )
}

# ------------------------------------------------------------------------------
# 6.6 Main robust MLE function
# ------------------------------------------------------------------------------

fit_lp_mle <- function(counts,
                       start_theta = NULL,
                       start_lambda = NULL,
                       lower_risk = 1,
                       upper_risk = 100,
                       compute_profile_ci = FALSE,
                       min_relative_se = 5e-4,
                       max_hessian_condition_number = 1e10,
                       verbose = FALSE) {
  
  if (nrow(counts) == 0) {
    return(list(
      estimates = tibble(
        theta = NA_real_,
        theta_se = NA_real_,
        theta_ci_low = NA_real_,
        theta_ci_high = NA_real_,
        theta_ci_status = NA_character_,
        lambda = NA_real_,
        lambda_se = NA_real_,
        lambda_ci_low = NA_real_,
        lambda_ci_high = NA_real_,
        lambda_ci_status = NA_character_,
        logLik = NA_real_,
        convergence = NA_integer_,
        message = "No complete-case crashes available",
        optimizer = NA_character_,
        se_method = NA_character_,
        se_status = "not_computed_no_data",
        n_cells = 0,
        df_used = NA_real_
      ),
      cell_prevalence = tibble()
    ))
  }
  
  # Simple starting values from pooled counts.
  DD_total <- sum(counts$DD)
  DS_total <- sum(counts$DS)
  SS_total <- sum(counts$SS)
  D_total  <- sum(counts$D)
  S_total  <- sum(counts$S)
  
  if (is.null(start_theta)) {
    R <- if_else(
      DD_total > 0 & SS_total > 0,
      (DS_total^2) / (DD_total * SS_total),
      NA_real_
    )
    
    start_theta <- case_when(
      is.na(R) ~ 3,
      R < 4 ~ 1.01,
      TRUE ~ ((R - 2) + sqrt(R^2 - 4 * R)) / 2
    )
  }
  
  if (is.null(start_lambda)) {
    # Rough starting value; does not need to be exact.
    one_vehicle_odds <- if_else(
      D_total > 0 & S_total > 0,
      D_total / S_total,
      NA_real_
    )
    
    start_lambda <- case_when(
      !is.na(one_vehicle_odds) & is.finite(one_vehicle_odds) ~ max(start_theta, 1.01),
      TRUE ~ max(start_theta, 1.01)
    )
  }
  
  start_theta <- min(max(start_theta, lower_risk + 1e-4), upper_risk)
  start_lambda <- min(max(start_lambda, lower_risk + 1e-4), upper_risk)
  
  start <- log(c(start_theta, start_lambda))
  lower <- log(c(lower_risk + 1e-6, lower_risk + 1e-6))
  upper <- log(c(upper_risk, upper_risk))
  
  objective <- function(par) {
    # Penalize out-of-bound values. This also helps the fallback optimizer.
    if (any(par < lower) || any(par > upper)) {
      return(1e100 + sum((pmin(0, par - lower))^2) + sum((pmax(0, par - upper))^2))
    }
    
    val <- -profile_loglik(par, counts)
    if (!is.finite(val)) return(1e100)
    val
  }
  
  # Main optimizer: nlminb tends to be stable for bounded likelihood problems.
  fit_nlminb <- tryCatch(
    nlminb(
      start = start,
      objective = objective,
      lower = lower,
      upper = upper,
      control = list(
        eval.max = 2000,
        iter.max = 1000,
        rel.tol = 1e-10,
        x.tol = 1e-8
      )
    ),
    error = function(e) e
  )
  
  if (!inherits(fit_nlminb, "error") && fit_nlminb$convergence == 0) {
    
    fit <- list(
      par = fit_nlminb$par,
      value = fit_nlminb$objective,
      convergence = fit_nlminb$convergence,
      message = fit_nlminb$message %||% NA_character_,
      optimizer = "nlminb"
    )
    
  } else {
    
    # Fallback optimizer: Nelder-Mead with boundary penalty in objective().
    fit_optim <- tryCatch(
      optim(
        par = start,
        fn = objective,
        method = "Nelder-Mead",
        control = list(
          maxit = 3000,
          reltol = 1e-10
        )
      ),
      error = function(e) e
    )
    
    if (inherits(fit_optim, "error")) {
      return(list(
        estimates = tibble(
          theta = NA_real_,
          theta_se = NA_real_,
          theta_ci_low = NA_real_,
          theta_ci_high = NA_real_,
          theta_ci_status = NA_character_,
          lambda = NA_real_,
          lambda_se = NA_real_,
          lambda_ci_low = NA_real_,
          lambda_ci_high = NA_real_,
          lambda_ci_status = NA_character_,
          logLik = NA_real_,
          convergence = NA_integer_,
          message = paste("optimizer error:", fit_optim$message),
          optimizer = NA_character_,
          se_method = NA_character_,
          se_status = "optimizer_error",
          n_cells = nrow(counts),
          df_used = 2 + nrow(counts)
        ),
        cell_prevalence = tibble()
      ))
    }
    
    fit <- list(
      par = pmin(pmax(fit_optim$par, lower), upper),
      value = fit_optim$value,
      convergence = fit_optim$convergence,
      message = fit_optim$message %||% NA_character_,
      optimizer = "optim_Nelder_Mead_fallback"
    )
  }
  
  prof <- profile_loglik(fit$par, counts, return_eta = TRUE)
  
  theta_hat <- prof$theta
  lambda_hat <- prof$lambda
  logLik_hat <- prof$loglik
  
  # --------------------------------------------------------------------------
  # Hessian / delta-method SEs
  # --------------------------------------------------------------------------
  
  theta_se <- NA_real_
  lambda_se <- NA_real_
  se_method <- NA_character_
  se_status <- "not_computed"
  hessian_condition_number <- NA_real_
  hessian_min_eigen <- NA_real_
  hessian_max_eigen <- NA_real_
  
  near_lower_boundary <- any(abs(fit$par - lower) < 1e-5)
  near_upper_boundary <- any(abs(fit$par - upper) < 1e-5)
  near_boundary <- near_lower_boundary | near_upper_boundary
  
  if (fit$convergence == 0 && !near_boundary) {
    
    hess <- tryCatch(
      optimHess(fit$par, fn = objective),
      error = function(e) NULL
    )
    
    if (!is.null(hess) && all(is.finite(hess))) {
      
      # Symmetrize Hessian to reduce numerical asymmetry.
      hess <- (hess + t(hess)) / 2
      
      eig <- tryCatch(
        eigen(hess, symmetric = TRUE, only.values = TRUE)$values,
        error = function(e) NA_real_
      )
      
      if (all(is.finite(eig))) {
        
        hessian_min_eigen <- min(eig)
        hessian_max_eigen <- max(eig)
        hessian_condition_number <- abs(hessian_max_eigen / hessian_min_eigen)
        
        if (all(eig > 1e-8) &&
            is.finite(hessian_condition_number) &&
            hessian_condition_number <= max_hessian_condition_number) {
          
          vcov_log <- tryCatch(
            solve(hess),
            error = function(e) NULL
          )
          
          if (!is.null(vcov_log) &&
              all(is.finite(vcov_log)) &&
              all(diag(vcov_log) > 0)) {
            
            se_log <- sqrt(diag(vcov_log))
            
            theta_se_candidate <- theta_hat * se_log[1]
            lambda_se_candidate <- lambda_hat * se_log[2]
            
            theta_rel_se <- theta_se_candidate / theta_hat
            lambda_rel_se <- lambda_se_candidate / lambda_hat
            
            # Guard against fake precision. If the relative SE is implausibly tiny,
            # do not report it as if it were reliable.
            if (theta_rel_se < min_relative_se || lambda_rel_se < min_relative_se) {
              se_status <- "suspiciously_tiny_relative_se"
            } else {
              theta_se <- theta_se_candidate
              lambda_se <- lambda_se_candidate
              se_method <- "profiled_hessian_delta"
              se_status <- "ok"
            }
            
          } else {
            se_status <- "bad_vcov"
          }
          
        } else if (any(eig <= 1e-8)) {
          se_status <- "hessian_not_positive_definite"
        } else {
          se_status <- "hessian_ill_conditioned"
        }
        
      } else {
        se_status <- "hessian_eigen_failed"
      }
      
    } else {
      se_status <- "hessian_failed"
    }
    
  } else if (fit$convergence != 0) {
    se_status <- "optimizer_not_converged"
  } else if (near_boundary) {
    se_status <- "estimate_near_boundary"
  }
  
  # --------------------------------------------------------------------------
  # Profile-likelihood CIs
  # --------------------------------------------------------------------------
  
  theta_ci_low <- NA_real_
  theta_ci_high <- NA_real_
  theta_ci_status <- "not_computed"
  lambda_ci_low <- NA_real_
  lambda_ci_high <- NA_real_
  lambda_ci_status <- "not_computed"
  
  if (compute_profile_ci && fit$convergence == 0) {
    
    theta_ci <- profile_ci_one_param(
      par_hat = fit$par,
      counts = counts,
      which_param = "theta",
      lower_risk = lower_risk,
      upper_risk = upper_risk
    )
    
    lambda_ci <- profile_ci_one_param(
      par_hat = fit$par,
      counts = counts,
      which_param = "lambda",
      lower_risk = lower_risk,
      upper_risk = upper_risk
    )
    
    theta_ci_low <- theta_ci$ci_low
    theta_ci_high <- theta_ci$ci_high
    theta_ci_status <- theta_ci$ci_status
    
    lambda_ci_low <- lambda_ci$ci_low
    lambda_ci_high <- lambda_ci$ci_high
    lambda_ci_status <- lambda_ci$ci_status
  }
  
  # --------------------------------------------------------------------------
  # Cell prevalence output and estimate table
  # --------------------------------------------------------------------------
  
  cell_prevalence <- counts %>%
    mutate(
      eta_hat = prof$eta,
      N_hat_drinking_to_sober = prof$N,
      frac_drinking_hat = N_hat_drinking_to_sober /
        (1 + N_hat_drinking_to_sober),
      theta = theta_hat,
      lambda = lambda_hat
    )
  
  estimates <- tibble(
    theta = theta_hat,
    theta_se = theta_se,
    theta_ci_low = theta_ci_low,
    theta_ci_high = theta_ci_high,
    theta_ci_status = theta_ci_status,
    
    lambda = lambda_hat,
    lambda_se = lambda_se,
    lambda_ci_low = lambda_ci_low,
    lambda_ci_high = lambda_ci_high,
    lambda_ci_status = lambda_ci_status,
    
    logLik = logLik_hat,
    convergence = fit$convergence,
    message = fit$message %||% NA_character_,
    optimizer = fit$optimizer,
    
    se_method = se_method,
    se_status = se_status,
    hessian_min_eigen = hessian_min_eigen,
    hessian_max_eigen = hessian_max_eigen,
    hessian_condition_number = hessian_condition_number,
    
    n_cells = nrow(counts),
    df_used = 2 + nrow(counts),
    
    n_two_vehicle_complete = sum(counts$n2),
    n_one_vehicle_complete = sum(counts$n1),
    n_two_vehicle_all = sum(counts$total_two_vehicle),
    n_one_vehicle_all = sum(counts$total_one_vehicle),
    
    pct_two_vehicle_missing_or_ambiguous =
      safe_divide(sum(counts$missing_two_vehicle), sum(counts$total_two_vehicle)),
    pct_one_vehicle_missing_or_ambiguous =
      safe_divide(sum(counts$missing_one_vehicle), sum(counts$total_one_vehicle)),
    
    implied_frac_drinking_weighted =
      weighted.mean(
        cell_prevalence$frac_drinking_hat,
        w = cell_prevalence$total_complete_case,
        na.rm = TRUE
      ),
    
    # Sparse-cell diagnostics
    n_sparse_cells_le_1 = sum(counts$total_complete_case <= 1),
    n_sparse_cells_le_2 = sum(counts$total_complete_case <= 2),
    n_sparse_cells_le_5 = sum(counts$total_complete_case <= 5),
    pct_sparse_cells_le_1 = mean(counts$total_complete_case <= 1),
    pct_sparse_cells_le_2 = mean(counts$total_complete_case <= 2),
    pct_sparse_cells_le_5 = mean(counts$total_complete_case <= 5),
    
    n_cells_no_two_vehicle_info = sum(counts$n2 == 0),
    n_cells_no_one_vehicle_info = sum(counts$n1 == 0),
    n_cells_no_drinking_info = sum((counts$DD + counts$DS + counts$D) == 0),
    n_cells_no_sober_info = sum((counts$SS + counts$DS + counts$S) == 0)
  )
  
  list(
    estimates = estimates,
    cell_prevalence = cell_prevalence
  )
}

# ------------------------------------------------------------------------------
# 7. Run one MLE specification safely
# ------------------------------------------------------------------------------

run_one_mle_spec <- function(driver_data,
                             definition_var,
                             definition_label,
                             sample_variant,
                             cell_spec,
                             cell_var,
                             cell_spec_label,
                             subgroup_label = "all",
                             subgroup_var = NA_character_,
                             subgroup_value = NA_character_,
                             compute_profile_ci = FALSE) {
  
  cat(
    "\nRunning:",
    definition_label,
    "| sample:", sample_variant,
    "| cell:", cell_spec_label,
    "| subgroup:", subgroup_label,
    "\n"
  )
  
  crash_types <- make_crash_types_for_definition(
    driver_data = driver_data,
    definition_var = definition_var,
    definition_label = definition_label,
    sample_variant = sample_variant
  )
  
  counts <- make_cell_counts(crash_types, cell_var = cell_var)
  
  fit <- fit_lp_mle(
    counts = counts,
    compute_profile_ci = compute_profile_ci
  )
  
  estimates <- fit$estimates %>%
    mutate(
      alcohol_definition = definition_var,
      alcohol_definition_label = definition_label,
      sample_variant = sample_variant,
      cell_spec = cell_spec,
      cell_spec_label = cell_spec_label,
      subgroup_label = subgroup_label,
      subgroup_var = subgroup_var,
      subgroup_value = subgroup_value,
      .before = 1
    )
  
  cell_prevalence <- fit$cell_prevalence %>%
    mutate(
      alcohol_definition = definition_var,
      alcohol_definition_label = definition_label,
      sample_variant = sample_variant,
      cell_spec = cell_spec,
      cell_spec_label = cell_spec_label,
      subgroup_label = subgroup_label,
      subgroup_var = subgroup_var,
      subgroup_value = subgroup_value,
      .before = 1
    )
  
  list(
    estimates = estimates,
    cell_prevalence = cell_prevalence
  )
}

# ------------------------------------------------------------------------------
# 8. Table 2-style estimates: increasingly fine equal-mixing cells
# ------------------------------------------------------------------------------

table2_results <- list()
table2_cell_prevalence <- list()

if (run_table2) {
  
  counter <- 1
  
  for (sv in sample_variants$sample_variant) {
    for (d in seq_len(nrow(alcohol_definitions))) {
      for (c in seq_len(nrow(cell_specs))) {
        
        res <- run_one_mle_spec(
          driver_data = driver_lp,
          definition_var = alcohol_definitions$alcohol_definition[d],
          definition_label = alcohol_definitions$alcohol_definition_label[d],
          sample_variant = sv,
          cell_spec = cell_specs$cell_spec[c],
          cell_var = cell_specs$cell_var[c],
          cell_spec_label = cell_specs$cell_spec_label[c],
          compute_profile_ci = compute_profile_ci_table2
        )
        
        table2_results[[counter]] <- res$estimates
        table2_cell_prevalence[[counter]] <- res$cell_prevalence
        counter <- counter + 1
      }
    }
  }
}

if (length(table2_results) > 0) {
  
  table2_long <- bind_rows(table2_results) %>%
    left_join(sample_variants, by = "sample_variant") %>%
    arrange(
      sample_variant,
      alcohol_definition,
      match(cell_spec, cell_specs$cell_spec)
    )
  
  table2_cells <- bind_rows(table2_cell_prevalence)
  
} else {
  
  table2_long <- tibble()
  table2_cells <- tibble()
}

write_csv(
  table2_long,
  file.path(out_dir, "lp_mle_table2_long.csv")
)

write_csv(
  table2_cells,
  file.path(out_dir, "lp_mle_cell_prevalence_table2.csv")
)

if (nrow(table2_long) > 0) {
  
  table2_wide_theta <- table2_long %>%
    mutate(
      theta_with_se = format_est_se(theta, theta_se, se_status),
      theta_with_ci = format_est_ci(theta, theta_ci_low, theta_ci_high)
    ) %>%
    select(
      sample_variant_label,
      alcohol_definition_label,
      cell_spec_label,
      theta_with_se,
      theta_with_ci
    ) %>%
    pivot_wider(
      names_from = cell_spec_label,
      values_from = c(theta_with_se, theta_with_ci)
    )
  
  table2_wide_lambda <- table2_long %>%
    mutate(
      lambda_with_se = format_est_se(lambda, lambda_se, se_status),
      lambda_with_ci = format_est_ci(lambda, lambda_ci_low, lambda_ci_high)
    ) %>%
    select(
      sample_variant_label,
      alcohol_definition_label,
      cell_spec_label,
      lambda_with_se,
      lambda_with_ci
    ) %>%
    pivot_wider(
      names_from = cell_spec_label,
      values_from = c(lambda_with_se, lambda_with_ci)
    )
  
} else {
  
  table2_wide_theta <- tibble()
  table2_wide_lambda <- tibble()
}

write_csv(
  table2_wide_theta,
  file.path(out_dir, "lp_mle_table2_wide_theta.csv")
)

write_csv(
  table2_wide_lambda,
  file.path(out_dir, "lp_mle_table2_wide_lambda.csv")
)

# ------------------------------------------------------------------------------
# 9. Table 3-style estimates: risk varies by year
# ------------------------------------------------------------------------------

table3_results <- list()

if (run_table3_by_year) {
  
  counter <- 1
  years_to_run <- sort(unique(driver_lp$YEAR))
  
  # For year-specific estimates, use the preferred LP-style cell structure
  # within each year: state x hour x weekend.
  for (sv in sample_variants$sample_variant) {
    for (d in seq_len(nrow(alcohol_definitions))) {
      for (yr in years_to_run) {
        
        driver_sub <- driver_lp %>%
          filter(YEAR == yr) %>%
          mutate(
            cell_table3_year = paste(STATE, HOUR, weekend_lp, sep = "_")
          )
        
        res <- run_one_mle_spec(
          driver_data = driver_sub,
          definition_var = alcohol_definitions$alcohol_definition[d],
          definition_label = alcohol_definitions$alcohol_definition_label[d],
          sample_variant = sv,
          cell_spec = "state_hour_weekend_within_year",
          cell_var = "cell_table3_year",
          cell_spec_label = "state x hour x weekend within year",
          subgroup_label = paste0("YEAR = ", yr),
          subgroup_var = "YEAR",
          subgroup_value = as.character(yr),
          compute_profile_ci = compute_profile_ci_table3
        )
        
        table3_results[[counter]] <- res$estimates
        counter <- counter + 1
      }
    }
  }
}

if (length(table3_results) > 0) {
  
  table3_by_year <- bind_rows(table3_results) %>%
    left_join(sample_variants, by = "sample_variant") %>%
    mutate(YEAR = as.integer(subgroup_value)) %>%
    arrange(sample_variant, alcohol_definition, YEAR)
  
} else {
  
  table3_by_year <- tibble()
}

write_csv(
  table3_by_year,
  file.path(out_dir, "lp_mle_table3_by_year.csv")
)

if (nrow(table3_by_year) > 0) {
  
  cat("\nTable 3-style by-year estimates:\n")
  print(
    table3_by_year %>%
      select(
        sample_variant_label,
        alcohol_definition_label,
        YEAR,
        theta,
        theta_se,
        theta_ci_low,
        theta_ci_high,
        lambda,
        lambda_se,
        lambda_ci_low,
        lambda_ci_high,
        se_status,
        convergence,
        optimizer,
        implied_frac_drinking_weighted,
        n_two_vehicle_complete,
        n_one_vehicle_complete
      ),
    n = Inf,
    width = Inf
  )
}

# ------------------------------------------------------------------------------
# 10. Table 4-style estimates: risk varies by hour
# ------------------------------------------------------------------------------

table4_results <- list()

if (run_table4_by_hour) {
  
  counter <- 1
  hours_to_run <- sort(unique(driver_lp$HOUR))
  
  # The LP nighttime sample should contain 20, 21, 22, 23, 0, 1, 2, 3, 4.
  # For hour-specific estimates, use the preferred LP-style cell structure
  # within each hour: state x year x weekend.
  for (sv in sample_variants$sample_variant) {
    for (d in seq_len(nrow(alcohol_definitions))) {
      for (hr in hours_to_run) {
        
        driver_sub <- driver_lp %>%
          filter(HOUR == hr) %>%
          mutate(
            cell_table4_hour = paste(STATE, YEAR, weekend_lp, sep = "_")
          )
        
        res <- run_one_mle_spec(
          driver_data = driver_sub,
          definition_var = alcohol_definitions$alcohol_definition[d],
          definition_label = alcohol_definitions$alcohol_definition_label[d],
          sample_variant = sv,
          cell_spec = "state_year_weekend_within_hour",
          cell_var = "cell_table4_hour",
          cell_spec_label = "state x year x weekend within hour",
          subgroup_label = paste0("HOUR = ", hr),
          subgroup_var = "HOUR",
          subgroup_value = as.character(hr),
          compute_profile_ci = compute_profile_ci_table4
        )
        
        table4_results[[counter]] <- res$estimates
        counter <- counter + 1
      }
    }
  }
}

if (length(table4_results) > 0) {
  
  table4_by_hour <- bind_rows(table4_results) %>%
    left_join(sample_variants, by = "sample_variant") %>%
    mutate(
      HOUR = as.integer(subgroup_value),
      hour_label = case_when(
        HOUR == 20 ~ "8:00 p.m.-9:00 p.m.",
        HOUR == 21 ~ "9:00 p.m.-10:00 p.m.",
        HOUR == 22 ~ "10:00 p.m.-11:00 p.m.",
        HOUR == 23 ~ "11:00 p.m.-midnight",
        HOUR == 0  ~ "midnight-1:00 a.m.",
        HOUR == 1  ~ "1:00 a.m.-2:00 a.m.",
        HOUR == 2  ~ "2:00 a.m.-3:00 a.m.",
        HOUR == 3  ~ "3:00 a.m.-4:00 a.m.",
        HOUR == 4  ~ "4:00 a.m.-5:00 a.m.",
        TRUE ~ as.character(HOUR)
      ),
      hour_order = case_when(
        HOUR == 20 ~ 1,
        HOUR == 21 ~ 2,
        HOUR == 22 ~ 3,
        HOUR == 23 ~ 4,
        HOUR == 0  ~ 5,
        HOUR == 1  ~ 6,
        HOUR == 2  ~ 7,
        HOUR == 3  ~ 8,
        HOUR == 4  ~ 9,
        TRUE ~ 99
      )
    ) %>%
    arrange(sample_variant, alcohol_definition, hour_order)
  
} else {
  
  table4_by_hour <- tibble()
}

write_csv(
  table4_by_hour,
  file.path(out_dir, "lp_mle_table4_by_hour.csv")
)

if (nrow(table4_by_hour) > 0) {
  
  cat("\nTable 4-style by-hour estimates:\n")
  print(
    table4_by_hour %>%
      select(
        sample_variant_label,
        alcohol_definition_label,
        hour_label,
        theta,
        theta_se,
        theta_ci_low,
        theta_ci_high,
        lambda,
        lambda_se,
        lambda_ci_low,
        lambda_ci_high,
        se_status,
        convergence,
        optimizer,
        implied_frac_drinking_weighted,
        n_two_vehicle_complete,
        n_one_vehicle_complete
      ),
    n = Inf,
    width = Inf
  )
}

# ------------------------------------------------------------------------------
# 11. Basic plots
# ------------------------------------------------------------------------------

if (nrow(table2_long) > 0) {
  
  plot_table2_theta <- table2_long %>%
    mutate(
      cell_spec_label = factor(
        cell_spec_label,
        levels = cell_specs$cell_spec_label
      )
    ) %>%
    ggplot(
      aes(
        x = cell_spec_label,
        y = theta,
        group = alcohol_definition_label
      )
    ) +
    geom_line() +
    geom_point(size = 2) +
    facet_wrap(~ sample_variant_label + alcohol_definition_label, scales = "free_y") +
    coord_flip() +
    labs(
      title = "LP MLE Table 2-Style Estimates: Two-Car Relative Risk",
      subtitle = "Common theta with cell-specific drinking/sober prevalence ratio N",
      x = "Equal-mixing cell specification",
      y = expression(theta)
    ) +
    theme_minimal(base_size = 12)
  
  ggsave(
    filename = file.path(fig_dir, "lp_mle_table2_theta_by_cell_spec.png"),
    plot = plot_table2_theta,
    width = 14,
    height = 10,
    dpi = 300
  )
  
  plot_table2_lambda <- table2_long %>%
    mutate(
      cell_spec_label = factor(
        cell_spec_label,
        levels = cell_specs$cell_spec_label
      )
    ) %>%
    ggplot(
      aes(
        x = cell_spec_label,
        y = lambda,
        group = alcohol_definition_label
      )
    ) +
    geom_line() +
    geom_point(size = 2) +
    facet_wrap(~ sample_variant_label + alcohol_definition_label, scales = "free_y") +
    coord_flip() +
    labs(
      title = "LP MLE Table 2-Style Estimates: One-Car Relative Risk",
      subtitle = "Common lambda with cell-specific drinking/sober prevalence ratio N",
      x = "Equal-mixing cell specification",
      y = expression(lambda)
    ) +
    theme_minimal(base_size = 12)
  
  ggsave(
    filename = file.path(fig_dir, "lp_mle_table2_lambda_by_cell_spec.png"),
    plot = plot_table2_lambda,
    width = 14,
    height = 10,
    dpi = 300
  )
}

if (nrow(table3_by_year) > 0) {
  
  plot_table3_year_theta <- table3_by_year %>%
    ggplot(
      aes(
        x = YEAR,
        y = theta,
        group = alcohol_definition_label
      )
    ) +
    geom_line() +
    geom_point(size = 2) +
    facet_wrap(~ sample_variant_label + alcohol_definition_label, scales = "free_y") +
    scale_x_continuous(breaks = sort(unique(table3_by_year$YEAR))) +
    labs(
      title = "LP MLE Table 3-Style Estimates by Year",
      subtitle = "Preferred within-year cells: state x hour x weekend",
      x = "Year",
      y = expression(theta)
    ) +
    theme_minimal(base_size = 12)
  
  ggsave(
    filename = file.path(fig_dir, "lp_mle_table3_theta_by_year.png"),
    plot = plot_table3_year_theta,
    width = 14,
    height = 10,
    dpi = 300
  )
}

if (nrow(table4_by_hour) > 0) {
  
  plot_table4_hour_theta <- table4_by_hour %>%
    mutate(
      hour_label = factor(hour_label, levels = unique(hour_label[order(hour_order)]))
    ) %>%
    ggplot(
      aes(
        x = hour_label,
        y = theta,
        group = alcohol_definition_label
      )
    ) +
    geom_line() +
    geom_point(size = 2) +
    facet_wrap(~ sample_variant_label + alcohol_definition_label, scales = "free_y") +
    coord_flip() +
    labs(
      title = "LP MLE Table 4-Style Estimates by Hour",
      subtitle = "Preferred within-hour cells: state x year x weekend",
      x = "Hour",
      y = expression(theta)
    ) +
    theme_minimal(base_size = 12)
  
  ggsave(
    filename = file.path(fig_dir, "lp_mle_table4_theta_by_hour.png"),
    plot = plot_table4_hour_theta,
    width = 14,
    height = 10,
    dpi = 300
  )
}

# ------------------------------------------------------------------------------
# 12. Final notes printed to console
# ------------------------------------------------------------------------------

cat("\nDone! Main outputs saved in:\n")
cat("  ", out_dir, "\n\n")

cat("Key files:\n")
cat("  lp_mle_table2_long.csv\n")
cat("  lp_mle_table2_wide_theta.csv\n")
cat("  lp_mle_table2_wide_lambda.csv\n")
cat("  lp_mle_table3_by_year.csv\n")
cat("  lp_mle_table4_by_hour.csv\n")
cat("  lp_mle_cell_prevalence_table2.csv\n\n")

cat("Interpretation note:\n")
cat("  theta is the estimated two-car relative fatal crash risk for drinking vs sober drivers.\n")
cat("  lambda is the estimated one-car relative fatal crash risk for drinking vs sober drivers.\n")
cat("  N_hat_drinking_to_sober is the estimated cell-specific ratio of drinking to sober drivers.\n")
cat("  frac_drinking_hat = N / (1 + N).\n\n")

cat("Robust SE note:\n")
cat("  Hessian/delta-method SEs are only reported when convergence and Hessian diagnostics pass.\n")
cat("  If se_status is not 'ok', the SE should be treated as unreliable.\n")
cat("  Profile-likelihood CIs, when computed, are shown in the wide tables as [lower, upper].\n")