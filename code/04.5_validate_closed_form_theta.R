################################################################################
# 04.5_validate_closed_form_theta.R
#
# Purpose:
#   Validate the two-group Levitt-Porter theta code in the simplest setting:
#
#     * one all-data cell;
#     * two-car crashes only (DD, DS, SS);
#     * the same alcohol definitions and sample variants used in Script 04.
#
#   For one cell, the numerical maximum-likelihood estimate should agree with
#   the Levitt-Porter closed-form estimate based on
#
#       R = DS^2 / (DD * SS)
#
#       theta_hat = [(R - 2) + sqrt(R^2 - 4R)] / 2,  when R >= 4,
#       theta_hat = 1,                                 when R < 4.
#
#   The script performs three comparisons:
#
#     1. Closed-form theta vs. an independently rerun profile-likelihood
#        optimizer using the same likelihood as Script 04.
#     2. Closed-form theta vs. Script 04's saved all-data estimate.
#     3. Independent optimizer vs. Script 04's saved all-data estimate.
#
# Required prior scripts:
#   00_read_data.R through 04_lp_mle_period.R
#
# Required inputs:
#   data/processed/driver_analysis_lp_2015_2024.csv
#   data/processed/script04_simple_mle_period/lp_period_table2_long.csv
#
# Outputs:
#   data/processed/script04.5_closed_form_validation/
#     theta_closed_form_validation.csv
#     theta_closed_form_validation_summary.csv
#
# Note:
#   Cases with DD = 0 or SS = 0 do not have a finite R formula and are marked
#   "not_comparable" rather than failed.
################################################################################

rm(list = ls())

library(tidyverse)
library(scales)

# ------------------------------------------------------------------------------
# 0. User options and paths
# ------------------------------------------------------------------------------

driver_path <- "data/processed/driver_analysis_lp_2015_2024.csv"
script04_results_path <- paste0(
  "data/processed/script04_simple_mle_period/",
  "lp_period_table2_long.csv"
)

out_dir <- "data/processed/script04.5_closed_form_validation"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Match Script 04's parameter bounds.
lower_risk <- 1
upper_risk <- 200

# Agreement rule:
#   absolute difference <= abs_tolerance + rel_tolerance * |reference value|
abs_tolerance <- 1e-6
rel_tolerance <- 1e-6

# Set TRUE if you want the script to stop after saving outputs when a
# comparable validation check fails. FALSE prints a warning and continues.
stop_on_failure <- FALSE

# ------------------------------------------------------------------------------
# 1. Helper functions
# ------------------------------------------------------------------------------

safe_divide <- function(num, den) {
  ifelse(is.na(den) | den == 0, NA_real_, num / den)
}

is_close <- function(x, y,
                     abs_tol = abs_tolerance,
                     rel_tol = rel_tolerance) {
  is.finite(x) & is.finite(y) &
    abs(x - y) <= abs_tol + rel_tol * abs(y)
}

# Return both the unconstrained LP closed-form estimate and the estimate expected
# after imposing Script 04's [lower_risk, upper_risk] optimization bounds.
closed_form_theta <- function(DD, DS, SS,
                              lower_bound = lower_risk,
                              upper_bound = upper_risk) {
  R_value <- safe_divide(DS^2, DD * SS)

  if (is.na(R_value) || DD <= 0 || SS <= 0) {
    return(tibble(
      R = R_value,
      theta_closed_form_unbounded = NA_real_,
      theta_expected_with_bounds = NA_real_,
      closed_form_status = "undefined_zero_DD_or_SS",
      comparable = FALSE
    ))
  }

  # Small numerical protection around R = 4.
  if (R_value < 4) {
    theta_unbounded <- 1
    status <- "lower_boundary_R_below_4"
  } else {
    discriminant <- max(R_value^2 - 4 * R_value, 0)
    theta_unbounded <- ((R_value - 2) + sqrt(discriminant)) / 2
    status <- "interior_closed_form"
  }

  theta_bounded <- min(max(theta_unbounded, lower_bound), upper_bound)

  if (theta_unbounded > upper_bound) {
    status <- "upper_optimizer_bound_binding"
  }

  tibble(
    R = R_value,
    theta_closed_form_unbounded = theta_unbounded,
    theta_expected_with_bounds = theta_bounded,
    closed_form_status = status,
    comparable = TRUE
  )
}

# This is the same cell-specific profiling step used in Script 04. For a fixed
# theta, it maximizes the two-car likelihood over eta = log(N), where
# N = N_D / N_S is the drinking-to-sober exposure ratio.
profile_eta_two_car <- function(theta, counts_two,
                                max_iter = 100,
                                tol = 1e-10) {
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

    gradient <- (2 * DD + DS) - n2 * (A / Z)
    hessian <- -n2 * ((Aprime * Z - A^2) / Z^2)

    bad_hessian <- is.na(hessian) |
      !is.finite(hessian) |
      abs(hessian) < 1e-10

    step <- gradient / hessian
    step[bad_hessian] <- 0
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

two_car_loglik <- function(log_theta, counts_two, return_details = FALSE) {
  theta <- exp(log_theta)

  if (!is.finite(theta) || theta <= 0) {
    if (return_details) {
      return(list(
        loglik = -Inf,
        theta = theta,
        eta = NA_real_,
        N = NA_real_
      ))
    }
    return(-Inf)
  }

  DD <- counts_two$DD
  DS <- counts_two$DS
  n2 <- counts_two$n2

  eta <- profile_eta_two_car(theta, counts_two)
  N <- exp(eta)
  Z <- theta * N^2 + (theta + 1) * N + 1

  # Multinomial log likelihood up to the count-only constant.
  loglik_cell <- DD * (log(theta) + 2 * eta) +
    DS * (log(theta + 1) + eta) -
    n2 * log(Z)

  total_loglik <- sum(loglik_cell, na.rm = TRUE)

  if (return_details) {
    return(list(
      loglik = total_loglik,
      theta = theta,
      eta = eta,
      N = N
    ))
  }

  total_loglik
}

fit_theta_optimizer <- function(DD, DS, SS,
                                lower_bound = lower_risk,
                                upper_bound = upper_risk) {
  counts_two <- tibble(
    cell_id = "all_data",
    DD = DD,
    DS = DS,
    SS = SS,
    n2 = DD + DS + SS
  )

  if (counts_two$n2 <= 0) {
    return(tibble(
      theta_optimizer = NA_real_,
      N_hat_optimizer = NA_real_,
      optimizer_logLik = NA_real_,
      optimizer_status = "no_two_car_data"
    ))
  }

  # Match Script 04 exactly at the lower boundary.
  lower_log <- log(lower_bound + 1e-8)
  upper_log <- log(upper_bound)

  objective <- function(log_theta) {
    value <- -two_car_loglik(log_theta, counts_two)
    if (!is.finite(value)) return(1e100)
    value
  }

  optimization <- optimize(
    f = objective,
    interval = c(lower_log, upper_log),
    tol = 1e-10
  )

  fit <- two_car_loglik(
    optimization$minimum,
    counts_two,
    return_details = TRUE
  )

  near_lower <- abs(optimization$minimum - lower_log) < 1e-4
  near_upper <- abs(optimization$minimum - upper_log) < 1e-4

  status <- case_when(
    near_lower ~ "near_lower_bound",
    near_upper ~ "near_upper_bound",
    TRUE ~ "interior"
  )

  tibble(
    theta_optimizer = fit$theta,
    N_hat_optimizer = as.numeric(fit$N),
    optimizer_logLik = fit$loglik,
    optimizer_status = status
  )
}

# ------------------------------------------------------------------------------
# 2. Load and prepare the Script 04 analysis data
# ------------------------------------------------------------------------------

if (!file.exists(driver_path)) {
  stop(
    "Missing input file: ", driver_path,
    "\nRun Scripts 00-01 before Script 04.5."
  )
}

if (!file.exists(script04_results_path)) {
  stop(
    "Missing Script 04 output: ", script04_results_path,
    "\nRun 04_lp_mle_period.R before Script 04.5."
  )
}

driver_lp <- read_csv(driver_path, show_col_types = FALSE) %>%
  rename_with(toupper) %>%
  mutate(
    period = case_when(
      YEAR %in% 2015:2019 ~ "pre_covid",
      YEAR %in% 2020:2021 ~ "covid",
      YEAR %in% 2022:2024 ~ "post_covid",
      TRUE ~ NA_character_
    ),

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

    bac_ge_05 = case_when(
      bac_valid & ALC_RES >= 50 ~ 1,
      bac_valid & ALC_RES < 50 ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(!is.na(period))

# Recreate the 95%-complete indicator exactly as Script 04 does if it is absent.
if (!"BAC_COMPLETE_95" %in% names(driver_lp)) {
  message("BAC_COMPLETE_95 not found; recreating it as in Script 04.")

  bac_completeness_state_year <- driver_lp %>%
    group_by(YEAR, STATE) %>%
    summarize(
      n_police_drinking_yes = sum(
        police_drinking_complete == 1,
        na.rm = TRUE
      ),
      n_police_drinking_yes_bac_valid = sum(
        police_drinking_complete == 1 & bac_valid,
        na.rm = TRUE
      ),
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
    BAC_COMPLETE_95 = if_else(
      is.na(BAC_COMPLETE_95),
      FALSE,
      as.logical(BAC_COMPLETE_95)
    )
  )

alcohol_definitions <- tribble(
  ~alcohol_definition,                 ~alcohol_definition_label,
  "police_drinking_complete",          "Police drinking, complete case",
  "police_drinking_unknown_as_no",     "Police drinking, unknown as no",
  "vehicle_dr_drink_complete",         "Vehicle DR_DRINK, complete case",
  "bac_ge_08",                         "Observed BAC >= 0.08",
  "bac_ge_10",                         "Observed BAC >= 0.10",
  "bac_ge_05",                         "Observed BAC >= 0.05"
)

sample_variants <- tribble(
  ~sample_variant,                    ~sample_variant_label,
  "all_state_years",                  "All state-years",
  "bac_complete_95_state_years",      "BAC-complete 95% state-years"
)

# ------------------------------------------------------------------------------
# 3. Build all-data DD/DS/SS counts exactly as in Script 04
# ------------------------------------------------------------------------------

make_all_data_counts <- function(driver_data,
                                 definition_var,
                                 definition_label,
                                 sample_variant,
                                 sample_variant_label) {
  analysis_data <- driver_data

  if (sample_variant == "bac_complete_95_state_years") {
    analysis_data <- analysis_data %>%
      filter(BAC_COMPLETE_95)
  }

  crash_types <- analysis_data %>%
    mutate(drinking_status = .data[[definition_var]]) %>%
    group_by(
      YEAR_FILE,
      YEAR,
      STATE,
      ST_CASE,
      VE_FORMS
    ) %>%
    summarize(
      n_drivers = n(),
      n_drinking = sum(drinking_status == 1, na.rm = TRUE),
      n_sober = sum(drinking_status == 0, na.rm = TRUE),
      n_missing = sum(is.na(drinking_status)),
      .groups = "drop"
    ) %>%
    mutate(
      crash_type = case_when(
        VE_FORMS == 2 & n_missing == 0 &
          n_drinking == 2 & n_sober == 0 ~ "DD",
        VE_FORMS == 2 & n_missing == 0 &
          n_drinking == 1 & n_sober == 1 ~ "DS",
        VE_FORMS == 2 & n_missing == 0 &
          n_drinking == 0 & n_sober == 2 ~ "SS",
        TRUE ~ "missing_or_ambiguous"
      )
    )

  crash_types %>%
    summarize(
      DD = sum(crash_type == "DD"),
      DS = sum(crash_type == "DS"),
      SS = sum(crash_type == "SS"),
      missing_or_ambiguous_two_car = sum(
        VE_FORMS == 2 & crash_type == "missing_or_ambiguous"
      ),
      total_two_car = sum(VE_FORMS == 2),
      .groups = "drop"
    ) %>%
    mutate(
      n2_complete = DD + DS + SS,
      alcohol_definition = definition_var,
      alcohol_definition_label = definition_label,
      sample_variant = sample_variant,
      sample_variant_label = sample_variant_label,
      .before = 1
    )
}

validation_grid <- expand_grid(
  sample_variant = sample_variants$sample_variant,
  alcohol_definition = alcohol_definitions$alcohol_definition
) %>%
  left_join(sample_variants, by = "sample_variant") %>%
  left_join(alcohol_definitions, by = "alcohol_definition") %>%
  select(
    sample_variant,
    sample_variant_label,
    alcohol_definition,
    alcohol_definition_label
  )

validation_counts <- validation_grid %>%
  pmap_dfr(
    function(sample_variant,
             sample_variant_label,
             alcohol_definition,
             alcohol_definition_label) {
      make_all_data_counts(
        driver_data = driver_lp,
        definition_var = alcohol_definition,
        definition_label = alcohol_definition_label,
        sample_variant = sample_variant,
        sample_variant_label = sample_variant_label
      )
    }
  )

# ------------------------------------------------------------------------------
# 4. Calculate the closed-form and independent optimizer estimates
# ------------------------------------------------------------------------------

validation_independent <- validation_counts %>%
  rowwise() %>%
  mutate(
    closed_form = list(closed_form_theta(DD, DS, SS)),
    optimizer = list(fit_theta_optimizer(DD, DS, SS))
  ) %>%
  ungroup() %>%
  unnest_wider(closed_form) %>%
  unnest_wider(optimizer)

# ------------------------------------------------------------------------------
# 5. Read Script 04's saved all-data estimates
# ------------------------------------------------------------------------------

script04_all_data <- read_csv(
  script04_results_path,
  show_col_types = FALSE
) %>%
  filter(cell_spec == "all_data") %>%
  select(
    sample_variant,
    alcohol_definition,
    min_n2,
    theta_script04 = theta,
    theta_script04_status = theta_status,
    theta_script04_logLik = theta_logLik,
    n_two_vehicle_complete_used
  ) %>%
  mutate(min_n2 = paste0("min_n2_", min_n2)) %>%
  pivot_wider(
    names_from = min_n2,
    values_from = c(
      theta_script04,
      theta_script04_status,
      theta_script04_logLik,
      n_two_vehicle_complete_used
    ),
    names_glue = "{.value}_{min_n2}"
  )

validation <- validation_independent %>%
  left_join(
    script04_all_data,
    by = c("sample_variant", "alcohol_definition")
  ) %>%
  mutate(
    # Closed-form vs independently rerun optimizer.
    abs_diff_formula_vs_optimizer = abs(
      theta_expected_with_bounds - theta_optimizer
    ),
    rel_diff_formula_vs_optimizer = safe_divide(
      abs_diff_formula_vs_optimizer,
      abs(theta_expected_with_bounds)
    ),
    pass_formula_vs_optimizer = case_when(
      !comparable ~ NA,
      TRUE ~ is_close(theta_optimizer, theta_expected_with_bounds)
    ),

    # Closed-form vs the min_n2 = 1 all-data result saved by Script 04.
    abs_diff_formula_vs_script04 = abs(
      theta_expected_with_bounds - theta_script04_min_n2_1
    ),
    rel_diff_formula_vs_script04 = safe_divide(
      abs_diff_formula_vs_script04,
      abs(theta_expected_with_bounds)
    ),
    pass_formula_vs_script04 = case_when(
      !comparable ~ NA,
      is.na(theta_script04_min_n2_1) ~ FALSE,
      TRUE ~ is_close(
        theta_script04_min_n2_1,
        theta_expected_with_bounds
      )
    ),

    # Independent optimizer vs the value saved by Script 04.
    abs_diff_optimizer_vs_script04 = abs(
      theta_optimizer - theta_script04_min_n2_1
    ),
    rel_diff_optimizer_vs_script04 = safe_divide(
      abs_diff_optimizer_vs_script04,
      abs(theta_optimizer)
    ),
    pass_optimizer_vs_script04 = case_when(
      is.na(theta_script04_min_n2_1) | is.na(theta_optimizer) ~ FALSE,
      TRUE ~ is_close(theta_script04_min_n2_1, theta_optimizer)
    ),

    # For one all-data cell, min_n2 = 1 and min_n2 = 5 should use the same
    # data whenever the total complete two-car count is at least five.
    pass_script04_min_n2_invariance = case_when(
      n2_complete < 5 ~ NA,
      is.na(theta_script04_min_n2_1) |
        is.na(theta_script04_min_n2_5) ~ FALSE,
      TRUE ~ is_close(
        theta_script04_min_n2_1,
        theta_script04_min_n2_5
      )
    ),

    overall_validation = case_when(
      !comparable ~ "not_comparable",
      pass_formula_vs_optimizer &
        pass_formula_vs_script04 &
        pass_optimizer_vs_script04 &
        (is.na(pass_script04_min_n2_invariance) |
           pass_script04_min_n2_invariance) ~ "PASS",
      TRUE ~ "FAIL"
    )
  ) %>%
  arrange(sample_variant, alcohol_definition) %>%
  select(
    sample_variant,
    sample_variant_label,
    alcohol_definition,
    alcohol_definition_label,
    DD,
    DS,
    SS,
    n2_complete,
    missing_or_ambiguous_two_car,
    total_two_car,
    R,
    closed_form_status,
    comparable,
    theta_closed_form_unbounded,
    theta_expected_with_bounds,
    theta_optimizer,
    N_hat_optimizer,
    optimizer_status,
    optimizer_logLik,
    theta_script04_min_n2_1,
    theta_script04_status_min_n2_1,
    theta_script04_logLik_min_n2_1,
    theta_script04_min_n2_5,
    theta_script04_status_min_n2_5,
    abs_diff_formula_vs_optimizer,
    rel_diff_formula_vs_optimizer,
    pass_formula_vs_optimizer,
    abs_diff_formula_vs_script04,
    rel_diff_formula_vs_script04,
    pass_formula_vs_script04,
    abs_diff_optimizer_vs_script04,
    rel_diff_optimizer_vs_script04,
    pass_optimizer_vs_script04,
    pass_script04_min_n2_invariance,
    overall_validation
  )

# ------------------------------------------------------------------------------
# 6. Save and print validation results
# ------------------------------------------------------------------------------

validation_path <- file.path(
  out_dir,
  "theta_closed_form_validation.csv"
)

write_csv(validation, validation_path)

validation_summary <- validation %>%
  summarize(
    n_checks = n(),
    n_comparable = sum(comparable, na.rm = TRUE),
    n_not_comparable = sum(!comparable, na.rm = TRUE),
    n_pass = sum(overall_validation == "PASS", na.rm = TRUE),
    n_fail = sum(overall_validation == "FAIL", na.rm = TRUE),
    all_comparable_checks_pass = n_fail == 0,
    max_abs_diff_formula_vs_optimizer = if_else(
      any(comparable),
      max(abs_diff_formula_vs_optimizer[comparable], na.rm = TRUE),
      NA_real_
    ),
    max_abs_diff_formula_vs_script04 = if_else(
      any(comparable),
      max(abs_diff_formula_vs_script04[comparable], na.rm = TRUE),
      NA_real_
    ),
    abs_tolerance = abs_tolerance,
    rel_tolerance = rel_tolerance,
    lower_risk_bound = lower_risk,
    upper_risk_bound = upper_risk
  )

summary_path <- file.path(
  out_dir,
  "theta_closed_form_validation_summary.csv"
)

write_csv(validation_summary, summary_path)

cat("\n", strrep("=", 78), "\n", sep = "")
cat("SCRIPT 04.5: CLOSED-FORM THETA VALIDATION\n")
cat(strrep("=", 78), "\n\n", sep = "")

print(
  validation %>%
    select(
      sample_variant_label,
      alcohol_definition_label,
      DD,
      DS,
      SS,
      R,
      theta_expected_with_bounds,
      theta_optimizer,
      theta_script04_min_n2_1,
      abs_diff_formula_vs_optimizer,
      abs_diff_formula_vs_script04,
      overall_validation
    ),
  n = Inf,
  width = Inf
)

cat("\nSummary:\n")
print(validation_summary, width = Inf)

cat("\nSaved detailed results to:\n  ", validation_path, "\n", sep = "")
cat("Saved summary to:\n  ", summary_path, "\n", sep = "")

failed_checks <- validation %>%
  filter(overall_validation == "FAIL")

if (nrow(failed_checks) > 0) {
  failure_message <- paste0(
    "Script 04.5 found ", nrow(failed_checks),
    " comparable specification(s) that failed the tolerance check. ",
    "Inspect theta_closed_form_validation.csv for details."
  )

  if (stop_on_failure) {
    stop(failure_message)
  } else {
    warning(failure_message)
  }
} else {
  cat("\nPASS: All comparable closed-form, independent-optimizer, and ",
      "Script-04 estimates agree within tolerance.\n", sep = "")
}
