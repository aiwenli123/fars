# Estimate two-group fatal crash risks using BAC cutoffs of 0.05 and 0.08.
# Refit the two-vehicle and one-vehicle models in 1,000 crash-level bootstrap
# samples, then combine sampling and imputation uncertainty using Rubin's rules.
#
# Input:  data/processed/driver_miper_long_2015_2024.csv
# Outputs:
#   results/two_group_bootstrap/two_group_bootstrap_combined_results.csv
#   results/two_group_bootstrap/two_group_bootstrap_replicate_estimates.csv
# The replicate file is retained so pooled results can be rebuilt without
# repeating the bootstrap. 
# If existing combined results are available, they will be reused by default.

# rm(list=ls())
library(tidyverse)

mi_driver_path <- "data/processed/driver_miper_long_2015_2024.csv"
out_dir <- "results/two_group_bootstrap"
combined_path <- file.path(out_dir, "two_group_bootstrap_combined_results.csv")
replicate_path <- file.path(out_dir, "two_group_bootstrap_replicate_estimates.csv")
force_recompute <- FALSE  # TRUE reruns all bootstrap fits and overwrites outputs.

bootstrap_reps <- 1000
bootstrap_seed <- 123
lower_risk <- 1
upper_risk <- 200

threshold_specs <- tibble(
  threshold = c(0.05, 0.08),
  threshold_label = c("BAC >= 0.05", "BAC >= 0.08"),
  mi_definition_var = c("MI_BAC_GE_05", "MI_BAC_GE_08")
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Reusing the final output does not read the driver data or rerun the bootstrap.
if (!force_recompute && file.exists(combined_path)) {
  message("Reusing existing results: ", combined_path)
  message("No bootstrap was run. Set force_recompute <- TRUE to rerun it.")
  print(read_csv(combined_path, show_col_types = FALSE), width = Inf)
  } else {

  # Classify one- and two-vehicle crashes at each BAC cutoff.
  make_crash_types <- function(driver_data, drinking_var) {
    driver_data %>%
      mutate(DRINKING_STATUS_15 = .data[[drinking_var]], CRASH_ID_15 = paste(YEAR_FILE, ST_CASE, sep = "_")) %>%
      group_by(YEAR_FILE, ST_CASE, VE_FORMS, CELL_RUR_URB_WEEKEND_15, CRASH_ID_15) %>%
      summarize(
        n_drivers = n(),
        n_drinking = sum(DRINKING_STATUS_15 == 1, na.rm = TRUE),
        n_sober = sum(DRINKING_STATUS_15 == 0, na.rm = TRUE),
        n_missing = sum(is.na(DRINKING_STATUS_15)),
        .groups = "drop"
      ) %>%
      mutate(
        cell_id = CELL_RUR_URB_WEEKEND_15,
        crash_type = case_when(
          VE_FORMS == 2 & n_drivers == 2 & n_missing == 0 & n_drinking == 2 & n_sober == 0 ~ "DD",
          VE_FORMS == 2 & n_drivers == 2 & n_missing == 0 & n_drinking == 1 & n_sober == 1 ~ "DS",
          VE_FORMS == 2 & n_drivers == 2 & n_missing == 0 & n_drinking == 0 & n_sober == 2 ~ "SS",
          VE_FORMS == 1 & n_drivers == 1 & n_missing == 0 & n_drinking == 1 & n_sober == 0 ~ "D",
          VE_FORMS == 1 & n_drivers == 1 & n_missing == 0 & n_drinking == 0 & n_sober == 1 ~ "S",
          TRUE ~ "missing_or_ambiguous"
        )
      ) %>%
      select(CRASH_ID_15, cell_id, VE_FORMS, crash_type)
  }

  # Tabulate the crash types in each rural/urban-by-weekend cell.
  # Missing/ambiguous crash types contribute no likelihood counts.
  make_cell_counts <- function(crash_types_df, bootstrap_weight = NULL) {
    if (is.null(bootstrap_weight)) bootstrap_weight <- rep(1L, nrow(crash_types_df))

    crash_types_df %>%
      mutate(bootstrap_weight = bootstrap_weight) %>%
      group_by(cell_id) %>%
      summarize(
        DD = sum(bootstrap_weight * (crash_type == "DD")),
        DS = sum(bootstrap_weight * (crash_type == "DS")),
        SS = sum(bootstrap_weight * (crash_type == "SS")),
        D = sum(bootstrap_weight * (crash_type == "D")),
        S = sum(bootstrap_weight * (crash_type == "S")),
        .groups = "drop"
      ) %>%
      mutate(n2 = DD + DS + SS, n1 = D + S) %>%
      filter(n1 + n2 > 0)
  }

  # Fit theta by profiling the cell exposure ratios; then fit lambda.
  # All original likelihoods, bounds, tolerances, and optimizers are unchanged.
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
    ll <- DD * (log(theta) + 2 * eta) + DS * (log(theta + 1) + eta) - n2 * log(Z)
    total_ll <- sum(ll, na.rm = TRUE)

    if (return_eta) return(list(loglik = total_ll, theta = theta, eta = eta, N = N))
    total_ll
  }

  fit_theta_simple <- function(counts, min_n2 = 1, lower_risk = 1, upper_risk = 100) {
    counts_two <- counts %>% filter(n2 >= min_n2)
    if (nrow(counts_two) == 0) return(list(theta = NA_real_, cell_prevalence = tibble()))

    lower_log <- log(lower_risk + 1e-8)
    upper_log <- log(upper_risk)
    objective <- function(log_theta) {
      val <- -two_car_loglik(log_theta, counts_two)
      if (!is.finite(val)) return(1e100)
      val
    }
    opt <- optimize(f = objective, interval = c(lower_log, upper_log), tol = 1e-10)
    prof <- two_car_loglik(opt$minimum, counts_two, return_eta = TRUE)
    cell_prevalence <- counts_two %>%
      mutate(theta = prof$theta, eta_hat = prof$eta, N_hat_drinking_to_sober = prof$N)
    list(theta = prof$theta, cell_prevalence = cell_prevalence)
  }

  lambda_plugin_loglik <- function(log_lambda, plugin_counts) {
    lambda <- exp(log_lambda)
    if (!is.finite(lambda) || lambda <= 0) return(-Inf)
    D <- plugin_counts$D
    n1 <- plugin_counts$n1
    Nhat <- plugin_counts$N_hat_drinking_to_sober
    ll <- D * (log(lambda) + log(Nhat)) - n1 * log(1 + lambda * Nhat)
    sum(ll, na.rm = TRUE)
  }

  fit_lambda_simple <- function(counts, cell_prevalence, lower_risk = 1, upper_risk = 100) {
    if (nrow(cell_prevalence) == 0) return(NA_real_)

    plugin_counts <- counts %>%
      select(cell_id, D, S, n1) %>%
      inner_join(cell_prevalence %>% select(cell_id, N_hat_drinking_to_sober), by = "cell_id") %>%
      filter(n1 > 0, is.finite(N_hat_drinking_to_sober), N_hat_drinking_to_sober > 0)
    if (nrow(plugin_counts) == 0) return(NA_real_)

    lower_log <- log(lower_risk + 1e-8)
    upper_log <- log(upper_risk)
    objective <- function(log_lambda) {
      val <- -lambda_plugin_loglik(log_lambda, plugin_counts)
      if (!is.finite(val)) return(1e100)
      val
    }
    opt <- optimize(f = objective, interval = c(lower_log, upper_log), tol = 1e-10)
    exp(opt$minimum)
  }

  fit_benchmark <- function(crash_types_df, bootstrap_weight = NULL) {
    counts <- make_cell_counts(crash_types_df, bootstrap_weight)
    theta_fit <- fit_theta_simple(counts, min_n2 = 1, lower_risk = lower_risk, upper_risk = upper_risk)
    lambda_fit <- fit_lambda_simple(counts, theta_fit$cell_prevalence, lower_risk = lower_risk, upper_risk = upper_risk)
    tibble(theta = theta_fit$theta, lambda = lambda_fit)
  }

  # Fit the unresampled data in each of the ten imputations and at each cutoff.
  driver_mi_long <- read_csv(mi_driver_path, show_col_types = FALSE) %>%
    rename_with(toupper) %>%
    mutate(
      WEEKEND_15 = case_when(
        "WEEKEND_LP" %in% names(.) ~ as.logical(WEEKEND_LP),
        DAY_WEEK %in% c(1, 7) ~ TRUE,
        TRUE ~ FALSE
      ),
      RUR_URB_GROUP_15 = case_when(RUR_URB == 1 ~ "rural", RUR_URB == 2 ~ "urban", TRUE ~ NA_character_),
      CELL_RUR_URB_WEEKEND_15 = paste(RUR_URB_GROUP_15, WEEKEND_15, sep = "_")
    ) %>%
    filter(YEAR %in% 2015:2024, VE_FORMS %in% c(1, 2), !is.na(RUR_URB_GROUP_15))

  imputations <- sort(unique(driver_mi_long$IMPUTATION))
  classified_samples <- list()
  original_results <- list()
  result_index <- 1

  for (imp in imputations) {
    imp_data <- driver_mi_long %>% filter(IMPUTATION == imp)
    for (t in seq_len(nrow(threshold_specs))) {
      key <- paste(imp, threshold_specs$threshold[t], sep = "_")
      classified_samples[[key]] <- make_crash_types(imp_data, threshold_specs$mi_definition_var[t])
      original_results[[result_index]] <- fit_benchmark(classified_samples[[key]]) %>%
        mutate(imputation = imp, threshold = threshold_specs$threshold[t], threshold_label = threshold_specs$threshold_label[t], .before = 1)
      result_index <- result_index + 1
    }
  }

  original_results <- bind_rows(original_results)

  # If original replicate estimates exist, reuse them rather than refitting
  # all 20,000 bootstrap models. Otherwise, draw exactly the original samples.
  if (!force_recompute && file.exists(replicate_path)) {
    message("Reusing existing bootstrap replicates: ", replicate_path)
    bootstrap_results <- read_csv(replicate_path, show_col_types = FALSE)
  } else {

    base_crashes <- classified_samples[[paste(imputations[1], threshold_specs$threshold[1], sep = "_")]] %>%
      select(CRASH_ID_15, cell_id, VE_FORMS)
    strata_indices <- split(seq_len(nrow(base_crashes)), interaction(base_crashes$cell_id, base_crashes$VE_FORMS, drop = TRUE))
    set.seed(bootstrap_seed)

    bootstrap_results <- vector("list", bootstrap_reps * length(imputations) * nrow(threshold_specs))
    result_index <- 1

    for (b in seq_len(bootstrap_reps)) {
      bootstrap_weight <- integer(nrow(base_crashes))
      for (idx in strata_indices) {
        bootstrap_weight[idx] <- tabulate(sample.int(length(idx), length(idx), replace = TRUE), nbins = length(idx))
      }

      for (imp in imputations) {
        for (t in seq_len(nrow(threshold_specs))) {
          key <- paste(imp, threshold_specs$threshold[t], sep = "_")
          sample_crashes <- classified_samples[[key]]
          sample_weight <- bootstrap_weight[match(sample_crashes$CRASH_ID_15, base_crashes$CRASH_ID_15)]

          bootstrap_results[[result_index]] <- tryCatch(
            fit_benchmark(sample_crashes, sample_weight) %>%
              mutate(replicate = b, imputation = imp, threshold = threshold_specs$threshold[t], threshold_label = threshold_specs$threshold_label[t], .before = 1),
            error = function(e) tibble(
              replicate = b, imputation = imp, threshold = threshold_specs$threshold[t],
              threshold_label = threshold_specs$threshold_label[t], theta = NA_real_, lambda = NA_real_
            )
          )
          result_index <- result_index + 1
        }
      }

      if (b %% 50 == 0) message("Completed ", b, " of ", bootstrap_reps, " bootstrap replicates.")
    }

    bootstrap_results <- bind_rows(bootstrap_results)

    # Keep the replicate estimates as a computational input for future reuse.
    write_csv(bootstrap_results, replicate_path)
  }

  # Compute within-imputation bootstrap variances and pool the estimates.
  within_variances <- bootstrap_results %>%
    pivot_longer(c(theta, lambda), names_to = "parameter", values_to = "bootstrap_estimate") %>%
    group_by(imputation, threshold, threshold_label, parameter) %>%
    summarize(
      bootstrap_reps_requested = bootstrap_reps,
      bootstrap_reps_used = sum(is.finite(bootstrap_estimate)),
      within_var = var(bootstrap_estimate, na.rm = TRUE),
      bootstrap_se = sqrt(within_var),
      .groups = "drop"
    )

  original_long <- original_results %>%
    pivot_longer(c(theta, lambda), names_to = "parameter", values_to = "imputation_estimate") %>%
    left_join(within_variances, by = c("imputation", "threshold", "threshold_label", "parameter"))

  combined_results <- original_long %>%
    group_by(threshold, threshold_label, parameter) %>%
    summarize(
      m = n(),
      estimate = mean(imputation_estimate),
      within_var = mean(within_var),
      between_var = var(imputation_estimate),
      total_var = within_var + (1 + 1 / m) * between_var,
      se_total = sqrt(total_var),
      ci_low = estimate - 1.96 * se_total,
      ci_high = estimate + 1.96 * se_total,
      min_bootstrap_reps_used = min(bootstrap_reps_used),
      .groups = "drop"
    )

  # These paired percentile intervals are useful as a distributional check. The
  # Rubin-bootstrap intervals above remain primary because they also add MI uncertainty.
  paired_percentile_results <- bootstrap_results %>%
    pivot_longer(c(theta, lambda), names_to = "parameter", values_to = "bootstrap_estimate") %>%
    group_by(replicate, threshold, threshold_label, parameter) %>%
    summarize(estimate = mean(bootstrap_estimate, na.rm = TRUE), .groups = "drop") %>%
    group_by(threshold, threshold_label, parameter) %>%
    summarize(
      percentile_ci_low = quantile(estimate, 0.025, na.rm = TRUE),
      percentile_ci_high = quantile(estimate, 0.975, na.rm = TRUE),
      .groups = "drop"
    )

  combined_results <- combined_results %>%
    left_join(paired_percentile_results, by = c("threshold", "threshold_label", "parameter")) %>%
    arrange(threshold, factor(parameter, levels = c("theta", "lambda")))

  write_csv(combined_results, combined_path)
  print(combined_results, width = Inf)
  message("Saved combined results: ", combined_path)
}
