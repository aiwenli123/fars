################################################################################
# 13_assumption2_sensitivity.R
#
# Purpose:
#   Extend the four-BAC-group model from Script 11 by relaxing Assumption 2
#   (equal mixing within cells) in two ways.
#
#   Model 1: sober versus positive-BAC clustering
#     * Keep the four risk groups S, L, M, H.
#     * Treat L/M/H as one positive-BAC unit only for the mixing mechanism.
#     * Fit Delta = 0, 0.05, 0.10, and 0.15.
#
#   Model 2: exact-BAC-group clustering
#     * Increase S-S, L-L, M-M, and H-H interactions.
#     * Fit a small rho grid concentrated near zero.
#     * Produce diagnostics showing how large the implied diagonal increases are,
#       especially for the relatively rare L and M groups.
#
# Estimation:
#   * Follow Script 11's two-stage procedure.
#   * First estimate theta_L/theta_M/theta_H and cell-specific r_L/r_M/r_H
#     from two-car crashes.
#   * Then use the new r estimates in the one-car model to estimate
#     lambda_L/lambda_M/lambda_H.
#   * Re-estimate both stages in every imputed data set and pool positive
#     parameters on the log scale using the same Rubin-style calculation as
#     Script 11.
#
# Main input:
#   data/processed/script10_four_group_prep/
#     four_group_crash_counts_by_cell_imp.rds
#
# Notes:
#   * By default, this script estimates only the preferred
#     rural/urban x weekend specification. Add other Script 10 specifications to
#     cell_specs_to_estimate if desired.
#   * Delta = 0 and rho = 0 are both the original equal-mixing model. Their
#     agreement is used as an internal numerical check.
#   * Positive rho scenarios are fitted sequentially within each imputation.
#     The previous rho solution is added as a warm start, while all original
#     fixed starts are retained.
#   * The script saves objective values from every optimizer start, plots
#     imputation-specific rho paths, and creates scenario-level fit summaries.
#   * The script also keeps both the global-best solution and a
#     baseline-connected solution that stays closest to the previous sensitivity
#     value on the log-parameter scale. This helps reveal near-tied likelihood
#     modes instead of averaging across qualitatively different solutions.
#   * For a road share p_g, the relative diagonal increase is
#     rho * (1 - p_g) / p_g, so even small rho values can matter for rare groups.
################################################################################

rm(list = ls())

library(tidyverse)

# ------------------------------------------------------------------------------
# 0. User options
# ------------------------------------------------------------------------------

prep_dir <- "data/processed/script10_four_group_prep"
base_out_dir <- "data/processed/script13_assumption2_sensitivity"
base_fig_dir <- "figures/script13_assumption2_sensitivity"

counts_path <- file.path(prep_dir, "four_group_crash_counts_by_cell_imp.rds")

# Turn either sensitivity model on or off.
# For the final complete run, set both to TRUE.
run_delta_models <- TRUE
run_rho_models <- TRUE

if (!run_delta_models && !run_rho_models) {
  stop("At least one of run_delta_models or run_rho_models must be TRUE.")
}

# Partial runs are saved in a subfolder so they do not overwrite the most recent
# full Delta-and-rho results.
save_partial_run_in_subfolder <- TRUE

run_label <- case_when(
  run_delta_models && run_rho_models ~ "delta_and_rho",
  run_delta_models ~ "delta_only",
  TRUE ~ "rho_only"
)

if (save_partial_run_in_subfolder && run_label != "delta_and_rho") {
  out_dir <- file.path(base_out_dir, run_label)
  fig_dir <- file.path(base_fig_dir, run_label)
} else {
  out_dir <- base_out_dir
  fig_dir <- base_fig_dir
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Preferred specification from Scripts 11-12.
cell_specs_to_estimate <- c("rur_urb_weekend")

# Same default as Script 11.
min_n2 <- 1

# Values that will actually be fitted.
delta_values <- c(0, 0.05, 0.10, 0.15)
rho_values_to_fit <- c(0, 0.0025, 0.005, 0.0075, 0.01)

# Values used only for formula-based diagnostic tables.
rho_diagnostic_grid <- c(0.0025, 0.005, 0.0075, 0.01, 0.025, 0.05)

# Reference increases for questions such as:
# "What rho would increase the L-L interaction probability by 25 percent?"
target_diagonal_increases <- c(0.10, 0.25, 0.50, 1.00)

# Bounds for every positive theta, lambda, and r parameter.
log_param_lower <- log(1e-4)
log_param_upper <- log(500)

optim_control <- list(maxit = 5000, factr = 1e7)

# A "distinct near-tied mode" warning is created when the global-best and
# baseline-connected two-car solutions are far apart in parameter space but
# have almost the same negative log-likelihood. These cutoffs are diagnostics,
# not formal hypothesis-test thresholds.
mode_distance_warning_threshold <- 1.00
mode_objective_gap_warning_threshold <- 0.10
same_solution_max_abs_difference <- 1e-4

# Multiple fixed starts reduce dependence on one arbitrary initialization.
theta_start_grid <- list(
  c(theta_L = 1.2, theta_M = 2.0, theta_H = 5.0),
  c(theta_L = 1.0, theta_M = 1.0, theta_H = 1.0),
  c(theta_L = 1.5, theta_M = 3.0, theta_H = 10.0),
  c(theta_L = 2.0, theta_M = 5.0, theta_H = 20.0)
)

lambda_start_grid <- list(
  c(lambda_L = 1.2, lambda_M = 2.0, lambda_H = 5.0),
  c(lambda_L = 1.0, lambda_M = 1.0, lambda_H = 1.0),
  c(lambda_L = 1.5, lambda_M = 3.0, lambda_H = 10.0),
  c(lambda_L = 2.0, lambda_M = 5.0, lambda_H = 20.0)
)

# ------------------------------------------------------------------------------
# 1. Sensitivity scenarios
# ------------------------------------------------------------------------------

scenario_parts <- list()

if (run_delta_models) {
  scenario_parts[["delta"]] <- tibble(
    mixing_model = "sober_positive_delta",
    sensitivity_name = "Delta",
    sensitivity_value = delta_values
  )
}

if (run_rho_models) {
  scenario_parts[["rho"]] <- tibble(
    mixing_model = "exact_group_rho",
    sensitivity_name = "rho",
    sensitivity_value = rho_values_to_fit
  )
}

scenario_grid <- bind_rows(scenario_parts) %>%
  mutate(
    scenario_id = paste0(
      if_else(mixing_model == "sober_positive_delta", "delta_", "rho_"),
      formatC(sensitivity_value, format = "f", digits = 4)
    ),
    scenario_label = case_when(
      mixing_model == "sober_positive_delta" ~
        paste0("S vs positive BAC: Delta = ", sensitivity_value),
      mixing_model == "exact_group_rho" ~
        paste0("Exact BAC groups: rho = ", sensitivity_value),
      TRUE ~ scenario_id
    )
  ) %>%
  arrange(mixing_model, sensitivity_value)

# Equal-mixing results are used as the reference for percent changes and the
# formula-based rho diagnostics. Prefer Delta = 0 when Delta was run; otherwise
# use rho = 0.
equal_mixing_reference_model <- if (
  any(
    scenario_grid$mixing_model == "sober_positive_delta" &
      scenario_grid$sensitivity_value == 0
  )
) {
  "sober_positive_delta"
} else {
  "exact_group_rho"
}

# ------------------------------------------------------------------------------
# 2. Small helper functions
# ------------------------------------------------------------------------------

safe_log <- function(x) {
  log(pmax(x, 1e-300))
}

safe_divide <- function(num, den) {
  ifelse(is.na(den) | den == 0, NA_real_, num / den)
}

clamp <- function(x, lower, upper) {
  pmin(pmax(x, lower), upper)
}

bad_likelihood_penalty <- 1e100

two_car_types <- c("SS", "SL", "SM", "SH", "LL", "LM", "LH", "MM", "MH", "HH")
same_group_types <- c("SS", "LL", "MM", "HH")
groups4 <- c("S", "L", "M", "H")

# Approximate r starting values from observed two-car driver endpoints.
make_r_starts <- function(dat) {
  dat %>%
    transmute(
      cell_id,
      endpoint_S = 2 * SS + SL + SM + SH,
      endpoint_L = SL + 2 * LL + LM + LH,
      endpoint_M = SM + LM + 2 * MM + MH,
      endpoint_H = SH + LH + MH + 2 * HH,
      r_L = safe_divide(endpoint_L + 0.5, endpoint_S + 0.5),
      r_M = safe_divide(endpoint_M + 0.5, endpoint_S + 0.5),
      r_H = safe_divide(endpoint_H + 0.5, endpoint_S + 0.5)
    ) %>%
    select(cell_id, r_L, r_M, r_H)
}

unpack_two_car_par <- function(par, cell_ids) {
  theta <- c(
    S = 1,
    L = unname(exp(par["log_theta_L"])),
    M = unname(exp(par["log_theta_M"])),
    H = unname(exp(par["log_theta_H"]))
  )

  r_mat <- tibble(cell_id = cell_ids) %>%
    mutate(
      r_S = 1,
      r_L = exp(par[paste0("log_r_L__", cell_id)]),
      r_M = exp(par[paste0("log_r_M__", cell_id)]),
      r_H = exp(par[paste0("log_r_H__", cell_id)])
    )

  list(theta = theta, r_mat = r_mat)
}

# ------------------------------------------------------------------------------
# 3. Two-car probabilities under the two Assumption-2 extensions
# ------------------------------------------------------------------------------

make_two_car_weights <- function(theta, r_mat, mixing_model, sensitivity_value) {
  theta_L <- unname(theta["L"])
  theta_M <- unname(theta["M"])
  theta_H <- unname(theta["H"])

  if (mixing_model == "sober_positive_delta") {
    Delta <- sensitivity_value

    out <- r_mat %>%
      mutate(
        R_positive = r_L + r_M + r_H,
        T = 1 + R_positive,
        validity_margin = 1 - Delta * R_positive,

        # S-S and S-positive weights.
        w_SS = 1 + Delta * R_positive^2,
        w_SL = (1 + theta_L) * r_L * validity_margin,
        w_SM = (1 + theta_M) * r_M * validity_margin,
        w_SH = (1 + theta_H) * r_H * validity_margin,

        # Positive-positive weights.
        w_LL = (1 + Delta) * theta_L * r_L^2,
        w_LM = (1 + Delta) * (theta_L + theta_M) * r_L * r_M,
        w_LH = (1 + Delta) * (theta_L + theta_H) * r_L * r_H,
        w_MM = (1 + Delta) * theta_M * r_M^2,
        w_MH = (1 + Delta) * (theta_M + theta_H) * r_M * r_H,
        w_HH = (1 + Delta) * theta_H * r_H^2
      )

  } else if (mixing_model == "exact_group_rho") {
    rho <- sensitivity_value

    if (rho < 0 || rho > 1) {
      stop("rho must be between 0 and 1.")
    }

    out <- r_mat %>%
      mutate(
        R_positive = r_L + r_M + r_H,
        T = 1 + R_positive,
        validity_margin = 1 - rho,

        # Exact-group diagonal weights.
        w_SS = (1 - rho) + rho * T,
        w_LL = theta_L * ((1 - rho) * r_L^2 + rho * r_L * T),
        w_MM = theta_M * ((1 - rho) * r_M^2 + rho * r_M * T),
        w_HH = theta_H * ((1 - rho) * r_H^2 + rho * r_H * T),

        # Every mixed pair is reduced by the common factor 1-rho.
        w_SL = (1 - rho) * (1 + theta_L) * r_L,
        w_SM = (1 - rho) * (1 + theta_M) * r_M,
        w_SH = (1 - rho) * (1 + theta_H) * r_H,
        w_LM = (1 - rho) * (theta_L + theta_M) * r_L * r_M,
        w_LH = (1 - rho) * (theta_L + theta_H) * r_L * r_H,
        w_MH = (1 - rho) * (theta_M + theta_H) * r_M * r_H
      )

  } else {
    stop("Unknown mixing_model: ", mixing_model)
  }

  out %>%
    mutate(
      D_weight_sum = w_SS + w_SL + w_SM + w_SH + w_LL +
        w_LM + w_LH + w_MM + w_MH + w_HH,

      # Both sensitivity models should have this same denominator.
      D_expected = T * (1 + theta_L * r_L + theta_M * r_M + theta_H * r_H),
      denominator_abs_diff = abs(D_weight_sum - D_expected)
    )
}

make_two_car_probabilities <- function(theta, r_mat, mixing_model, sensitivity_value) {
  weights <- make_two_car_weights(
    theta = theta,
    r_mat = r_mat,
    mixing_model = mixing_model,
    sensitivity_value = sensitivity_value
  )

  weight_cols <- paste0("w_", two_car_types)

  # Delta requires 1 - Delta * R_positive > 0. The rho model is valid for
  # 0 <= rho <= 1, but rho = 1 would make all mixed probabilities zero.
  invalid <-
    any(!is.finite(as.matrix(weights[, weight_cols]))) ||
    any(!is.finite(weights$D_weight_sum)) ||
    any(weights$D_weight_sum <= 0) ||
    any(weights$validity_margin <= 0) ||
    any(as.matrix(weights[, weight_cols]) < 0)

  if (invalid) {
    return(NULL)
  }

  probabilities <- weights
  for (type in two_car_types) {
    probabilities[[paste0("p_", type)]] <-
      probabilities[[paste0("w_", type)]] / probabilities$D_weight_sum
  }

  probabilities
}

negloglik_two_car <- function(par, dat, cell_ids, mixing_model, sensitivity_value) {
  unpacked <- unpack_two_car_par(par, cell_ids)

  probs <- make_two_car_probabilities(
    theta = unpacked$theta,
    r_mat = unpacked$r_mat,
    mixing_model = mixing_model,
    sensitivity_value = sensitivity_value
  )

  if (is.null(probs)) {
    return(bad_likelihood_penalty)
  }

  ll_dat <- dat %>%
    left_join(
      probs %>% select(cell_id, starts_with("p_")),
      by = "cell_id"
    )

  prob_cols <- paste0("p_", two_car_types)
  if (any(!is.finite(as.matrix(ll_dat[, prob_cols])))) {
    return(bad_likelihood_penalty)
  }

  loglik_terms <- rep(0, nrow(ll_dat))
  for (type in two_car_types) {
    loglik_terms <- loglik_terms +
      ll_dat[[type]] * safe_log(ll_dat[[paste0("p_", type)]])
  }

  if (any(!is.finite(loglik_terms))) {
    return(bad_likelihood_penalty)
  }

  -sum(loglik_terms)
}

make_two_car_start <- function(dat, theta_start) {
  r_starts <- make_r_starts(dat)

  par <- c(
    log_theta_L = log(unname(theta_start["theta_L"])),
    log_theta_M = log(unname(theta_start["theta_M"])),
    log_theta_H = log(unname(theta_start["theta_H"]))
  )

  for (i in seq_len(nrow(r_starts))) {
    cell_i <- r_starts$cell_id[i]
    par[paste0("log_r_L__", cell_i)] <- log(clamp(r_starts$r_L[i], 1e-4, 500))
    par[paste0("log_r_M__", cell_i)] <- log(clamp(r_starts$r_M[i], 1e-4, 500))
    par[paste0("log_r_H__", cell_i)] <- log(clamp(r_starts$r_H[i], 1e-4, 500))
  }

  par
}


# Add one named warm start if its parameter names match the fixed starts.
add_one_warm_start <- function(starts, warm_start_par, start_label) {
  if (is.null(warm_start_par)) {
    return(starts)
  }

  template_names <- names(starts[[1]])

  if (!setequal(names(warm_start_par), template_names)) {
    warning(start_label, " parameter names did not match. This warm start was skipped.")
    return(starts)
  }

  warm_start_par <- warm_start_par[template_names]
  warm_list <- list(warm_start_par)
  names(warm_list) <- start_label

  c(warm_list, starts)
}

# Keep both previous branches as possible starts. The fixed starts are always
# retained, so branch tracking does not force the optimizer to follow one path.
add_warm_starts <- function(
    fixed_starts,
    warm_start_global = NULL,
    warm_start_connected = NULL
) {
  starts <- fixed_starts

  starts <- add_one_warm_start(
    starts = starts,
    warm_start_par = warm_start_connected,
    start_label = "warm_start_connected_previous"
  )

  starts <- add_one_warm_start(
    starts = starts,
    warm_start_par = warm_start_global,
    start_label = "warm_start_global_previous"
  )

  starts
}

parameter_distance <- function(par_1, par_2) {
  if (is.null(par_1) || is.null(par_2)) {
    return(NA_real_)
  }

  if (!setequal(names(par_1), names(par_2))) {
    return(NA_real_)
  }

  par_2 <- par_2[names(par_1)]
  sqrt(sum((par_1 - par_2)^2))
}

maximum_absolute_parameter_difference <- function(par_1, par_2) {
  if (is.null(par_1) || is.null(par_2)) {
    return(NA_real_)
  }

  if (!setequal(names(par_1), names(par_2))) {
    return(NA_real_)
  }

  par_2 <- par_2[names(par_1)]
  max(abs(par_1 - par_2))
}

# Summarize all optimizer runs from one stage/scenario.
#
# Two selections are retained:
#   1. global_best: the smallest negative log-likelihood;
#   2. baseline_connected: the successful solution closest to the preceding
#      sensitivity value's connected solution on the log-parameter scale.
#
# At sensitivity value zero, the connected solution is defined to equal the
# global-best solution.
summarize_optimizer_starts <- function(
    starts,
    fits,
    initial_objectives,
    connected_reference_par = NULL
) {
  start_labels <- names(starts)

  final_objectives <- map_dbl(
    fits,
    ~ if (is.null(.x)) NA_real_ else .x$value
  )

  successful <- is.finite(final_objectives) &
    final_objectives < bad_likelihood_penalty

  if (!any(successful)) {
    stop("All optimizer starts failed.")
  }

  successful_ids <- which(successful)
  ordered_ids <- successful_ids[order(final_objectives[successful_ids])]

  global_best_id <- ordered_ids[1]
  second_id <- if (length(ordered_ids) >= 2) ordered_ids[2] else NA_integer_
  global_best_fit <- fits[[global_best_id]]
  global_best_par <- global_best_fit$par

  distance_to_connected_reference <- rep(NA_real_, length(starts))

  if (!is.null(connected_reference_par)) {
    for (i in successful_ids) {
      distance_to_connected_reference[i] <- parameter_distance(
        fits[[i]]$par,
        connected_reference_par
      )
    }
  }

  if (
    is.null(connected_reference_par) ||
      all(!is.finite(distance_to_connected_reference[successful_ids]))
  ) {
    connected_id <- global_best_id
  } else {
    candidate_table <- tibble(
      start_id = successful_ids,
      distance = distance_to_connected_reference[successful_ids],
      objective = final_objectives[successful_ids]
    ) %>%
      filter(is.finite(distance)) %>%
      arrange(distance, objective)

    connected_id <- candidate_table$start_id[1]
  }

  connected_fit <- fits[[connected_id]]
  connected_par <- connected_fit$par

  start_table <- map_dfr(seq_along(starts), function(i) {
    fit_i <- fits[[i]]

    if (is.null(fit_i)) {
      return(tibble(
        start_id = i,
        start_label = start_labels[i],
        used_warm_start = str_starts(start_labels[i], "warm_start_"),
        used_global_warm_start = start_labels[i] == "warm_start_global_previous",
        used_connected_warm_start = start_labels[i] == "warm_start_connected_previous",
        initial_objective = initial_objectives[i],
        final_objective = NA_real_,
        optim_convergence = NA_integer_,
        optimizer_message = NA_character_,
        optimizer_failed = TRUE,
        is_best_start = FALSE,
        is_global_best_start = FALSE,
        is_connected_start = FALSE,
        distance_to_best_solution = NA_real_,
        max_abs_diff_to_best_solution = NA_real_,
        distance_to_connected_reference = NA_real_,
        distance_to_connected_solution = NA_real_,
        max_abs_diff_to_connected_solution = NA_real_
      ))
    }

    message_i <- if (is.null(fit_i$message)) NA_character_ else as.character(fit_i$message)

    tibble(
      start_id = i,
      start_label = start_labels[i],
      used_warm_start = str_starts(start_labels[i], "warm_start_"),
      used_global_warm_start = start_labels[i] == "warm_start_global_previous",
      used_connected_warm_start = start_labels[i] == "warm_start_connected_previous",
      initial_objective = initial_objectives[i],
      final_objective = fit_i$value,
      optim_convergence = fit_i$convergence,
      optimizer_message = message_i,
      optimizer_failed = FALSE,
      is_best_start = i == global_best_id,
      is_global_best_start = i == global_best_id,
      is_connected_start = i == connected_id,
      distance_to_best_solution = parameter_distance(fit_i$par, global_best_par),
      max_abs_diff_to_best_solution = maximum_absolute_parameter_difference(
        fit_i$par,
        global_best_par
      ),
      distance_to_connected_reference = distance_to_connected_reference[i],
      distance_to_connected_solution = parameter_distance(fit_i$par, connected_par),
      max_abs_diff_to_connected_solution = maximum_absolute_parameter_difference(
        fit_i$par,
        connected_par
      )
    )
  })

  if (is.na(second_id)) {
    second_best_objective <- NA_real_
    objective_gap <- NA_real_
    distance_best_to_second <- NA_real_
    max_abs_diff_best_to_second <- NA_real_
  } else {
    second_fit <- fits[[second_id]]
    second_best_objective <- second_fit$value
    objective_gap <- second_fit$value - global_best_fit$value
    distance_best_to_second <- parameter_distance(
      second_fit$par,
      global_best_par
    )
    max_abs_diff_best_to_second <- maximum_absolute_parameter_difference(
      second_fit$par,
      global_best_par
    )
  }

  successful_objectives <- final_objectives[successful]
  global_connected_distance <- parameter_distance(global_best_par, connected_par)
  global_connected_max_abs <- maximum_absolute_parameter_difference(
    global_best_par,
    connected_par
  )
  connected_objective_gap <- connected_fit$value - global_best_fit$value

  start_summary <- tibble(
    n_optimizer_starts = length(starts),
    n_successful_starts = sum(successful),
    best_start_id = global_best_id,
    best_start_label = start_labels[global_best_id],
    warm_start_available = any(str_starts(start_labels, "warm_start_")),
    warm_start_was_best = str_starts(start_labels[global_best_id], "warm_start_"),
    global_best_start_id = global_best_id,
    global_best_start_label = start_labels[global_best_id],
    connected_start_id = connected_id,
    connected_start_label = start_labels[connected_id],
    connected_reference_available = !is.null(connected_reference_par),
    best_objective = global_best_fit$value,
    connected_objective = connected_fit$value,
    connected_minus_global_objective = connected_objective_gap,
    global_connected_parameter_distance = global_connected_distance,
    global_connected_max_abs_parameter_difference = global_connected_max_abs,
    global_and_connected_same_solution =
      is.finite(global_connected_max_abs) &&
      global_connected_max_abs <= same_solution_max_abs_difference,
    second_best_objective = second_best_objective,
    objective_gap_second_minus_best = objective_gap,
    objective_range_across_successful_starts =
      max(successful_objectives) - min(successful_objectives),
    parameter_distance_best_to_second = distance_best_to_second,
    max_abs_parameter_diff_best_to_second = max_abs_diff_best_to_second
  )

  parameter_solutions <- map_dfr(successful_ids, function(i) {
    fit_i <- fits[[i]]

    tibble(
      start_id = i,
      start_label = start_labels[i],
      final_objective = fit_i$value,
      optim_convergence = fit_i$convergence,
      is_global_best_start = i == global_best_id,
      is_connected_start = i == connected_id,
      parameter = names(fit_i$par),
      estimate_log = as.numeric(fit_i$par),
      estimate = exp(as.numeric(fit_i$par))
    )
  })

  list(
    global_best_id = global_best_id,
    global_best_fit = global_best_fit,
    connected_id = connected_id,
    connected_fit = connected_fit,
    fits = fits,
    start_table = start_table,
    start_summary = start_summary,
    parameter_solutions = parameter_solutions
  )
}

build_two_car_result <- function(
    fit,
    selected_start_id,
    selected_start_label,
    selection_name,
    optimizer_runs,
    dat,
    cell_ids,
    mixing_model,
    sensitivity_value
) {
  if (!is.finite(fit$value) || fit$value >= bad_likelihood_penalty) {
    stop(
      "Selected two-car optimizer solution failed for ", mixing_model,
      " = ", sensitivity_value
    )
  }

  hess <- tryCatch(
    optimHess(
      par = fit$par,
      fn = negloglik_two_car,
      dat = dat,
      cell_ids = cell_ids,
      mixing_model = mixing_model,
      sensitivity_value = sensitivity_value
    ),
    error = function(e) NULL
  )

  se_log <- rep(NA_real_, length(fit$par))
  names(se_log) <- names(fit$par)
  hessian_ok <- FALSE

  if (!is.null(hess) && all(is.finite(hess))) {
    eigenvalues <- tryCatch(
      eigen(hess, symmetric = TRUE, only.values = TRUE)$values,
      error = function(e) NA_real_
    )

    if (all(is.finite(eigenvalues)) && min(eigenvalues) > 1e-8) {
      vcov_log <- tryCatch(solve(hess), error = function(e) NULL)

      if (!is.null(vcov_log) && all(is.finite(vcov_log))) {
        se_log <- sqrt(pmax(diag(vcov_log), 0))
        names(se_log) <- names(fit$par)
        hessian_ok <- TRUE
      }
    }
  }

  unpacked <- unpack_two_car_par(fit$par, cell_ids)

  probabilities <- make_two_car_probabilities(
    theta = unpacked$theta,
    r_mat = unpacked$r_mat,
    mixing_model = mixing_model,
    sensitivity_value = sensitivity_value
  )

  if (is.null(probabilities)) {
    stop("Selected two-car fit produced invalid probabilities.")
  }

  list(
    selection_name = selection_name,
    fit = fit,
    best_start_id = selected_start_id,
    best_start_label = selected_start_label,
    start_diagnostics = optimizer_runs$start_table,
    start_summary = optimizer_runs$start_summary,
    hessian_ok = hessian_ok,
    se_log = se_log,
    theta = unpacked$theta,
    r_mat = unpacked$r_mat,
    probabilities = probabilities,
    n_cells = nrow(dat),
    n2_total = sum(dat$n2),
    max_denominator_abs_diff = max(probabilities$denominator_abs_diff),
    min_validity_margin = min(probabilities$validity_margin),
    n_params_near_lower = sum(abs(fit$par - log_param_lower) < 1e-4),
    n_params_near_upper = sum(abs(fit$par - log_param_upper) < 1e-4)
  )
}

estimate_two_car_stage <- function(
    dat,
    mixing_model,
    sensitivity_value,
    warm_start_global = NULL,
    warm_start_connected = NULL,
    connected_reference_par = NULL
) {
  dat <- dat %>%
    filter(n2 >= min_n2) %>%
    arrange(cell_id)

  if (nrow(dat) == 0) {
    stop("No two-car cells remain after applying min_n2 = ", min_n2)
  }

  cell_ids <- dat$cell_id

  fixed_starts <- map(theta_start_grid, ~ make_two_car_start(dat, .x))
  names(fixed_starts) <- paste0("fixed_start_", seq_along(fixed_starts))

  starts <- add_warm_starts(
    fixed_starts = fixed_starts,
    warm_start_global = warm_start_global,
    warm_start_connected = warm_start_connected
  )

  initial_objectives <- map_dbl(starts, function(start_par) {
    negloglik_two_car(
      par = start_par,
      dat = dat,
      cell_ids = cell_ids,
      mixing_model = mixing_model,
      sensitivity_value = sensitivity_value
    )
  })

  fits <- map(starts, function(start_par) {
    tryCatch(
      optim(
        par = start_par,
        fn = negloglik_two_car,
        dat = dat,
        cell_ids = cell_ids,
        mixing_model = mixing_model,
        sensitivity_value = sensitivity_value,
        method = "L-BFGS-B",
        lower = rep(log_param_lower, length(start_par)),
        upper = rep(log_param_upper, length(start_par)),
        control = optim_control
      ),
      error = function(e) NULL
    )
  })

  optimizer_runs <- summarize_optimizer_starts(
    starts = starts,
    fits = fits,
    initial_objectives = initial_objectives,
    connected_reference_par = connected_reference_par
  )

  global_result <- build_two_car_result(
    fit = optimizer_runs$global_best_fit,
    selected_start_id = optimizer_runs$global_best_id,
    selected_start_label = names(starts)[optimizer_runs$global_best_id],
    selection_name = "global_best",
    optimizer_runs = optimizer_runs,
    dat = dat,
    cell_ids = cell_ids,
    mixing_model = mixing_model,
    sensitivity_value = sensitivity_value
  )

  connected_result <- build_two_car_result(
    fit = optimizer_runs$connected_fit,
    selected_start_id = optimizer_runs$connected_id,
    selected_start_label = names(starts)[optimizer_runs$connected_id],
    selection_name = "baseline_connected",
    optimizer_runs = optimizer_runs,
    dat = dat,
    cell_ids = cell_ids,
    mixing_model = mixing_model,
    sensitivity_value = sensitivity_value
  )

  list(
    global = global_result,
    connected = connected_result,
    optimizer_runs = optimizer_runs
  )
}

# ------------------------------------------------------------------------------
# 4. One-car stage
# ------------------------------------------------------------------------------

make_one_car_probabilities <- function(lambda, r_mat) {
  r_mat %>%
    transmute(
      cell_id,
      w_S = lambda["S"] * r_S,
      w_L = lambda["L"] * r_L,
      w_M = lambda["M"] * r_M,
      w_H = lambda["H"] * r_H
    ) %>%
    mutate(D = w_S + w_L + w_M + w_H) %>%
    transmute(
      cell_id,
      q_S = w_S / D,
      q_L = w_L / D,
      q_M = w_M / D,
      q_H = w_H / D
    )
}

negloglik_one_car <- function(par, dat, r_mat) {
  lambda <- c(
    S = 1,
    L = unname(exp(par["log_lambda_L"])),
    M = unname(exp(par["log_lambda_M"])),
    H = unname(exp(par["log_lambda_H"]))
  )

  probs <- make_one_car_probabilities(lambda, r_mat)
  ll_dat <- dat %>% left_join(probs, by = "cell_id")

  prob_cols <- c("q_S", "q_L", "q_M", "q_H")
  if (any(!is.finite(as.matrix(ll_dat[, prob_cols])))) {
    return(bad_likelihood_penalty)
  }

  loglik_terms <-
    ll_dat$one_S * safe_log(ll_dat$q_S) +
    ll_dat$one_L * safe_log(ll_dat$q_L) +
    ll_dat$one_M * safe_log(ll_dat$q_M) +
    ll_dat$one_H * safe_log(ll_dat$q_H)

  if (any(!is.finite(loglik_terms))) {
    return(bad_likelihood_penalty)
  }

  -sum(loglik_terms)
}

make_one_car_start <- function(lambda_start) {
  c(
    log_lambda_L = log(unname(lambda_start["lambda_L"])),
    log_lambda_M = log(unname(lambda_start["lambda_M"])),
    log_lambda_H = log(unname(lambda_start["lambda_H"]))
  )
}

build_one_car_result <- function(
    fit,
    selected_start_id,
    selected_start_label,
    selection_name,
    optimizer_runs,
    dat,
    r_mat
) {
  hess <- tryCatch(
    optimHess(
      par = fit$par,
      fn = negloglik_one_car,
      dat = dat,
      r_mat = r_mat
    ),
    error = function(e) NULL
  )

  se_log <- rep(NA_real_, length(fit$par))
  names(se_log) <- names(fit$par)
  hessian_ok <- FALSE

  if (!is.null(hess) && all(is.finite(hess))) {
    eigenvalues <- tryCatch(
      eigen(hess, symmetric = TRUE, only.values = TRUE)$values,
      error = function(e) NA_real_
    )

    if (all(is.finite(eigenvalues)) && min(eigenvalues) > 1e-8) {
      vcov_log <- tryCatch(solve(hess), error = function(e) NULL)

      if (!is.null(vcov_log) && all(is.finite(vcov_log))) {
        se_log <- sqrt(pmax(diag(vcov_log), 0))
        names(se_log) <- names(fit$par)
        hessian_ok <- TRUE
      }
    }
  }

  lambda <- c(
    S = 1,
    L = unname(exp(fit$par["log_lambda_L"])),
    M = unname(exp(fit$par["log_lambda_M"])),
    H = unname(exp(fit$par["log_lambda_H"]))
  )

  list(
    selection_name = selection_name,
    fit = fit,
    best_start_id = selected_start_id,
    best_start_label = selected_start_label,
    start_diagnostics = optimizer_runs$start_table,
    start_summary = optimizer_runs$start_summary,
    hessian_ok = hessian_ok,
    se_log = se_log,
    lambda = lambda,
    n_cells = nrow(dat),
    n1_total = sum(dat$n1),
    n_params_near_lower = sum(abs(fit$par - log_param_lower) < 1e-4),
    n_params_near_upper = sum(abs(fit$par - log_param_upper) < 1e-4)
  )
}

estimate_one_car_stage <- function(
    dat,
    r_mat,
    warm_start_global = NULL,
    warm_start_connected = NULL,
    connected_reference_par = NULL
) {
  dat <- dat %>%
    filter(n1 > 0, cell_id %in% r_mat$cell_id) %>%
    arrange(cell_id)

  r_mat <- r_mat %>%
    filter(cell_id %in% dat$cell_id) %>%
    arrange(cell_id)

  fixed_starts <- map(lambda_start_grid, make_one_car_start)
  names(fixed_starts) <- paste0("fixed_start_", seq_along(fixed_starts))

  starts <- add_warm_starts(
    fixed_starts = fixed_starts,
    warm_start_global = warm_start_global,
    warm_start_connected = warm_start_connected
  )

  initial_objectives <- map_dbl(starts, function(start_par) {
    negloglik_one_car(
      par = start_par,
      dat = dat,
      r_mat = r_mat
    )
  })

  fits <- map(starts, function(start_par) {
    tryCatch(
      optim(
        par = start_par,
        fn = negloglik_one_car,
        dat = dat,
        r_mat = r_mat,
        method = "L-BFGS-B",
        lower = rep(log_param_lower, length(start_par)),
        upper = rep(log_param_upper, length(start_par)),
        control = optim_control
      ),
      error = function(e) NULL
    )
  })

  optimizer_runs <- summarize_optimizer_starts(
    starts = starts,
    fits = fits,
    initial_objectives = initial_objectives,
    connected_reference_par = connected_reference_par
  )

  global_result <- build_one_car_result(
    fit = optimizer_runs$global_best_fit,
    selected_start_id = optimizer_runs$global_best_id,
    selected_start_label = names(starts)[optimizer_runs$global_best_id],
    selection_name = "global_best",
    optimizer_runs = optimizer_runs,
    dat = dat,
    r_mat = r_mat
  )

  connected_result <- build_one_car_result(
    fit = optimizer_runs$connected_fit,
    selected_start_id = optimizer_runs$connected_id,
    selected_start_label = names(starts)[optimizer_runs$connected_id],
    selection_name = "baseline_connected",
    optimizer_runs = optimizer_runs,
    dat = dat,
    r_mat = r_mat
  )

  list(
    global = global_result,
    connected = connected_result,
    optimizer_runs = optimizer_runs
  )
}

estimate_one_imp_scenario <- function(
    dat,
    mixing_model,
    sensitivity_value,
    warm_start_two_car_global = NULL,
    warm_start_two_car_connected = NULL,
    warm_start_one_car_global = NULL,
    warm_start_one_car_connected = NULL
) {
  two_stage <- estimate_two_car_stage(
    dat = dat,
    mixing_model = mixing_model,
    sensitivity_value = sensitivity_value,
    warm_start_global = warm_start_two_car_global,
    warm_start_connected = warm_start_two_car_connected,
    connected_reference_par = warm_start_two_car_connected
  )

  # The global branch uses the global-best exposure ratios and the global-best
  # conditional one-car fit.
  one_stage_global_r <- estimate_one_car_stage(
    dat = dat,
    r_mat = two_stage$global$r_mat,
    warm_start_global = warm_start_one_car_global,
    warm_start_connected = warm_start_one_car_connected,
    connected_reference_par = warm_start_one_car_global
  )

  # The connected branch uses the connected exposure ratios and follows the
  # previous connected one-car solution when multiple one-car solutions exist.
  one_stage_connected_r <- estimate_one_car_stage(
    dat = dat,
    r_mat = two_stage$connected$r_mat,
    warm_start_global = warm_start_one_car_global,
    warm_start_connected = warm_start_one_car_connected,
    connected_reference_par = warm_start_one_car_connected
  )

  global_fit <- list(
    two_fit = two_stage$global,
    one_fit = one_stage_global_r$global
  )

  connected_fit <- list(
    two_fit = two_stage$connected,
    one_fit = one_stage_connected_r$connected
  )

  list(
    # These aliases preserve all of the original global-best outputs below.
    two_fit = global_fit$two_fit,
    one_fit = global_fit$one_fit,
    global = global_fit,
    connected = connected_fit,
    two_optimizer = two_stage$optimizer_runs,
    one_global_optimizer = one_stage_global_r$optimizer_runs,
    one_connected_optimizer = one_stage_connected_r$optimizer_runs
  )
}

# ------------------------------------------------------------------------------
# 5. Pool positive estimates on the log scale
# ------------------------------------------------------------------------------

pool_log_estimates <- function(df, group_vars) {
  df %>%
    group_by(across(all_of(group_vars))) %>%
    summarize(
      m_total = n(),
      m_est = sum(is.finite(estimate_log)),
      m_se = sum(is.finite(se_log)),
      estimate_log_bar = mean(estimate_log[is.finite(estimate_log)]),
      within_var_log = if_else(
        m_se == 0,
        NA_real_,
        mean(se_log[is.finite(se_log)]^2)
      ),
      between_var_log = if_else(
        m_est <= 1,
        0,
        var(estimate_log[is.finite(estimate_log)])
      ),
      total_var_log = if_else(
        is.na(within_var_log),
        NA_real_,
        within_var_log + (1 + 1 / m_total) * between_var_log
      ),
      se_log_total = sqrt(total_var_log),
      estimate = exp(estimate_log_bar),
      ci_low = if_else(
        is.na(se_log_total),
        NA_real_,
        exp(estimate_log_bar - 1.96 * se_log_total)
      ),
      ci_high = if_else(
        is.na(se_log_total),
        NA_real_,
        exp(estimate_log_bar + 1.96 * se_log_total)
      ),
      mean_n_crashes = mean(n_crashes, na.rm = TRUE),
      .groups = "drop"
    )
}

# ------------------------------------------------------------------------------
# 6. Load Script 10 crash-count data
# ------------------------------------------------------------------------------

if (!file.exists(counts_path)) {
  stop(
    "Missing input file: ", counts_path,
    "\nPlease run 10_four_group_prep.R first."
  )
}

counts_all <- readRDS(counts_path)

if ("imputation" %in% names(counts_all) && !"IMPUTATION" %in% names(counts_all)) {
  counts_all <- counts_all %>% rename(IMPUTATION = imputation)
}

counts_all <- counts_all %>%
  mutate(cell_id = as.character(cell_id)) %>%
  filter(cell_spec %in% cell_specs_to_estimate)

required_count_cols <- c(
  "IMPUTATION", "cell_spec", "cell_id",
  two_car_types,
  "one_S", "one_L", "one_M", "one_H", "n2", "n1"
)

missing_count_cols <- setdiff(required_count_cols, names(counts_all))
if (length(missing_count_cols) > 0) {
  stop("Missing required count columns: ", paste(missing_count_cols, collapse = ", "))
}

cat("Loaded Script 10 crash-count data:\n")
cat("  Rows:", nrow(counts_all), "\n")
cat("  Cell specs:", paste(sort(unique(counts_all$cell_spec)), collapse = ", "), "\n")
cat("  Imputations:", paste(sort(unique(counts_all$IMPUTATION)), collapse = ", "), "\n")
cat("  Sensitivity scenarios:", nrow(scenario_grid), "\n\n")

# ------------------------------------------------------------------------------
# 7. Run every scenario separately in every imputed data set
# ------------------------------------------------------------------------------

all_results <- list()
result_index <- 0

for (cell_spec_i in sort(unique(counts_all$cell_spec))) {
  for (imputation_i in sort(unique(counts_all$IMPUTATION))) {

    dat_i <- counts_all %>%
      filter(
        IMPUTATION == imputation_i,
        cell_spec == cell_spec_i
      ) %>%
      arrange(cell_id)

    for (mixing_model_i in unique(scenario_grid$mixing_model)) {

      scenarios_i <- scenario_grid %>%
        filter(mixing_model == mixing_model_i) %>%
        arrange(sensitivity_value)

      # Keep separate warm starts for the global-best and baseline-connected
      # branches. Both are offered to the optimizer at the next value.
      previous_two_car_global_par <- NULL
      previous_two_car_connected_par <- NULL
      previous_one_car_global_par <- NULL
      previous_one_car_connected_par <- NULL

      for (j in seq_len(nrow(scenarios_i))) {
        scenario_i <- scenarios_i[j, ]

        cat(
          "Estimating:", scenario_i$scenario_label,
          "| cell_spec =", cell_spec_i,
          "| imputation =", imputation_i,
          "| previous branches available =",
          !is.null(previous_two_car_global_par), "\n"
        )

        fit_i <- estimate_one_imp_scenario(
          dat = dat_i,
          mixing_model = scenario_i$mixing_model,
          sensitivity_value = scenario_i$sensitivity_value,
          warm_start_two_car_global = previous_two_car_global_par,
          warm_start_two_car_connected = previous_two_car_connected_par,
          warm_start_one_car_global = previous_one_car_global_par,
          warm_start_one_car_connected = previous_one_car_connected_par
        )

        result_index <- result_index + 1

        all_results[[result_index]] <- list(
          IMPUTATION = imputation_i,
          cell_spec = cell_spec_i,
          mixing_model = scenario_i$mixing_model,
          sensitivity_name = scenario_i$sensitivity_name,
          sensitivity_value = scenario_i$sensitivity_value,
          scenario_id = scenario_i$scenario_id,
          scenario_label = scenario_i$scenario_label,
          data = dat_i,
          fit = fit_i
        )

        previous_two_car_global_par <- fit_i$global$two_fit$fit$par
        previous_two_car_connected_par <- fit_i$connected$two_fit$fit$par
        previous_one_car_global_par <- fit_i$global$one_fit$fit$par
        previous_one_car_connected_par <- fit_i$connected$one_fit$fit$par
      }
    }
  }
}

# ------------------------------------------------------------------------------
# 8. Tidy theta and lambda estimates
# ------------------------------------------------------------------------------

parameter_estimates_by_imp <- map_dfr(all_results, function(res) {
  two_fit <- res$fit$two_fit
  one_fit <- res$fit$one_fit

  common <- tibble(
    IMPUTATION = res$IMPUTATION,
    cell_spec = res$cell_spec,
    mixing_model = res$mixing_model,
    sensitivity_name = res$sensitivity_name,
    sensitivity_value = res$sensitivity_value,
    scenario_id = res$scenario_id,
    scenario_label = res$scenario_label
  )

  bind_rows(
    common %>%
      slice(rep(1, 3)) %>%
      mutate(
        stage = "two_car",
        parameter = c("theta_L", "theta_M", "theta_H"),
        estimate = as.numeric(two_fit$theta[c("L", "M", "H")]),
        estimate_log = log(estimate),
        se_log = as.numeric(
          two_fit$se_log[c("log_theta_L", "log_theta_M", "log_theta_H")]
        ),
        se_delta = estimate * se_log,
        n_cells = two_fit$n_cells,
        n_crashes = two_fit$n2_total
      ),
    common %>%
      slice(rep(1, 3)) %>%
      mutate(
        stage = "one_car",
        parameter = c("lambda_L", "lambda_M", "lambda_H"),
        estimate = as.numeric(one_fit$lambda[c("L", "M", "H")]),
        estimate_log = log(estimate),
        se_log = as.numeric(
          one_fit$se_log[c("log_lambda_L", "log_lambda_M", "log_lambda_H")]
        ),
        se_delta = estimate * se_log,
        n_cells = one_fit$n_cells,
        n_crashes = one_fit$n1_total
      )
  )
})

parameter_estimates_pooled <- pool_log_estimates(
  parameter_estimates_by_imp,
  group_vars = c(
    "mixing_model", "sensitivity_name", "sensitivity_value",
    "scenario_id", "scenario_label", "cell_spec", "stage", "parameter"
  )
)

# Delta = 0 is the common equal-mixing reference.
parameter_baseline <- parameter_estimates_pooled %>%
  filter(mixing_model == equal_mixing_reference_model, sensitivity_value == 0) %>%
  select(
    cell_spec, stage, parameter,
    baseline_estimate = estimate,
    baseline_ci_low = ci_low,
    baseline_ci_high = ci_high
  )

parameter_percent_changes <- parameter_estimates_pooled %>%
  left_join(parameter_baseline, by = c("cell_spec", "stage", "parameter")) %>%
  mutate(
    percent_change_from_equal_mixing =
      100 * (estimate / baseline_estimate - 1)
  ) %>%
  arrange(mixing_model, sensitivity_value, cell_spec, stage, parameter)

# ------------------------------------------------------------------------------
# 9. Tidy and pool cell-specific exposure ratios r_g,c
# ------------------------------------------------------------------------------

r_estimates_by_imp <- map_dfr(all_results, function(res) {
  two_fit <- res$fit$two_fit

  map_dfr(c("L", "M", "H"), function(group_i) {
    parameter_i <- paste0("r_", group_i)

    tibble(
      IMPUTATION = res$IMPUTATION,
      cell_spec = res$cell_spec,
      mixing_model = res$mixing_model,
      sensitivity_name = res$sensitivity_name,
      sensitivity_value = res$sensitivity_value,
      scenario_id = res$scenario_id,
      scenario_label = res$scenario_label,
      cell_id = two_fit$r_mat$cell_id,
      BAC_group = group_i,
      parameter = parameter_i,
      estimate = two_fit$r_mat[[parameter_i]],
      estimate_log = log(estimate),
      se_log = as.numeric(
        two_fit$se_log[paste0("log_", parameter_i, "__", two_fit$r_mat$cell_id)]
      ),
      se_delta = estimate * se_log,
      n_cells = two_fit$n_cells,
      n_crashes = two_fit$n2_total
    )
  })
})

r_estimates_pooled <- pool_log_estimates(
  r_estimates_by_imp,
  group_vars = c(
    "mixing_model", "sensitivity_name", "sensitivity_value",
    "scenario_id", "scenario_label", "cell_spec", "cell_id",
    "BAC_group", "parameter"
  )
)

r_baseline <- r_estimates_pooled %>%
  filter(mixing_model == equal_mixing_reference_model, sensitivity_value == 0) %>%
  select(
    cell_spec, cell_id, BAC_group, parameter,
    baseline_estimate = estimate
  )

r_percent_changes <- r_estimates_pooled %>%
  left_join(
    r_baseline,
    by = c("cell_spec", "cell_id", "BAC_group", "parameter")
  ) %>%
  mutate(
    percent_change_from_equal_mixing =
      100 * (estimate / baseline_estimate - 1)
  ) %>%
  arrange(mixing_model, sensitivity_value, cell_spec, cell_id, BAC_group)

# ------------------------------------------------------------------------------
# 10. Road-share summaries implied by r
# ------------------------------------------------------------------------------

road_shares_by_imp <- map_dfr(all_results, function(res) {
  res$fit$two_fit$r_mat %>%
    mutate(
      R_positive = r_L + r_M + r_H,
      T = 1 + R_positive,
      p_S = 1 / T,
      p_L = r_L / T,
      p_M = r_M / T,
      p_H = r_H / T,
      IMPUTATION = res$IMPUTATION,
      cell_spec = res$cell_spec,
      mixing_model = res$mixing_model,
      sensitivity_name = res$sensitivity_name,
      sensitivity_value = res$sensitivity_value,
      scenario_id = res$scenario_id,
      scenario_label = res$scenario_label,
      .before = 1
    ) %>%
    select(
      IMPUTATION, cell_spec, mixing_model, sensitivity_name,
      sensitivity_value, scenario_id, scenario_label, cell_id,
      R_positive, T, p_S, p_L, p_M, p_H
    ) %>%
    pivot_longer(
      cols = c(p_S, p_L, p_M, p_H),
      names_to = "BAC_group",
      values_to = "road_share"
    ) %>%
    mutate(BAC_group = str_remove(BAC_group, "^p_"))
})

road_share_summary <- road_shares_by_imp %>%
  group_by(
    mixing_model, sensitivity_name, sensitivity_value,
    scenario_id, scenario_label, cell_spec, cell_id, BAC_group
  ) %>%
  summarize(
    mean_road_share = mean(road_share),
    sd_across_imputations = sd(road_share),
    min_road_share = min(road_share),
    max_road_share = max(road_share),
    .groups = "drop"
  )

road_share_baseline <- road_share_summary %>%
  filter(mixing_model == equal_mixing_reference_model, sensitivity_value == 0) %>%
  select(
    cell_spec, cell_id, BAC_group,
    baseline_mean_road_share = mean_road_share
  )

road_share_changes <- road_share_summary %>%
  left_join(
    road_share_baseline,
    by = c("cell_spec", "cell_id", "BAC_group")
  ) %>%
  mutate(
    percentage_point_change =
      100 * (mean_road_share - baseline_mean_road_share),
    percent_change_from_equal_mixing =
      100 * (mean_road_share / baseline_mean_road_share - 1)
  ) %>%
  arrange(mixing_model, sensitivity_value, cell_spec, cell_id, BAC_group)

# ------------------------------------------------------------------------------
# 11. Optimizer and numerical diagnostics
# ------------------------------------------------------------------------------

fit_diagnostics <- map_dfr(all_results, function(res) {
  two_fit <- res$fit$two_fit
  one_fit <- res$fit$one_fit

  two_start <- two_fit$start_summary
  one_start <- one_fit$start_summary

  bind_rows(
    tibble(
      IMPUTATION = res$IMPUTATION,
      cell_spec = res$cell_spec,
      mixing_model = res$mixing_model,
      sensitivity_name = res$sensitivity_name,
      sensitivity_value = res$sensitivity_value,
      scenario_id = res$scenario_id,
      scenario_label = res$scenario_label,
      stage = "two_car",
      n_cells = two_fit$n_cells,
      n_crashes = two_fit$n2_total,
      negloglik = two_fit$fit$value,
      optim_convergence = two_fit$fit$convergence,
      hessian_ok = two_fit$hessian_ok,
      best_start_id = two_fit$best_start_id,
      best_start_label = two_fit$best_start_label,
      warm_start_available = two_start$warm_start_available,
      warm_start_was_best = two_start$warm_start_was_best,
      n_optimizer_starts = two_start$n_optimizer_starts,
      n_successful_starts = two_start$n_successful_starts,
      second_best_negloglik = two_start$second_best_objective,
      objective_gap_second_minus_best =
        two_start$objective_gap_second_minus_best,
      objective_range_across_successful_starts =
        two_start$objective_range_across_successful_starts,
      parameter_distance_best_to_second =
        two_start$parameter_distance_best_to_second,
      max_abs_parameter_diff_best_to_second =
        two_start$max_abs_parameter_diff_best_to_second,
      n_params_near_lower = two_fit$n_params_near_lower,
      n_params_near_upper = two_fit$n_params_near_upper,
      max_denominator_abs_diff = two_fit$max_denominator_abs_diff,
      min_validity_margin = two_fit$min_validity_margin
    ),
    tibble(
      IMPUTATION = res$IMPUTATION,
      cell_spec = res$cell_spec,
      mixing_model = res$mixing_model,
      sensitivity_name = res$sensitivity_name,
      sensitivity_value = res$sensitivity_value,
      scenario_id = res$scenario_id,
      scenario_label = res$scenario_label,
      stage = "one_car",
      n_cells = one_fit$n_cells,
      n_crashes = one_fit$n1_total,
      negloglik = one_fit$fit$value,
      optim_convergence = one_fit$fit$convergence,
      hessian_ok = one_fit$hessian_ok,
      best_start_id = one_fit$best_start_id,
      best_start_label = one_fit$best_start_label,
      warm_start_available = one_start$warm_start_available,
      warm_start_was_best = one_start$warm_start_was_best,
      n_optimizer_starts = one_start$n_optimizer_starts,
      n_successful_starts = one_start$n_successful_starts,
      second_best_negloglik = one_start$second_best_objective,
      objective_gap_second_minus_best =
        one_start$objective_gap_second_minus_best,
      objective_range_across_successful_starts =
        one_start$objective_range_across_successful_starts,
      parameter_distance_best_to_second =
        one_start$parameter_distance_best_to_second,
      max_abs_parameter_diff_best_to_second =
        one_start$max_abs_parameter_diff_best_to_second,
      n_params_near_lower = one_fit$n_params_near_lower,
      n_params_near_upper = one_fit$n_params_near_upper,
      max_denominator_abs_diff = NA_real_,
      min_validity_margin = NA_real_
    )
  )
})

# Detailed results for every fixed and warm optimizer start.
optimizer_start_diagnostics <- map_dfr(all_results, function(res) {
  common <- tibble(
    IMPUTATION = res$IMPUTATION,
    cell_spec = res$cell_spec,
    mixing_model = res$mixing_model,
    sensitivity_name = res$sensitivity_name,
    sensitivity_value = res$sensitivity_value,
    scenario_id = res$scenario_id,
    scenario_label = res$scenario_label
  )

  bind_rows(
    res$fit$two_fit$start_diagnostics %>%
      mutate(stage = "two_car", .before = 1),
    res$fit$one_fit$start_diagnostics %>%
      mutate(stage = "one_car", .before = 1)
  ) %>%
    crossing(common) %>%
    select(
      IMPUTATION, cell_spec, mixing_model, sensitivity_name,
      sensitivity_value, scenario_id, scenario_label, stage,
      everything()
    )
})

# Delta validity is 1 - Delta * (r_L + r_M + r_H) > 0.
delta_feasibility_diagnostics <- road_shares_by_imp %>%
  filter(mixing_model == "sober_positive_delta") %>%
  distinct(
    IMPUTATION, cell_spec, sensitivity_value, scenario_id,
    scenario_label, cell_id, R_positive
  ) %>%
  mutate(
    validity_margin = 1 - sensitivity_value * R_positive,
    maximum_valid_delta = safe_divide(1, R_positive),
    valid = validity_margin > 0
  ) %>%
  arrange(sensitivity_value, cell_spec, IMPUTATION, cell_id)

# ------------------------------------------------------------------------------
# 12. Diagnostics for choosing rho
# ------------------------------------------------------------------------------

# 12a. Direct count sparsity by two-car crash type.
rho_count_sparsity <- counts_all %>%
  select(IMPUTATION, cell_spec, cell_id, all_of(two_car_types)) %>%
  pivot_longer(
    cols = all_of(two_car_types),
    names_to = "crash_type",
    values_to = "count"
  ) %>%
  mutate(same_group = crash_type %in% same_group_types) %>%
  group_by(cell_spec, crash_type, same_group) %>%
  summarize(
    n_cell_imputation_rows = n(),
    average_count_per_cell = mean(count),
    median_count_per_cell = median(count),
    minimum_count_per_cell = min(count),
    maximum_count_per_cell = max(count),
    percent_cell_imputation_rows_zero = 100 * mean(count == 0),
    average_total_count_per_imputation =
      mean(tapply(count, IMPUTATION, sum)),
    .groups = "drop"
  ) %>%
  arrange(desc(same_group), crash_type)

# 12b. Formula-based effect of candidate rho values, evaluated at the equal-
# mixing road shares. This does not require refitting each candidate rho.
equal_mixing_shares <- road_shares_by_imp %>%
  filter(mixing_model == equal_mixing_reference_model, sensitivity_value == 0) %>%
  select(IMPUTATION, cell_spec, cell_id, BAC_group, road_share)

rho_effect_by_cell_group <- crossing(
  equal_mixing_shares,
  rho = rho_diagnostic_grid
) %>%
  mutate(
    independent_diagonal_probability = road_share^2,
    rho_diagonal_probability =
      (1 - rho) * road_share^2 + rho * road_share,
    absolute_diagonal_increase =
      rho_diagonal_probability - independent_diagonal_probability,
    relative_diagonal_increase =
      safe_divide(
        absolute_diagonal_increase,
        independent_diagonal_probability
      ),
    percent_diagonal_increase = 100 * relative_diagonal_increase,
    percent_reduction_in_each_mixed_probability = 100 * rho
  ) %>%
  arrange(rho, BAC_group, cell_spec, IMPUTATION, cell_id)

rho_effect_summary <- rho_effect_by_cell_group %>%
  group_by(rho, cell_spec, BAC_group) %>%
  summarize(
    minimum_road_share = min(road_share),
    median_road_share = median(road_share),
    maximum_road_share = max(road_share),
    minimum_percent_diagonal_increase = min(percent_diagonal_increase),
    median_percent_diagonal_increase = median(percent_diagonal_increase),
    p90_percent_diagonal_increase =
      as.numeric(quantile(percent_diagonal_increase, 0.90)),
    maximum_percent_diagonal_increase = max(percent_diagonal_increase),
    median_absolute_diagonal_increase = median(absolute_diagonal_increase),
    percent_reduction_in_each_mixed_probability = first(
      percent_reduction_in_each_mixed_probability
    ),
    .groups = "drop"
  ) %>%
  arrange(rho, BAC_group, cell_spec)

# 12c. rho values corresponding to selected target relative diagonal increases.
# From target = rho * (1-p)/p, solve rho = target * p/(1-p).
rho_reference_for_target_increases <- crossing(
  equal_mixing_shares,
  target_relative_increase = target_diagonal_increases
) %>%
  mutate(
    target_percent_increase = 100 * target_relative_increase,
    rho_needed = target_relative_increase * road_share / (1 - road_share)
  ) %>%
  group_by(cell_spec, BAC_group, target_percent_increase) %>%
  summarize(
    minimum_rho_needed = min(rho_needed),
    median_rho_needed = median(rho_needed),
    maximum_rho_needed = max(rho_needed),
    .groups = "drop"
  ) %>%
  arrange(BAC_group, target_percent_increase, cell_spec)

# ------------------------------------------------------------------------------
# 13. Fitted two-car counts and residual diagnostics
# ------------------------------------------------------------------------------

fitted_two_car_counts_by_type <- map_dfr(all_results, function(res) {
  observed_long <- res$data %>%
    filter(n2 >= min_n2) %>%
    select(cell_id, n2, all_of(two_car_types)) %>%
    pivot_longer(
      cols = all_of(two_car_types),
      names_to = "crash_type",
      values_to = "observed_count"
    )

  probability_long <- res$fit$two_fit$probabilities %>%
    select(cell_id, starts_with("p_")) %>%
    pivot_longer(
      cols = starts_with("p_"),
      names_to = "crash_type",
      values_to = "fitted_probability"
    ) %>%
    mutate(crash_type = str_remove(crash_type, "^p_"))

  observed_long %>%
    left_join(probability_long, by = c("cell_id", "crash_type")) %>%
    mutate(
      expected_count = n2 * fitted_probability,
      observed_to_expected = safe_divide(observed_count, expected_count),
      pearson_residual = safe_divide(
        observed_count - expected_count,
        sqrt(expected_count)
      ),
      same_group = crash_type %in% same_group_types,
      IMPUTATION = res$IMPUTATION,
      cell_spec = res$cell_spec,
      mixing_model = res$mixing_model,
      sensitivity_name = res$sensitivity_name,
      sensitivity_value = res$sensitivity_value,
      scenario_id = res$scenario_id,
      scenario_label = res$scenario_label,
      .before = 1
    ) %>%
    select(
      IMPUTATION, cell_spec, mixing_model, sensitivity_name,
      sensitivity_value, scenario_id, scenario_label, cell_id,
      crash_type, same_group, n2, observed_count, fitted_probability,
      expected_count, observed_to_expected, pearson_residual
    )
})

# Summary by crash type across all cells and imputations. This is useful for
# tables such as the LL and MM observed-to-expected comparisons.
fitted_two_car_summary_by_type <- fitted_two_car_counts_by_type %>%
  group_by(
    mixing_model, sensitivity_name, sensitivity_value,
    scenario_id, scenario_label, cell_spec, crash_type, same_group
  ) %>%
  summarize(
    mean_observed_count = mean(observed_count),
    mean_expected_count = mean(expected_count),
    ratio_of_total_observed_to_expected =
      safe_divide(sum(observed_count), sum(expected_count)),
    mean_pearson_residual = mean(pearson_residual, na.rm = TRUE),
    mean_absolute_pearson_residual =
      mean(abs(pearson_residual), na.rm = TRUE),
    rms_pearson_residual =
      sqrt(mean(pearson_residual^2, na.rm = TRUE)),
    maximum_absolute_pearson_residual =
      max(abs(pearson_residual), na.rm = TRUE),
    percent_cells_with_zero_observed = 100 * mean(observed_count == 0),
    .groups = "drop"
  ) %>%
  arrange(mixing_model, sensitivity_value, desc(same_group), crash_type)

# First summarize observed/expected counts separately within every imputation.
# This avoids treating the 10 imputed data sets as independent copies.
two_car_fit_by_type_by_imputation <- fitted_two_car_counts_by_type %>%
  group_by(
    IMPUTATION, cell_spec, mixing_model, sensitivity_name,
    sensitivity_value, scenario_id, scenario_label,
    crash_type, same_group
  ) %>%
  summarize(
    total_observed_count = sum(observed_count),
    total_expected_count = sum(expected_count),
    observed_to_expected =
      safe_divide(total_observed_count, total_expected_count),
    mean_absolute_pearson_residual =
      mean(abs(pearson_residual), na.rm = TRUE),
    rms_pearson_residual =
      sqrt(mean(pearson_residual^2, na.rm = TRUE)),
    maximum_absolute_pearson_residual =
      max(abs(pearson_residual), na.rm = TRUE),
    .groups = "drop"
  )

# Residual summaries over every cell-by-crash-type contribution.
residual_summary_by_imputation <- fitted_two_car_counts_by_type %>%
  group_by(
    IMPUTATION, cell_spec, mixing_model, sensitivity_name,
    sensitivity_value, scenario_id, scenario_label
  ) %>%
  summarize(
    mean_absolute_pearson_residual =
      mean(abs(pearson_residual), na.rm = TRUE),
    rms_pearson_residual =
      sqrt(mean(pearson_residual^2, na.rm = TRUE)),
    maximum_absolute_pearson_residual =
      max(abs(pearson_residual), na.rm = TRUE),
    .groups = "drop"
  )

# RMS log observed/expected discrepancy across the 10 crash types.
log_oe_summary_by_imputation <- two_car_fit_by_type_by_imputation %>%
  group_by(
    IMPUTATION, cell_spec, mixing_model, sensitivity_name,
    sensitivity_value, scenario_id, scenario_label
  ) %>%
  summarize(
    rms_log_observed_to_expected = {
      ratios <- observed_to_expected[
        is.finite(observed_to_expected) & observed_to_expected > 0
      ]
      if (length(ratios) == 0) NA_real_ else sqrt(mean(log(ratios)^2))
    },
    .groups = "drop"
  )

# Put the crash-type observed/expected ratios into separate columns.
oe_ratios_wide <- two_car_fit_by_type_by_imputation %>%
  select(
    IMPUTATION, cell_spec, mixing_model, sensitivity_name,
    sensitivity_value, scenario_id, scenario_label,
    crash_type, observed_to_expected
  ) %>%
  mutate(metric_name = paste0("observed_to_expected_", crash_type)) %>%
  select(-crash_type) %>%
  pivot_wider(
    names_from = metric_name,
    values_from = observed_to_expected
  )

# Add each scenario's two-car likelihood and compare it with the zero-sensitivity
# scenario from the same mixing model and imputation.
two_car_nll <- fit_diagnostics %>%
  filter(stage == "two_car") %>%
  select(
    IMPUTATION, cell_spec, mixing_model, sensitivity_name,
    sensitivity_value, scenario_id, scenario_label,
    negloglik
  )

baseline_nll <- two_car_nll %>%
  filter(sensitivity_value == 0) %>%
  select(
    IMPUTATION, cell_spec, mixing_model,
    baseline_negloglik = negloglik
  )

two_car_fit_summary_by_imputation <- residual_summary_by_imputation %>%
  left_join(
    log_oe_summary_by_imputation,
    by = c(
      "IMPUTATION", "cell_spec", "mixing_model", "sensitivity_name",
      "sensitivity_value", "scenario_id", "scenario_label"
    )
  ) %>%
  left_join(
    oe_ratios_wide,
    by = c(
      "IMPUTATION", "cell_spec", "mixing_model", "sensitivity_name",
      "sensitivity_value", "scenario_id", "scenario_label"
    )
  ) %>%
  left_join(
    two_car_nll,
    by = c(
      "IMPUTATION", "cell_spec", "mixing_model", "sensitivity_name",
      "sensitivity_value", "scenario_id", "scenario_label"
    )
  ) %>%
  left_join(
    baseline_nll,
    by = c("IMPUTATION", "cell_spec", "mixing_model")
  ) %>%
  mutate(
    negloglik_change_from_model_zero = negloglik - baseline_negloglik
  ) %>%
  arrange(mixing_model, sensitivity_value, cell_spec, IMPUTATION)

# Compact scenario-level table requested for comparing candidate rho values.
scenario_level_fit_summary <- two_car_fit_summary_by_imputation %>%
  group_by(
    mixing_model, sensitivity_name, sensitivity_value,
    scenario_id, scenario_label, cell_spec
  ) %>%
  summarize(
    n_imputations = n(),
    mean_two_car_negloglik = mean(negloglik, na.rm = TRUE),
    sd_two_car_negloglik = sd(negloglik, na.rm = TRUE),
    mean_negloglik_change_from_model_zero =
      mean(negloglik_change_from_model_zero, na.rm = TRUE),
    mean_rms_log_observed_to_expected =
      mean(rms_log_observed_to_expected, na.rm = TRUE),
    mean_absolute_pearson_residual =
      mean(mean_absolute_pearson_residual, na.rm = TRUE),
    mean_rms_pearson_residual =
      mean(rms_pearson_residual, na.rm = TRUE),
    mean_maximum_absolute_pearson_residual =
      mean(maximum_absolute_pearson_residual, na.rm = TRUE),
    maximum_absolute_pearson_residual_across_imputations =
      max(maximum_absolute_pearson_residual, na.rm = TRUE),
    mean_observed_to_expected_LL =
      mean(observed_to_expected_LL, na.rm = TRUE),
    mean_observed_to_expected_LM =
      mean(observed_to_expected_LM, na.rm = TRUE),
    mean_observed_to_expected_MM =
      mean(observed_to_expected_MM, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(mixing_model, sensitivity_value, cell_spec)

# ------------------------------------------------------------------------------
# 14. Baseline validation: Delta = 0 should equal rho = 0
# ------------------------------------------------------------------------------

if (run_delta_models && run_rho_models) {
  baseline_delta_parameters <- parameter_estimates_pooled %>%
    filter(mixing_model == "sober_positive_delta", sensitivity_value == 0) %>%
    select(
      cell_spec, stage, parameter,
      estimate_delta0 = estimate
    )

  baseline_rho_parameters <- parameter_estimates_pooled %>%
    filter(mixing_model == "exact_group_rho", sensitivity_value == 0) %>%
    select(
      cell_spec, stage, parameter,
      estimate_rho0 = estimate
    )

  baseline_parameter_validation <- baseline_delta_parameters %>%
    inner_join(
      baseline_rho_parameters,
      by = c("cell_spec", "stage", "parameter")
    ) %>%
    mutate(
      comparison_type = "theta_lambda",
      cell_id = NA_character_,
      absolute_difference = abs(estimate_delta0 - estimate_rho0),
      relative_difference = safe_divide(
        absolute_difference,
        abs(estimate_delta0)
      ),
      validation_note = "Delta = 0 compared with rho = 0"
    ) %>%
    select(
      comparison_type, cell_spec, cell_id, stage, parameter,
      estimate_delta0, estimate_rho0,
      absolute_difference, relative_difference, validation_note
    )

  baseline_delta_r <- r_estimates_pooled %>%
    filter(mixing_model == "sober_positive_delta", sensitivity_value == 0) %>%
    select(
      cell_spec, cell_id, parameter,
      estimate_delta0 = estimate
    )

  baseline_rho_r <- r_estimates_pooled %>%
    filter(mixing_model == "exact_group_rho", sensitivity_value == 0) %>%
    select(
      cell_spec, cell_id, parameter,
      estimate_rho0 = estimate
    )

  baseline_r_validation <- baseline_delta_r %>%
    inner_join(
      baseline_rho_r,
      by = c("cell_spec", "cell_id", "parameter")
    ) %>%
    mutate(
      comparison_type = "r",
      stage = "two_car",
      absolute_difference = abs(estimate_delta0 - estimate_rho0),
      relative_difference = safe_divide(
        absolute_difference,
        abs(estimate_delta0)
      ),
      validation_note = "Delta = 0 compared with rho = 0"
    ) %>%
    select(
      comparison_type, cell_spec, cell_id, stage, parameter,
      estimate_delta0, estimate_rho0,
      absolute_difference, relative_difference, validation_note
    )

  baseline_validation <- bind_rows(
    baseline_parameter_validation,
    baseline_r_validation
  ) %>%
    arrange(comparison_type, cell_spec, cell_id, stage, parameter)

} else {
  baseline_validation <- tibble(
    comparison_type = "not_run",
    cell_spec = NA_character_,
    cell_id = NA_character_,
    stage = NA_character_,
    parameter = NA_character_,
    estimate_delta0 = NA_real_,
    estimate_rho0 = NA_real_,
    absolute_difference = NA_real_,
    relative_difference = NA_real_,
    validation_note = paste0(
      "Baseline comparison requires both models. Current run: ", run_label
    )
  )
}

# ------------------------------------------------------------------------------
# 15. Imputation-specific rho paths
# ------------------------------------------------------------------------------

rho_parameter_paths_by_imp <- parameter_estimates_by_imp %>%
  filter(mixing_model == "exact_group_rho") %>%
  arrange(cell_spec, stage, parameter, IMPUTATION, sensitivity_value)

rho_r_paths_by_imp <- r_estimates_by_imp %>%
  filter(mixing_model == "exact_group_rho") %>%
  arrange(cell_spec, cell_id, BAC_group, IMPUTATION, sensitivity_value)

rho_parameter_paths_pooled <- parameter_estimates_pooled %>%
  filter(mixing_model == "exact_group_rho") %>%
  arrange(cell_spec, stage, parameter, sensitivity_value)

rho_r_paths_pooled <- r_estimates_pooled %>%
  filter(mixing_model == "exact_group_rho") %>%
  arrange(cell_spec, cell_id, BAC_group, sensitivity_value)

if (nrow(rho_parameter_paths_by_imp) > 0) {
  rho_parameter_path_plot <- ggplot(
    rho_parameter_paths_by_imp,
    aes(
      x = sensitivity_value,
      y = estimate,
      group = factor(IMPUTATION)
    )
  ) +
    geom_line(alpha = 0.45, linewidth = 0.45) +
    geom_point(alpha = 0.65, size = 1.1) +
    geom_line(
      data = rho_parameter_paths_pooled,
      aes(x = sensitivity_value, y = estimate, group = 1),
      inherit.aes = FALSE,
      linewidth = 1.1
    ) +
    geom_point(
      data = rho_parameter_paths_pooled,
      aes(x = sensitivity_value, y = estimate),
      inherit.aes = FALSE,
      size = 2
    ) +
    facet_wrap(~ parameter, scales = "free_y", ncol = 3) +
    scale_y_log10() +
    labs(
      title = "Risk-parameter paths across rho values",
      subtitle = paste0(
        "Thin lines are individual imputations; thick line is the pooled estimate. ",
        "Cell specification: ", paste(cell_specs_to_estimate, collapse = ", ")
      ),
      x = "rho",
      y = "Estimate on log scale"
    ) +
    theme_bw() +
    theme(panel.grid.minor = element_blank())

  ggsave(
    file.path(fig_dir, "fig13_01_rho_parameter_paths_by_imputation.png"),
    rho_parameter_path_plot,
    width = 10,
    height = 7,
    dpi = 300
  )
}

if (nrow(rho_r_paths_by_imp) > 0) {
  rho_r_path_plot <- ggplot(
    rho_r_paths_by_imp,
    aes(
      x = sensitivity_value,
      y = estimate,
      group = factor(IMPUTATION)
    )
  ) +
    geom_line(alpha = 0.45, linewidth = 0.45) +
    geom_point(alpha = 0.65, size = 1.0) +
    geom_line(
      data = rho_r_paths_pooled,
      aes(x = sensitivity_value, y = estimate, group = 1),
      inherit.aes = FALSE,
      linewidth = 1.1
    ) +
    geom_point(
      data = rho_r_paths_pooled,
      aes(x = sensitivity_value, y = estimate),
      inherit.aes = FALSE,
      size = 1.8
    ) +
    facet_grid(BAC_group ~ cell_id, scales = "free_y") +
    scale_y_log10() +
    labs(
      title = "Cell-specific exposure-ratio paths across rho values",
      subtitle = "Thin lines are individual imputations; thick line is the pooled estimate",
      x = "rho",
      y = "Exposure ratio relative to S, on log scale"
    ) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 30, hjust = 1)
    )

  ggsave(
    file.path(fig_dir, "fig13_02_rho_r_paths_by_imputation.png"),
    rho_r_path_plot,
    width = 13,
    height = 8,
    dpi = 300
  )
}

# ------------------------------------------------------------------------------
# 16. Compare the global-best and baseline-connected branches
# ------------------------------------------------------------------------------

get_branch_fits <- function(res) {
  list(
    global_best = res$fit$global,
    baseline_connected = res$fit$connected
  )
}

branch_parameter_estimates_by_imp <- map_dfr(all_results, function(res) {
  imap_dfr(get_branch_fits(res), function(branch_fit, branch_name) {
    two_fit <- branch_fit$two_fit
    one_fit <- branch_fit$one_fit

    common <- tibble(
      branch = branch_name,
      IMPUTATION = res$IMPUTATION,
      cell_spec = res$cell_spec,
      mixing_model = res$mixing_model,
      sensitivity_name = res$sensitivity_name,
      sensitivity_value = res$sensitivity_value,
      scenario_id = res$scenario_id,
      scenario_label = res$scenario_label
    )

    bind_rows(
      common %>%
        slice(rep(1, 3)) %>%
        mutate(
          stage = "two_car",
          parameter = c("theta_L", "theta_M", "theta_H"),
          estimate = as.numeric(two_fit$theta[c("L", "M", "H")]),
          estimate_log = log(estimate),
          se_log = as.numeric(
            two_fit$se_log[c("log_theta_L", "log_theta_M", "log_theta_H")]
          ),
          se_delta = estimate * se_log,
          n_cells = two_fit$n_cells,
          n_crashes = two_fit$n2_total
        ),
      common %>%
        slice(rep(1, 3)) %>%
        mutate(
          stage = "one_car",
          parameter = c("lambda_L", "lambda_M", "lambda_H"),
          estimate = as.numeric(one_fit$lambda[c("L", "M", "H")]),
          estimate_log = log(estimate),
          se_log = as.numeric(
            one_fit$se_log[c("log_lambda_L", "log_lambda_M", "log_lambda_H")]
          ),
          se_delta = estimate * se_log,
          n_cells = one_fit$n_cells,
          n_crashes = one_fit$n1_total
        )
    )
  })
})

branch_parameter_estimates_pooled <- pool_log_estimates(
  branch_parameter_estimates_by_imp,
  group_vars = c(
    "branch", "mixing_model", "sensitivity_name", "sensitivity_value",
    "scenario_id", "scenario_label", "cell_spec", "stage", "parameter"
  )
)

branch_parameter_baseline <- branch_parameter_estimates_pooled %>%
  filter(
    mixing_model == equal_mixing_reference_model,
    sensitivity_value == 0
  ) %>%
  select(
    branch, cell_spec, stage, parameter,
    baseline_estimate = estimate
  )

branch_parameter_percent_changes <- branch_parameter_estimates_pooled %>%
  left_join(
    branch_parameter_baseline,
    by = c("branch", "cell_spec", "stage", "parameter")
  ) %>%
  mutate(
    percent_change_from_branch_zero =
      100 * (estimate / baseline_estimate - 1)
  ) %>%
  arrange(branch, mixing_model, sensitivity_value, stage, parameter)

branch_r_estimates_by_imp <- map_dfr(all_results, function(res) {
  imap_dfr(get_branch_fits(res), function(branch_fit, branch_name) {
    two_fit <- branch_fit$two_fit

    map_dfr(c("L", "M", "H"), function(group_i) {
      parameter_i <- paste0("r_", group_i)

      tibble(
        branch = branch_name,
        IMPUTATION = res$IMPUTATION,
        cell_spec = res$cell_spec,
        mixing_model = res$mixing_model,
        sensitivity_name = res$sensitivity_name,
        sensitivity_value = res$sensitivity_value,
        scenario_id = res$scenario_id,
        scenario_label = res$scenario_label,
        cell_id = two_fit$r_mat$cell_id,
        BAC_group = group_i,
        parameter = parameter_i,
        estimate = two_fit$r_mat[[parameter_i]],
        estimate_log = log(estimate),
        se_log = as.numeric(
          two_fit$se_log[paste0("log_", parameter_i, "__", two_fit$r_mat$cell_id)]
        ),
        se_delta = estimate * se_log,
        n_cells = two_fit$n_cells,
        n_crashes = two_fit$n2_total
      )
    })
  })
})

branch_r_estimates_pooled <- pool_log_estimates(
  branch_r_estimates_by_imp,
  group_vars = c(
    "branch", "mixing_model", "sensitivity_name", "sensitivity_value",
    "scenario_id", "scenario_label", "cell_spec", "cell_id",
    "BAC_group", "parameter"
  )
)

branch_objectives_by_imp <- map_dfr(all_results, function(res) {
  imap_dfr(get_branch_fits(res), function(branch_fit, branch_name) {
    tibble(
      branch = branch_name,
      IMPUTATION = res$IMPUTATION,
      cell_spec = res$cell_spec,
      mixing_model = res$mixing_model,
      sensitivity_name = res$sensitivity_name,
      sensitivity_value = res$sensitivity_value,
      scenario_id = res$scenario_id,
      scenario_label = res$scenario_label,
      two_car_negloglik = branch_fit$two_fit$fit$value,
      one_car_negloglik = branch_fit$one_fit$fit$value
    )
  })
})

branch_objective_baseline <- branch_objectives_by_imp %>%
  filter(sensitivity_value == 0) %>%
  select(
    branch, IMPUTATION, cell_spec, mixing_model,
    baseline_two_car_negloglik = two_car_negloglik,
    baseline_one_car_negloglik = one_car_negloglik
  )

branch_objectives_by_imp <- branch_objectives_by_imp %>%
  left_join(
    branch_objective_baseline,
    by = c("branch", "IMPUTATION", "cell_spec", "mixing_model")
  ) %>%
  mutate(
    two_car_negloglik_change_from_branch_zero =
      two_car_negloglik - baseline_two_car_negloglik,
    one_car_negloglik_change_from_branch_zero =
      one_car_negloglik - baseline_one_car_negloglik
  )

branch_objective_summary <- branch_objectives_by_imp %>%
  group_by(
    branch, mixing_model, sensitivity_name, sensitivity_value,
    scenario_id, scenario_label, cell_spec
  ) %>%
  summarize(
    n_imputations = n(),
    mean_two_car_negloglik = mean(two_car_negloglik),
    mean_two_car_negloglik_change_from_branch_zero =
      mean(two_car_negloglik_change_from_branch_zero),
    mean_one_car_negloglik = mean(one_car_negloglik),
    mean_one_car_negloglik_change_from_branch_zero =
      mean(one_car_negloglik_change_from_branch_zero),
    .groups = "drop"
  )

# Directly compare the two selected branches within every imputation/scenario.
branch_comparison_by_imp <- map_dfr(all_results, function(res) {
  global_two <- res$fit$global$two_fit
  connected_two <- res$fit$connected$two_fit
  global_one <- res$fit$global$one_fit
  connected_one <- res$fit$connected$one_fit

  two_distance <- parameter_distance(
    global_two$fit$par,
    connected_two$fit$par
  )
  two_max_abs <- maximum_absolute_parameter_difference(
    global_two$fit$par,
    connected_two$fit$par
  )

  one_distance <- parameter_distance(
    global_one$fit$par,
    connected_one$fit$par
  )
  one_max_abs <- maximum_absolute_parameter_difference(
    global_one$fit$par,
    connected_one$fit$par
  )

  tibble(
    IMPUTATION = res$IMPUTATION,
    cell_spec = res$cell_spec,
    mixing_model = res$mixing_model,
    sensitivity_name = res$sensitivity_name,
    sensitivity_value = res$sensitivity_value,
    scenario_id = res$scenario_id,
    scenario_label = res$scenario_label,
    global_two_car_negloglik = global_two$fit$value,
    connected_two_car_negloglik = connected_two$fit$value,
    connected_minus_global_two_car_negloglik =
      connected_two$fit$value - global_two$fit$value,
    two_car_parameter_distance = two_distance,
    two_car_max_abs_parameter_difference = two_max_abs,
    same_two_car_solution =
      is.finite(two_max_abs) &&
      two_max_abs <= same_solution_max_abs_difference,
    near_tied_distinct_two_car_modes =
      is.finite(two_distance) &&
      is.finite(connected_two$fit$value - global_two$fit$value) &&
      two_distance >= mode_distance_warning_threshold &&
      (connected_two$fit$value - global_two$fit$value) <=
        mode_objective_gap_warning_threshold,
    global_one_car_negloglik = global_one$fit$value,
    connected_one_car_negloglik = connected_one$fit$value,
    connected_minus_global_one_car_negloglik =
      connected_one$fit$value - global_one$fit$value,
    one_car_parameter_distance = one_distance,
    one_car_max_abs_parameter_difference = one_max_abs,
    same_one_car_solution =
      is.finite(one_max_abs) &&
      one_max_abs <= same_solution_max_abs_difference
  )
})

branch_warning_summary <- branch_comparison_by_imp %>%
  group_by(
    mixing_model, sensitivity_name, sensitivity_value,
    scenario_id, scenario_label, cell_spec
  ) %>%
  summarize(
    n_imputations = n(),
    n_imputations_same_two_car_solution = sum(same_two_car_solution),
    n_imputations_different_two_car_solution = sum(!same_two_car_solution),
    n_imputations_near_tied_distinct_modes =
      sum(near_tied_distinct_two_car_modes),
    mean_connected_minus_global_two_car_negloglik =
      mean(connected_minus_global_two_car_negloglik),
    maximum_connected_minus_global_two_car_negloglik =
      max(connected_minus_global_two_car_negloglik),
    mean_two_car_parameter_distance = mean(two_car_parameter_distance),
    maximum_two_car_parameter_distance = max(two_car_parameter_distance),
    .groups = "drop"
  ) %>%
  arrange(mixing_model, sensitivity_value, cell_spec)

# Save every successful optimizer solution's parameter vector. This makes it
# possible to inspect the actual modes instead of relying only on distances.
optimizer_parameter_solutions <- map_dfr(all_results, function(res) {
  common <- tibble(
    IMPUTATION = res$IMPUTATION,
    cell_spec = res$cell_spec,
    mixing_model = res$mixing_model,
    sensitivity_name = res$sensitivity_name,
    sensitivity_value = res$sensitivity_value,
    scenario_id = res$scenario_id,
    scenario_label = res$scenario_label
  )

  bind_rows(
    res$fit$two_optimizer$parameter_solutions %>%
      mutate(
        stage = "two_car",
        optimizer_context = "shared_two_car"
      ),
    res$fit$one_global_optimizer$parameter_solutions %>%
      mutate(
        stage = "one_car",
        optimizer_context = "conditional_on_global_two_car_r"
      ),
    res$fit$one_connected_optimizer$parameter_solutions %>%
      mutate(
        stage = "one_car",
        optimizer_context = "conditional_on_connected_two_car_r"
      )
  ) %>%
    crossing(common) %>%
    select(
      IMPUTATION, cell_spec, mixing_model, sensitivity_name,
      sensitivity_value, scenario_id, scenario_label,
      stage, optimizer_context, everything()
    )
})

# Compact rho table for model reporting. A later policy-simulation script can
# join its policy-effect estimate to this table by branch and rho.
rho_reporting_parameters <- branch_parameter_percent_changes %>%
  filter(
    mixing_model == "exact_group_rho",
    parameter %in% c("theta_M", "theta_H", "lambda_M", "lambda_H")
  ) %>%
  select(
    branch, sensitivity_value, cell_spec, parameter,
    estimate, ci_low, ci_high, percent_change_from_branch_zero
  ) %>%
  pivot_wider(
    names_from = parameter,
    values_from = c(
      estimate, ci_low, ci_high, percent_change_from_branch_zero
    ),
    names_glue = "{parameter}_{.value}"
  )

rho_reporting_effects <- rho_effect_summary %>%
  filter(BAC_group %in% c("L", "M", "H")) %>%
  select(
    rho, cell_spec, BAC_group,
    median_percent_diagonal_increase,
    percent_reduction_in_each_mixed_probability
  ) %>%
  pivot_wider(
    names_from = BAC_group,
    values_from = median_percent_diagonal_increase,
    names_glue = "median_percent_diagonal_increase_{BAC_group}"
  ) %>%
  distinct()

rho_sensitivity_reporting_table <- rho_reporting_parameters %>%
  left_join(
    branch_objective_summary %>%
      filter(mixing_model == "exact_group_rho") %>%
      select(
        branch, sensitivity_value, cell_spec,
        mean_two_car_negloglik,
        mean_two_car_negloglik_change_from_branch_zero
      ),
    by = c("branch", "sensitivity_value", "cell_spec")
  ) %>%
  left_join(
    branch_warning_summary %>%
      filter(mixing_model == "exact_group_rho") %>%
      select(
        sensitivity_value, cell_spec,
        n_imputations_different_two_car_solution,
        n_imputations_near_tied_distinct_modes,
        mean_connected_minus_global_two_car_negloglik,
        maximum_two_car_parameter_distance
      ),
    by = c("sensitivity_value", "cell_spec")
  ) %>%
  left_join(
    rho_reporting_effects,
    by = c("sensitivity_value" = "rho", "cell_spec")
  ) %>%
  mutate(
    percent_reduction_in_each_mixed_probability = if_else(
      sensitivity_value == 0,
      0,
      percent_reduction_in_each_mixed_probability
    ),
    median_percent_diagonal_increase_L = if_else(
      sensitivity_value == 0,
      0,
      median_percent_diagonal_increase_L
    ),
    median_percent_diagonal_increase_M = if_else(
      sensitivity_value == 0,
      0,
      median_percent_diagonal_increase_M
    ),
    median_percent_diagonal_increase_H = if_else(
      sensitivity_value == 0,
      0,
      median_percent_diagonal_increase_H
    )
  ) %>%
  arrange(branch, sensitivity_value, cell_spec)

# Branch-specific path plots. The original figures above still show the
# unrestricted global-best branch for backward compatibility.
rho_branch_parameters_by_imp <- branch_parameter_estimates_by_imp %>%
  filter(mixing_model == "exact_group_rho") %>%
  mutate(
    branch_label = recode(
      branch,
      global_best = "Global-best solution",
      baseline_connected = "Baseline-connected solution"
    ),
    panel_label = paste(branch_label, parameter, sep = ": ")
  )

rho_branch_parameters_pooled <- branch_parameter_estimates_pooled %>%
  filter(mixing_model == "exact_group_rho") %>%
  mutate(
    branch_label = recode(
      branch,
      global_best = "Global-best solution",
      baseline_connected = "Baseline-connected solution"
    ),
    panel_label = paste(branch_label, parameter, sep = ": ")
  )

if (nrow(rho_branch_parameters_by_imp) > 0) {
  rho_branch_parameter_plot <- ggplot(
    rho_branch_parameters_by_imp,
    aes(
      x = sensitivity_value,
      y = estimate,
      group = factor(IMPUTATION)
    )
  ) +
    geom_line(alpha = 0.40, linewidth = 0.40) +
    geom_point(alpha = 0.60, size = 1.0) +
    geom_line(
      data = rho_branch_parameters_pooled,
      aes(x = sensitivity_value, y = estimate, group = 1),
      inherit.aes = FALSE,
      linewidth = 1.1
    ) +
    geom_point(
      data = rho_branch_parameters_pooled,
      aes(x = sensitivity_value, y = estimate),
      inherit.aes = FALSE,
      size = 1.8
    ) +
    facet_wrap(~ panel_label, scales = "free_y", ncol = 3) +
    scale_y_log10() +
    labs(
      title = "Global-best and baseline-connected risk-parameter paths",
      subtitle = paste0(
        "Thin lines are imputations; thick lines are pooled within each branch. ",
        "Cell specification: ", paste(cell_specs_to_estimate, collapse = ", ")
      ),
      x = "rho",
      y = "Estimate on log scale"
    ) +
    theme_bw() +
    theme(panel.grid.minor = element_blank())

  ggsave(
    file.path(fig_dir, "fig13_03_rho_parameter_paths_global_vs_connected.png"),
    rho_branch_parameter_plot,
    width = 11,
    height = 12,
    dpi = 300
  )
}

rho_branch_r_by_imp <- branch_r_estimates_by_imp %>%
  filter(mixing_model == "exact_group_rho") %>%
  mutate(
    branch_label = recode(
      branch,
      global_best = "Global-best",
      baseline_connected = "Baseline-connected"
    ),
    branch_group = paste(branch_label, BAC_group, sep = ": ")
  )

rho_branch_r_pooled <- branch_r_estimates_pooled %>%
  filter(mixing_model == "exact_group_rho") %>%
  mutate(
    branch_label = recode(
      branch,
      global_best = "Global-best",
      baseline_connected = "Baseline-connected"
    ),
    branch_group = paste(branch_label, BAC_group, sep = ": ")
  )

if (nrow(rho_branch_r_by_imp) > 0) {
  rho_branch_r_plot <- ggplot(
    rho_branch_r_by_imp,
    aes(
      x = sensitivity_value,
      y = estimate,
      group = factor(IMPUTATION)
    )
  ) +
    geom_line(alpha = 0.40, linewidth = 0.40) +
    geom_point(alpha = 0.60, size = 0.9) +
    geom_line(
      data = rho_branch_r_pooled,
      aes(x = sensitivity_value, y = estimate, group = 1),
      inherit.aes = FALSE,
      linewidth = 1.0
    ) +
    geom_point(
      data = rho_branch_r_pooled,
      aes(x = sensitivity_value, y = estimate),
      inherit.aes = FALSE,
      size = 1.6
    ) +
    facet_grid(branch_group ~ cell_id, scales = "free_y") +
    scale_y_log10() +
    labs(
      title = "Global-best and baseline-connected exposure-ratio paths",
      subtitle = "Thin lines are imputations; thick lines are pooled within each branch",
      x = "rho",
      y = "Exposure ratio relative to S, on log scale"
    ) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 30, hjust = 1)
    )

  ggsave(
    file.path(fig_dir, "fig13_04_rho_r_paths_global_vs_connected.png"),
    rho_branch_r_plot,
    width = 13,
    height = 13,
    dpi = 300
  )
}


# ------------------------------------------------------------------------------
# 17. Save outputs
# ------------------------------------------------------------------------------

write_csv(
  parameter_estimates_by_imp,
  file.path(out_dir, "01_parameter_estimates_by_imp.csv")
)

write_csv(
  parameter_estimates_pooled,
  file.path(out_dir, "02_parameter_estimates_pooled.csv")
)

write_csv(
  parameter_percent_changes,
  file.path(out_dir, "03_parameter_percent_changes.csv")
)

write_csv(
  r_estimates_by_imp,
  file.path(out_dir, "04_r_estimates_by_imp.csv")
)

write_csv(
  r_estimates_pooled,
  file.path(out_dir, "05_r_estimates_pooled.csv")
)

write_csv(
  r_percent_changes,
  file.path(out_dir, "06_r_percent_changes.csv")
)

write_csv(
  road_shares_by_imp,
  file.path(out_dir, "07_road_shares_by_imp.csv")
)

write_csv(
  road_share_changes,
  file.path(out_dir, "08_road_share_changes.csv")
)

write_csv(
  fit_diagnostics,
  file.path(out_dir, "09_fit_diagnostics.csv")
)

write_csv(
  delta_feasibility_diagnostics,
  file.path(out_dir, "10_delta_feasibility_diagnostics.csv")
)

write_csv(
  rho_count_sparsity,
  file.path(out_dir, "11_rho_count_sparsity.csv")
)

write_csv(
  rho_effect_by_cell_group,
  file.path(out_dir, "12_rho_effect_by_cell_group.csv")
)

write_csv(
  rho_effect_summary,
  file.path(out_dir, "13_rho_effect_summary.csv")
)

write_csv(
  rho_reference_for_target_increases,
  file.path(out_dir, "14_rho_reference_for_target_diagonal_increases.csv")
)

write_csv(
  fitted_two_car_counts_by_type,
  file.path(out_dir, "15_two_car_fitted_counts_by_type.csv")
)

write_csv(
  fitted_two_car_summary_by_type,
  file.path(out_dir, "16_two_car_fit_summary_by_type.csv")
)

write_csv(
  baseline_validation,
  file.path(out_dir, "17_baseline_validation.csv")
)

write_csv(
  rho_parameter_paths_by_imp,
  file.path(out_dir, "18_rho_parameter_paths_by_imputation.csv")
)

write_csv(
  rho_r_paths_by_imp,
  file.path(out_dir, "19_rho_r_paths_by_imputation.csv")
)

write_csv(
  optimizer_start_diagnostics,
  file.path(out_dir, "20_optimizer_start_diagnostics.csv")
)

write_csv(
  two_car_fit_by_type_by_imputation,
  file.path(out_dir, "21_two_car_fit_by_type_by_imputation.csv")
)

write_csv(
  two_car_fit_summary_by_imputation,
  file.path(out_dir, "22_two_car_fit_summary_by_imputation.csv")
)

write_csv(
  scenario_level_fit_summary,
  file.path(out_dir, "23_scenario_level_fit_summary.csv")
)

write_csv(
  optimizer_parameter_solutions,
  file.path(out_dir, "24_optimizer_parameter_solutions.csv")
)

write_csv(
  branch_parameter_estimates_by_imp,
  file.path(out_dir, "25_branch_parameter_estimates_by_imp.csv")
)

write_csv(
  branch_parameter_estimates_pooled,
  file.path(out_dir, "26_branch_parameter_estimates_pooled.csv")
)

write_csv(
  branch_parameter_percent_changes,
  file.path(out_dir, "27_branch_parameter_percent_changes.csv")
)

write_csv(
  branch_r_estimates_by_imp,
  file.path(out_dir, "28_branch_r_estimates_by_imp.csv")
)

write_csv(
  branch_r_estimates_pooled,
  file.path(out_dir, "29_branch_r_estimates_pooled.csv")
)

write_csv(
  branch_objectives_by_imp,
  file.path(out_dir, "30_branch_objectives_by_imp.csv")
)

write_csv(
  branch_objective_summary,
  file.path(out_dir, "31_branch_objective_summary.csv")
)

write_csv(
  branch_comparison_by_imp,
  file.path(out_dir, "32_global_vs_connected_by_imputation.csv")
)

write_csv(
  branch_warning_summary,
  file.path(out_dir, "33_branch_warning_summary.csv")
)

write_csv(
  rho_sensitivity_reporting_table,
  file.path(out_dir, "34_rho_sensitivity_reporting_table.csv")
)

saveRDS(
  list(
    counts_all = counts_all,
    scenario_grid = scenario_grid,
    all_results = all_results,
    parameter_estimates_by_imp = parameter_estimates_by_imp,
    parameter_estimates_pooled = parameter_estimates_pooled,
    parameter_percent_changes = parameter_percent_changes,
    r_estimates_by_imp = r_estimates_by_imp,
    r_estimates_pooled = r_estimates_pooled,
    r_percent_changes = r_percent_changes,
    road_shares_by_imp = road_shares_by_imp,
    road_share_changes = road_share_changes,
    fit_diagnostics = fit_diagnostics,
    optimizer_start_diagnostics = optimizer_start_diagnostics,
    delta_feasibility_diagnostics = delta_feasibility_diagnostics,
    rho_count_sparsity = rho_count_sparsity,
    rho_effect_by_cell_group = rho_effect_by_cell_group,
    rho_effect_summary = rho_effect_summary,
    rho_reference_for_target_increases = rho_reference_for_target_increases,
    fitted_two_car_counts_by_type = fitted_two_car_counts_by_type,
    fitted_two_car_summary_by_type = fitted_two_car_summary_by_type,
    two_car_fit_by_type_by_imputation = two_car_fit_by_type_by_imputation,
    two_car_fit_summary_by_imputation = two_car_fit_summary_by_imputation,
    scenario_level_fit_summary = scenario_level_fit_summary,
    optimizer_parameter_solutions = optimizer_parameter_solutions,
    branch_parameter_estimates_by_imp = branch_parameter_estimates_by_imp,
    branch_parameter_estimates_pooled = branch_parameter_estimates_pooled,
    branch_parameter_percent_changes = branch_parameter_percent_changes,
    branch_r_estimates_by_imp = branch_r_estimates_by_imp,
    branch_r_estimates_pooled = branch_r_estimates_pooled,
    branch_objectives_by_imp = branch_objectives_by_imp,
    branch_objective_summary = branch_objective_summary,
    branch_comparison_by_imp = branch_comparison_by_imp,
    branch_warning_summary = branch_warning_summary,
    rho_sensitivity_reporting_table = rho_sensitivity_reporting_table,
    baseline_validation = baseline_validation,
    rho_parameter_paths_by_imp = rho_parameter_paths_by_imp,
    rho_r_paths_by_imp = rho_r_paths_by_imp,
    options = list(
      run_delta_models = run_delta_models,
      run_rho_models = run_rho_models,
      run_label = run_label,
      cell_specs_to_estimate = cell_specs_to_estimate,
      min_n2 = min_n2,
      delta_values = delta_values,
      rho_values_to_fit = rho_values_to_fit,
      rho_diagnostic_grid = rho_diagnostic_grid,
      target_diagonal_increases = target_diagonal_increases,
      log_param_lower = log_param_lower,
      log_param_upper = log_param_upper,
      mode_distance_warning_threshold = mode_distance_warning_threshold,
      mode_objective_gap_warning_threshold = mode_objective_gap_warning_threshold,
      same_solution_max_abs_difference = same_solution_max_abs_difference
    )
  ),
  file.path(out_dir, "script13_branch_tracking_results_bundle.rds")
)

cat("\nScript 13 branch-tracking version finished.\n")
cat("  Run type:", run_label, "\n")
cat("  Output folder:", out_dir, "\n")
cat("  Figure folder:", fig_dir, "\n")



############################################
### code to make tables
############################################

library(tidyverse)

# ------------------------------------------------------------------------------
# 1. File locations
# ------------------------------------------------------------------------------

out_dir <- "data/processed/script13_assumption2_sensitivity"

parameter_file <- file.path(
  out_dir,
  "26_branch_parameter_estimates_pooled.csv"
)

r_file <- file.path(
  out_dir,
  "29_branch_r_estimates_pooled.csv"
)

parameters <- read_csv(parameter_file, show_col_types = FALSE)
r_estimates <- read_csv(r_file, show_col_types = FALSE)

# Use character keys for sensitivity values to avoid numerical matching issues.
parameters <- parameters %>%
  mutate(value_key = sprintf("%.4f", sensitivity_value))

r_estimates <- r_estimates %>%
  mutate(value_key = sprintf("%.4f", sensitivity_value))

# ------------------------------------------------------------------------------
# 2. Choose the scenarios to report
# ------------------------------------------------------------------------------

scenarios_to_report <- tribble(
  ~scenario_order, ~scenario, ~mixing_model,          ~sensitivity_value, ~branch,
  1, "Equal mixing", "sober_positive_delta",          0.000, "global_best",
  2, "S/positive clustering", "sober_positive_delta", 0.100, "global_best",
  3, "S/positive clustering, stronger",
  "sober_positive_delta",                           0.150, "global_best",
  4, "Exact-group clustering", "exact_group_rho",     0.005, "baseline_connected",
  5, "Exact-group clustering, stronger",
  "exact_group_rho",                                0.010, "baseline_connected"
) %>%
  mutate(value_key = sprintf("%.4f", sensitivity_value))

# To make a shorter table, delete Scenario 3 and/or Scenario 5 above.

# ------------------------------------------------------------------------------
# 3. Obtain the common equal-mixing baseline
# ------------------------------------------------------------------------------

baseline_parameters <- parameters %>%
  filter(
    branch == "global_best",
    mixing_model == "sober_positive_delta",
    value_key == "0.0000",
    parameter %in% c(
      "theta_L", "theta_M", "theta_H",
      "lambda_L", "lambda_M", "lambda_H"
    )
  ) %>%
  mutate(
    BAC_group = str_remove(parameter, "^(theta|lambda)_"),
    parameter_type = if_else(
      str_starts(parameter, "theta_"),
      "theta",
      "lambda"
    )
  ) %>%
  select(BAC_group, parameter_type, baseline_estimate = estimate) %>%
  pivot_wider(
    names_from = parameter_type,
    values_from = baseline_estimate,
    names_glue = "baseline_{parameter_type}"
  )

baseline_r <- r_estimates %>%
  filter(
    branch == "global_best",
    mixing_model == "sober_positive_delta",
    value_key == "0.0000"
  ) %>%
  select(
    BAC_group,
    cell_id,
    baseline_r = estimate
  )

# ------------------------------------------------------------------------------
# 4. Select theta and lambda estimates for the desired scenarios
# ------------------------------------------------------------------------------

selected_parameters <- scenarios_to_report %>%
  inner_join(
    parameters,
    by = c("mixing_model", "branch", "value_key")
  ) %>%
  filter(
    parameter %in% c(
      "theta_L", "theta_M", "theta_H",
      "lambda_L", "lambda_M", "lambda_H"
    )
  ) %>%
  mutate(
    BAC_group = str_remove(parameter, "^(theta|lambda)_"),
    parameter_type = if_else(
      str_starts(parameter, "theta_"),
      "theta",
      "lambda"
    )
  ) %>%
  select(
    scenario_order,
    scenario,
    mixing_model,
    sensitivity_value = sensitivity_value.x,
    branch,
    BAC_group,
    parameter_type,
    estimate,
    ci_low,
    ci_high
  ) %>%
  pivot_wider(
    names_from = parameter_type,
    values_from = c(estimate, ci_low, ci_high),
    names_glue = "{parameter_type}_{.value}"
  ) %>%
  left_join(
    baseline_parameters,
    by = "BAC_group"
  ) %>%
  mutate(
    theta_transition = paste0(
      sprintf("%.3f", baseline_theta),
      " $\\rightarrow$ ",
      sprintf("%.3f", theta_estimate)
    ),
    lambda_transition = paste0(
      sprintf("%.3f", baseline_lambda),
      " $\\rightarrow$ ",
      sprintf("%.3f", lambda_estimate)
    )
  )

# ------------------------------------------------------------------------------
# 5. Select the four cell-specific exposure ratios
# ------------------------------------------------------------------------------

selected_r <- scenarios_to_report %>%
  inner_join(
    r_estimates,
    by = c("mixing_model", "branch", "value_key")
  ) %>%
  left_join(
    baseline_r,
    by = c("BAC_group", "cell_id")
  ) %>%
  mutate(
    sensitivity_value = sensitivity_value.x,
    
    # Names corresponding to rural/urban x weekend.
    cell_column = case_when(
      cell_id == "rural_FALSE" ~ "r_rural_weekday",
      cell_id == "rural_TRUE"  ~ "r_rural_weekend",
      cell_id == "urban_FALSE" ~ "r_urban_weekday",
      cell_id == "urban_TRUE"  ~ "r_urban_weekend",
      TRUE ~ paste0("r_", make.names(cell_id))
    ),
    
    r_transition = paste0(
      sprintf("%.3f", baseline_r),
      " $\\rightarrow$ ",
      sprintf("%.3f", estimate)
    )
  )

# Numeric long-form file. This is useful for checking values or making
# alternative table layouts later.
table_numeric_long <- selected_r %>%
  select(
    scenario_order,
    scenario,
    mixing_model,
    sensitivity_value,
    branch,
    BAC_group,
    cell_id,
    baseline_r,
    scenario_r = estimate,
    r_ci_low = ci_low,
    r_ci_high = ci_high
  ) %>%
  left_join(
    selected_parameters %>%
      select(
        scenario_order,
        BAC_group,
        baseline_theta,
        theta_estimate,
        theta_ci_low,
        theta_ci_high,
        baseline_lambda,
        lambda_estimate,
        lambda_ci_low,
        lambda_ci_high
      ),
    by = c("scenario_order", "BAC_group")
  ) %>%
  arrange(
    scenario_order,
    factor(BAC_group, levels = c("L", "M", "H")),
    cell_id
  )

# Put the four exposure transitions into separate columns.
r_formatted_wide <- selected_r %>%
  select(
    scenario_order,
    scenario,
    mixing_model,
    sensitivity_value,
    branch,
    BAC_group,
    cell_column,
    r_transition
  ) %>%
  pivot_wider(
    names_from = cell_column,
    values_from = r_transition
  )

# ------------------------------------------------------------------------------
# 6. Final formatted table
# ------------------------------------------------------------------------------

table_formatted <- selected_parameters %>%
  select(
    scenario_order,
    scenario,
    mixing_model,
    sensitivity_value,
    branch,
    BAC_group,
    theta_transition,
    lambda_transition
  ) %>%
  left_join(
    r_formatted_wide,
    by = c(
      "scenario_order",
      "scenario",
      "mixing_model",
      "sensitivity_value",
      "branch",
      "BAC_group"
    )
  ) %>%
  mutate(
    BAC_group = factor(BAC_group, levels = c("L", "M", "H"))
  ) %>%
  arrange(scenario_order, BAC_group)

# ------------------------------------------------------------------------------
# 7. Save the table-ready outputs
# ------------------------------------------------------------------------------

write_csv(
  table_numeric_long,
  file.path(out_dir, "35_assumption2_table_numeric_long.csv")
)

write_csv(
  table_formatted,
  file.path(out_dir, "36_assumption2_table_formatted.csv")
)

print(table_formatted, n = Inf, width = Inf)