################################################################################
# Validate four-group bootstrap CIs using simulated fatal crashes.
# Run as 100 PARCC Slurm array tasks, each with 10 datasets and 200 bootstraps.
################################################################################

library(tidyverse)
# setwd("your_hpc_working_directory")

# Identify the Slurm batch and keep its random stream distinct.
task_id <- suppressWarnings(as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")))
if (is.na(task_id) || !task_id %in% 1:100) stop("Submit this script using a Slurm array (tasks 1-100).")
set.seed(123 + 10000 * task_id + 2)

out_dir <- file.path("results/simulation_validation_full", sprintf("batch_%03d", task_id))
output_files <- c("bootstrap_original_point_estimates.csv", "bootstrap_diagnostics.csv",
                  "bootstrap_confidence_intervals.csv")
if (any(file.exists(file.path(out_dir, output_files)))) stop("Batch outputs already exist: ", out_dir)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

two_types <- c("SS", "SL", "SM", "SH", "LL", "LM", "LH", "MM", "MH", "HH")
one_types <- c("S", "L", "M", "H")
parameters <- c("theta_L", "theta_M", "theta_H", "lambda_L", "lambda_M", "lambda_H")
theta_true <- c(S = 1, L = 1.2, M = 1.6, H = 3.5)
lambda_true <- c(S = 1, L = 1.3, M = 2.0, H = 4.2)
e_two <- e_one <- 1e-4
log_param_lower <- log(1e-4)
log_param_upper <- log(500)
optim_control <- list(maxit = 5000, factr = 1e7)

# Use fitted FARS exposure ratios for data generation and starting values.
risk_starts <- read_csv("results/four_group_bootstrap/four_group_original_imputation_estimates.csv", show_col_types = FALSE) %>%
  filter(cell_spec == "rur_urb_weekend") %>%
  group_by(parameter) %>%
  summarize(estimate = mean(estimate), .groups = "drop") %>%
  deframe()
theta_start <- c(S = 1, L = risk_starts[["theta_L"]], M = risk_starts[["theta_M"]], H = risk_starts[["theta_H"]])
lambda_start <- c(S = 1, L = risk_starts[["lambda_L"]], M = risk_starts[["lambda_M"]], H = risk_starts[["lambda_H"]])
r_true <- read_csv("results/four_group_bootstrap/four_group_original_cell_exposure_ratios.csv", show_col_types = FALSE) %>%
  filter(cell_spec == "rur_urb_weekend") %>%
  group_by(cell_id) %>%
  summarize(r_L = mean(r_L), r_M = mean(r_M), r_H = mean(r_H), .groups = "drop")
fars_counts <- read_csv("data/processed/script03_four_group_prep/four_group_crash_counts_by_cell_imp.csv", show_col_types = FALSE) %>%
  filter(cell_spec == "rur_urb_weekend") %>%
  group_by(cell_id) %>%
  summarize(target_n2 = mean(n2), target_n1 = mean(n1), .groups = "drop")
simulation_inputs <- left_join(fars_counts, r_true, by = "cell_id")

# Simulate fatal crashes from encounters and crash opportunities in each cell.
simulate_cell <- function(cell, target_n2, target_n1, r_L, r_M, r_H) {
  r <- c(S = 1, L = r_L, M = r_M, H = r_H)
  road_share <- r / sum(r)
  encounter_prob <- c(
    SS = road_share["S"]^2, SL = 2 * road_share["S"] * road_share["L"],
    SM = 2 * road_share["S"] * road_share["M"], SH = 2 * road_share["S"] * road_share["H"],
    LL = road_share["L"]^2, LM = 2 * road_share["L"] * road_share["M"],
    LH = 2 * road_share["L"] * road_share["H"], MM = road_share["M"]^2,
    MH = 2 * road_share["M"] * road_share["H"], HH = road_share["H"]^2
  )
  two_fatal_prob <- e_two * c(
    SS = 2 * theta_true["S"], SL = theta_true["S"] + theta_true["L"],
    SM = theta_true["S"] + theta_true["M"], SH = theta_true["S"] + theta_true["H"],
    LL = 2 * theta_true["L"], LM = theta_true["L"] + theta_true["M"],
    LH = theta_true["L"] + theta_true["H"], MM = 2 * theta_true["M"],
    MH = theta_true["M"] + theta_true["H"], HH = 2 * theta_true["H"]
  )
  n_two_opportunities <- round(target_n2 / sum(encounter_prob * two_fatal_prob))
  two_opportunities <- as.vector(rmultinom(1, n_two_opportunities, encounter_prob))
  two_fatal_counts <- rbinom(length(two_types), two_opportunities, two_fatal_prob)
  names(two_fatal_counts) <- two_types
  
  one_fatal_prob <- e_one * lambda_true
  n_one_opportunities <- round(target_n1 / sum(road_share * one_fatal_prob))
  one_opportunities <- as.vector(rmultinom(1, n_one_opportunities, road_share))
  one_fatal_counts <- rbinom(length(one_types), one_opportunities, one_fatal_prob)
  names(one_fatal_counts) <- one_types
  
  counts <- tibble(cell_id = cell, !!!as.list(two_fatal_counts),
                   !!!set_names(as.list(one_fatal_counts), paste0("one_", one_types)),
                   n2 = sum(two_fatal_counts), n1 = sum(one_fatal_counts))
  crashes <- bind_rows(
    map_dfr(two_types, ~ tibble(cell_id = cell, VE_FORMS = 2L, crash_type_4 = rep(.x, two_fatal_counts[.x]))),
    map_dfr(one_types, ~ tibble(cell_id = cell, VE_FORMS = 1L, crash_type_4 = rep(.x, one_fatal_counts[.x])))
  )
  list(counts = counts, crashes = crashes)
}

# Reconstruct crash composition counts after each stratified bootstrap sample.
crashes_to_counts <- function(crashes) {
  two_counts <- crashes %>%
    filter(VE_FORMS == 2) %>% count(cell_id, crash_type_4) %>%
    complete(cell_id = simulation_inputs$cell_id, crash_type_4 = two_types, fill = list(n = 0)) %>%
    pivot_wider(names_from = crash_type_4, values_from = n, values_fill = 0)
  one_counts <- crashes %>%
    filter(VE_FORMS == 1) %>% count(cell_id, crash_type_4) %>%
    complete(cell_id = simulation_inputs$cell_id, crash_type_4 = one_types, fill = list(n = 0)) %>%
    pivot_wider(names_from = crash_type_4, values_from = n, values_fill = 0, names_prefix = "one_")
  simulation_inputs %>% select(cell_id) %>%
    left_join(two_counts, by = "cell_id") %>% left_join(one_counts, by = "cell_id") %>%
    mutate(n2 = rowSums(across(all_of(two_types))),
           n1 = rowSums(across(all_of(paste0("one_", one_types)))))
}

# Fit the two-vehicle likelihood, then the one-vehicle likelihood with fitted exposure.
safe_log <- function(x) log(pmax(x, 1e-300))
bad_likelihood_penalty <- 1e100

unpack_two_car_par <- function(par, cell_ids) {
  theta <- c(S = 1, L = unname(exp(par["log_theta_L"])),
             M = unname(exp(par["log_theta_M"])), H = unname(exp(par["log_theta_H"])))
  r_mat <- tibble(cell_id = cell_ids) %>%
    mutate(r_S = 1, r_L = exp(par[paste0("log_r_L__", cell_id)]),
           r_M = exp(par[paste0("log_r_M__", cell_id)]),
           r_H = exp(par[paste0("log_r_H__", cell_id)]))
  list(theta = theta, r_mat = r_mat)
}

make_two_car_probabilities <- function(theta, r_mat) {
  r_mat %>% transmute(
    cell_id, w_SS = theta["S"] * r_S^2,
    w_SL = (theta["S"] + theta["L"]) * r_S * r_L,
    w_SM = (theta["S"] + theta["M"]) * r_S * r_M,
    w_SH = (theta["S"] + theta["H"]) * r_S * r_H,
    w_LL = theta["L"] * r_L^2,
    w_LM = (theta["L"] + theta["M"]) * r_L * r_M,
    w_LH = (theta["L"] + theta["H"]) * r_L * r_H,
    w_MM = theta["M"] * r_M^2,
    w_MH = (theta["M"] + theta["H"]) * r_M * r_H,
    w_HH = theta["H"] * r_H^2
  ) %>%
    mutate(D = w_SS + w_SL + w_SM + w_SH + w_LL + w_LM + w_LH + w_MM + w_MH + w_HH) %>%
    transmute(cell_id, p_SS = w_SS / D, p_SL = w_SL / D, p_SM = w_SM / D,
              p_SH = w_SH / D, p_LL = w_LL / D, p_LM = w_LM / D,
              p_LH = w_LH / D, p_MM = w_MM / D, p_MH = w_MH / D, p_HH = w_HH / D)
}

negloglik_two_car <- function(par, dat, cell_ids) {
  unpacked <- unpack_two_car_par(par, cell_ids)
  probs <- make_two_car_probabilities(unpacked$theta, unpacked$r_mat)
  ll_dat <- left_join(dat, probs, by = "cell_id")
  if (any(!is.finite(as.matrix(ll_dat[, paste0("p_", two_types)])))) return(bad_likelihood_penalty)
  loglik_terms <- rep(0, nrow(ll_dat))
  for (type in two_types) loglik_terms <- loglik_terms + ll_dat[[type]] * safe_log(ll_dat[[paste0("p_", type)]])
  if (any(!is.finite(loglik_terms))) return(bad_likelihood_penalty)
  -sum(loglik_terms)
}

make_two_car_start <- function(cell_ids) {
  r_start <- r_true %>% filter(cell_id %in% cell_ids) %>% arrange(match(cell_id, cell_ids))
  par <- c(log_theta_L = log(unname(theta_start["L"])), log_theta_M = log(unname(theta_start["M"])),
           log_theta_H = log(unname(theta_start["H"])))
  for (i in seq_along(cell_ids)) {
    cell_i <- cell_ids[i]
    par[paste0("log_r_L__", cell_i)] <- log(r_start$r_L[i])
    par[paste0("log_r_M__", cell_i)] <- log(r_start$r_M[i])
    par[paste0("log_r_H__", cell_i)] <- log(r_start$r_H[i])
  }
  par
}

estimate_two_car_stage <- function(dat) {
  dat <- dat %>% filter(n2 >= 1) %>% arrange(cell_id)
  if (nrow(dat) == 0) stop("No two-vehicle crash cells available.")
  cell_ids <- dat$cell_id
  start_par <- make_two_car_start(cell_ids)
  fit <- optim(par = start_par, fn = negloglik_two_car, dat = dat, cell_ids = cell_ids,
               method = "L-BFGS-B", lower = rep(log_param_lower, length(start_par)),
               upper = rep(log_param_upper, length(start_par)), control = optim_control)
  unpacked <- unpack_two_car_par(fit$par, cell_ids)
  list(fit = fit, theta = unpacked$theta, r_mat = unpacked$r_mat)
}

make_one_car_probabilities <- function(lambda, r_mat) {
  r_mat %>% transmute(cell_id, w_S = lambda["S"] * r_S, w_L = lambda["L"] * r_L,
                      w_M = lambda["M"] * r_M, w_H = lambda["H"] * r_H) %>%
    mutate(D = w_S + w_L + w_M + w_H) %>%
    transmute(cell_id, q_S = w_S / D, q_L = w_L / D, q_M = w_M / D, q_H = w_H / D)
}

negloglik_one_car <- function(par, dat, r_mat) {
  lambda <- c(S = 1, L = unname(exp(par["log_lambda_L"])),
              M = unname(exp(par["log_lambda_M"])), H = unname(exp(par["log_lambda_H"])))
  ll_dat <- left_join(dat, make_one_car_probabilities(lambda, r_mat), by = "cell_id")
  if (any(!is.finite(as.matrix(ll_dat[, paste0("q_", one_types)])))) return(bad_likelihood_penalty)
  loglik_terms <- ll_dat$one_S * safe_log(ll_dat$q_S) + ll_dat$one_L * safe_log(ll_dat$q_L) +
    ll_dat$one_M * safe_log(ll_dat$q_M) + ll_dat$one_H * safe_log(ll_dat$q_H)
  if (any(!is.finite(loglik_terms))) return(bad_likelihood_penalty)
  -sum(loglik_terms)
}

estimate_one_car_stage <- function(dat, r_mat) {
  dat <- dat %>% filter(n1 > 0) %>% arrange(cell_id)
  r_mat <- r_mat %>% filter(cell_id %in% dat$cell_id) %>% arrange(cell_id)
  start_par <- c(log_lambda_L = log(unname(lambda_start["L"])), log_lambda_M = log(unname(lambda_start["M"])),
                 log_lambda_H = log(unname(lambda_start["H"])))
  fit <- optim(par = start_par, fn = negloglik_one_car, dat = dat, r_mat = r_mat,
               method = "L-BFGS-B", lower = rep(log_param_lower, length(start_par)),
               upper = rep(log_param_upper, length(start_par)), control = optim_control)
  lambda <- c(S = 1, L = unname(exp(fit$par["log_lambda_L"])),
              M = unname(exp(fit$par["log_lambda_M"])), H = unname(exp(fit$par["log_lambda_H"])))
  list(fit = fit, lambda = lambda)
}

fit_four_group <- function(count_data) {
  two_fit <- estimate_two_car_stage(count_data)
  one_fit <- estimate_one_car_stage(count_data, two_fit$r_mat)
  two_car_convergence <- two_fit$fit$convergence
  one_car_convergence <- one_fit$fit$convergence
  convergence <- ifelse(two_car_convergence == 0 && one_car_convergence == 0, 0L, 1L)
  risk_estimates <- tibble(
    parameter = parameters,
    true_value = c(theta_true[c("L", "M", "H")], lambda_true[c("L", "M", "H")]),
    estimate = c(two_fit$theta[c("L", "M", "H")], one_fit$lambda[c("L", "M", "H")])
  )
  list(risk_estimates = risk_estimates, convergence = convergence,
       two_car_convergence = two_car_convergence, one_car_convergence = one_car_convergence)
}

# Generate 10 datasets and refit the two-stage model for 200 bootstraps each.
n_simulations <- 10
n_bootstrap <- 200
original <- vector("list", n_simulations)
boot <- vector("list", n_simulations * n_bootstrap)
diagnostics <- vector("list", n_simulations * n_bootstrap)

for (i in seq_len(n_simulations)) {
  simulated <- pmap(simulation_inputs, function(cell_id, target_n2, target_n1, r_L, r_M, r_H) {
    simulate_cell(cell_id, target_n2, target_n1, r_L, r_M, r_H)
  })
  counts_i <- map_dfr(simulated, "counts")
  crashes_i <- map_dfr(simulated, "crashes")
  fit_i <- tryCatch(fit_four_group(counts_i), error = function(e) e)
  
  if (inherits(fit_i, "error")) {
    original[[i]] <- tibble(simulation = i, parameter = NA_character_, true_value = NA_real_,
                            estimate = NA_real_, convergence = NA_integer_)
  } else {
    original[[i]] <- fit_i$risk_estimates %>%
      transmute(simulation = i, parameter, true_value, estimate, convergence = fit_i$convergence)
  }
  
  for (j in seq_len(n_bootstrap)) {
    k <- (i - 1L) * n_bootstrap + j
    crashes_boot <- crashes_i %>% group_by(cell_id, VE_FORMS) %>%
      slice_sample(prop = 1, replace = TRUE) %>% ungroup()
    fit_boot <- tryCatch(fit_four_group(crashes_to_counts(crashes_boot)), error = function(e) e)
    
    if (inherits(fit_boot, "error")) {
      diagnostics[[k]] <- tibble(simulation = i, bootstrap = j, convergence = NA_integer_,
                                 two_car_convergence = NA_integer_, one_car_convergence = NA_integer_,
                                 message = conditionMessage(fit_boot))
    } else {
      boot[[k]] <- fit_boot$risk_estimates %>%
        transmute(simulation = i, bootstrap = j, parameter, true_value,
                  estimate, convergence = fit_boot$convergence)
      diagnostics[[k]] <- tibble(simulation = i, bootstrap = j, convergence = fit_boot$convergence,
                                 two_car_convergence = fit_boot$two_car_convergence,
                                 one_car_convergence = fit_boot$one_car_convergence, message = "")
    }
  }
}

original <- bind_rows(original)
boot <- bind_rows(boot)
diagnostics <- bind_rows(diagnostics)

# Form log-scale Wald intervals using the successful bootstrap fits.
intervals <- boot %>% filter(convergence == 0) %>%
  mutate(log_estimate = log(estimate)) %>%
  group_by(simulation, parameter) %>%
  summarize(true_value = first(true_value), n_bootstrap_successful = n(),
            bootstrap_se_log = sd(log_estimate), .groups = "drop") %>%
  left_join(original %>% filter(convergence == 0) %>%
              select(simulation, parameter, original_estimate = estimate),
            by = c("simulation", "parameter")) %>%
  mutate(lower = exp(log(original_estimate) - 1.96 * bootstrap_se_log),
         upper = exp(log(original_estimate) + 1.96 * bootstrap_se_log),
         covers_truth = lower <= true_value & upper >= true_value)

# Script 09b reads these three files from every batch directory.
write_csv(original, file.path(out_dir, output_files[1]))
write_csv(diagnostics, file.path(out_dir, output_files[2]))
write_csv(intervals, file.path(out_dir, output_files[3]))
