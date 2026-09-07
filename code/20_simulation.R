################################################################################
# Choose one set of true four-group parameters
# Simulate fake datasets from crash opportunities, and compare their crash counts with FARS
# Then fit the paper's model to fake datasets + bootstrap samples, and compare with true parameters
################################################################################

rm(list = ls())

library(tidyverse)

# Settings
cell_counts_path <- "data/processed/script10_four_group_prep/four_group_crash_counts_by_cell_imp.csv"
estimate_path <- "results/four_group_bootstrap/four_group_original_imputation_estimates.csv"
exposure_path <- "results/four_group_bootstrap/four_group_original_cell_exposure_ratios.csv"
out_dir <- "results/simulation_validation"
simulation_seed <- 123
e_two <- 1e-4
e_one <- 1e-4

# Model settings matching script 16.
min_n2 <- 1
log_param_lower <- log(1e-4)
log_param_upper <- log(500)
optim_control <- list(maxit = 5000, factr = 1e7)

# True relative-risk parameters.
theta_true <- c(S = 1, L = 1.2, M = 1.6, H = 3.5)
lambda_true <- c(S = 1, L = 1.3, M = 2.0, H = 4.2)

groups <- c("S", "L", "M", "H")
two_types <- c("SS", "SL", "SM", "SH", "LL", "LM", "LH", "MM", "MH", "HH")
one_types <- groups

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(simulation_seed)

# Use the mean fitted estimates from script 16 as optimization starting values.
risk_starts <- read_csv(estimate_path, show_col_types = FALSE) %>%
  filter(cell_spec == "rur_urb_weekend") %>%
  group_by(parameter) %>%
  summarize(estimate = mean(estimate), .groups = "drop") %>%
  deframe()

theta_start <- c(
  S = 1,
  L = unname(risk_starts["theta_L"]),
  M = unname(risk_starts["theta_M"]),
  H = unname(risk_starts["theta_H"])
)

lambda_start <- c(
  S = 1,
  L = unname(risk_starts["lambda_L"]),
  M = unname(risk_starts["lambda_M"]),
  H = unname(risk_starts["lambda_H"])
)

# Use the mean fitted exposure ratios across imputations as both the true
# exposure ratios and the optimization starting values.
r_true <- read_csv(exposure_path, show_col_types = FALSE) %>%
  filter(cell_spec == "rur_urb_weekend") %>%
  group_by(cell_id) %>%
  summarize(r_L = mean(r_L), r_M = mean(r_M), r_H = mean(r_H), .groups = "drop")

# Use FARS average one- and two-vehicle total counts by cell as target expected crash counts.
fars_counts <- read_csv(cell_counts_path, show_col_types = FALSE) %>%
  filter(cell_spec == "rur_urb_weekend") %>%
  group_by(cell_id) %>%
  summarize(
    across(all_of(two_types), mean),
    across(all_of(paste0("one_", one_types)), mean),
    target_n2 = mean(n2),
    target_n1 = mean(n1),
    .groups = "drop"
  )

true_parameters <- r_true %>%
  mutate(
    theta_L = theta_true["L"], theta_M = theta_true["M"], theta_H = theta_true["H"],
    lambda_L = lambda_true["L"], lambda_M = lambda_true["M"], lambda_H = lambda_true["H"]
  )

##############################################
# Generate one fake dataset
##############################################

simulate_cell <- function(cell, target_n2, target_n1, r_L, r_M, r_H) {
  r <- c(S = 1, L = r_L, M = r_M, H = r_H)
  road_share <- r / sum(r)
  
  encounter_prob <- c(
    SS = road_share["S"]^2,
    SL = 2 * road_share["S"] * road_share["L"],
    SM = 2 * road_share["S"] * road_share["M"],
    SH = 2 * road_share["S"] * road_share["H"],
    LL = road_share["L"]^2,
    LM = 2 * road_share["L"] * road_share["M"],
    LH = 2 * road_share["L"] * road_share["H"],
    MM = road_share["M"]^2,
    MH = 2 * road_share["M"] * road_share["H"],
    HH = road_share["H"]^2
  )
  
  two_fatal_prob <- e_two * c(
    SS = 2 * theta_true["S"],
    SL = theta_true["S"] + theta_true["L"],
    SM = theta_true["S"] + theta_true["M"],
    SH = theta_true["S"] + theta_true["H"],
    LL = 2 * theta_true["L"],
    LM = theta_true["L"] + theta_true["M"],
    LH = theta_true["L"] + theta_true["H"],
    MM = 2 * theta_true["M"],
    MH = theta_true["M"] + theta_true["H"],
    HH = 2 * theta_true["H"]
  )
  
  expected_two_fatal_prob <- sum(encounter_prob * two_fatal_prob)
  n_two_opportunities <- round(target_n2 / expected_two_fatal_prob)
  two_opportunities <- as.vector(rmultinom(1, n_two_opportunities, encounter_prob))
  names(two_opportunities) <- two_types
  two_fatal_counts <- rbinom(length(two_types), two_opportunities, two_fatal_prob)
  names(two_fatal_counts) <- two_types
  
  one_fatal_prob <- e_one * lambda_true
  expected_one_fatal_prob <- sum(road_share * one_fatal_prob)
  n_one_opportunities <- round(target_n1 / expected_one_fatal_prob)
  one_opportunities <- as.vector(rmultinom(1, n_one_opportunities, road_share))
  names(one_opportunities) <- one_types
  one_fatal_counts <- rbinom(length(one_types), one_opportunities, one_fatal_prob)
  names(one_fatal_counts) <- one_types
  
  counts <- tibble(
    cell_id = cell,
    !!!as.list(two_fatal_counts),
    !!!set_names(as.list(one_fatal_counts), paste0("one_", one_types)),
    n2 = sum(two_fatal_counts),
    n1 = sum(one_fatal_counts),
    n_two_opportunities = n_two_opportunities,
    n_one_opportunities = n_one_opportunities
  )
  
  crash_level <- bind_rows(
    map_dfr(two_types, ~ tibble(cell_id = cell, VE_FORMS = 2L, crash_type_4 = rep(.x, two_fatal_counts[.x]))),
    map_dfr(one_types, ~ tibble(cell_id = cell, VE_FORMS = 1L, crash_type_4 = rep(.x, one_fatal_counts[.x])))
  )
  
  list(counts = counts, crash_level = crash_level)
}

simulation_inputs <- fars_counts %>% left_join(r_true, by = "cell_id")

simulated <- pmap(simulation_inputs, function(cell_id, ..., target_n2, target_n1, r_L, r_M, r_H) {
  simulate_cell(cell_id, target_n2, target_n1, r_L, r_M, r_H)
})

fake_counts <- map_dfr(simulated, "counts")
fake_crashes <- map_dfr(simulated, "crash_level") %>%
  mutate(CRASH_ID = row_number(), .before = 1)

count_check <- fars_counts %>%
  select(cell_id, all_of(two_types), all_of(paste0("one_", one_types)), target_n2, target_n1) %>%
  pivot_longer(-c(cell_id, target_n2, target_n1), names_to = "type", values_to = "fars_mean_count") %>%
  left_join(
    fake_counts %>%
      select(cell_id, all_of(two_types), all_of(paste0("one_", one_types))) %>%
      pivot_longer(-cell_id, names_to = "type", values_to = "simulated_count"),
    by = c("cell_id", "type")
  ) %>%
  mutate(
    crash_size = if_else(str_starts(type, "one_"), "one_vehicle", "two_vehicle"),
    crash_type = str_remove(type, "^one_"),
    difference = simulated_count - fars_mean_count,
    percent_difference = 100 * difference / fars_mean_count
  ) %>%
  select(cell_id, crash_size, crash_type, fars_mean_count, simulated_count, difference, percent_difference)

write_csv(true_parameters, file.path(out_dir, "simulation_true_parameters.csv"))
write_csv(fake_counts, file.path(out_dir, "simulation_fake_counts_by_cell.csv"))
write_csv(fake_crashes, file.path(out_dir, "simulation_fake_crash_level.csv"))
write_csv(count_check, file.path(out_dir, "simulation_count_check.csv"))

cat("\nSaved outputs to:", out_dir, "\n")

##############################################
# Fit four-group model to the fake dataset
##############################################

safe_log <- function(x) log(pmax(x, 1e-300))
bad_likelihood_penalty <- 1e100

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

make_two_car_probabilities <- function(theta, r_mat) {
  r_mat %>%
    transmute(
      cell_id,
      w_SS = theta["S"] * r_S^2,
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
    transmute(
      cell_id,
      p_SS = w_SS / D,
      p_SL = w_SL / D,
      p_SM = w_SM / D,
      p_SH = w_SH / D,
      p_LL = w_LL / D,
      p_LM = w_LM / D,
      p_LH = w_LH / D,
      p_MM = w_MM / D,
      p_MH = w_MH / D,
      p_HH = w_HH / D
    )
}

negloglik_two_car <- function(par, dat, cell_ids) {
  unpacked <- unpack_two_car_par(par, cell_ids)
  probs <- make_two_car_probabilities(unpacked$theta, unpacked$r_mat)
  ll_dat <- dat %>% left_join(probs, by = "cell_id")
  
  prob_cols <- paste0("p_", two_types)
  if (any(!is.finite(as.matrix(ll_dat[, prob_cols])))) return(bad_likelihood_penalty)
  
  loglik_terms <- rep(0, nrow(ll_dat))
  for (type in two_types) {
    loglik_terms <- loglik_terms + ll_dat[[type]] * safe_log(ll_dat[[paste0("p_", type)]])
  }
  
  if (any(!is.finite(loglik_terms))) return(bad_likelihood_penalty)
  -sum(loglik_terms)
}

# Start the two-car model at the estimates from script 16.
make_two_car_start <- function(cell_ids) {
  r_start <- r_true %>%
    filter(cell_id %in% cell_ids) %>%
    arrange(match(cell_id, cell_ids))
  
  par <- c(
    log_theta_L = log(unname(theta_start["L"])),
    log_theta_M = log(unname(theta_start["M"])),
    log_theta_H = log(unname(theta_start["H"]))
  )
  
  for (i in seq_along(cell_ids)) {
    cell_i <- cell_ids[i]
    par[paste0("log_r_L__", cell_i)] <- log(r_start$r_L[i])
    par[paste0("log_r_M__", cell_i)] <- log(r_start$r_M[i])
    par[paste0("log_r_H__", cell_i)] <- log(r_start$r_H[i])
  }
  
  par
}

estimate_two_car_stage <- function(dat) {
  dat <- dat %>% filter(n2 >= min_n2) %>% arrange(cell_id)
  if (nrow(dat) == 0) stop("No two-car cells remain after applying min_n2.")
  
  cell_ids <- dat$cell_id
  start_par <- make_two_car_start(cell_ids)
  
  fit <- optim(
    par = start_par,
    fn = negloglik_two_car,
    dat = dat,
    cell_ids = cell_ids,
    method = "L-BFGS-B",
    lower = rep(log_param_lower, length(start_par)),
    upper = rep(log_param_upper, length(start_par)),
    control = optim_control
  )
  
  unpacked <- unpack_two_car_par(fit$par, cell_ids)
  
  list(
    fit = fit,
    theta = unpacked$theta,
    r_mat = unpacked$r_mat,
    n_cells = nrow(dat),
    n2_total = sum(dat$n2)
  )
}

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
  if (any(!is.finite(as.matrix(ll_dat[, prob_cols])))) return(bad_likelihood_penalty)
  
  loglik_terms <-
    ll_dat$one_S * safe_log(ll_dat$q_S) +
    ll_dat$one_L * safe_log(ll_dat$q_L) +
    ll_dat$one_M * safe_log(ll_dat$q_M) +
    ll_dat$one_H * safe_log(ll_dat$q_H)
  
  if (any(!is.finite(loglik_terms))) return(bad_likelihood_penalty)
  -sum(loglik_terms)
}

make_one_car_start <- function() {
  c(
    log_lambda_L = log(unname(lambda_start["L"])),
    log_lambda_M = log(unname(lambda_start["M"])),
    log_lambda_H = log(unname(lambda_start["H"]))
  )
}

estimate_one_car_stage <- function(dat, r_mat) {
  dat <- dat %>% filter(n1 > 0) %>% arrange(cell_id)
  r_mat <- r_mat %>% filter(cell_id %in% dat$cell_id) %>% arrange(cell_id)
  start_par <- make_one_car_start()
  
  fit <- optim(
    par = start_par,
    fn = negloglik_one_car,
    dat = dat,
    r_mat = r_mat,
    method = "L-BFGS-B",
    lower = rep(log_param_lower, length(start_par)),
    upper = rep(log_param_upper, length(start_par)),
    control = optim_control
  )
  
  lambda <- c(
    S = 1,
    L = unname(exp(fit$par["log_lambda_L"])),
    M = unname(exp(fit$par["log_lambda_M"])),
    H = unname(exp(fit$par["log_lambda_H"]))
  )
  
  list(
    fit = fit,
    lambda = lambda,
    n_cells = nrow(dat),
    n1_total = sum(dat$n1)
  )
}

# Fit the two-car stage first, then use its exposure estimates in the one-car stage.
fit_four_group <- function(count_data) {
  two_fit <- estimate_two_car_stage(count_data)
  one_fit <- estimate_one_car_stage(count_data, two_fit$r_mat)
  
  two_convergence <- two_fit$fit$convergence
  one_convergence <- one_fit$fit$convergence
  convergence <- ifelse(two_convergence == 0 && one_convergence == 0, 0L, 1L)
  
  risk_estimates <- tibble(
    parameter = c("theta_L", "theta_M", "theta_H", "lambda_L", "lambda_M", "lambda_H"),
    true_value = c(theta_true[c("L", "M", "H")], lambda_true[c("L", "M", "H")]),
    estimate = c(
      two_fit$theta[c("L", "M", "H")],
      one_fit$lambda[c("L", "M", "H")]
    )
  ) %>%
    mutate(
      error = estimate - true_value,
      percent_error = 100 * error / true_value
    )
  
  exposure_estimates <- two_fit$r_mat %>%
    select(cell_id, r_L, r_M, r_H)
  
  list(
    risk_estimates = risk_estimates,
    exposure_estimates = exposure_estimates,
    two_car_fit = two_fit,
    one_car_fit = one_fit,
    two_car_convergence = two_convergence,
    one_car_convergence = one_convergence,
    convergence = convergence
  )
}

fake_fit <- fit_four_group(fake_counts)

cat("\nTwo-car convergence code:", fake_fit$two_car_convergence, "\n")
cat("One-car convergence code:", fake_fit$one_car_convergence, "\n")

cat("\nTrue and estimated risk parameters:\n")
print(fake_fit$risk_estimates)

cat("\nEstimated exposure ratios by cell:\n")
print(fake_fit$exposure_estimates)

##############################################
# Simulate many more datasets, fit the model to each, and summarize the results
##############################################

simulation_reps <- 100
set.seed(simulation_seed + 1)

simulation_estimates <- vector("list", simulation_reps)
simulation_diagnostics <- vector("list", simulation_reps)

for (b in seq_len(simulation_reps)) {
  
  # Generate a fresh dataset using the same true parameters and target sample sizes.
  simulated_b <- pmap(simulation_inputs, function(cell_id, ..., target_n2, target_n1, r_L, r_M, r_H) {
    simulate_cell(cell_id, target_n2, target_n1, r_L, r_M, r_H)
  })
  counts_b <- map_dfr(simulated_b, "counts")
  
  # Record fitting errors so one failed fit does not stop the simulation.
  fit_b <- tryCatch(fit_four_group(counts_b), error = function(e) e)
  
  if (inherits(fit_b, "error")) {
    simulation_diagnostics[[b]] <- tibble(
      simulation = b,
      convergence = NA_integer_,
      two_car_convergence = NA_integer_,
      one_car_convergence = NA_integer_,
      message = conditionMessage(fit_b)
    )
  } else {
    simulation_estimates[[b]] <- fit_b$risk_estimates %>%
      mutate(simulation = b, convergence = fit_b$convergence, .before = 1)
    
    simulation_diagnostics[[b]] <- tibble(
      simulation = b,
      convergence = fit_b$convergence,
      two_car_convergence = fit_b$two_car_convergence,
      one_car_convergence = fit_b$one_car_convergence,
      message = ""
    )
  }
  
  if (b %% 10 == 0 || b == simulation_reps) {
    cat("Completed", b, "of", simulation_reps, "simulated datasets\n")
  }
}

simulation_estimates <- bind_rows(simulation_estimates)
simulation_diagnostics <- bind_rows(simulation_diagnostics)

write_csv(simulation_estimates, file.path(out_dir, "simulation_point_estimates.csv"))
write_csv(simulation_diagnostics, file.path(out_dir, "simulation_fit_diagnostics.csv"))

# Summarize converged fits and record how many estimates approach the optimizer bounds.
n_converged <- sum(simulation_diagnostics$convergence == 0, na.rm = TRUE)
n_converged

simulation_summary <- simulation_estimates %>%
  filter(convergence == 0) %>%
  group_by(parameter) %>%
  summarize(
    true_value = first(true_value),
    n = n(),
    mean_estimate = mean(estimate),
    median_estimate = median(estimate),
    bias = mean(error),
    percent_bias = 100 * mean(error) / first(true_value),
    empirical_sd = sd(estimate),
    bias_mcse = sd(estimate) / sqrt(n()),
    rmse = sqrt(mean(error^2)),
    fraction_at_lower_bound = mean(estimate <= exp(log_param_lower) * (1 + 1e-6)),
    fraction_at_upper_bound = mean(estimate >= exp(log_param_upper) * (1 - 1e-6)),
    .groups = "drop"
  )

print(simulation_summary, width = Inf)
write_csv(simulation_summary, file.path(out_dir, "simulation_point_estimate_summary.csv"))

# Show the distribution of estimates relative to the imposed truth.
estimate_plot <- simulation_estimates %>%
  filter(convergence == 0) %>%
  ggplot(aes(x = estimate)) +
  geom_histogram(bins = 25, fill = "steelblue", color = "white") +
  geom_vline(
    data = simulation_summary,
    aes(xintercept = true_value),
    color = "firebrick",
    linetype = "dashed",
    linewidth = 0.7
  ) +
  facet_wrap(~ parameter, scales = "free", ncol = 3) +
  labs(
    x = "Estimated parameter",
    y = "Number of simulated datasets",
    caption = "Dashed red lines indicate the true parameter values."
  ) +
  theme_minimal()

print(estimate_plot)

ggsave(
  file.path(out_dir, "simulation_point_estimate_distributions.png"),
  estimate_plot,
  width = 10,
  height = 6,
  dpi = 300
)

##############################################
# Bootstrap simulation
##############################################

n_simulations_boot <- 100
n_bootstrap <- 100
checkpoint_every <- 10
checkpoint_dir <- file.path(out_dir, "bootstrap_checkpoints")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(simulation_seed + 2)

# Convert crash-level bootstrap samples back to counts for model fitting.
crashes_to_counts <- function(crash_data) {
  two_counts <- crash_data %>%
    filter(VE_FORMS == 2) %>%
    count(cell_id, crash_type_4) %>%
    complete(cell_id = simulation_inputs$cell_id, crash_type_4 = two_types, fill = list(n = 0)) %>%
    pivot_wider(names_from = crash_type_4, values_from = n, values_fill = 0)
  
  one_counts <- crash_data %>%
    filter(VE_FORMS == 1) %>%
    count(cell_id, crash_type_4) %>%
    complete(cell_id = simulation_inputs$cell_id, crash_type_4 = one_types, fill = list(n = 0)) %>%
    pivot_wider(
      names_from = crash_type_4,
      values_from = n,
      values_fill = 0,
      names_prefix = "one_"
    )
  
  simulation_inputs %>%
    select(cell_id) %>%
    left_join(two_counts, by = "cell_id") %>%
    left_join(one_counts, by = "cell_id") %>%
    mutate(
      n2 = rowSums(across(all_of(two_types))),
      n1 = rowSums(across(all_of(paste0("one_", one_types))))
    )
}

bootstrap_estimates <- vector("list", n_simulations_boot * n_bootstrap)
bootstrap_diagnostics <- vector("list", n_simulations_boot * n_bootstrap)
original_point_estimates <- vector("list", n_simulations_boot)

index <- 1
start_simulation <- 1
checkpoint_model <- "script16_two_stage"

# Resume from the most recent checkpoint, if one exists.
progress_files <- list.files(
  checkpoint_dir,
  pattern = "^checkpoint_[0-9]+_progress\\.csv$",
  full.names = TRUE
)

if (length(progress_files) > 0) {
  checkpoint_numbers <- as.integer(str_extract(basename(progress_files), "(?<=checkpoint_)\\d+"))
  latest_checkpoint <- max(checkpoint_numbers)
  prefix <- file.path(checkpoint_dir, paste0("checkpoint_", sprintf("%03d", latest_checkpoint)))
  
  progress <- read_csv(paste0(prefix, "_progress.csv"), show_col_types = FALSE)
  
  if (!"model_version" %in% names(progress) || progress$model_version != checkpoint_model) {
    stop("Existing bootstrap checkpoints were created with an older model. Delete bootstrap_checkpoints and rerun.")
  }
  
  original_saved <- read_csv(paste0(prefix, "_original_point_estimates.csv"), show_col_types = FALSE)
  bootstrap_saved <- read_csv(paste0(prefix, "_bootstrap_estimates.csv"), show_col_types = FALSE)
  diagnostics_saved <- read_csv(paste0(prefix, "_diagnostics.csv"), show_col_types = FALSE)
  random_seed_saved <- read_csv(paste0(prefix, "_random_seed.csv"), show_col_types = FALSE)
  
  for (i in seq_len(latest_checkpoint)) {
    original_point_estimates[[i]] <- original_saved %>% filter(simulation == i)
  }
  
  for (k in seq_len(nrow(diagnostics_saved))) {
    bootstrap_diagnostics[[k]] <- diagnostics_saved[k, ]
  }
  
  estimate_groups <- bootstrap_saved %>%
    group_by(simulation, bootstrap) %>%
    group_split()
  
  for (k in seq_along(estimate_groups)) {
    bootstrap_estimates[[k]] <- estimate_groups[[k]]
  }
  
  index <- progress$index
  start_simulation <- progress$completed_simulation + 1
  .Random.seed <- as.integer(random_seed_saved$random_seed)
  
  cat("Resuming after simulated dataset", progress$completed_simulation, "\n")
}

for (i in seq.int(start_simulation, n_simulations_boot)) {
  
  # Generate one fresh fake dataset.
  simulated_i <- pmap(simulation_inputs, function(cell_id, ..., target_n2, target_n1, r_L, r_M, r_H) {
    simulate_cell(cell_id, target_n2, target_n1, r_L, r_M, r_H)
  })
  
  counts_i <- map_dfr(simulated_i, "counts")
  crashes_i <- map_dfr(simulated_i, "crash_level")
  
  # Fit the model to the original simulated dataset.
  fit_i <- tryCatch(fit_four_group(counts_i), error = function(e) e)
  
  if (inherits(fit_i, "error")) {
    original_point_estimates[[i]] <- tibble(
      simulation = i,
      parameter = NA_character_,
      true_value = NA_real_,
      estimate = NA_real_,
      convergence = NA_integer_
    )
  } else {
    original_point_estimates[[i]] <- fit_i$risk_estimates %>%
      transmute(
        simulation = i,
        parameter,
        true_value,
        estimate,
        convergence = fit_i$convergence
      )
  }
  
  # Bootstrap fatal crashes within each cell and vehicle type.
  for (j in seq_len(n_bootstrap)) {
    crashes_boot <- crashes_i %>%
      group_by(cell_id, VE_FORMS) %>%
      slice_sample(prop = 1, replace = TRUE) %>%
      ungroup()
    
    counts_boot <- crashes_to_counts(crashes_boot)
    fit_boot <- tryCatch(fit_four_group(counts_boot), error = function(e) e)
    
    if (inherits(fit_boot, "error")) {
      bootstrap_diagnostics[[index]] <- tibble(
        simulation = i,
        bootstrap = j,
        convergence = NA_integer_,
        two_car_convergence = NA_integer_,
        one_car_convergence = NA_integer_,
        message = conditionMessage(fit_boot)
      )
    } else {
      bootstrap_estimates[[index]] <- fit_boot$risk_estimates %>%
        transmute(
          simulation = i,
          bootstrap = j,
          parameter,
          true_value,
          estimate,
          convergence = fit_boot$convergence
        )
      
      bootstrap_diagnostics[[index]] <- tibble(
        simulation = i,
        bootstrap = j,
        convergence = fit_boot$convergence,
        two_car_convergence = fit_boot$two_car_convergence,
        one_car_convergence = fit_boot$one_car_convergence,
        message = ""
      )
    }
    
    index <- index + 1
  }
  
  cat("Completed simulated dataset", i, "of", n_simulations_boot, "\n")
  
  # Save a cumulative checkpoint every 10 simulated datasets.
  if (i %% checkpoint_every == 0 || i == n_simulations_boot) {
    prefix <- file.path(checkpoint_dir, paste0("checkpoint_", sprintf("%03d", i)))
    
    write_csv(
      bind_rows(original_point_estimates),
      paste0(prefix, "_original_point_estimates.csv")
    )
    write_csv(
      bind_rows(bootstrap_estimates),
      paste0(prefix, "_bootstrap_estimates.csv")
    )
    write_csv(
      bind_rows(bootstrap_diagnostics),
      paste0(prefix, "_diagnostics.csv")
    )
    write_csv(
      tibble(
        completed_simulation = i,
        index = index,
        model_version = checkpoint_model
      ),
      paste0(prefix, "_progress.csv")
    )
    write_csv(
      tibble(random_seed = .Random.seed),
      paste0(prefix, "_random_seed.csv")
    )
    
    cat("Saved checkpoint after simulated dataset", i, "\n")
  }
}

original_point_estimates <- bind_rows(original_point_estimates)
bootstrap_estimates <- bind_rows(bootstrap_estimates)
bootstrap_diagnostics <- bind_rows(bootstrap_diagnostics)

##############################################
# Compute bootstrap confidence intervals and coverage
##############################################

# Estimate bootstrap variance on the log scale, 
# then form a 95% interval around the original point estimate.
bootstrap_ci <- bootstrap_estimates %>%
  filter(convergence == 0) %>%
  mutate(estimate_log = log(estimate)) %>%
  group_by(simulation, parameter) %>%
  summarize(
    true_value = first(true_value),
    n_bootstrap_successful = n(),
    bootstrap_se_log = sd(estimate_log),
    .groups = "drop"
  ) %>%
  left_join(
    original_point_estimates %>%
      filter(convergence == 0) %>%
      select(simulation, parameter, original_estimate = estimate),
    by = c("simulation", "parameter")
  ) %>%
  mutate(
    lower = exp(log(original_estimate) - 1.96 * bootstrap_se_log),
    upper = exp(log(original_estimate) + 1.96 * bootstrap_se_log),
    covers_truth = lower <= true_value & upper >= true_value
  )

coverage_summary <- bootstrap_ci %>%
  filter(!is.na(covers_truth)) %>%
  group_by(parameter) %>%
  summarize(
    true_value = first(true_value),
    n_simulations = n(),
    coverage = mean(covers_truth),
    .groups = "drop"
  )

print(bootstrap_ci, width = Inf)
print(coverage_summary, width = Inf)

write_csv(original_point_estimates, file.path(out_dir, "bootstrap_original_point_estimates.csv"))
write_csv(bootstrap_estimates, file.path(out_dir, "bootstrap_sampled_estimates.csv"))
write_csv(bootstrap_diagnostics, file.path(out_dir, "bootstrap_diagnostics.csv"))
write_csv(bootstrap_ci, file.path(out_dir, "bootstrap_confidence_intervals.csv"))
write_csv(coverage_summary, file.path(out_dir, "bootstrap_coverage_summary.csv"))
