################################################################################
# 11_four_group_mle.R
#
# Purpose:
#   Estimate the four-BAC-group Levitt-Porter-style MLE using the crash-count
#   data prepared by 10_four_group_prep.R.
#
# Driver groups:
#   S = BAC == 0                         baseline group
#   L = 0 < BAC < 0.05
#   M = 0.05 <= BAC < 0.08
#   H = BAC >= 0.08
#
# Parameters estimated in this script:
#   Two-car stage:
#     theta_L, theta_M, theta_H       global relative two-car fatal crash risks
#     r_Lc, r_Mc, r_Hc               cell-specific road exposure ratios vs S
#
#   One-car stage:
#     lambda_L, lambda_M, lambda_H    global relative one-car fatal crash risks
#
# Model details:
#   This script implements the probabilities and log likelihood from the
#   Overleaf four-group model notes. The main difference from scripts 04/07 is
#   that the binary drunk/sober categories are replaced by four BAC groups.
#
# Main input:
#   data/processed/script10_four_group_prep/four_group_crash_counts_by_cell_imp.rds
#
# Main outputs:
#   data/processed/script11_four_group_mle/
#     four_group_mle_estimates_by_imp.csv
#     four_group_mle_estimates_pooled.csv
#     four_group_mle_cell_exposure_ratios.csv
#     four_group_mle_fit_diagnostics.csv
#     four_group_mle_results_bundle.rds
#
# Important implementation choice:
#   All positive parameters are estimated on the log scale.
#   Example: eta_L = log(theta_L), so theta_L = exp(eta_L).
#   This guarantees all theta/lambda/r estimates are positive.
################################################################################

rm(list = ls())

library(tidyverse)
library(scales)

# ------------------------------------------------------------------------------
# 0. User options
# ------------------------------------------------------------------------------

prep_dir <- "data/processed/script10_four_group_prep"
out_dir  <- "data/processed/script11_four_group_mle"
fig_dir  <- "figures/script11_four_group_mle"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

counts_path <- file.path(prep_dir, "four_group_crash_counts_by_cell_imp.rds")

# The four coarse specifications produced by script 10.
cell_specs_to_estimate <- c("all_data", "weekend", "rur_urb", "rur_urb_weekend")

# Keep cells with at least this many complete two-car crashes for the two-car MLE.
# For the current coarse specs, this should not remove much/anything, but it is
# useful to keep the option explicit.
min_n2 <- 1

# Parameter bounds on the log scale. These are intentionally wide but finite.
# Natural-scale interpretation: exp(log_param_lower) to exp(log_param_upper).
# exp(log(1e-4)) = 0.0001 and exp(log(500)) = 500.
log_param_lower <- log(1e-4)
log_param_upper <- log(500)

# Optimizer controls.
optim_control <- list(maxit = 5000, factr = 1e7)

# Multiple starting values make the numerical optimizer less dependent on one
# arbitrary initialization.
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
# 1. Helper functions
# ------------------------------------------------------------------------------

safe_log <- function(x) log(pmax(x, 1e-300))

safe_divide <- function(num, den) {
  ifelse(is.na(den) | den == 0, NA_real_, num / den)
}

clamp <- function(x, lower, upper) {
  pmin(pmax(x, lower), upper)
}

# If a probability calculation accidentally produces NA/Inf, return a large
# penalty rather than letting sum(..., na.rm = TRUE) hide the problem.
# This makes optimizer/debugging failures much easier to see.
bad_likelihood_penalty <- 1e100

# Approximate starting values for r_Lc, r_Mc, r_Hc using two-car crash endpoints.
# If all theta values were 1, the two-car driver endpoint composition would be
# roughly proportional to the road exposure composition. This is only a starting
# value, not a final estimate.
make_r_starts_from_two_car_counts <- function(dat) {
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

# Convert a named parameter vector on the log scale into theta and r objects.
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

# Compute the 10 two-car probabilities for every cell.
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

# Negative two-car log likelihood. Constants from the multinomial likelihood are
# ignored because they do not affect the MLE.
negloglik_two_car <- function(par, dat, cell_ids) {
  unpacked <- unpack_two_car_par(par, cell_ids)
  probs <- make_two_car_probabilities(unpacked$theta, unpacked$r_mat)

  ll_dat <- dat %>%
    left_join(probs, by = "cell_id")

  prob_cols <- c("p_SS", "p_SL", "p_SM", "p_SH", "p_LL", "p_LM", "p_LH", "p_MM", "p_MH", "p_HH")
  if (any(!is.finite(as.matrix(ll_dat[, prob_cols])))) return(bad_likelihood_penalty)

  loglik_terms <-
    ll_dat$SS * safe_log(ll_dat$p_SS) +
    ll_dat$SL * safe_log(ll_dat$p_SL) +
    ll_dat$SM * safe_log(ll_dat$p_SM) +
    ll_dat$SH * safe_log(ll_dat$p_SH) +
    ll_dat$LL * safe_log(ll_dat$p_LL) +
    ll_dat$LM * safe_log(ll_dat$p_LM) +
    ll_dat$LH * safe_log(ll_dat$p_LH) +
    ll_dat$MM * safe_log(ll_dat$p_MM) +
    ll_dat$MH * safe_log(ll_dat$p_MH) +
    ll_dat$HH * safe_log(ll_dat$p_HH)

  if (any(!is.finite(loglik_terms))) return(bad_likelihood_penalty)

  -sum(loglik_terms)
}

# Build one starting parameter vector for the two-car optimizer.
make_two_car_start <- function(dat, theta_start) {
  r_starts <- make_r_starts_from_two_car_counts(dat)

  par <- c(
    log_theta_L = log(unname(theta_start["theta_L"])),
    log_theta_M = log(unname(theta_start["theta_M"])),
    log_theta_H = log(unname(theta_start["theta_H"]))
  )

  for (i in seq_len(nrow(r_starts))) {
    cc <- r_starts$cell_id[i]
    par[paste0("log_r_L__", cc)] <- log(clamp(r_starts$r_L[i], 1e-4, 500))
    par[paste0("log_r_M__", cc)] <- log(clamp(r_starts$r_M[i], 1e-4, 500))
    par[paste0("log_r_H__", cc)] <- log(clamp(r_starts$r_H[i], 1e-4, 500))
  }

  par
}

# Run the two-car MLE for one imputation and one cell specification.
estimate_two_car_stage <- function(dat) {
  dat <- dat %>% filter(n2 >= min_n2) %>% arrange(cell_id)
  cell_ids <- dat$cell_id

  starts <- map(theta_start_grid, ~ make_two_car_start(dat, .x))

  fits <- map(starts, function(start_par) {
    optim(
      par = start_par,
      fn = negloglik_two_car,
      dat = dat,
      cell_ids = cell_ids,
      method = "L-BFGS-B",
      lower = rep(log_param_lower, length(start_par)),
      upper = rep(log_param_upper, length(start_par)),
      control = optim_control
    )
  })

  best_id <- which.min(map_dbl(fits, "value"))
  fit <- fits[[best_id]]

  # Hessian on log-parameter scale. This is used only for approximate standard
  # errors and MI pooling, not for the point estimates.
  hess <- tryCatch(
    optimHess(
      par = fit$par,
      fn = negloglik_two_car,
      dat = dat,
      cell_ids = cell_ids
    ),
    error = function(e) NULL
  )

  vcov_log <- NULL
  se_log <- rep(NA_real_, length(fit$par))
  names(se_log) <- names(fit$par)
  hessian_ok <- FALSE

  if (!is.null(hess) && all(is.finite(hess))) {
    eig <- tryCatch(eigen(hess, symmetric = TRUE, only.values = TRUE)$values, error = function(e) NA_real_)
    if (all(is.finite(eig)) && min(eig) > 1e-8) {
      vcov_log <- solve(hess)
      se_log <- sqrt(pmax(diag(vcov_log), 0))
      names(se_log) <- names(fit$par)
      hessian_ok <- TRUE
    }
  }

  unpacked <- unpack_two_car_par(fit$par, cell_ids)

  list(
    fit = fit,
    best_start_id = best_id,
    hessian_ok = hessian_ok,
    se_log = se_log,
    theta = unpacked$theta,
    r_mat = unpacked$r_mat,
    n_cells = nrow(dat),
    n2_total = sum(dat$n2)
  )
}

# One-car probabilities q_g in each cell, given lambda and the two-car r estimates.
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

  ll_dat <- dat %>%
    left_join(probs, by = "cell_id")

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

make_one_car_start <- function(lambda_start) {
  c(
    log_lambda_L = log(unname(lambda_start["lambda_L"])),
    log_lambda_M = log(unname(lambda_start["lambda_M"])),
    log_lambda_H = log(unname(lambda_start["lambda_H"]))
  )
}

estimate_one_car_stage <- function(dat, r_mat) {
  dat <- dat %>% filter(n1 > 0) %>% arrange(cell_id)
  r_mat <- r_mat %>% filter(cell_id %in% dat$cell_id) %>% arrange(cell_id)

  starts <- map(lambda_start_grid, make_one_car_start)

  fits <- map(starts, function(start_par) {
    optim(
      par = start_par,
      fn = negloglik_one_car,
      dat = dat,
      r_mat = r_mat,
      method = "L-BFGS-B",
      lower = rep(log_param_lower, length(start_par)),
      upper = rep(log_param_upper, length(start_par)),
      control = optim_control
    )
  })

  best_id <- which.min(map_dbl(fits, "value"))
  fit <- fits[[best_id]]

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
    eig <- tryCatch(eigen(hess, symmetric = TRUE, only.values = TRUE)$values, error = function(e) NA_real_)
    if (all(is.finite(eig)) && min(eig) > 1e-8) {
      vcov_log <- solve(hess)
      se_log <- sqrt(pmax(diag(vcov_log), 0))
      names(se_log) <- names(fit$par)
      hessian_ok <- TRUE
    }
  }

  lambda <- c(
    S = 1,
    L = unname(exp(fit$par["log_lambda_L"])),
    M = unname(exp(fit$par["log_lambda_M"])),
    H = unname(exp(fit$par["log_lambda_H"]))
  )

  list(
    fit = fit,
    best_start_id = best_id,
    hessian_ok = hessian_ok,
    se_log = se_log,
    lambda = lambda,
    n_cells = nrow(dat),
    n1_total = sum(dat$n1)
  )
}

# Estimate both stages for one imputation/spec combination.
estimate_one_imp_spec <- function(dat_one) {
  two_fit <- estimate_two_car_stage(dat_one)
  one_fit <- estimate_one_car_stage(dat_one, two_fit$r_mat)

  list(two_fit = two_fit, one_fit = one_fit)
}

# Rubin-style pooling on the log scale. This is most natural because all point
# estimates are positive, and the optimizer works on log parameters.
pool_log_estimates <- function(est_by_imp) {
  est_by_imp %>%
    group_by(cell_spec, stage, parameter) %>%
    summarize(
      m = n(),
      estimate_log_bar = mean(estimate_log, na.rm = TRUE),
      # Hessian SEs can fail for sparse/flat likelihoods. Keep point estimates
      # pooled even if se_log is missing, and report interval columns only when
      # we have usable SEs.
      within_var_log = ifelse(all(is.na(se_log)), NA_real_, mean(se_log^2, na.rm = TRUE)),
      between_var_log = var(estimate_log, na.rm = TRUE),
      total_var_log = ifelse(
        is.na(within_var_log),
        NA_real_,
        within_var_log + (1 + 1 / m) * between_var_log
      ),
      se_log_total = sqrt(total_var_log),
      estimate = exp(estimate_log_bar),
      ci_low = ifelse(is.na(se_log_total), NA_real_, exp(estimate_log_bar - 1.96 * se_log_total)),
      ci_high = ifelse(is.na(se_log_total), NA_real_, exp(estimate_log_bar + 1.96 * se_log_total)),
      mean_n_crashes = mean(n_crashes, na.rm = TRUE),
      .groups = "drop"
    )
}

# ------------------------------------------------------------------------------
# 2. Load prepared crash-count data
# ------------------------------------------------------------------------------

if (!file.exists(counts_path)) {
  stop(
    "Missing input file: ", counts_path,
    "\nPlease run 10_four_group_prep.R first."
  )
}

counts_all <- read_rds(counts_path) %>%
  rename_with(toupper, .cols = "IMPUTATION") %>%
  mutate(cell_id = as.character(cell_id)) %>%
  filter(cell_spec %in% cell_specs_to_estimate)

cat("Loaded script 10 crash-count data:\n")
cat("  Rows:", comma(nrow(counts_all)), "\n")
cat("  Cell specs:", paste(sort(unique(counts_all$cell_spec)), collapse = ", "), "\n")
cat("  Imputations:", paste(sort(unique(counts_all$IMPUTATION)), collapse = ", "), "\n\n")

required_count_cols <- c(
  "IMPUTATION", "cell_spec", "cell_id",
  "SS", "SL", "SM", "SH", "LL", "LM", "LH", "MM", "MH", "HH",
  "one_S", "one_L", "one_M", "one_H", "n2", "n1"
)
missing_count_cols <- setdiff(required_count_cols, names(counts_all))
if (length(missing_count_cols) > 0) {
  stop("Missing required count columns: ", paste(missing_count_cols, collapse = ", "))
}

# ------------------------------------------------------------------------------
# 3. Run MLE separately for each imputation and cell specification
# ------------------------------------------------------------------------------

analysis_keys <- counts_all %>%
  distinct(IMPUTATION, cell_spec) %>%
  arrange(cell_spec, IMPUTATION)

all_results <- vector("list", nrow(analysis_keys))

for (ii in seq_len(nrow(analysis_keys))) {
  imp_i <- analysis_keys$IMPUTATION[ii]
  spec_i <- analysis_keys$cell_spec[ii]

  cat("Estimating:", "cell_spec =", spec_i, ", IMPUTATION =", imp_i, "\n")

  dat_i <- counts_all %>%
    filter(IMPUTATION == imp_i, cell_spec == spec_i) %>%
    arrange(cell_id)

  fit_i <- estimate_one_imp_spec(dat_i)

  all_results[[ii]] <- list(
    IMPUTATION = imp_i,
    cell_spec = spec_i,
    fit = fit_i
  )
}

# ------------------------------------------------------------------------------
# 4. Tidy estimates and diagnostics
# ------------------------------------------------------------------------------

estimates_by_imp <- map_dfr(all_results, function(res) {
  imp_i <- res$IMPUTATION
  spec_i <- res$cell_spec
  two_fit <- res$fit$two_fit
  one_fit <- res$fit$one_fit

  bind_rows(
    tibble(
      IMPUTATION = imp_i,
      cell_spec = spec_i,
      stage = "two_car",
      parameter = c("theta_L", "theta_M", "theta_H"),
      estimate = as.numeric(two_fit$theta[c("L", "M", "H")]),
      estimate_log = log(estimate),
      se_log = as.numeric(two_fit$se_log[c("log_theta_L", "log_theta_M", "log_theta_H")]),
      se_delta = estimate * se_log,
      n_cells = two_fit$n_cells,
      n_crashes = two_fit$n2_total,
      optim_value = two_fit$fit$value,
      optim_convergence = two_fit$fit$convergence,
      hessian_ok = two_fit$hessian_ok,
      best_start_id = two_fit$best_start_id
    ),
    tibble(
      IMPUTATION = imp_i,
      cell_spec = spec_i,
      stage = "one_car",
      parameter = c("lambda_L", "lambda_M", "lambda_H"),
      estimate = as.numeric(one_fit$lambda[c("L", "M", "H")]),
      estimate_log = log(estimate),
      se_log = as.numeric(one_fit$se_log[c("log_lambda_L", "log_lambda_M", "log_lambda_H")]),
      se_delta = estimate * se_log,
      n_cells = one_fit$n_cells,
      n_crashes = one_fit$n1_total,
      optim_value = one_fit$fit$value,
      optim_convergence = one_fit$fit$convergence,
      hessian_ok = one_fit$hessian_ok,
      best_start_id = one_fit$best_start_id
    )
  )
})

cell_exposure_ratios <- map_dfr(all_results, function(res) {
  res$fit$two_fit$r_mat %>%
    mutate(
      IMPUTATION = res$IMPUTATION,
      cell_spec = res$cell_spec,
      .before = 1
    )
})

fit_diagnostics <- estimates_by_imp %>%
  distinct(
    IMPUTATION, cell_spec, stage, n_cells, n_crashes,
    optim_value, optim_convergence, hessian_ok, best_start_id
  ) %>%
  arrange(cell_spec, IMPUTATION, stage)

estimates_pooled <- pool_log_estimates(estimates_by_imp)

# ------------------------------------------------------------------------------
# 5. Save outputs
# ------------------------------------------------------------------------------

write_csv(estimates_by_imp, file.path(out_dir, "four_group_mle_estimates_by_imp.csv"))
write_csv(estimates_pooled, file.path(out_dir, "four_group_mle_estimates_pooled.csv"))
write_csv(cell_exposure_ratios, file.path(out_dir, "four_group_mle_cell_exposure_ratios.csv"))
write_csv(fit_diagnostics, file.path(out_dir, "four_group_mle_fit_diagnostics.csv"))

saveRDS(
  list(
    counts_all = counts_all,
    all_results = all_results,
    estimates_by_imp = estimates_by_imp,
    estimates_pooled = estimates_pooled,
    cell_exposure_ratios = cell_exposure_ratios,
    fit_diagnostics = fit_diagnostics,
    options = list(
      min_n2 = min_n2,
      cell_specs_to_estimate = cell_specs_to_estimate,
      log_param_lower = log_param_lower,
      log_param_upper = log_param_upper,
      theta_start_grid = theta_start_grid,
      lambda_start_grid = lambda_start_grid
    )
  ),
  file.path(out_dir, "four_group_mle_results_bundle.rds")
)
