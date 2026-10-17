# Fit the four-group model to each of the ten completed BAC datasets, then
# quantify uncertainty with 1,000 stratified crash-level bootstrap samples.
# Refit both likelihood stages for every sample and pool on the log scale.
# Four parallel workers to speed up the bootstrapping (note: increases memory use).
#
# Input:  data/processed/script03_four_group_prep/four_group_crash_level_by_imp.csv
# Outputs: results/four_group_bootstrap/
#   four_group_bootstrap_combined_results.csv
#   four_group_bootstrap_replicate_estimates.csv
#   four_group_original_imputation_estimates.csv
#   four_group_original_cell_exposure_ratios.csv
#   checkpoints/bootstrap_replicates_XXXX_XXXX.csv
# Existing results and checkpoints are reused by default; no need to rerun.

# rm(list=ls())
library(tidyverse)

crash_path <- "data/processed/script03_four_group_prep/four_group_crash_level_by_imp.csv"
out_dir <- "results/four_group_bootstrap"
checkpoint_dir <- file.path(out_dir, "checkpoints")
combined_path <- file.path(out_dir, "four_group_bootstrap_combined_results.csv")
replicate_path <- file.path(out_dir, "four_group_bootstrap_replicate_estimates.csv")
original_path <- file.path(out_dir, "four_group_original_imputation_estimates.csv")
exposure_path <- file.path(out_dir, "four_group_original_cell_exposure_ratios.csv")
force_recompute <- FALSE  # TRUE discards reuse and reruns all fits and bootstrap samples.

bootstrap_reps <- 1000
bootstrap_seed <- 123
parallel_workers <- 5
checkpoint_every <- 5
min_n2 <- 1
log_param_lower <- log(1e-4)
log_param_upper <- log(500)
optim_control <- list(maxit = 5000, factr = 1e7)

cell_specs <- tibble(
  cell_spec = "rur_urb_weekend",
  cell_var = "CELL_RUR_URB_WEEKEND_10",
  cell_spec_label = "rural/urban x weekend"
)

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

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# An existing final result takes priority: skip file loading, fitting, and resampling.
if (!force_recompute && file.exists(combined_path)) {
  message("Reusing existing results: ", combined_path)
  message("No bootstrap was run. Set force_recompute <- TRUE to refit.")
  print(read_csv(combined_path, show_col_types = FALSE), width = Inf)
} else {
  safe_log <- function(x) log(pmax(x, 1e-300))

  safe_divide <- function(num, den) {
    ifelse(is.na(den) | den == 0, NA_real_, num / den)
  }

  clamp <- function(x, lower, upper) {
    pmin(pmax(x, lower), upper)
  }

  bad_likelihood_penalty <- 1e100

  make_cell_counts <- function(crash_data, cell_var, bootstrap_weight = NULL) {
    if (is.null(bootstrap_weight)) bootstrap_weight <- rep(1L, nrow(crash_data))

    crash_data %>%
      mutate(cell_id = as.character(.data[[cell_var]]), bootstrap_weight = bootstrap_weight) %>%
      group_by(cell_id) %>%
      summarize(
        SS = sum(bootstrap_weight * (crash_type_4 == "SS"), na.rm = TRUE),
        SL = sum(bootstrap_weight * (crash_type_4 == "SL"), na.rm = TRUE),
        SM = sum(bootstrap_weight * (crash_type_4 == "SM"), na.rm = TRUE),
        SH = sum(bootstrap_weight * (crash_type_4 == "SH"), na.rm = TRUE),
        LL = sum(bootstrap_weight * (crash_type_4 == "LL"), na.rm = TRUE),
        LM = sum(bootstrap_weight * (crash_type_4 == "LM"), na.rm = TRUE),
        LH = sum(bootstrap_weight * (crash_type_4 == "LH"), na.rm = TRUE),
        MM = sum(bootstrap_weight * (crash_type_4 == "MM"), na.rm = TRUE),
        MH = sum(bootstrap_weight * (crash_type_4 == "MH"), na.rm = TRUE),
        HH = sum(bootstrap_weight * (crash_type_4 == "HH"), na.rm = TRUE),
        one_S = sum(bootstrap_weight * (crash_type_4 == "S"), na.rm = TRUE),
        one_L = sum(bootstrap_weight * (crash_type_4 == "L"), na.rm = TRUE),
        one_M = sum(bootstrap_weight * (crash_type_4 == "M"), na.rm = TRUE),
        one_H = sum(bootstrap_weight * (crash_type_4 == "H"), na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        n2 = SS + SL + SM + SH + LL + LM + LH + MM + MH + HH,
        n1 = one_S + one_L + one_M + one_H
      ) %>%
      filter(n1 + n2 > 0)
  }

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
    unpacked <- unpack_two_car_par(fit$par, cell_ids)

    list(
      fit = fit,
      best_start_id = best_id,
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
      transmute(cell_id, q_S = w_S / D, q_L = w_L / D, q_M = w_M / D, q_H = w_H / D)
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
    lambda <- c(
      S = 1,
      L = unname(exp(fit$par["log_lambda_L"])),
      M = unname(exp(fit$par["log_lambda_M"])),
      H = unname(exp(fit$par["log_lambda_H"]))
    )

    list(
      fit = fit,
      best_start_id = best_id,
      lambda = lambda,
      n_cells = nrow(dat),
      n1_total = sum(dat$n1)
    )
  }

  estimate_one_imp_spec <- function(dat) {
    two_fit <- estimate_two_car_stage(dat)
    one_fit <- estimate_one_car_stage(dat, two_fit$r_mat)

    estimates <- bind_rows(
      tibble(
        stage = "two_car",
        parameter = c("theta_L", "theta_M", "theta_H"),
        estimate = as.numeric(two_fit$theta[c("L", "M", "H")]),
        estimate_log = log(estimate),
        n_cells = two_fit$n_cells,
        n_crashes = two_fit$n2_total,
        optim_value = two_fit$fit$value,
        optim_convergence = two_fit$fit$convergence,
        best_start_id = two_fit$best_start_id
      ),
      tibble(
        stage = "one_car",
        parameter = c("lambda_L", "lambda_M", "lambda_H"),
        estimate = as.numeric(one_fit$lambda[c("L", "M", "H")]),
        estimate_log = log(estimate),
        n_cells = one_fit$n_cells,
        n_crashes = one_fit$n1_total,
        optim_value = one_fit$fit$value,
        optim_convergence = one_fit$fit$convergence,
        best_start_id = one_fit$best_start_id
      )
    )

    list(estimates = estimates, r_mat = two_fit$r_mat)
  }

  replicate_batches <- split(seq_len(bootstrap_reps),
                             ceiling(seq_len(bootstrap_reps) / checkpoint_every))
  checkpoint_path <- function(batch) {
    file.path(checkpoint_dir,
              sprintf("bootstrap_replicates_%04d_%04d.csv", min(batch), max(batch)))
  }
  checkpoint_paths <- map_chr(replicate_batches, checkpoint_path)

  # Only read the large crash-level dataset if original fits or bootstrap
  # samples actually need computing. Existing outputs can otherwise be pooled.
  need_bootstrap <- force_recompute ||
    (!file.exists(replicate_path) && !all(file.exists(checkpoint_paths)))
  need_original <- force_recompute || !file.exists(original_path)

  if (need_original || need_bootstrap) {
    crash_all <- read_csv(crash_path, show_col_types = FALSE) %>%
      mutate(
        IMPUTATION = as.integer(IMPUTATION),
        CRASH_ID_10 = as.character(CRASH_ID_10),
        crash_type_4 = as.character(crash_type_4)
      )
    imputations <- sort(unique(crash_all$IMPUTATION))
    crashes_by_imp <- split(crash_all, crash_all$IMPUTATION)
  }

  # The original (unresampled) estimates supply the between-imputation
  # variance. Reuse the original CSV whenever possible.
  if (!need_original) {
    message("Reusing original-data estimates: ", original_path)
    original_results <- read_csv(original_path, show_col_types = FALSE) %>%
      select(IMPUTATION, cell_spec, cell_spec_label, stage, parameter,
             estimate, estimate_log, n_cells, n_crashes, optim_value,
             optim_convergence, best_start_id)
  } else {
    original_results <- list()
    original_exposures <- list()
    result_index <- 1

    for (s in seq_len(nrow(cell_specs))) {
      for (imp in imputations) {
        crash_imp <- crashes_by_imp[[as.character(imp)]]
        counts_imp <- make_cell_counts(crash_imp, cell_specs$cell_var[s])
        fit_imp <- estimate_one_imp_spec(counts_imp)
        original_results[[result_index]] <- fit_imp$estimates %>%
          mutate(IMPUTATION = imp, cell_spec = cell_specs$cell_spec[s],
                 cell_spec_label = cell_specs$cell_spec_label[s], .before = 1)
        original_exposures[[result_index]] <- fit_imp$r_mat %>%
          mutate(IMPUTATION = imp, cell_spec = cell_specs$cell_spec[s],
                 cell_spec_label = cell_specs$cell_spec_label[s], .before = 1)
        result_index <- result_index + 1
      }
    }
    original_results <- bind_rows(original_results)
    original_exposures <- bind_rows(original_exposures)
    # These files allow later pooling without refitting and provide
    # the cell exposure estimates that may be useful for simulation validation.
    write_csv(original_results, original_path)
    write_csv(original_exposures, exposure_path)
  }

  # Reuse the aggregate replicate CSV first, then complete checkpoints, before
  # ever generating new bootstrap samples.
  if (!force_recompute && file.exists(replicate_path)) {
    message("Reusing existing bootstrap estimates: ", replicate_path)
    bootstrap_results <- read_csv(replicate_path, show_col_types = FALSE)
  } else {
    if (!force_recompute && all(file.exists(checkpoint_paths))) {
      message("All bootstrap checkpoints found; combining without refitting.")
    } else {
      dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

      # Sample crashes separately by cell and crash size. The same sampled
      # crash weights are applied to all ten BAC imputations.
      bootstrap_crashes_by_imp <- map(crashes_by_imp, ~ .x %>%
        select(CRASH_ID_10, VE_FORMS, crash_type_4,
               all_of(unique(cell_specs$cell_var))))
      base_crashes <- bootstrap_crashes_by_imp[[1]] %>% arrange(CRASH_ID_10)
      bootstrap_specs <- map(seq_len(nrow(cell_specs)), function(s) {
        base_spec <- base_crashes %>%
          transmute(CRASH_ID_10, VE_FORMS,
                    cell_id = as.character(.data[[cell_specs$cell_var[s]]]))
        list(
          base_spec = base_spec,
          strata_indices = split(seq_len(nrow(base_spec)),
                                  interaction(base_spec$cell_id, base_spec$VE_FORMS, drop = TRUE)),
          match_indices = map(bootstrap_crashes_by_imp,
                              ~ match(.x$CRASH_ID_10, base_spec$CRASH_ID_10))
        )
      })

      run_bootstrap_replicate <- function(b) {
        set.seed(bootstrap_seed + b)
        replicate_results <- vector("list", length(imputations) * nrow(cell_specs))
        result_index <- 1

        for (s in seq_len(nrow(cell_specs))) {
          base_spec <- bootstrap_specs[[s]]$base_spec
          strata_indices <- bootstrap_specs[[s]]$strata_indices
          bootstrap_weight <- integer(nrow(base_spec))

          for (idx in strata_indices) {
            bootstrap_weight[idx] <- tabulate(
              sample.int(length(idx), length(idx), replace = TRUE), nbins = length(idx))
          }

          for (imp in imputations) {
            crash_imp <- bootstrap_crashes_by_imp[[as.character(imp)]]
            sample_weight <- bootstrap_weight[
              bootstrap_specs[[s]]$match_indices[[as.character(imp)]]]

            replicate_results[[result_index]] <- tryCatch(
              {
                counts_imp <- make_cell_counts(crash_imp, cell_specs$cell_var[s], sample_weight)
                estimate_one_imp_spec(counts_imp)$estimates %>%
                  mutate(replicate = b, IMPUTATION = imp,
                         cell_spec = cell_specs$cell_spec[s],
                         cell_spec_label = cell_specs$cell_spec_label[s], .before = 1)
              },
              error = function(e) tibble(
                replicate = b, IMPUTATION = imp,
                cell_spec = cell_specs$cell_spec[s],
                cell_spec_label = cell_specs$cell_spec_label[s],
                stage = rep(c("two_car", "one_car"), each = 3),
                parameter = c("theta_L", "theta_M", "theta_H",
                              "lambda_L", "lambda_M", "lambda_H"),
                estimate = NA_real_, estimate_log = NA_real_,
                n_cells = NA_integer_, n_crashes = NA_real_,
                optim_value = NA_real_, optim_convergence = NA_integer_,
                best_start_id = NA_integer_
              )
            )
            result_index <- result_index + 1
          }
        }
        bind_rows(replicate_results)
      }

      # Export everything the parallel workers need
      worker_objects <- c(
        "safe_log", "safe_divide", "clamp", "bad_likelihood_penalty", "make_cell_counts",
        "make_r_starts_from_two_car_counts", "unpack_two_car_par", "make_two_car_probabilities",
        "negloglik_two_car", "make_two_car_start", "estimate_two_car_stage",
        "make_one_car_probabilities", "negloglik_one_car", "make_one_car_start",
        "estimate_one_car_stage", "estimate_one_imp_spec", "run_bootstrap_replicate",
        "bootstrap_seed", "bootstrap_specs", "bootstrap_crashes_by_imp", "imputations",
        "cell_specs", "theta_start_grid", "lambda_start_grid", "min_n2",
        "log_param_lower", "log_param_upper", "optim_control"
      )

      pending_batches <- if (force_recompute) replicate_batches else
        replicate_batches[!file.exists(checkpoint_paths)]
      message("Running ", sum(lengths(pending_batches)),
              " missing/recomputed bootstrap replicates with ", parallel_workers, " workers.")

      cluster <- parallel::makeCluster(parallel_workers)
      tryCatch(
        {
          parallel::clusterEvalQ(cluster, library(tidyverse))
          parallel::clusterExport(cluster, worker_objects, envir = environment())
          for (batch in pending_batches) {
            batch_path <- checkpoint_path(batch)
            batch_results <- parallel::parLapplyLB(cluster, batch, run_bootstrap_replicate) %>%
              bind_rows()
            temporary_path <- paste0(batch_path, ".tmp")
            write_csv(batch_results, temporary_path)
            if (file.exists(batch_path)) file.remove(batch_path)
            if (!file.rename(temporary_path, batch_path)) {
              stop("Cannot save checkpoint: ", batch_path)
            }
            message("Completed through replicate ", max(batch), " of ", bootstrap_reps)
          }
        },
        finally = parallel::stopCluster(cluster)
      )
    }

    # Store replicate estimates so future pooling never needs resampling.
    if (!all(file.exists(checkpoint_paths))) stop("Missing bootstrap checkpoints.")
    bootstrap_results <- map_dfr(checkpoint_paths,
                                 ~ read_csv(.x, show_col_types = FALSE))
    write_csv(bootstrap_results, replicate_path)
  }

  # Rubin's rules on the log scale, using bootstrap sampling variance within
  # each completed dataset and original-fit variability between datasets.
  within_variances <- bootstrap_results %>%
    group_by(IMPUTATION, cell_spec, cell_spec_label, stage, parameter) %>%
    summarize(
      bootstrap_reps_requested = bootstrap_reps,
      bootstrap_reps_used = sum(is.finite(estimate_log)),
      within_var_log = var(estimate_log, na.rm = TRUE),
      bootstrap_se_log = sqrt(within_var_log),
      .groups = "drop"
    )

  original_results <- original_results %>%
    left_join(within_variances,
              by = c("IMPUTATION", "cell_spec", "cell_spec_label", "stage", "parameter"))

  combined_results <- original_results %>%
    group_by(cell_spec, cell_spec_label, stage, parameter) %>%
    summarize(
      m = n(),
      estimate_log_bar = mean(estimate_log),
      within_var_log = mean(within_var_log),
      between_var_log = var(estimate_log),
      total_var_log = within_var_log + (1 + 1 / m) * between_var_log,
      se_log_total = sqrt(total_var_log),
      estimate = exp(estimate_log_bar),
      ci_low = exp(estimate_log_bar - 1.96 * se_log_total),
      ci_high = exp(estimate_log_bar + 1.96 * se_log_total),
      min_bootstrap_reps_used = min(bootstrap_reps_used),
      mean_n_crashes = mean(n_crashes),
      .groups = "drop"
    )

  paired_percentile_results <- bootstrap_results %>%
    group_by(replicate, cell_spec, cell_spec_label, stage, parameter) %>%
    summarize(estimate = exp(mean(estimate_log, na.rm = TRUE)), .groups = "drop") %>%
    group_by(cell_spec, cell_spec_label, stage, parameter) %>%
    summarize(
      percentile_ci_low = quantile(estimate, 0.025, na.rm = TRUE),
      percentile_ci_high = quantile(estimate, 0.975, na.rm = TRUE),
      .groups = "drop"
    )

  combined_results <- combined_results %>%
    left_join(paired_percentile_results,
              by = c("cell_spec", "cell_spec_label", "stage", "parameter")) %>%
    arrange(stage, parameter)

  write_csv(combined_results, combined_path)
  print(combined_results, width = Inf)
  message("Saved combined results: ", combined_path)
}
