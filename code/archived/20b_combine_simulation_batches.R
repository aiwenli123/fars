################################################################################
# Combine and summarize the 100 simulation batches (10 datasets per batch)
# 1,000 total simulated datasets, with 200 bootstrap samples per dataset
################################################################################

library(tidyverse)

setwd("~/projects/aiwen/fars")

input_dir <- "results/simulation_validation_full"
output_dir <- file.path(input_dir, "combined")

parameters <- c("theta_L", "theta_M", "theta_H", "lambda_L", "lambda_M", "lambda_H")
true_values <- c(theta_L = 1.2, theta_M = 1.6, theta_H = 3.5,
                 lambda_L = 1.3, lambda_M = 2.0, lambda_H = 4.2)

# Read each batch and assign unique dataset IDs (1 through 1,000).
read_batches <- function(filename) {
  map_dfr(1:100, function(b) {
    path <- file.path(input_dir, sprintf("batch_%03d", b), filename)
    if (!file.exists(path)) stop("Missing file: ", path)

    read_csv(path, show_col_types = FALSE, progress = FALSE) %>%
      mutate(dataset_id = 10 * (b - 1) + simulation)
  })
}

original <- read_batches("bootstrap_original_point_estimates.csv")
diagnostics <- read_batches("bootstrap_diagnostics.csv")
intervals <- read_batches("bootstrap_confidence_intervals.csv")

# Check convergence and completeness of all original and bootstrap fits.
original_ok <- nrow(original) == 6000 &&
  n_distinct(original$dataset_id, original$parameter) == 6000 &&
  all(original$dataset_id %in% 1:1000, original$parameter %in% parameters) &&
  isTRUE(all(original$convergence == 0)) &&
  all(is.finite(original$estimate) & original$estimate > 0)

bootstrap_ok <- nrow(diagnostics) == 200000 &&
  n_distinct(diagnostics$dataset_id, diagnostics$bootstrap) == 200000 &&
  all(diagnostics$dataset_id %in% 1:1000, diagnostics$bootstrap %in% 1:200) &&
  isTRUE(all(diagnostics$convergence == 0)) &&
  isTRUE(all(diagnostics$two_car_convergence == 0)) &&
  isTRUE(all(diagnostics$one_car_convergence == 0))

intervals_ok <- nrow(intervals) == 6000 &&
  n_distinct(intervals$dataset_id, intervals$parameter) == 6000 &&
  all(intervals$dataset_id %in% 1:1000, intervals$parameter %in% parameters) &&
  isTRUE(all(intervals$n_bootstrap_successful == 200)) &&
  all(is.finite(intervals$lower) & is.finite(intervals$upper) &
        intervals$lower > 0 & intervals$lower <= intervals$upper)

cat("All runs converged and outputs complete:",
    original_ok && bootstrap_ok && intervals_ok, "\n")

# Build one row per simulated dataset and parameter, including missing intervals.
if (anyDuplicated(original[c("dataset_id", "parameter")]) > 0 ||
    anyDuplicated(intervals[c("dataset_id", "parameter")]) > 0) {
  stop("Duplicate dataset-parameter rows found.")
}

results <- expand_grid(dataset_id = 1:1000, parameter = parameters) %>%
  left_join(original %>% select(dataset_id, parameter, estimate, convergence),
            by = c("dataset_id", "parameter")) %>%
  left_join(intervals %>% select(dataset_id, parameter, lower, upper,
                                  n_bootstrap_successful),
            by = c("dataset_id", "parameter")) %>%
  mutate(
    true_value = unname(true_values[parameter]),
    valid_estimate = !is.na(convergence) & convergence == 0 &
      is.finite(estimate) & estimate > 0,
    valid_ci = valid_estimate & is.finite(lower) & is.finite(upper) &
      lower > 0 & lower <= upper &
      !is.na(n_bootstrap_successful) & n_bootstrap_successful >= 2,
    covers_truth = if_else(valid_ci, lower <= true_value & upper >= true_value, NA)
  )

# Summarize empirical coverage, Monte Carlo uncertainty, and estimation accuracy.
coverage_summary <- results %>%
  group_by(parameter) %>%
  summarize(
    true_value = first(true_value),
    n_simulations = n(),
    n_valid_ci = sum(valid_ci),
    n_covering = sum(covers_truth, na.rm = TRUE),
    coverage = if_else(n_valid_ci > 0, n_covering / n_valid_ci, NA_real_),
    mc_se = if_else(n_valid_ci > 0,
                    sqrt(coverage * (1 - coverage) / n_valid_ci), NA_real_),
    mean_bias = mean(if_else(valid_estimate, estimate - true_value, NA_real_),
                     na.rm = TRUE),
    mean_ci_width = mean(if_else(valid_ci, upper - lower, NA_real_),
                         na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(match(parameter, parameters))

print(coverage_summary, width = Inf)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write_csv(coverage_summary, file.path(output_dir, "simulation_coverage_summary.csv"))
