################################################################################
# 12_four_group_results.R
#
# Purpose:
#   Visualize the four-BAC-group MLE results from script 11 and compare the
#   preferred four-group specification (rural/urban x weekend) with the binary
#   2-BAC-group LP/MI model results from script 07.
#
# Main inputs:
#   data/processed/script11_four_group_mle/four_group_mle_estimates_pooled.csv
#   data/processed/script11_four_group_mle/four_group_mle_estimates_by_imp.csv
#   data/processed/script11_four_group_mle/four_group_mle_fit_diagnostics.csv
#   data/processed/script07_mi_lp_mle/location_weekend/lp_mle_raw_vs_mi_comparison.csv
#
# Main outputs:
#   data/processed/script12_four_group_visuals/
#     script12_four_group_estimates_for_plot.csv
#     script12_binary_vs_four_group_preferred.csv
#     script12_binary_vs_four_group_preferred_display.csv
#
#   figures/script12_four_group_visuals/
#     fig12_01_four_group_pooled_estimates.png
#     fig12_02_four_group_by_imputation.png
#     fig12_03_four_group_preferred_estimates.png
#     fig12_04_binary_vs_four_group_preferred.png
#     fig12_05_binary_vs_four_group_ratio.png
#     fig12_06_four_group_fit_diagnostics.png
################################################################################

rm(list = ls())

library(tidyverse)
library(scales)

# ------------------------------------------------------------------------------
# 0. User options
# ------------------------------------------------------------------------------

script11_dir <- "data/processed/script11_four_group_mle"
script07_dir <- "data/processed/script07_mi_lp_mle/location_weekend"

out_dir <- "data/processed/script12_four_group_visuals"
fig_dir <- "figures/script12_four_group_visuals"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Preferred specification for the main comparison requested.
preferred_cell_spec <- "rur_urb_weekend"

# Script 07 comparison choice. For comparability with script 11, the default is
# MI BAC, all state-years, no n2 filter, and the same coarse cell specification.
preferred_script07_analysis_sample <- "mi_all_state_years"
preferred_min_n2 <- 1
thresholds_to_compare <- c(0.05, 0.08)

# Plot dimensions.
plot_width <- 11
plot_height <- 6.5
plot_dpi <- 300

# ------------------------------------------------------------------------------
# 1. Helper functions
# ------------------------------------------------------------------------------

assert_file_exists <- function(path, hint = NULL) {
  if (!file.exists(path)) {
    msg <- paste0("Missing file: ", path)
    if (!is.null(hint)) msg <- paste0(msg, "\n", hint)
    stop(msg, call. = FALSE)
  }
}

clean_cell_spec_label <- function(x) {
  recode(
    x,
    "all_data" = "All data",
    "weekend" = "Weekend",
    "rur_urb" = "Rural/urban",
    "rur_urb_weekend" = "Rural/urban x weekend",
    .default = str_replace_all(x, "_", " ") |> str_to_title()
  )
}

clean_stage_label <- function(x) {
  recode(
    x,
    "two_car" = "Two-car risk (theta)",
    "one_car" = "One-car risk (lambda)",
    .default = x
  )
}

clean_group_label <- function(parameter) {
  case_when(
    str_detect(parameter, "_L$") ~ "L: 0 < BAC < 0.05",
    str_detect(parameter, "_M$") ~ "M: 0.05 <= BAC < 0.08",
    str_detect(parameter, "_H$") ~ "H: BAC >= 0.08",
    TRUE ~ parameter
  )
}

clean_binary_threshold_label <- function(threshold) {
  case_when(
    threshold == 0.05 ~ "2-group: BAC >= 0.05\n(M/H combined vs S/L)",
    threshold == 0.08 ~ "2-group: BAC >= 0.08\n(H vs S/L/M)",
    TRUE ~ paste0("2-group: BAC >= ", threshold)
  )
}

make_ratio_ci <- function(est_four, low_four, high_four, est_binary, low_binary, high_binary) {
  # Approximate CI for a ratio using the endpoints on the log scale. This is only
  # used for plotting/diagnostics, not for re-estimating the models.
  tibble(
    ratio = est_four / est_binary,
    ratio_low = low_four / high_binary,
    ratio_high = high_four / low_binary
  )
}

save_plot <- function(plot, filename, width = plot_width, height = plot_height) {
  ggsave(
    filename = file.path(fig_dir, filename),
    plot = plot,
    width = width,
    height = height,
    dpi = plot_dpi
  )
}

# ------------------------------------------------------------------------------
# 2. Load script 11 four-group outputs
# ------------------------------------------------------------------------------

four_pooled_path <- file.path(script11_dir, "four_group_mle_estimates_pooled.csv")
four_by_imp_path <- file.path(script11_dir, "four_group_mle_estimates_by_imp.csv")
four_diag_path <- file.path(script11_dir, "four_group_mle_fit_diagnostics.csv")

assert_file_exists(
  four_pooled_path,
  "Please run 11_four_group_mle.R first."
)
assert_file_exists(
  four_by_imp_path,
  "Please run 11_four_group_mle.R first."
)
assert_file_exists(
  four_diag_path,
  "Please run 11_four_group_mle.R first."
)

four_pooled_raw <- read_csv(four_pooled_path, show_col_types = FALSE)
four_by_imp_raw <- read_csv(four_by_imp_path, show_col_types = FALSE)
four_diag_raw <- read_csv(four_diag_path, show_col_types = FALSE)

cell_spec_order <- c("all_data", "weekend", "rur_urb", "rur_urb_weekend")
group_order <- c("L: 0 < BAC < 0.05", "M: 0.05 <= BAC < 0.08", "H: BAC >= 0.08")
stage_order <- c("Two-car risk (theta)", "One-car risk (lambda)")

four_pooled <- four_pooled_raw %>%
  mutate(
    cell_spec = factor(cell_spec, levels = cell_spec_order),
    cell_spec_label = factor(clean_cell_spec_label(as.character(cell_spec)),
                             levels = clean_cell_spec_label(cell_spec_order)),
    stage_label = factor(clean_stage_label(stage), levels = stage_order),
    bac_group = factor(clean_group_label(parameter), levels = group_order),
    param_family = if_else(stage == "two_car", "theta", "lambda"),
    ci_low_plot = if_else(is.na(ci_low), estimate, ci_low),
    ci_high_plot = if_else(is.na(ci_high), estimate, ci_high)
  )

four_by_imp <- four_by_imp_raw %>%
  mutate(
    cell_spec = factor(cell_spec, levels = cell_spec_order),
    cell_spec_label = factor(clean_cell_spec_label(as.character(cell_spec)),
                             levels = clean_cell_spec_label(cell_spec_order)),
    stage_label = factor(clean_stage_label(stage), levels = stage_order),
    bac_group = factor(clean_group_label(parameter), levels = group_order),
    param_family = if_else(stage == "two_car", "theta", "lambda")
  )

write_csv(four_pooled, file.path(out_dir, "script12_four_group_estimates_for_plot.csv"))

cat("Loaded script 11 results:\n")
cat("  Pooled rows:", nrow(four_pooled), "\n")
cat("  Imputation-level rows:", nrow(four_by_imp), "\n")
cat("  Cell specs:", paste(levels(four_pooled$cell_spec_label), collapse = ", "), "\n\n")

# ------------------------------------------------------------------------------
# 3. Figure 1: pooled four-group estimates across all four cell specifications
# ------------------------------------------------------------------------------

fig12_01 <- ggplot(
  four_pooled,
  aes(x = bac_group, y = estimate, ymin = ci_low_plot, ymax = ci_high_plot)
) +
  geom_hline(yintercept = 1, linewidth = 0.3, linetype = "dashed") +
  geom_pointrange(position = position_dodge(width = 0.35)) +
  facet_grid(stage_label ~ cell_spec_label) +
  scale_y_log10(labels = label_number(accuracy = 0.1)) +
  labs(
    title = "Four-BAC-group MLE estimates across cell specifications",
    subtitle = "Pooled over imputations on the log scale; vertical intervals are Rubin-style 95% CIs when available",
    x = NULL,
    y = "Relative fatal crash risk estimate (log scale)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    panel.grid.minor = element_blank()
  )

save_plot(fig12_01, "fig12_01_four_group_pooled_estimates.png", width = 13, height = 7)

# ------------------------------------------------------------------------------
# 4. Figure 2: imputation-level variability around pooled estimates
# ------------------------------------------------------------------------------

fig12_02 <- ggplot(four_by_imp, aes(x = factor(IMPUTATION), y = estimate)) +
  geom_hline(yintercept = 1, linewidth = 0.3, linetype = "dashed") +
  geom_point(alpha = 0.75) +
  geom_hline(
    data = four_pooled,
    aes(yintercept = estimate),
    linewidth = 0.35
  ) +
  facet_grid(stage_label + bac_group ~ cell_spec_label, scales = "free_y") +
  scale_y_log10(labels = label_number(accuracy = 0.1)) +
  labs(
    title = "Four-group estimates by imputation",
    subtitle = "Points are imputation-specific estimates; horizontal line is pooled estimate",
    x = "Imputation",
    y = "Relative fatal crash risk estimate (log scale)"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    strip.text.y = element_text(angle = 0)
  )

save_plot(fig12_02, "fig12_02_four_group_by_imputation.png", width = 13, height = 11)

# ------------------------------------------------------------------------------
# 5. Figure 3: preferred four-group specification only
# ------------------------------------------------------------------------------

preferred_four <- four_pooled %>%
  filter(as.character(cell_spec) == preferred_cell_spec)

if (nrow(preferred_four) == 0) {
  stop("No script 11 rows found for preferred_cell_spec = ", preferred_cell_spec, call. = FALSE)
}

fig12_03 <- ggplot(
  preferred_four,
  aes(x = bac_group, y = estimate, ymin = ci_low_plot, ymax = ci_high_plot)
) +
  geom_hline(yintercept = 1, linewidth = 0.3, linetype = "dashed") +
  geom_pointrange() +
  facet_wrap(~ stage_label, nrow = 1) +
  scale_y_log10(labels = label_number(accuracy = 0.1)) +
  labs(
    title = "Preferred four-group specification: rural/urban x weekend",
    subtitle = "This is the specification used for the binary-vs-four-group comparison below",
    x = NULL,
    y = "Relative fatal crash risk estimate (log scale)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid.minor = element_blank()
  )

save_plot(fig12_03, "fig12_03_four_group_preferred_estimates.png", width = 10, height = 5.5)

# ------------------------------------------------------------------------------
# 6. Load script 07 binary LP/MI comparison output
# ------------------------------------------------------------------------------

script07_compare_path <- file.path(script07_dir, "lp_mle_raw_vs_mi_comparison.csv")

assert_file_exists(
  script07_compare_path,
  paste0(
    "Please run 07_mi_lp_mle.R with location_spec = 'LOCATION_WEEKEND'. ",
    "Expected file under: ", script07_dir
  )
)

binary_raw <- read_csv(script07_compare_path, show_col_types = FALSE)

# Make this robust to slightly different analysis_sample naming from earlier script
# versions. The current script 07 uses mi_all_state_years.
available_analysis_samples <- sort(unique(binary_raw$analysis_sample))
if (!preferred_script07_analysis_sample %in% available_analysis_samples) {
  stop(
    "Could not find analysis_sample = '", preferred_script07_analysis_sample, "' in script 07 output.\n",
    "Available analysis_sample values are:\n  ",
    paste(available_analysis_samples, collapse = "\n  "),
    call. = FALSE
  )
}

binary_pref <- binary_raw %>%
  filter(
    analysis_sample == preferred_script07_analysis_sample,
    threshold %in% thresholds_to_compare,
    cell_spec == preferred_cell_spec,
    min_n2 == preferred_min_n2,
    parameter %in% c("theta", "lambda")
  ) %>%
  mutate(
    stage = if_else(parameter == "theta", "two_car", "one_car"),
    stage_label = factor(clean_stage_label(stage), levels = stage_order),
    binary_model_label = clean_binary_threshold_label(threshold),
    binary_model_label = factor(binary_model_label, levels = clean_binary_threshold_label(thresholds_to_compare)),
    ci_low_plot = if_else(is.na(ci_low), estimate, ci_low),
    ci_high_plot = if_else(is.na(ci_high), estimate, ci_high)
  )

if (nrow(binary_pref) == 0) {
  stop(
    "No matching script 07 rows found for the preferred comparison.\n",
    "Check preferred_cell_spec, preferred_min_n2, thresholds_to_compare, and analysis_sample.",
    call. = FALSE
  )
}

# ------------------------------------------------------------------------------
# 7. Build binary-vs-four-group comparison table
# ------------------------------------------------------------------------------

four_for_compare <- preferred_four %>%
  mutate(
    threshold = case_when(
      str_detect(as.character(bac_group), "0.05 <=") ~ 0.05,
      str_detect(as.character(bac_group), "BAC >= 0.08") ~ 0.08,
      TRUE ~ NA_real_
    ),
    four_group_comparison_label = case_when(
      threshold == 0.05 ~ "4-group: M group only\n0.05 <= BAC < 0.08",
      threshold == 0.08 ~ "4-group: H group only\nBAC >= 0.08",
      TRUE ~ as.character(bac_group)
    )
  ) %>%
  filter(!is.na(threshold)) %>%
  select(
    threshold,
    stage,
    stage_label,
    parameter_four = parameter,
    four_group_comparison_label,
    four_estimate = estimate,
    four_ci_low = ci_low,
    four_ci_high = ci_high,
    four_mean_n_crashes = mean_n_crashes
  )

binary_vs_four <- binary_pref %>%
  select(
    threshold,
    threshold_label,
    analysis_sample,
    analysis_label,
    cell_spec,
    cell_spec_label,
    min_n2,
    min_n2_label,
    stage,
    stage_label,
    parameter_binary = parameter,
    binary_model_label,
    binary_estimate = estimate,
    binary_ci_low = ci_low,
    binary_ci_high = ci_high,
    binary_status = status,
    binary_m_est = m_est,
    binary_m_total = m_total
  ) %>%
  inner_join(four_for_compare, by = c("threshold", "stage", "stage_label")) %>%
  bind_cols(
    make_ratio_ci(
      est_four = .$four_estimate,
      low_four = .$four_ci_low,
      high_four = .$four_ci_high,
      est_binary = .$binary_estimate,
      low_binary = .$binary_ci_low,
      high_binary = .$binary_ci_high
    )
  ) %>%
  mutate(
    stage_label = factor(stage_label, levels = stage_order),
    threshold_label_short = paste0("BAC >= ", number(threshold, accuracy = 0.01)),
    comparison_note = case_when(
      threshold == 0.05 ~ paste(
        "The binary threshold model estimates one combined risk for BAC >= 0.05,",
        "while the four-group model separates the medium group from the high group."
      ),
      threshold == 0.08 ~ paste(
        "The binary threshold model's treated group corresponds most directly to H in the four-group model,",
        "but the reference group differs because the four-group model keeps L and M separate from S."
      ),
      TRUE ~ NA_character_
    )
  )

binary_vs_four_display <- binary_vs_four %>%
  transmute(
    stage = as.character(stage_label),
    threshold = threshold_label_short,
    binary_model = as.character(binary_model_label),
    four_group_model = four_group_comparison_label,
    binary_estimate = round(binary_estimate, 3),
    binary_ci = paste0("(", round(binary_ci_low, 3), ", ", round(binary_ci_high, 3), ")"),
    four_group_estimate = round(four_estimate, 3),
    four_group_ci = paste0("(", round(four_ci_low, 3), ", ", round(four_ci_high, 3), ")"),
    four_div_binary = round(ratio, 3),
    comparison_note
  )

write_csv(binary_vs_four, file.path(out_dir, "script12_binary_vs_four_group_preferred.csv"))
write_csv(binary_vs_four_display, file.path(out_dir, "script12_binary_vs_four_group_preferred_display.csv"))

cat("Script 07 comparison rows used:\n")
cat("  Rows:", nrow(binary_pref), "\n")
cat("  Analysis sample:", preferred_script07_analysis_sample, "\n")
cat("  Cell spec:", preferred_cell_spec, "\n")
cat("  min_n2:", preferred_min_n2, "\n\n")

# ------------------------------------------------------------------------------
# 8. Figure 4: direct estimate comparison, binary vs four-group
# ------------------------------------------------------------------------------

comparison_plot_dat <- bind_rows(
  binary_vs_four %>%
    transmute(
      stage_label,
      threshold_label_short,
      model = as.character(binary_model_label),
      estimate = binary_estimate,
      ci_low = binary_ci_low,
      ci_high = binary_ci_high
    ),
  binary_vs_four %>%
    transmute(
      stage_label,
      threshold_label_short,
      model = four_group_comparison_label,
      estimate = four_estimate,
      ci_low = four_ci_low,
      ci_high = four_ci_high
    )
) %>%
  mutate(
    model = factor(model, levels = unique(model)),
    ci_low_plot = if_else(is.na(ci_low), estimate, ci_low),
    ci_high_plot = if_else(is.na(ci_high), estimate, ci_high)
  )

fig12_04 <- ggplot(
  comparison_plot_dat,
  aes(x = model, y = estimate, ymin = ci_low_plot, ymax = ci_high_plot)
) +
  geom_hline(yintercept = 1, linewidth = 0.3, linetype = "dashed") +
  geom_pointrange() +
  facet_grid(stage_label ~ threshold_label_short, scales = "free_x", space = "free_x") +
  scale_y_log10(labels = label_number(accuracy = 0.1)) +
  labs(
    title = "Preferred specification: 2-BAC-group model vs four-BAC-group model",
    subtitle = "Binary model from script 07; four-group model from script 11; both use rural/urban x weekend cells",
    x = NULL,
    y = "Relative fatal crash risk estimate (log scale)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    panel.grid.minor = element_blank()
  )

save_plot(fig12_04, "fig12_04_binary_vs_four_group_preferred.png", width = 12, height = 7)

# ------------------------------------------------------------------------------
# 9. Figure 5: ratio of four-group estimate to binary estimate
# ------------------------------------------------------------------------------

ratio_plot_dat <- binary_vs_four %>%
  mutate(
    four_group_comparison_label = factor(
      four_group_comparison_label,
      levels = unique(four_group_comparison_label)
    ),
    ratio_low_plot = if_else(is.na(ratio_low) | !is.finite(ratio_low), ratio, ratio_low),
    ratio_high_plot = if_else(is.na(ratio_high) | !is.finite(ratio_high), ratio, ratio_high)
  )

fig12_05 <- ggplot(
  ratio_plot_dat,
  aes(x = four_group_comparison_label, y = ratio, ymin = ratio_low_plot, ymax = ratio_high_plot)
) +
  geom_hline(yintercept = 1, linewidth = 0.3, linetype = "dashed") +
  geom_pointrange() +
  facet_grid(stage_label ~ threshold_label_short, scales = "free_x", space = "free_x") +
  scale_y_log10(labels = label_number(accuracy = 0.1)) +
  labs(
    title = "How different is the four-group estimate from the binary threshold estimate?",
    subtitle = "Values above 1 mean the four-group estimate is larger than the corresponding binary threshold estimate",
    x = NULL,
    y = "Four-group estimate / binary estimate (log scale)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    panel.grid.minor = element_blank()
  )

save_plot(fig12_05, "fig12_05_binary_vs_four_group_ratio.png", width = 11, height = 6.5)

# ------------------------------------------------------------------------------
# 10. Figure 6: fit diagnostics across specifications
# ------------------------------------------------------------------------------

fit_diag_summary <- four_diag_raw %>%
  mutate(
    cell_spec = factor(cell_spec, levels = cell_spec_order),
    cell_spec_label = factor(clean_cell_spec_label(as.character(cell_spec)),
                             levels = clean_cell_spec_label(cell_spec_order)),
    stage_label = factor(clean_stage_label(stage), levels = stage_order),
    convergence_ok = optim_convergence == 0,
    hessian_ok = if_else(is.na(hessian_ok), FALSE, hessian_ok)
  ) %>%
  group_by(cell_spec_label, stage_label) %>%
  summarize(
    n_fits = n(),
    convergence_ok_n = sum(convergence_ok, na.rm = TRUE),
    hessian_ok_n = sum(hessian_ok, na.rm = TRUE),
    convergence_ok_share = convergence_ok_n / n_fits,
    hessian_ok_share = hessian_ok_n / n_fits,
    mean_n_cells = mean(n_cells, na.rm = TRUE),
    mean_n_crashes = mean(n_crashes, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(convergence_ok_share, hessian_ok_share),
    names_to = "diagnostic",
    values_to = "share"
  ) %>%
  mutate(
    diagnostic = recode(
      diagnostic,
      "convergence_ok_share" = "Optimizer converged",
      "hessian_ok_share" = "Hessian SE available"
    )
  )

fig12_06 <- ggplot(
  fit_diag_summary,
  aes(x = cell_spec_label, y = share)
) +
  geom_col() +
  facet_grid(stage_label ~ diagnostic) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = "Script 11 fit diagnostics across specifications",
    subtitle = "Each bar summarizes the share of imputation-level fits satisfying the diagnostic",
    x = NULL,
    y = "Share of fits"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    panel.grid.minor = element_blank()
  )

save_plot(fig12_06, "fig12_06_four_group_fit_diagnostics.png", width = 11, height = 6.5)

write_csv(fit_diag_summary, file.path(out_dir, "script12_four_group_fit_diagnostics_summary.csv"))
