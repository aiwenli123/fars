################################################################################
# 14_better_figures.R
#
# Purpose: make cleaner figures, e.g., for papers or presentations
#
# Main input:
#   data/processed/script11_four_group_mle/
#     four_group_mle_estimates_pooled.csv
#
# Main outputs:
#   data/processed/script14_.../...
#
#   figures/script14_.../...
################################################################################

rm(list = ls())

library(tidyverse)
library(scales)

# ------------------------------------------------------------------------------
# 0. Paths and user options
# ------------------------------------------------------------------------------

script11_dir <- "data/processed/script11_four_group_mle"
input_path <- file.path(
  script11_dir,
  "four_group_mle_estimates_pooled.csv"
)

out_dir <- "data/processed/script14_four_group_figure"
fig_dir <- "figures/script14_four_group_figure"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

cell_specs_to_plot <- c("all_data", "rur_urb_weekend")

plot_width <- 11
plot_height <- 6.4
plot_dpi <- 400

cell_spec_offset <- c(
  "All data" = -0.11,
  "Rural/urban × weekend" = 0.11
)

# ------------------------------------------------------------------------------
# 1. Helper functions
# ------------------------------------------------------------------------------

assert_file_exists <- function(path) {
  if (!file.exists(path)) {
    stop(
      "Missing input file: ", path,
      "\nPlease run 11_four_group_mle.R first.",
      call. = FALSE
    )
  }
}

clean_cell_spec_label <- function(x) {
  recode(
    x,
    "all_data" = "All data",
    "rur_urb_weekend" = "Rural/urban × weekend",
    .default = str_replace_all(x, "_", " ") |> str_to_title()
  )
}

clean_stage_label <- function(x) {
  recode(
    x,
    "two_car" = "Two-car crashes",
    "one_car" = "One-car crashes",
    .default = x
  )
}

clean_group_label <- function(parameter) {
  case_when(
    str_detect(parameter, "_L$") ~ "Low BAC\n0 < BAC < 0.05",
    str_detect(parameter, "_M$") ~ "Middle BAC\n0.05 ≤ BAC < 0.08",
    str_detect(parameter, "_H$") ~ "High BAC\nBAC ≥ 0.08",
    TRUE ~ parameter
  )
}

# ------------------------------------------------------------------------------
# 2. Load and prepare pooled estimates
# ------------------------------------------------------------------------------

assert_file_exists(input_path)

pooled_raw <- read_csv(input_path, show_col_types = FALSE)

required_cols <- c(
  "cell_spec", "stage", "parameter",
  "estimate", "ci_low", "ci_high"
)

missing_cols <- setdiff(required_cols, names(pooled_raw))

if (length(missing_cols) > 0) {
  stop(
    "The pooled estimates file is missing required columns: ",
    paste(missing_cols, collapse = ", "),
    call. = FALSE
  )
}

bac_group_levels <- c(
  "Low BAC\n0 < BAC < 0.05",
  "Middle BAC\n0.05 ≤ BAC < 0.08",
  "High BAC\nBAC ≥ 0.08"
)

stage_levels <- c("Two-car crashes", "One-car crashes")

cell_spec_levels <- c(
  "All data",
  "Rural/urban × weekend"
)

plot_data <- pooled_raw %>%
  filter(
    cell_spec %in% cell_specs_to_plot,
    stage %in% c("two_car", "one_car"),
    parameter %in% c(
      "theta_L", "theta_M", "theta_H",
      "lambda_L", "lambda_M", "lambda_H"
    )
  ) %>%
  mutate(
    cell_spec_label = clean_cell_spec_label(cell_spec),
    stage_label = clean_stage_label(stage),
    bac_group = clean_group_label(parameter),
    cell_spec_label = factor(
      cell_spec_label,
      levels = cell_spec_levels
    ),
    stage_label = factor(
      stage_label,
      levels = stage_levels
    ),
    bac_group = factor(
      bac_group,
      levels = bac_group_levels
    ),
    bac_group_number = as.numeric(bac_group),
    y_position = bac_group_number +
      unname(cell_spec_offset[as.character(cell_spec_label)])
  ) %>%
  arrange(stage_label, bac_group, cell_spec_label)

expected_combinations <- crossing(
  cell_spec = cell_specs_to_plot,
  stage = c("two_car", "one_car"),
  group = c("L", "M", "H")
)

observed_combinations <- plot_data %>%
  transmute(
    cell_spec,
    stage,
    group = str_extract(parameter, "[LMH]$")
  ) %>%
  distinct()

missing_combinations <- anti_join(
  expected_combinations,
  observed_combinations,
  by = c("cell_spec", "stage", "group")
)

if (nrow(missing_combinations) > 0) {
  stop(
    "Some requested estimates are missing from the Script 11 output:\n",
    paste(
      apply(missing_combinations, 1, paste, collapse = " | "),
      collapse = "\n"
    ),
    call. = FALSE
  )
}

if (any(!is.finite(plot_data$estimate) | plot_data$estimate <= 0)) {
  stop(
    "All plotted risk estimates must be finite and greater than zero.",
    call. = FALSE
  )
}

if (any(is.na(plot_data$ci_low) | is.na(plot_data$ci_high))) {
  warning(
    "Some confidence intervals are missing. Point estimates will still be shown."
  )
}

write_csv(
  plot_data %>%
    select(
      cell_spec, cell_spec_label,
      stage, stage_label,
      parameter, bac_group,
      estimate, ci_low, ci_high,
      any_of(c(
        "m", "m_total", "m_est", "se_log_total",
        "mean_n_crashes"
      ))
    ),
  file.path(out_dir, "four_group_estimates_selected_specs.csv")
)

cat("Prepared four-group estimates:\n")
print(
  plot_data %>%
    select(
      cell_spec_label, stage_label, bac_group,
      estimate, ci_low, ci_high
    ),
  n = Inf
)
cat("\n")

# ------------------------------------------------------------------------------
# 3. Relative-risk figure
# ------------------------------------------------------------------------------

interval_data <- plot_data %>%
  filter(
    is.finite(ci_low),
    is.finite(ci_high),
    ci_low > 0,
    ci_high > 0
  )

fig14_01 <- ggplot() +
  geom_vline(
    xintercept = 1,
    linewidth = 0.55,
    linetype = "dashed",
    color = "grey45"
  ) +
  geom_segment(
    data = interval_data,
    aes(
      x = ci_low,
      xend = ci_high,
      y = y_position,
      yend = y_position,
      color = cell_spec_label
    ),
    linewidth = 0.85,
    lineend = "round"
  ) +
  geom_point(
    data = plot_data,
    aes(
      x = estimate,
      y = y_position,
      color = cell_spec_label,
      shape = cell_spec_label
    ),
    size = 3.4,
    stroke = 1
  ) +
  facet_wrap(
    ~ stage_label,
    nrow = 1
  ) +
  scale_x_log10(
    breaks = c(0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500),
    labels = label_number(accuracy = 0.1),
    expand = expansion(mult = c(0.04, 0.08))
  ) +
  scale_y_continuous(
    breaks = seq_along(bac_group_levels),
    labels = bac_group_levels,
    limits = c(0.65, length(bac_group_levels) + 0.35),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_color_manual(
    values = c(
      "All data" = "#4D4D4D",
      "Rural/urban × weekend" = "#0072B2"
    ),
    drop = FALSE
  ) +
  scale_shape_manual(
    values = c(
      "All data" = 16,
      "Rural/urban × weekend" = 17
    ),
    drop = FALSE
  ) +
  labs(
    title = "Relative fatal-crash risk compared with sober drivers",
    subtitle = "Pooled estimates across 10 imputations with 95% confidence intervals",
    x = "Relative fatal-crash risk",
    y = NULL,
    color = "Cell specification",
    shape = "Cell specification",
    caption = paste("Points show relative-risk estimates on a logarithmically scaled axis;",
    "the dashed line at 1 indicates risk equal to that of sober drivers.")
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        shape = c(16, 17),
        size = 3.2,
        linewidth = 0.9
      )
    ),
    shape = "none"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16,
      margin = margin(b = 5)
    ),
    plot.subtitle = element_text(
      size = 11.5,
      color = "grey30",
      margin = margin(b = 12)
    ),
    plot.caption = element_text(
      size = 9.5,
      color = "grey35",
      hjust = 0,
      margin = margin(t = 10)
    ),
    axis.title.x = element_text(
      face = "bold",
      margin = margin(t = 9)
    ),
    axis.text.x = element_text(
      color = "grey20"
    ),
    axis.text.y = element_text(
      color = "grey10",
      lineheight = 1.05,
      margin = margin(r = 8)
    ),
    strip.text = element_text(
      face = "bold",
      size = 12.5,
      margin = margin(b = 6)
    ),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(
      color = "grey88",
      linewidth = 0.35
    ),
    panel.spacing.x = unit(1.4, "lines"),
    legend.position = "top",
    legend.justification = "left",
    legend.box.just = "left",
    legend.title = element_text(face = "bold"),
    plot.margin = margin(12, 18, 10, 12)
  )

ggsave(
  filename = file.path(
    fig_dir,
    "fig14_01_four_group_relative_risk.png"
  ),
  plot = fig14_01,
  width = plot_width,
  height = plot_height,
  dpi = plot_dpi,
  bg = "white"
)

ggsave(
  filename = file.path(
    fig_dir,
    "fig14_01_four_group_relative_risk.pdf"
  ),
  plot = fig14_01,
  width = plot_width,
  height = plot_height,
  device = "pdf"
)

cat("Saved:\n")
cat(
  "  ",
  file.path(fig_dir, "fig14_01_four_group_relative_risk.png"),
  "\n",
  sep = ""
)
cat(
  "  ",
  file.path(fig_dir, "fig14_01_four_group_relative_risk.pdf"),
  "\n",
  sep = ""
)




# ------------------------------------------------------------------------------
# National fatal-crash alcohol involvement
# ------------------------------------------------------------------------------

library(tidyverse)
library(scales)

years_stats <- 2015:2024
raw_dir <- "data/raw"
accident_path <- "data/processed/accident_2015_2024_core.csv"

stats_out_dir <- "data/processed/script14_policy_motivation"
stats_fig_dir <- "figures/script14_policy_motivation"

dir.create(stats_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(stats_fig_dir, recursive = TRUE, showWarnings = FALSE)

find_midrvacc_path <- function(year, raw_dir = "data/raw") {
  year_dir <- file.path(raw_dir, paste0("FARS", year, "NationalCSV"))
  
  if (!dir.exists(year_dir)) {
    stop("Could not find FARS directory for ", year, ": ", year_dir)
  }
  
  paths <- list.files(
    year_dir,
    pattern = "^midrvacc\\.csv$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  if (length(paths) == 0) {
    stop("Could not find MIDRVACC.csv for ", year, " in ", year_dir)
  }
  
  if (length(paths) > 1) {
    warning("Multiple MIDRVACC files found for ", year, "; using ", paths[[1]])
  }
  
  paths[[1]]
}

clean_midrvacc_bac <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  
  case_when(
    is.na(x) ~ NA_real_,
    x >= 0 & x <= 0.94 ~ x,       # already on the decimal BAC scale
    x >= 0 & x <= 94 ~ x / 100,   # MIDRVACC hundredths scale: 8 = 0.08
    TRUE ~ NA_real_                # missing/unknown codes such as 95-99
  )
}

read_midrvacc_year <- function(year) {
  path <- find_midrvacc_path(year, raw_dir)
  
  dat <- read_csv(path, show_col_types = FALSE) %>%
    rename_with(toupper)
  
  a_cols <- paste0("A", 1:10)
  required_cols <- c("ST_CASE", a_cols)
  missing_cols <- setdiff(required_cols, names(dat))
  
  if (length(missing_cols) > 0) {
    stop(
      "MIDRVACC file for ", year, " is missing: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  dat %>%
    select(ST_CASE, all_of(a_cols)) %>%
    mutate(
      YEAR_FILE = as.integer(year),
      ST_CASE = as.integer(ST_CASE)
    ) %>%
    pivot_longer(
      cols = all_of(a_cols),
      names_to = "IMPUTATION",
      names_pattern = "A(\\d+)",
      values_to = "MAX_DRIVER_BAC_RAW"
    ) %>%
    mutate(
      IMPUTATION = as.integer(IMPUTATION),
      MAX_DRIVER_BAC = clean_midrvacc_bac(MAX_DRIVER_BAC_RAW)
    ) %>%
    select(YEAR_FILE, ST_CASE, IMPUTATION, MAX_DRIVER_BAC)
}

# One row per fatal crash from the full FARS accident files.
accidents <- read_csv(accident_path, show_col_types = FALSE) %>%
  rename_with(toupper) %>%
  transmute(
    YEAR_FILE = as.integer(YEAR_FILE),
    ST_CASE = as.integer(ST_CASE)
  ) %>%
  filter(YEAR_FILE %in% years_stats) %>%
  distinct()

# One row per fatal crash per imputation.
midrvacc_long <- map_dfr(years_stats, read_midrvacc_year)

duplicate_keys <- midrvacc_long %>%
  count(YEAR_FILE, ST_CASE, IMPUTATION) %>%
  filter(n > 1)

if (nrow(duplicate_keys) > 0) {
  stop("Duplicate crash-imputation rows were found in MIDRVACC.")
}

# Construct the full crash x imputation grid so the denominator includes
# every fatal crash in every imputation.
crash_bac <- crossing(
  accidents,
  IMPUTATION = 1:10
) %>%
  left_join(
    midrvacc_long,
    by = c("YEAR_FILE", "ST_CASE", "IMPUTATION")
  )

coverage_summary <- crash_bac %>%
  group_by(YEAR_FILE) %>%
  summarise(
    n_fatal_crashes = n_distinct(ST_CASE),
    n_crash_imputation_rows = n(),
    n_missing_max_driver_bac = sum(is.na(MAX_DRIVER_BAC)),
    pct_missing_max_driver_bac = mean(is.na(MAX_DRIVER_BAC)),
    .groups = "drop"
  )

if (any(coverage_summary$n_missing_max_driver_bac > 0)) {
  warning(
    "Some crash-imputation rows have missing maximum driver BAC. ",
    "See midrvacc_coverage_by_year.csv."
  )
}

# Calculate statistics separately in each imputed dataset.
by_year_imp <- crash_bac %>%
  group_by(YEAR_FILE, IMPUTATION) %>%
  summarise(
    n_fatal_crashes = n(),
    n_bac_available = sum(!is.na(MAX_DRIVER_BAC)),
    n_any_positive = sum(MAX_DRIVER_BAC > 0, na.rm = TRUE),
    n_max_bac_05_079 = sum(
      MAX_DRIVER_BAC >= 0.05 & MAX_DRIVER_BAC < 0.08,
      na.rm = TRUE
    ),
    n_bac_05plus = sum(MAX_DRIVER_BAC >= 0.05, na.rm = TRUE),
    n_bac_08plus = sum(MAX_DRIVER_BAC >= 0.08, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    pct_bac_available = n_bac_available / n_fatal_crashes,
    pct_any_positive = n_any_positive / n_fatal_crashes,
    pct_max_bac_05_079 = n_max_bac_05_079 / n_fatal_crashes,
    pct_bac_05plus = n_bac_05plus / n_fatal_crashes,
    pct_bac_08plus = n_bac_08plus / n_fatal_crashes
  )

# Average the descriptive quantities over the 10 imputations.
by_year <- by_year_imp %>%
  group_by(YEAR_FILE) %>%
  summarise(
    n_fatal_crashes = first(n_fatal_crashes),
    avg_n_any_positive = mean(n_any_positive),
    avg_n_max_bac_05_079 = mean(n_max_bac_05_079),
    avg_n_bac_05plus = mean(n_bac_05plus),
    avg_n_bac_08plus = mean(n_bac_08plus),
    pct_any_positive = mean(pct_any_positive),
    pct_max_bac_05_079 = mean(pct_max_bac_05_079),
    pct_bac_05plus = mean(pct_bac_05plus),
    pct_bac_08plus = mean(pct_bac_08plus),
    .groups = "drop"
  )

# Pool all years before calculating percentages, so this is weighted by the
# number of fatal crashes rather than giving every year equal weight.
pooled_imp <- crash_bac %>%
  group_by(IMPUTATION) %>%
  summarise(
    n_fatal_crashes = n(),
    n_any_positive = sum(MAX_DRIVER_BAC > 0, na.rm = TRUE),
    n_max_bac_05_079 = sum(
      MAX_DRIVER_BAC >= 0.05 & MAX_DRIVER_BAC < 0.08,
      na.rm = TRUE
    ),
    n_bac_05plus = sum(MAX_DRIVER_BAC >= 0.05, na.rm = TRUE),
    n_bac_08plus = sum(MAX_DRIVER_BAC >= 0.08, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    pct_any_positive = n_any_positive / n_fatal_crashes,
    pct_max_bac_05_079 = n_max_bac_05_079 / n_fatal_crashes,
    pct_bac_05plus = n_bac_05plus / n_fatal_crashes,
    pct_bac_08plus = n_bac_08plus / n_fatal_crashes
  )

pooled <- pooled_imp %>%
  summarise(
    first_year = min(years_stats),
    last_year = max(years_stats),
    n_fatal_crashes = first(n_fatal_crashes),
    avg_n_any_positive = mean(n_any_positive),
    avg_n_max_bac_05_079 = mean(n_max_bac_05_079),
    avg_n_bac_05plus = mean(n_bac_05plus),
    avg_n_bac_08plus = mean(n_bac_08plus),
    pct_any_positive = mean(pct_any_positive),
    pct_max_bac_05_079 = mean(pct_max_bac_05_079),
    pct_bac_05plus = mean(pct_bac_05plus),
    pct_bac_08plus = mean(pct_bac_08plus)
  )

write_csv(
  coverage_summary,
  file.path(stats_out_dir, "midrvacc_coverage_by_year.csv")
)

write_csv(
  by_year_imp,
  file.path(stats_out_dir, "fatal_crash_alcohol_shares_by_year_imputation.csv")
)

write_csv(
  by_year,
  file.path(stats_out_dir, "fatal_crash_alcohol_shares_by_year.csv")
)

write_csv(
  pooled,
  file.path(stats_out_dir, "fatal_crash_alcohol_shares_2015_2024.csv")
)

pooled_display <- pooled %>%
  transmute(
    period = paste0(first_year, "-", last_year),
    total_fatal_crashes = comma(n_fatal_crashes),
    any_positive_bac = percent(pct_any_positive, accuracy = 0.1),
    maximum_bac_05_to_079 = percent(pct_max_bac_05_079, accuracy = 0.1),
    at_least_one_driver_bac_05plus = percent(pct_bac_05plus, accuracy = 0.1),
    at_least_one_driver_bac_08plus = percent(pct_bac_08plus, accuracy = 0.1)
  )

print(pooled_display, width = Inf)

plot_data <- by_year %>%
  transmute(
    year = YEAR_FILE,
    `Any positive driver BAC` = pct_any_positive,
    `Maximum driver BAC 0.05-0.079` = pct_max_bac_05_079,
    `Driver BAC at least 0.05` = pct_bac_05plus,
    `Driver BAC at least 0.08` = pct_bac_08plus
  ) %>%
  pivot_longer(
    cols = -year,
    names_to = "measure",
    values_to = "share"
  )

fatal_crash_share_plot <- ggplot(
  plot_data,
  aes(x = year, y = share, group = measure, color = measure)
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  scale_color_manual(
    values = c(
      "Any positive driver BAC" = "#0072B2",
      "Driver BAC at least 0.05" = "#E69F00",
      "Driver BAC at least 0.08" = "#009E73",
      "Maximum driver BAC 0.05-0.079" = "#CC79A7"
    ),
    breaks = c(
      "Any positive driver BAC",
      "Driver BAC at least 0.05",
      "Driver BAC at least 0.08",
      "Maximum driver BAC 0.05-0.079"
    )
  ) +
  scale_x_continuous(breaks = years_stats) +
  scale_y_continuous(labels = label_percent(accuracy = 1)) +
  labs(
    title = "Alcohol involvement in U.S. fatal crashes",
    subtitle = "Maximum driver BAC from FARS multiple-imputation data",
    x = NULL,
    y = "Share of fatal crashes",
    color = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(stats_fig_dir, "fatal_crash_alcohol_shares_by_year.png"),
  fatal_crash_share_plot,
  width = 10,
  height = 6,
  dpi = 300
)
