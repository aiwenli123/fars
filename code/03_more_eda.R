################################################################################
# 03_more_eda.R
#
# Purpose:
#   Create cleaner descriptive tables and plots for the
#   FARS / Levitt-Porter replication project.
#
#
# Main inputs (from scripts 01 and 02):
#   data/processed/driver_analysis_lp_2015_2024.csv
#   data/processed/lp_crash_types_2015_2024.csv
#   data/processed/lp_theta_simple_2015_2024.csv
#   data/processed/lp_one_vehicle_counts_2015_2024.csv
#
# Main outputs:
#   CSV tables:
#     data/processed/script03_tables/table_01_sample_construction.csv
#     data/processed/script03_tables/table_02_missingness_by_definition_period.csv
#     data/processed/script03_tables/table_03_two_vehicle_composition_all_crashes.csv
#     data/processed/script03_tables/table_04_two_vehicle_composition_classifiable_only.csv
#     data/processed/script03_tables/table_05_theta_by_definition_period.csv
#     data/processed/script03_tables/table_06_one_vehicle_composition.csv
#
#   Plots:
#     figures/script03_plots/plot_A_two_vehicle_shares_with_missing.png
#     figures/script03_plots/plot_B_two_vehicle_complete_case_shares.png
#     figures/script03_plots/plot_C_theta_with_denominators.png
#     figures/script03_plots/plot_D_missingness_over_time_by_year.png
################################################################################

rm(list = ls())

# ------------------------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------------------------

library(tidyverse)
library(scales)
# install.packages("ggrepel") # run once if needed
library(ggrepel)

# ------------------------------------------------------------------------------
# 1. File paths
# ------------------------------------------------------------------------------

driver_path <- "data/processed/driver_analysis_lp_2015_2024.csv"

crash_types_path <- "data/processed/lp_crash_types_2015_2024.csv"
theta_path       <- "data/processed/lp_theta_simple_2015_2024.csv"
one_vehicle_path <- "data/processed/lp_one_vehicle_counts_2015_2024.csv"

table_dir <- "data/processed/script03_tables"
plot_dir <- "figures/script03_plots"

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 2. Helper functions
# ------------------------------------------------------------------------------

check_file_exists <- function(path) {
  if (!file.exists(path)) {
    stop(
      "Missing file: ", path, "\n",
      "Please run 02_crash_type_counts.R before running this script."
    )
  }
}

period_levels <- c(
  "2015-2024 all years",
  "2015-2019 pre-COVID",
  "2020-2021 COVID",
  "2022-2024 post-COVID"
)

definition_levels <- c(
  "Police drinking, complete case",
  "Police drinking, unknown as no",
  "Vehicle DR_DRINK, complete case",
  "Observed BAC >= 0.08",
  "Observed BAC >= 0.10"
)

# ------------------------------------------------------------------------------
# 3. Load data
# ------------------------------------------------------------------------------

walk(
  c(driver_path, crash_types_path, theta_path, one_vehicle_path),
  check_file_exists
)

driver_lp <- read_csv(driver_path, show_col_types = FALSE) %>%
  rename_with(toupper)

# PER_NO is the person number within the vehicle/crash record.
# It is not the driver/passenger indicator.
# We identify drivers using PER_TYP == 1, not PER_NO == 1.
# Therefore, rows with PER_NO > 1 can still be valid driver rows.

# driver_lp %>%
#   filter(PER_NO > 1) %>%
#   count(PER_TYP, PER_TYPNAME) # should be all PER_TYP == 1 and PER_TYPNAME is driver.

crash_types <- read_csv(crash_types_path, show_col_types = FALSE) # 5 rows per accident, one for each alcohol definition

theta_simple <- read_csv(theta_path, show_col_types = FALSE) # estimates from script 02

one_vehicle_counts <- read_csv(one_vehicle_path, show_col_types = FALSE)

# ------------------------------------------------------------------------------
# 4. Recreate pooled + COVID period version of crash_types
# ------------------------------------------------------------------------------

crash_types_for_summary <- crash_types %>%
  mutate(
    all_years_group = "2015-2024 all years"
  ) %>%
  pivot_longer(
    cols = c(all_years_group, covid_period),
    names_to = "summary_level",
    values_to = "period"
  ) %>%
  filter(!is.na(period)) %>%
  mutate(
    summary_level = case_when(
      summary_level == "all_years_group" ~ "All years pooled",
      summary_level == "covid_period" ~ "COVID period",
      TRUE ~ summary_level
    ),
    period = factor(period, levels = period_levels),
    alcohol_definition_label = factor(
      alcohol_definition_label,
      levels = definition_levels
    )
  )

# ------------------------------------------------------------------------------
# 5. Table 1: Sample construction
# ------------------------------------------------------------------------------

driver_rows_by_crash <- driver_lp %>%
  group_by(YEAR_FILE, ST_CASE) %>%
  summarize(
    YEAR = first(YEAR),
    VE_FORMS = first(VE_FORMS),
    n_driver_rows = n(),
    .groups = "drop"
  )

clean_crashes <- driver_rows_by_crash %>%
  filter(
    VE_FORMS %in% c(1, 2),
    n_driver_rows == VE_FORMS
  )

table_01 <- tibble(
  Step = c(
    "Driver rows in cleaned driver-level LP file",
    "Distinct crashes in cleaned driver-level LP file",
    "Driver rows in one- or two-vehicle crashes",
    "Distinct one- or two-vehicle crashes",
    "Clean crashes with driver row count equal to VE_FORMS",
    "Driver rows in clean one- or two-vehicle crashes",
    "Clean one-vehicle crashes",
    "Clean two-vehicle crashes"
  ),
  Count = c(
    nrow(driver_lp),
    n_distinct(paste(driver_lp$YEAR_FILE, driver_lp$ST_CASE)),
    driver_lp %>% filter(VE_FORMS %in% c(1, 2)) %>% nrow(),
    driver_lp %>%
      filter(VE_FORMS %in% c(1, 2)) %>%
      distinct(YEAR_FILE, ST_CASE) %>%
      nrow(),
    nrow(clean_crashes),
    driver_lp %>%
      semi_join(clean_crashes %>% select(YEAR_FILE, ST_CASE),
                by = c("YEAR_FILE", "ST_CASE")) %>%
      nrow(),
    clean_crashes %>% filter(VE_FORMS == 1) %>% nrow(),
    clean_crashes %>% filter(VE_FORMS == 2) %>% nrow()
  )
)

write_csv(
  table_01,
  file.path(table_dir, "table_01_sample_construction.csv")
)

# ------------------------------------------------------------------------------
# 6. Table 2: Missingness by alcohol definition and period
# ------------------------------------------------------------------------------

table_02 <- crash_types_for_summary %>%
  mutate(
    vehicle_sample = case_when(
      VE_FORMS == 1 ~ "One-vehicle crashes",
      VE_FORMS == 2 ~ "Two-vehicle crashes",
      TRUE ~ "Other"
    ),
    missing_or_ambiguous = crash_type == "missing_or_ambiguous"
  ) %>%
  group_by(period, vehicle_sample, alcohol_definition_label) %>%
  summarize(
    total_crashes = n(),
    classifiable_crashes = sum(!missing_or_ambiguous),
    missing_or_ambiguous_crashes = sum(missing_or_ambiguous),
    missing_or_ambiguous_share = missing_or_ambiguous_crashes / total_crashes,
    .groups = "drop"
  ) %>%
  arrange(period, vehicle_sample, alcohol_definition_label)

write_csv(
  table_02,
  file.path(table_dir, "table_02_missingness_by_definition_period.csv")
)

# ------------------------------------------------------------------------------
# 7. Table 3: DD/DS/SS composition among all two-vehicle crashes
# ------------------------------------------------------------------------------

table_03_long <- crash_types_for_summary %>%
  filter(VE_FORMS == 2) %>%
  mutate(
    two_vehicle_class = case_when(
      crash_type %in% c("DD", "DS", "SS") ~ crash_type,
      TRUE ~ "Missing/ambiguous"
    ),
    two_vehicle_class = factor(
      two_vehicle_class,
      levels = c("DD", "DS", "SS", "Missing/ambiguous")
    )
  ) %>%
  group_by(period, alcohol_definition_label, two_vehicle_class) %>%
  summarize(n = n(), .groups = "drop_last") %>%
  mutate(
    total_two_vehicle_crashes = sum(n),
    share_among_all_two_vehicle = n / total_two_vehicle_crashes
  ) %>%
  ungroup()

table_03 <- table_03_long %>%
  select(
    period,
    alcohol_definition_label,
    two_vehicle_class,
    n,
    share_among_all_two_vehicle
  ) %>%
  pivot_wider(
    names_from = two_vehicle_class,
    values_from = c(n, share_among_all_two_vehicle),
    values_fill = 0
  ) %>%
  arrange(period, alcohol_definition_label)

write_csv(
  table_03,
  file.path(table_dir, "table_03_two_vehicle_composition_all_crashes.csv")
)

# ------------------------------------------------------------------------------
# 8. Table 4: DD/DS/SS composition among classifiable two-vehicle crashes only
# ------------------------------------------------------------------------------

table_04_long <- crash_types_for_summary %>%
  filter(
    VE_FORMS == 2,
    crash_type %in% c("DD", "DS", "SS")
  ) %>%
  mutate(
    two_vehicle_class = factor(crash_type, levels = c("DD", "DS", "SS"))
  ) %>%
  group_by(period, alcohol_definition_label, two_vehicle_class) %>%
  summarize(n = n(), .groups = "drop_last") %>%
  mutate(
    classifiable_two_vehicle_crashes = sum(n),
    share_among_classifiable = n / classifiable_two_vehicle_crashes
  ) %>%
  ungroup()

table_04 <- table_04_long %>%
  select(
    period,
    alcohol_definition_label,
    two_vehicle_class,
    n,
    share_among_classifiable
  ) %>%
  pivot_wider(
    names_from = two_vehicle_class,
    values_from = c(n, share_among_classifiable),
    values_fill = 0
  ) %>%
  arrange(period, alcohol_definition_label)

write_csv(
  table_04,
  file.path(table_dir, "table_04_two_vehicle_composition_classifiable_only.csv")
)

# ------------------------------------------------------------------------------
# 9. Table 5: Simple theta estimates by definition and period
# ------------------------------------------------------------------------------

table_05 <- theta_simple %>%
  mutate(
    period = factor(period, levels = period_levels),
    alcohol_definition_label = factor(
      alcohol_definition_label,
      levels = definition_levels
    )
  ) %>%
  select(
    period,
    alcohol_definition_label,
    DD,
    DS,
    SS,
    missing_or_ambiguous,
    total_two_vehicle_crashes,
    complete_case_two_vehicle_crashes,
    pct_missing_or_ambiguous,
    R,
    theta_two_vehicle,
    theta_note
  ) %>%
  arrange(period, alcohol_definition_label)

write_csv(
  table_05,
  file.path(table_dir, "table_05_theta_by_definition_period.csv")
)

# ------------------------------------------------------------------------------
# 10. Table 6: One-vehicle D/S composition by definition and period
# ------------------------------------------------------------------------------

table_06 <- one_vehicle_counts %>%
  mutate(
    period = factor(period, levels = period_levels),
    alcohol_definition_label = factor(
      alcohol_definition_label,
      levels = definition_levels
    )
  ) %>%
  select(
    period,
    alcohol_definition_label,
    D,
    S,
    missing_or_ambiguous,
    total_one_vehicle_crashes,
    complete_case_one_vehicle_crashes,
    pct_missing_or_ambiguous,
    share_drinking_among_complete
  ) %>%
  arrange(period, alcohol_definition_label)

write_csv(
  table_06,
  file.path(table_dir, "table_06_one_vehicle_composition.csv")
)

# ------------------------------------------------------------------------------
# 11. Plot A: Two-vehicle DD/DS/SS/missing shares, missing as own segment
# ------------------------------------------------------------------------------

plot_A_data <- crash_types_for_summary %>%
  filter(VE_FORMS == 2) %>%
  mutate(
    two_vehicle_class = case_when(
      crash_type %in% c("DD", "DS", "SS") ~ crash_type,
      TRUE ~ "Missing/ambiguous"
    ),
    two_vehicle_class = factor(
      two_vehicle_class,
      levels = c("SS", "DS", "DD", "Missing/ambiguous")
    )
  ) %>%
  group_by(period, alcohol_definition_label, two_vehicle_class) %>%
  summarize(n = n(), .groups = "drop_last") %>%
  mutate(
    total = sum(n),
    share = n / total
  ) %>%
  ungroup()

plot_A <- ggplot(
  plot_A_data,
  aes(
    x = alcohol_definition_label,
    y = share,
    fill = two_vehicle_class
  )
) +
  geom_col(width = 0.75) +
  facet_wrap(~ period, ncol = 2) +
  coord_flip() +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.03))
  ) +
  labs(
    title = "Two-Vehicle Fatal Crash Composition Including Missing Classifications",
    subtitle = "Missing/ambiguous classifications are shown as their own segment",
    x = NULL,
    y = "Share of all two-vehicle crashes",
    fill = "Crash type"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.major.y = element_blank()
  )

ggsave(
  filename = file.path(plot_dir, "plot_A_two_vehicle_shares_with_missing.png"),
  plot = plot_A,
  width = 13,
  height = 8,
  dpi = 300
)

# ------------------------------------------------------------------------------
# 12. Plot B: Complete-case two-vehicle DD/DS/SS shares only
# ------------------------------------------------------------------------------

plot_B_data <- crash_types_for_summary %>%
  filter(
    VE_FORMS == 2,
    crash_type %in% c("DD", "DS", "SS")
  ) %>%
  mutate(
    two_vehicle_class = factor(crash_type, levels = c("SS", "DS", "DD"))
  ) %>%
  group_by(period, alcohol_definition_label, two_vehicle_class) %>%
  summarize(n = n(), .groups = "drop_last") %>%
  mutate(
    complete_case_total = sum(n),
    share = n / complete_case_total
  ) %>%
  ungroup()

plot_B <- ggplot(
  plot_B_data,
  aes(
    x = alcohol_definition_label,
    y = share,
    fill = two_vehicle_class
  )
) +
  geom_col(width = 0.75) +
  facet_wrap(~ period, ncol = 2) +
  coord_flip() +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.03))
  ) +
  labs(
    title = "Complete-Case Two-Vehicle Fatal Crash Composition",
    subtitle = "Shares are renormalized among DD, DS, and SS only",
    x = NULL,
    y = "Share among classifiable two-vehicle crashes",
    fill = "Crash type"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.major.y = element_blank()
  )

ggsave(
  filename = file.path(plot_dir, "plot_B_two_vehicle_complete_case_shares.png"),
  plot = plot_B,
  width = 13,
  height = 8,
  dpi = 300
)

# ------------------------------------------------------------------------------
# 13. Plot C: Theta by definition with count labels / denominators
# ------------------------------------------------------------------------------

plot_C_data <- theta_simple %>%
  mutate(
    period = factor(period, levels = period_levels),
    alcohol_definition_label = factor(
      alcohol_definition_label,
      levels = definition_levels
    ),
    denom_label = paste0(
      "n = ",
      comma(complete_case_two_vehicle_crashes),
      " / ",
      comma(total_two_vehicle_crashes)
    )
  )

plot_C <- ggplot(
  plot_C_data,
  aes(
    x = alcohol_definition_label,
    y = theta_two_vehicle
  )
) +
  geom_point(size = 2.6) +
  geom_text(
    aes(label = denom_label),
    hjust = -0.05,
    size = 3
  ) +
  facet_wrap(~ period, ncol = 2) +
  coord_flip(clip = "off") +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.35))
  ) +
  labs(
    title = "Simple Levitt-Porter Two-Car Relative Risk Estimate",
    subtitle = "Label shows complete-case two-vehicle crashes / all two-vehicle crashes",
    x = NULL,
    y = expression(theta)
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    plot.margin = margin(10, 40, 10, 10)
  )

ggsave(
  filename = file.path(plot_dir, "plot_C_theta_with_denominators.png"),
  plot = plot_C,
  width = 13,
  height = 8,
  dpi = 300
)

# ------------------------------------------------------------------------------
# 14. Plot D: Missingness over time by year, not just period
# ------------------------------------------------------------------------------

plot_D_data <- crash_types %>%
  filter(VE_FORMS == 2) %>%
  mutate(
    alcohol_definition_label = factor(
      alcohol_definition_label,
      levels = definition_levels
    ),
    missing_or_ambiguous = crash_type == "missing_or_ambiguous"
  ) %>%
  group_by(YEAR, alcohol_definition_label) %>%
  summarize(
    total_two_vehicle_crashes = n(),
    missing_or_ambiguous_crashes = sum(missing_or_ambiguous),
    missing_or_ambiguous_share = missing_or_ambiguous_crashes /
      total_two_vehicle_crashes,
    .groups = "drop"
  )

# Labels at the final year for each line
plot_D_labels <- plot_D_data %>%
  group_by(alcohol_definition_label) %>%
  filter(YEAR == max(YEAR, na.rm = TRUE)) %>%
  ungroup()

plot_D <- ggplot(
  plot_D_data,
  aes(
    x = YEAR,
    y = missing_or_ambiguous_share,
    group = alcohol_definition_label,
    color = alcohol_definition_label
  )
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  geom_text_repel(
    data = plot_D_labels,
    aes(label = alcohol_definition_label),
    nudge_x = 0.7,
    direction = "y",
    hjust = 0,
    size = 3.4,
    segment.size = 0.25,
    segment.alpha = 0.7,
    show.legend = FALSE,
    max.overlaps = Inf
  ) +
  scale_x_continuous(
    breaks = sort(unique(plot_D_data$YEAR)),
    limits = c(
      min(plot_D_data$YEAR, na.rm = TRUE),
      max(plot_D_data$YEAR, na.rm = TRUE) + 2
    )
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  labs(
    title = "Missing or Ambiguous Two-Vehicle Classifications Over Time",
    subtitle = "Year-by-year missingness by alcohol definition, 2015-2024",
    x = "Year",
    y = "Missing / ambiguous share",
    color = "Alcohol definition"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(10, 90, 10, 10)
  )

ggsave(
  filename = file.path(plot_dir, "plot_D_missingness_over_time_by_year.png"),
  plot = plot_D,
  width = 13.5,
  height = 7,
  dpi = 300
)
