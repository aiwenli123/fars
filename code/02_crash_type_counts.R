################################################################################
# 02_crash_type_counts.R
#
# Purpose:
#   Collapse the cleaned driver-level FARS sample into crash-level
#   drunk/sober crash types.
#
#   For two-vehicle crashes:
#     DD = two drinking drivers
#     DS = one drinking driver and one sober driver
#     SS = two sober drivers
#
#   For one-vehicle crashes:
#     D = one drinking driver
#     S = one sober driver
#
#   Then compute simple Levitt-Porter two-car relative risk estimates using:
#
#     R = A_DS^2 / (A_DD * A_SS)
#
#   and
#
#     theta = ((R - 2) + sqrt(R^2 - 4R)) / 2, if R >= 4
#
#   If R < 4, set theta = 1 as the boundary estimate.
#
# Inputs:
#   data/processed/driver_analysis_lp_2015_2024.csv
#   data/processed/vehicle_2015_2024_core.csv   # used to get DR_DRINK column
#
# Main outputs:
#   data/processed/lp_crash_types_2015_2024.csv
#   data/processed/lp_crash_type_counts_2015_2024.csv
#   data/processed/lp_theta_simple_2015_2024.csv
#   data/processed/lp_one_vehicle_counts_2015_2024.csv
#
#   figures/lp_crash_type_shares_by_definition_all_years.png
#   figures/lp_theta_by_definition_and_period.png
#   figures/lp_missing_by_definition_and_period.png
################################################################################

rm(list = ls())

library(tidyverse)

# -------------------------------------------------------------------------
# 0. File paths
# -------------------------------------------------------------------------

driver_path  <- "data/processed/driver_analysis_lp_2015_2024.csv"
vehicle_path <- "data/processed/vehicle_2015_2024_core.csv"

out_dir <- "data/processed"
fig_dir <- "figures"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# 1. Load cleaned driver-level LP sample
# -------------------------------------------------------------------------

driver_lp <- read_csv(driver_path) %>%
  rename_with(toupper) # convert col names to uppercase

cat("\nLoaded driver-level LP sample:\n")
cat("Rows:", nrow(driver_lp), "\n")
cat("Columns:", ncol(driver_lp), "\n\n")

cat("Variables in driver_lp:\n")
print(names(driver_lp))

# -------------------------------------------------------------------------
# 2. Add DR_DRINK from vehicle file if not already present
# -------------------------------------------------------------------------

if (!"DR_DRINK" %in% names(driver_lp)) {
  
  cat("\nDR_DRINK not found in driver-level file. Joining from vehicle file...\n")
  
  vehicle_core <- read_csv(vehicle_path, show_col_types = FALSE) %>%
    rename_with(toupper)
  
  vehicle_drink <- vehicle_core %>%
    select(
      YEAR_FILE,
      ST_CASE,
      VEH_NO,
      any_of(c("DR_DRINK", "DR_DRINKNAME"))
    )
  
  driver_lp <- driver_lp %>%
    left_join(
      vehicle_drink,
      by = c("YEAR_FILE", "ST_CASE", "VEH_NO")
    )
  
} else {
  cat("\nDR_DRINK already found in driver-level file.\n")
}

# -------------------------------------------------------------------------
# 3. Basic sanity checks
# -------------------------------------------------------------------------

required_vars <- c(
  "YEAR_FILE", "YEAR", "ST_CASE", "VEH_NO", "VE_FORMS",
  "DRINKING", "ALC_RES"
)

missing_required <- setdiff(required_vars, names(driver_lp))

if (length(missing_required) > 0) {
  stop(
    "The following required variables are missing from driver_lp: ",
    paste(missing_required, collapse = ", ")
  )
}

cat("\nCounts by year and VE_FORMS:\n")
print(
  driver_lp %>%
    count(YEAR, VE_FORMS) %>%
    arrange(YEAR, VE_FORMS),
  n = Inf
)

# Check number of driver rows per crash
driver_rows_by_crash <- driver_lp %>%
  group_by(YEAR_FILE, ST_CASE) %>%
  summarize(
    YEAR = first(YEAR),
    VE_FORMS = first(VE_FORMS),
    n_driver_rows = n(),
    .groups = "drop"
  )

cat("\nDriver rows by crash vs VE_FORMS:\n")
print(
  driver_rows_by_crash %>%
    count(VE_FORMS, n_driver_rows, sort = TRUE),
  n = Inf
)

bad_driver_row_counts <- driver_rows_by_crash %>%
  filter(VE_FORMS %in% c(1, 2), n_driver_rows != VE_FORMS)

if (nrow(bad_driver_row_counts) > 0) {
  warning(
    "Some crashes have n_driver_rows != VE_FORMS. ",
    "These should have been filtered in script 01, but will be excluded here."
  )
}

# Keep only clean one- and two-vehicle crashes with expected driver row counts
clean_crashes <- driver_rows_by_crash %>%
  filter(
    VE_FORMS %in% c(1, 2),
    n_driver_rows == VE_FORMS
  )

driver_lp_clean <- driver_lp %>%
  semi_join(
    clean_crashes %>% select(YEAR_FILE, ST_CASE),
    by = c("YEAR_FILE", "ST_CASE")
  )

cat("\nClean driver-level sample after sanity filtering:\n")
cat("Rows:", nrow(driver_lp_clean), "\n")
cat("Crashes:", n_distinct(paste(driver_lp_clean$YEAR_FILE, driver_lp_clean$ST_CASE)), "\n\n")

# -------------------------------------------------------------------------
# 4. Create / refresh alcohol variables
# -------------------------------------------------------------------------

# Notes:
#   - DRINKING:
#       1 = drinking
#       0 = not drinking
#       other values treated as missing in complete-case version
#
#   - DR_DRINK:
#       1 = drinking
#       0 = not drinking
#       other values treated as missing (although there's no "other" values for 2015-2024)
#
#   - ALC_RES:
#       valid observed BAC is assumed to be 0 <= ALC_RES < 995
#       BAC value is ALC_RES / 1000
#       e.g., ALC_RES = 80 means BAC = 0.080

driver_lp_clean <- driver_lp_clean %>%
  mutate(
    # Police drinking definition: complete-case version
    police_drinking_complete = case_when(
      DRINKING == 1 ~ 1,
      DRINKING == 0 ~ 0,
      TRUE ~ NA_real_
    ),
    
    # Police drinking definition: unknown / non-informative as no
    police_drinking_unknown_as_no = case_when(
      DRINKING == 1 ~ 1,
      TRUE ~ 0
    ),
    
    # Vehicle-level drinking definition
    vehicle_dr_drink_complete = case_when(
      DR_DRINK == 1 ~ 1,
      DR_DRINK == 0 ~ 0,
      TRUE ~ NA_real_
    ),
    
    # Valid observed BAC
    bac_valid = !is.na(ALC_RES) & ALC_RES >= 0 & ALC_RES < 995,
    bac_value = if_else(bac_valid, ALC_RES / 1000, NA_real_),
    
    # BAC threshold definitions
    bac_ge_08 = case_when(
      bac_valid & ALC_RES >= 80 ~ 1,
      bac_valid & ALC_RES < 80 ~ 0,
      TRUE ~ NA_real_
    ),
    
    bac_ge_10 = case_when(
      bac_valid & ALC_RES >= 100 ~ 1,
      bac_valid & ALC_RES < 100 ~ 0,
      TRUE ~ NA_real_
    ),
    
    # Period variable
    covid_period = case_when(
      YEAR %in% 2015:2019 ~ "2015-2019 pre-COVID",
      YEAR %in% 2020:2021 ~ "2020-2021 COVID",
      YEAR %in% 2022:2024 ~ "2022-2024 post-COVID",
      TRUE ~ NA_character_
    )
  )

# -------------------------------------------------------------------------
# 5. Put alcohol definitions into long format
# -------------------------------------------------------------------------

alcohol_definition_vars <- c(
  "police_drinking_complete",
  "police_drinking_unknown_as_no",
  "vehicle_dr_drink_complete",
  "bac_ge_08",
  "bac_ge_10"
)

driver_long <- driver_lp_clean %>%
  pivot_longer(
    cols = all_of(alcohol_definition_vars),
    names_to = "alcohol_definition",
    values_to = "drinking_status"
  ) %>%
  mutate(
    alcohol_definition_label = case_when(
      alcohol_definition == "police_drinking_complete" ~
        "Police drinking, complete case",
      alcohol_definition == "police_drinking_unknown_as_no" ~
        "Police drinking, unknown as no",
      alcohol_definition == "vehicle_dr_drink_complete" ~
        "Vehicle DR_DRINK, complete case",
      alcohol_definition == "bac_ge_08" ~
        "Observed BAC >= 0.08",
      alcohol_definition == "bac_ge_10" ~
        "Observed BAC >= 0.10",
      TRUE ~ alcohol_definition
    ),
    driver_status_label = case_when(
      drinking_status == 1 ~ "D",
      drinking_status == 0 ~ "S",
      is.na(drinking_status) ~ "missing_or_noninformative",
      TRUE ~ "missing_or_noninformative"
    )
  )

cat("\nDriver-level alcohol definition summary:\n")
print(
  driver_long %>%
    count(alcohol_definition_label, driver_status_label) %>%
    group_by(alcohol_definition_label) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup() %>%
    arrange(alcohol_definition_label, driver_status_label),
  n = Inf
)

# -------------------------------------------------------------------------
# 6. Collapse driver rows to crash-level types
# -------------------------------------------------------------------------

crash_types <- driver_long %>%
  group_by(
    YEAR_FILE,
    YEAR,
    STATE,
    ST_CASE,
    VE_FORMS,
    HOUR,
    DAY_WEEK,
    covid_period,
    alcohol_definition,
    alcohol_definition_label
  ) %>%
  summarize(
    n_drivers = n(),
    n_drinking = sum(drinking_status == 1, na.rm = TRUE),
    n_sober = sum(drinking_status == 0, na.rm = TRUE),
    n_missing = sum(is.na(drinking_status)),
    .groups = "drop"
  ) %>%
  mutate(
    complete_case_crash = n_missing == 0,
    
    crash_type = case_when(
      VE_FORMS == 2 & n_missing == 0 & n_drinking == 2 & n_sober == 0 ~ "DD",
      VE_FORMS == 2 & n_missing == 0 & n_drinking == 1 & n_sober == 1 ~ "DS",
      VE_FORMS == 2 & n_missing == 0 & n_drinking == 0 & n_sober == 2 ~ "SS",
      
      VE_FORMS == 1 & n_missing == 0 & n_drinking == 1 & n_sober == 0 ~ "D",
      VE_FORMS == 1 & n_missing == 0 & n_drinking == 0 & n_sober == 1 ~ "S",
      
      TRUE ~ "missing_or_ambiguous"
    )
  )

cat("\nCrash-level type summary:\n")
print(
  crash_types %>%
    count(alcohol_definition_label, VE_FORMS, crash_type) %>%
    arrange(alcohol_definition_label, VE_FORMS, crash_type),
  n = Inf
)

# Save crash-level file
write_csv(
  crash_types,
  file.path(out_dir, "lp_crash_types_2015_2024.csv")
)

# -------------------------------------------------------------------------
# 7. Create pooled and COVID-period grouping variable
# -------------------------------------------------------------------------

crash_types_for_summary <- crash_types %>%
  mutate(all_years_group = "2015-2024 all years") %>%
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
    period = factor(
      period,
      levels = c(
        "2015-2024 all years",
        "2015-2019 pre-COVID",
        "2020-2021 COVID",
        "2022-2024 post-COVID"
      )
    )
  )

# -------------------------------------------------------------------------
# 8. Crash type count tables
# -------------------------------------------------------------------------

crash_type_counts <- crash_types_for_summary %>%
  count(
    summary_level,
    period,
    alcohol_definition,
    alcohol_definition_label,
    VE_FORMS,
    crash_type,
    name = "n_crashes"
  ) %>%
  group_by(
    summary_level,
    period,
    alcohol_definition,
    alcohol_definition_label,
    VE_FORMS
  ) %>%
  mutate(
    total_crashes = sum(n_crashes),
    share = n_crashes / total_crashes
  ) %>%
  ungroup() %>%
  arrange(
    summary_level,
    period,
    alcohol_definition,
    VE_FORMS,
    crash_type
  )

write_csv(
  crash_type_counts,
  file.path(out_dir, "lp_crash_type_counts_2015_2024.csv")
)

cat("\nTwo-vehicle crash type counts:\n")
print(
  crash_type_counts %>%
    filter(VE_FORMS == 2) %>%
    select(
      summary_level,
      period,
      alcohol_definition_label,
      crash_type,
      n_crashes,
      total_crashes,
      share
    ),
  n = Inf
)

cat("\nOne-vehicle crash type counts:\n")
print(
  crash_type_counts %>%
    filter(VE_FORMS == 1) %>%
    select(
      summary_level,
      period,
      alcohol_definition_label,
      crash_type,
      n_crashes,
      total_crashes,
      share
    ),
  n = Inf
)

# -------------------------------------------------------------------------
# 9. Simple Levitt-Porter theta estimates from two-vehicle crash counts
# -------------------------------------------------------------------------

compute_theta <- function(A_DD, A_DS, A_SS) {
  
  # Need positive DD and SS counts for denominator
  if (is.na(A_DD) | is.na(A_DS) | is.na(A_SS)) {
    return(NA_real_)
  }
  
  if (A_DD <= 0 | A_SS <= 0) {
    return(NA_real_)
  }
  
  R <- (A_DS^2) / (A_DD * A_SS)
  
  if (is.na(R)) {
    return(NA_real_)
  }
  
  # Boundary case: R < 4 is not consistent with the simple binomial model.
  # Levitt-Porter MLE places theta at 1 in this case.
  if (R < 4) {
    return(1)
  }
  
  theta <- ((R - 2) + sqrt(R^2 - 4 * R)) / 2
  
  return(theta)
}

theta_simple <- crash_types_for_summary %>%
  filter(VE_FORMS == 2) %>%
  mutate(
    two_vehicle_class = case_when(
      crash_type %in% c("DD", "DS", "SS") ~ crash_type,
      TRUE ~ "missing_or_ambiguous"
    )
  ) %>%
  count(
    summary_level,
    period,
    alcohol_definition,
    alcohol_definition_label,
    two_vehicle_class,
    name = "n"
  ) %>%
  pivot_wider(
    names_from = two_vehicle_class,
    values_from = n,
    values_fill = 0
  )

# Make sure all expected columns exist, even if a category is absent
for (col in c("DD", "DS", "SS", "missing_or_ambiguous")) {
  if (!col %in% names(theta_simple)) {
    theta_simple[[col]] <- 0L
  }
}

theta_simple <- theta_simple %>%
  mutate(
    total_two_vehicle_crashes = DD + DS + SS + missing_or_ambiguous,
    complete_case_two_vehicle_crashes = DD + DS + SS,
    pct_missing_or_ambiguous = missing_or_ambiguous / total_two_vehicle_crashes,
    
    R = case_when(
      DD > 0 & SS > 0 ~ (DS^2) / (DD * SS),
      TRUE ~ NA_real_
    ),
    
    theta_two_vehicle = pmap_dbl(
      list(DD, DS, SS),
      ~ compute_theta(..1, ..2, ..3)
    ),
    
    theta_note = case_when(
      is.na(theta_two_vehicle) ~ "Cannot estimate: DD or SS count is zero/missing",
      !is.na(R) & R < 4 ~ "Boundary estimate: R < 4, theta set to 1",
      TRUE ~ "Interior estimate"
    )
  ) %>%
  arrange(
    summary_level,
    period,
    alcohol_definition
  )

write_csv(
  theta_simple,
  file.path(out_dir, "lp_theta_simple_2015_2024.csv")
)

cat("\nSimple Levitt-Porter theta estimates:\n")
print(theta_simple, n = Inf)

# -------------------------------------------------------------------------
# 10. One-vehicle count table
# -------------------------------------------------------------------------

one_vehicle_counts <- crash_types_for_summary %>%
  filter(VE_FORMS == 1) %>%
  mutate(
    one_vehicle_class = case_when(
      crash_type %in% c("D", "S") ~ crash_type,
      TRUE ~ "missing_or_ambiguous"
    )
  ) %>%
  count(
    summary_level,
    period,
    alcohol_definition,
    alcohol_definition_label,
    one_vehicle_class,
    name = "n"
  ) %>%
  pivot_wider(
    names_from = one_vehicle_class,
    values_from = n,
    values_fill = 0
  )

# Make sure all expected columns exist
for (col in c("D", "S", "missing_or_ambiguous")) {
  if (!col %in% names(one_vehicle_counts)) {
    one_vehicle_counts[[col]] <- 0L
  }
}

one_vehicle_counts <- one_vehicle_counts %>%
  mutate(
    total_one_vehicle_crashes = D + S + missing_or_ambiguous,
    complete_case_one_vehicle_crashes = D + S,
    pct_missing_or_ambiguous = missing_or_ambiguous / total_one_vehicle_crashes,
    share_drinking_among_complete = D / complete_case_one_vehicle_crashes
  ) %>%
  arrange(
    summary_level,
    period,
    alcohol_definition
  )

write_csv(
  one_vehicle_counts,
  file.path(out_dir, "lp_one_vehicle_counts_2015_2024.csv")
)

cat("\nOne-vehicle crash counts:\n")
print(one_vehicle_counts, n = Inf)

# -------------------------------------------------------------------------
# 11. Helpful wide-format summary tables for quick inspection
# -------------------------------------------------------------------------

theta_all_years_wide <- theta_simple %>%
  filter(summary_level == "All years pooled") %>%
  select(
    alcohol_definition_label,
    DD,
    DS,
    SS,
    missing_or_ambiguous,
    pct_missing_or_ambiguous,
    R,
    theta_two_vehicle,
    theta_note
  ) %>%
  arrange(alcohol_definition_label)

cat("\nAll-years pooled theta table:\n")
print(theta_all_years_wide, n = Inf)

theta_by_covid_period_wide <- theta_simple %>%
  filter(summary_level == "COVID period") %>%
  select(
    period,
    alcohol_definition_label,
    DD,
    DS,
    SS,
    missing_or_ambiguous,
    pct_missing_or_ambiguous,
    R,
    theta_two_vehicle,
    theta_note
  ) %>%
  arrange(period, alcohol_definition_label)

cat("\nTheta table by COVID period:\n")
print(theta_by_covid_period_wide, n = Inf)

# Optional: save these quick-view tables too
write_csv(
  theta_all_years_wide,
  file.path(out_dir, "lp_theta_simple_all_years_quick_table.csv")
)

write_csv(
  theta_by_covid_period_wide,
  file.path(out_dir, "lp_theta_simple_by_covid_period_quick_table.csv")
)

# -------------------------------------------------------------------------
# 12. Plots
# -------------------------------------------------------------------------

# Plot 1: DD / DS / SS shares by alcohol definition, all years pooled
plot_two_vehicle_shares_all <- crash_type_counts %>%
  filter(
    summary_level == "All years pooled",
    VE_FORMS == 2,
    crash_type %in% c("DD", "DS", "SS")
  ) %>%
  mutate(
    crash_type = factor(crash_type, levels = c("SS", "DS", "DD"))
  ) %>%
  ggplot(
    aes(
      x = alcohol_definition_label,
      y = share,
      fill = crash_type
    )
  ) +
  geom_col() +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Two-Vehicle Fatal Crash Types by Alcohol Definition",
    subtitle = "2015-2024 nighttime one- and two-vehicle LP-style sample",
    x = NULL,
    y = "Share of two-vehicle crashes",
    fill = "Crash type"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  filename = file.path(fig_dir, "lp_crash_type_shares_by_definition_all_years.png"),
  plot = plot_two_vehicle_shares_all,
  width = 10,
  height = 6,
  dpi = 300
)

# Plot 2: theta by alcohol definition and period
plot_theta_period <- theta_simple %>%
  mutate(
    period = factor(
      period,
      levels = c(
        "2015-2024 all years",
        "2015-2019 pre-COVID",
        "2020-2021 COVID",
        "2022-2024 post-COVID"
      )
    )
  ) %>%
  ggplot(
    aes(
      x = alcohol_definition_label,
      y = theta_two_vehicle
    )
  ) +
  geom_point(size = 2.5) +
  facet_wrap(~ period) +
  coord_flip() +
  labs(
    title = "Simple Levitt-Porter Two-Car Relative Risk Estimate",
    subtitle = "Theta estimated from DD, DS, and SS two-vehicle fatal crash counts",
    x = NULL,
    y = expression(theta)
  ) +
  theme_minimal(base_size = 12)

ggsave(
  filename = file.path(fig_dir, "lp_theta_by_definition_and_period.png"),
  plot = plot_theta_period,
  width = 12,
  height = 7,
  dpi = 300
)

# Plot 3: missing / ambiguous share by definition and period
plot_missing_period <- theta_simple %>%
  mutate(
    period = factor(
      period,
      levels = c(
        "2015-2024 all years",
        "2015-2019 pre-COVID",
        "2020-2021 COVID",
        "2022-2024 post-COVID"
      )
    )
  ) %>%
  ggplot(
    aes(
      x = alcohol_definition_label,
      y = pct_missing_or_ambiguous
    )
  ) +
  geom_col() +
  facet_wrap(~ period) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Missing or Ambiguous Two-Vehicle Crash Classifications",
    subtitle = "Share of two-vehicle crashes not classifiable as DD, DS, or SS",
    x = NULL,
    y = "Missing / ambiguous share"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  filename = file.path(fig_dir, "lp_missing_by_definition_and_period.png"),
  plot = plot_missing_period,
  width = 12,
  height = 7,
  dpi = 300
)
