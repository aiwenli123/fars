################################################################################
# 10_four_group_prep.R
#
# Purpose:
#   Prep crash-count data for the four-BAC-group Levitt-Porter-style MLE.
#
#   Driver groups:
#     S = BAC == 0
#     L = 0 < BAC < 0.05
#     M = 0.05 <= BAC < 0.08
#     H = BAC >= 0.08
#
#   Two-car crash categories:
#     SS, SL, SM, SH, LL, LM, LH, MM, MH, HH
#
#   One-car crash categories:
#     S, L, M, H
#
# Main input:
#   data/processed/driver_miper_long_2015_2024.csv
#
# Main outputs:
#   data/processed/script10_four_group_prep/
#     four_group_crash_counts_by_cell_imp.csv
#     four_group_crash_counts_by_cell_imp.rds
#     four_group_crash_level_by_imp.csv
#     four_group_cell_sparsity_summary.csv
#     four_group_overall_counts_by_imp.csv
#     four_group_driver_bucket_summary_by_imp.csv
#     four_group_mle_prep_bundle.rds
#
# Notes:
#   - This script prepares the data only. It does not estimate theta/lambda.
#   - It uses the variable ANALYSIS_BAC_VALUE, which is measured BAC when available and
#     MIPER BAC when measured BAC is missing/invalid.
#   - The main analysis sample is MI all state-years. From script 09, we saw that
#     the data is sparse for some crash types, so better to keep as much data as possible.
################################################################################

rm(list = ls())

library(tidyverse)
library(scales)

# ------------------------------------------------------------------------------
# 0. User options
# ------------------------------------------------------------------------------

mi_driver_path <- "data/processed/driver_miper_long_2015_2024.csv"

out_dir <- "data/processed/script10_four_group_prep"
fig_dir <- "figures/script10_four_group_prep"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# use only broad cell specifications for now to prevent too much sparsity
cell_specs <- tibble(
  cell_spec = c(
    "all_data",
    "weekend",
    "rur_urb",
    "rur_urb_weekend"
  ),
  cell_var = c(
    "CELL_ALL_10",
    "CELL_WEEKEND_10",
    "CELL_RUR_URB_10",
    "CELL_RUR_URB_WEEKEND_10"
  ),
  cell_spec_label = c(
    "all data",
    "weekend",
    "rural/urban",
    "rural/urban x weekend"
  )
)

# Keep only clean rural/urban values for cell specs that use rural/urban.
rur_urb_keep_clear_only <- TRUE

# Sparse-cell cutoffs to summarize. Script 11 can decide which cutoff to use.
min_n2_cutoffs <- c(1, 5, 10, 25)

# Canonical order for categories.
groups4 <- c("S", "L", "M", "H")
two_car_types4 <- c("SS", "SL", "SM", "SH", "LL", "LM", "LH", "MM", "MH", "HH")
one_car_types4 <- c("S", "L", "M", "H")

# ------------------------------------------------------------------------------
# 1. Helper functions
# ------------------------------------------------------------------------------

safe_divide <- function(num, den) {
  ifelse(is.na(den) | den == 0, NA_real_, num / den)
}

make_crash_id <- function(year_file, st_case) {
  paste(year_file, st_case, sep = "_")
}

make_bucket4 <- function(bac) {
  case_when(
    is.na(bac) ~ NA_character_,
    bac == 0 ~ "S",
    bac > 0 & bac < 0.05 ~ "L",
    bac >= 0.05 & bac < 0.08 ~ "M",
    bac >= 0.08 ~ "H",
    TRUE ~ NA_character_
  )
}

sort_pair4 <- function(g1, g2) {
  # Force the pair order S < L < M < H, so e.g. LS becomes SL.
  map2_chr(g1, g2, function(x, y) {
    paste(sort(c(x, y), decreasing = FALSE, method = "radix"), collapse = "")
  }) %>%
    str_replace_all("LS", "SL") %>%
    str_replace_all("MS", "SM") %>%
    str_replace_all("HS", "SH") %>%
    str_replace_all("ML", "LM") %>%
    str_replace_all("HL", "LH") %>%
    str_replace_all("HM", "MH")
}

complete_count_columns <- function(df) {
  count_cols <- c(two_car_types4, paste0("one_", one_car_types4))
  missing_cols <- setdiff(count_cols, names(df))
  for (cc in missing_cols) df[[cc]] <- 0L
  df %>% select(any_of(names(df)[!names(df) %in% count_cols]), all_of(count_cols))
}

make_cell_counts_from_crashes <- function(crash_level_df, cell_var) {
  crash_level_df %>%
    mutate(cell_id = .data[[cell_var]]) %>%
    group_by(IMPUTATION, cell_id) %>%
    summarize(
      SS = sum(crash_type_4 == "SS", na.rm = TRUE),
      SL = sum(crash_type_4 == "SL", na.rm = TRUE),
      SM = sum(crash_type_4 == "SM", na.rm = TRUE),
      SH = sum(crash_type_4 == "SH", na.rm = TRUE),
      LL = sum(crash_type_4 == "LL", na.rm = TRUE),
      LM = sum(crash_type_4 == "LM", na.rm = TRUE),
      LH = sum(crash_type_4 == "LH", na.rm = TRUE),
      MM = sum(crash_type_4 == "MM", na.rm = TRUE),
      MH = sum(crash_type_4 == "MH", na.rm = TRUE),
      HH = sum(crash_type_4 == "HH", na.rm = TRUE),
      one_S = sum(crash_type_4 == "S", na.rm = TRUE),
      one_L = sum(crash_type_4 == "L", na.rm = TRUE),
      one_M = sum(crash_type_4 == "M", na.rm = TRUE),
      one_H = sum(crash_type_4 == "H", na.rm = TRUE),
      missing_two_vehicle = sum(VE_FORMS == 2 & crash_type_4 == "missing_or_ambiguous", na.rm = TRUE),
      missing_one_vehicle = sum(VE_FORMS == 1 & crash_type_4 == "missing_or_ambiguous", na.rm = TRUE),
      total_two_vehicle = sum(VE_FORMS == 2, na.rm = TRUE),
      total_one_vehicle = sum(VE_FORMS == 1, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      n2 = SS + SL + SM + SH + LL + LM + LH + MM + MH + HH,
      n1 = one_S + one_L + one_M + one_H,
      total_complete_case = n1 + n2,
      pct_two_vehicle_missing_or_ambiguous = safe_divide(missing_two_vehicle, total_two_vehicle),
      pct_one_vehicle_missing_or_ambiguous = safe_divide(missing_one_vehicle, total_one_vehicle)
    ) %>%
    filter(total_complete_case > 0)
}

summarize_cell_sparsity <- function(counts_df, min_n2_value) {
  counts_two <- counts_df %>% filter(n2 > 0)
  counts_kept <- counts_two %>% filter(n2 >= min_n2_value)
  pair_cols <- two_car_types4
  one_cols <- paste0("one_", one_car_types4)
  
  tibble(
    min_n2 = min_n2_value,
    n_cells_total = nrow(counts_df),
    n_theta_cells_available = nrow(counts_two),
    n_theta_cells_kept = nrow(counts_kept),
    pct_theta_cells_kept = safe_divide(nrow(counts_kept), nrow(counts_two)),
    two_car_crashes_available = sum(counts_two$n2, na.rm = TRUE),
    two_car_crashes_kept = sum(counts_kept$n2, na.rm = TRUE),
    pct_two_car_crashes_kept = safe_divide(
      sum(counts_kept$n2, na.rm = TRUE),
      sum(counts_two$n2, na.rm = TRUE)
    ),
    one_car_crashes_available = sum(counts_df$n1, na.rm = TRUE),
    n_cells_any_zero_pair = sum(rowSums(counts_two[pair_cols] == 0) > 0, na.rm = TRUE),
    pct_cells_any_zero_pair = safe_divide(
      sum(rowSums(counts_two[pair_cols] == 0) > 0, na.rm = TRUE),
      nrow(counts_two)
    ),
    n_cells_zero_LM_related_pair = sum(
      counts_two$SL == 0 | counts_two$SM == 0 | counts_two$LL == 0 |
        counts_two$LM == 0 | counts_two$MM == 0,
      na.rm = TRUE
    ),
    pct_cells_zero_LM_related_pair = safe_divide(
      sum(
        counts_two$SL == 0 | counts_two$SM == 0 | counts_two$LL == 0 |
          counts_two$LM == 0 | counts_two$MM == 0,
        na.rm = TRUE
      ),
      nrow(counts_two)
    ),
    min_pair_count_total = min(colSums(counts_df[pair_cols], na.rm = TRUE)),
    min_one_count_total = min(colSums(counts_df[one_cols], na.rm = TRUE))
  )
}

make_overall_counts <- function(crash_level_df) {
  crash_level_df %>%
    group_by(IMPUTATION, VE_FORMS, crash_type_4) %>%
    summarize(n_crashes = n(), .groups = "drop") %>%
    mutate(
      crash_size = case_when(
        VE_FORMS == 1 ~ "one-car",
        VE_FORMS == 2 ~ "two-car",
        TRUE ~ "other"
      )
    ) %>%
    select(IMPUTATION, crash_size, crash_type_4, n_crashes) %>%
    complete(
      IMPUTATION,
      nesting(crash_size),
      crash_type_4 = c(two_car_types4, one_car_types4, "missing_or_ambiguous"),
      fill = list(n_crashes = 0)
    ) %>%
    filter(
      (crash_size == "two-car" & crash_type_4 %in% c(two_car_types4, "missing_or_ambiguous")) |
        (crash_size == "one-car" & crash_type_4 %in% c(one_car_types4, "missing_or_ambiguous"))
    ) %>%
    arrange(IMPUTATION, crash_size, factor(crash_type_4, levels = c(one_car_types4, two_car_types4, "missing_or_ambiguous")))
}

# ------------------------------------------------------------------------------
# 2. Load MI long file and create analysis variables
# ------------------------------------------------------------------------------

if (!file.exists(mi_driver_path)) {
  stop("Missing input file: ", mi_driver_path, "\nPlease run 06_prep_mi_data.R first.")
}

driver_mi_long <- read_csv(mi_driver_path, show_col_types = FALSE) %>%
  rename_with(toupper)

cat("Loaded MI driver-level long file:\n")
cat("  Rows:", comma(nrow(driver_mi_long)), "\n")
cat("  Distinct imputations:", paste(sort(unique(driver_mi_long$IMPUTATION)), collapse = ", "), "\n\n")

required_cols <- c(
  "YEAR_FILE", "YEAR", "STATE", "ST_CASE", "VEH_NO", "PER_NO",
  "VE_FORMS", "HOUR", "IMPUTATION", "RUR_URB", "ANALYSIS_BAC_VALUE"
)

missing_required <- setdiff(required_cols, names(driver_mi_long))
if (length(missing_required) > 0) {
  stop("Missing required columns in driver_mi_long: ", paste(missing_required, collapse = ", "))
}

driver_mi_long <- driver_mi_long %>%
  mutate(
    PERIOD_10 = case_when(
      YEAR %in% 2015:2019 ~ "pre_covid",
      YEAR %in% 2020:2021 ~ "covid",
      YEAR %in% 2022:2024 ~ "post_covid",
      TRUE ~ NA_character_
    ),
    PERIOD_10 = factor(PERIOD_10, levels = c("pre_covid", "covid", "post_covid")),
    
    WEEKEND_10 = case_when(
      "WEEKEND_LP" %in% names(.) ~ as.logical(WEEKEND_LP),
      "DAY_WEEK" %in% names(.) & DAY_WEEK %in% c(1, 7) ~ TRUE,
      "DAY_WEEK" %in% names(.) & !is.na(DAY_WEEK) ~ FALSE,
      TRUE ~ NA
    ),
    
    RUR_URB_GROUP_10 = case_when(
      RUR_URB == 1 ~ "rural",
      RUR_URB == 2 ~ "urban",
      RUR_URB %in% c(6, 8, 9) ~ "other_unknown",
      TRUE ~ NA_character_
    ),
    RUR_URB_GROUP_10 = factor(
      RUR_URB_GROUP_10,
      levels = c("rural", "urban", "other_unknown")
    ),
    
    CRASH_ID_10 = make_crash_id(YEAR_FILE, ST_CASE),
    DRIVER_ID_10 = paste(YEAR_FILE, ST_CASE, VEH_NO, PER_NO, sep = "_"),
    
    BAC_GROUP_4 = make_bucket4(ANALYSIS_BAC_VALUE),
    BAC_GROUP_4 = factor(BAC_GROUP_4, levels = groups4),
    
    CELL_ALL_10 = "all_data",
    CELL_WEEKEND_10 = paste0("weekend_", WEEKEND_10),
    CELL_RUR_URB_10 = as.character(RUR_URB_GROUP_10),
    CELL_RUR_URB_WEEKEND_10 = paste(RUR_URB_GROUP_10, WEEKEND_10, sep = "_"),
    CELL_STATE_10 = as.character(STATE),
    CELL_STATE_WEEKEND_10 = paste(STATE, WEEKEND_10, sep = "_"),
    CELL_RUR_URB_STATE_10 = paste(RUR_URB_GROUP_10, STATE, sep = "_"),
    CELL_RUR_URB_STATE_WEEKEND_10 = paste(RUR_URB_GROUP_10, STATE, WEEKEND_10, sep = "_")
  ) %>%
  filter(!is.na(PERIOD_10), VE_FORMS %in% c(1, 2))

# Main four-group model uses MI all state-years. We keep rows with missing
# ANALYSIS_BAC_VALUE here so that incomplete crashes can be marked as
# missing_or_ambiguous after collapsing to the crash level.

if (isTRUE(rur_urb_keep_clear_only)) {
  # Keep a copy before filtering, but use clean rural/urban for main outputs because
  # the recommended starting cells include rural/urban.
  driver_mi_long_unfiltered_location <- driver_mi_long
  driver_mi_long <- driver_mi_long %>%
    filter(RUR_URB_GROUP_10 %in% c("rural", "urban"))
} else {
  driver_mi_long_unfiltered_location <- driver_mi_long
}

cat("Rows after analysis filters:\n")
cat("  Main clean-location rows:", comma(nrow(driver_mi_long)), "\n")
cat("  Distinct crashes:", comma(n_distinct(driver_mi_long$CRASH_ID_10)), "\n\n")

# ------------------------------------------------------------------------------
# 3. Driver-level bucket summaries
# ------------------------------------------------------------------------------

driver_bucket_summary_by_imp <- driver_mi_long %>%
  filter(!is.na(BAC_GROUP_4)) %>%
  group_by(IMPUTATION, BAC_GROUP_4) %>%
  summarize(n_driver_rows = n(), .groups = "drop_last") %>%
  mutate(share_driver_rows = n_driver_rows / sum(n_driver_rows)) %>%
  ungroup() %>%
  complete(
    IMPUTATION,
    BAC_GROUP_4 = factor(groups4, levels = groups4),
    fill = list(n_driver_rows = 0, share_driver_rows = 0)
  ) %>%
  arrange(IMPUTATION, BAC_GROUP_4)

write_csv(
  driver_bucket_summary_by_imp,
  file.path(out_dir, "four_group_driver_bucket_summary_by_imp.csv")
)

# ------------------------------------------------------------------------------
# 4. Collapse driver rows to one crash row per imputation
# ------------------------------------------------------------------------------

# Important: a two-car crash must have exactly two driver rows with non-missing BAC
# groups to be used in the complete-case two-car MLE counts. A one-car crash must
# have exactly one driver row with a non-missing BAC group.
crash_level_by_imp <- driver_mi_long %>%
  arrange(IMPUTATION, CRASH_ID_10, VEH_NO, PER_NO) %>%
  group_by(
    IMPUTATION, CRASH_ID_10,
    YEAR_FILE, YEAR, PERIOD_10, STATE, RUR_URB_GROUP_10, VE_FORMS,
    HOUR, WEEKEND_10,
    CELL_ALL_10, CELL_WEEKEND_10, CELL_RUR_URB_10, CELL_RUR_URB_WEEKEND_10,
    CELL_STATE_10, CELL_STATE_WEEKEND_10, CELL_RUR_URB_STATE_10,
    CELL_RUR_URB_STATE_WEEKEND_10
  ) %>%
  summarize(
    n_driver_rows = n(),
    n_bac_group_missing = sum(is.na(BAC_GROUP_4)),
    g1 = as.character(BAC_GROUP_4[1]),
    g2 = if_else(n() >= 2, as.character(BAC_GROUP_4[2]), NA_character_),
    .groups = "drop"
  ) %>%
  mutate(
    crash_type_4 = case_when(
      VE_FORMS == 1 & n_driver_rows == 1 & n_bac_group_missing == 0 ~ g1,
      VE_FORMS == 2 & n_driver_rows == 2 & n_bac_group_missing == 0 ~ sort_pair4(g1, g2),
      TRUE ~ "missing_or_ambiguous"
    ),
    crash_type_4 = factor(crash_type_4, levels = c(one_car_types4, two_car_types4, "missing_or_ambiguous"))
  )

write_csv(
  crash_level_by_imp,
  file.path(out_dir, "four_group_crash_level_by_imp.csv")
)

# ------------------------------------------------------------------------------
# 5. Overall four-group crash counts by imputation
# ------------------------------------------------------------------------------

overall_counts_by_imp <- make_overall_counts(crash_level_by_imp)

overall_counts_avg <- overall_counts_by_imp %>%
  group_by(crash_size, crash_type_4) %>%
  summarize(
    n_crashes_avg = mean(n_crashes, na.rm = TRUE),
    n_crashes_min = min(n_crashes, na.rm = TRUE),
    n_crashes_max = max(n_crashes, na.rm = TRUE),
    n_imputations = n_distinct(IMPUTATION),
    .groups = "drop"
  ) %>%
  arrange(crash_size, factor(crash_type_4, levels = c(one_car_types4, two_car_types4, "missing_or_ambiguous")))

write_csv(
  overall_counts_by_imp,
  file.path(out_dir, "four_group_overall_counts_by_imp.csv")
)

write_csv(
  overall_counts_avg,
  file.path(out_dir, "four_group_overall_counts_avg.csv")
)

# ------------------------------------------------------------------------------
# 6. Cell-level count tables for script 11 MLE
# ------------------------------------------------------------------------------

cell_count_list <- list()
cell_sparsity_list <- list()
counter <- 1

for (s in seq_len(nrow(cell_specs))) {
  this_cell_spec <- cell_specs$cell_spec[s]
  this_cell_var <- cell_specs$cell_var[s]
  this_cell_label <- cell_specs$cell_spec_label[s]
  
  cat("Creating cell counts for:", this_cell_spec, "\n")
  
  counts_s <- make_cell_counts_from_crashes(
    crash_level_df = crash_level_by_imp,
    cell_var = this_cell_var
  ) %>%
    mutate(
      cell_spec = this_cell_spec,
      cell_var = this_cell_var,
      cell_spec_label = this_cell_label,
      .before = cell_id
    )
  
  cell_count_list[[counter]] <- counts_s
  
  sparsity_s <- map_dfr(min_n2_cutoffs, ~ summarize_cell_sparsity(counts_s, .x)) %>%
    mutate(
      cell_spec = this_cell_spec,
      cell_var = this_cell_var,
      cell_spec_label = this_cell_label,
      .before = 1
    )
  
  cell_sparsity_list[[counter]] <- sparsity_s
  counter <- counter + 1
}

four_group_counts_by_cell_imp <- bind_rows(cell_count_list) %>%
  arrange(cell_spec, IMPUTATION, cell_id)

four_group_sparsity_summary <- bind_rows(cell_sparsity_list) %>%
  arrange(cell_spec, min_n2)

write_csv(
  four_group_counts_by_cell_imp,
  file.path(out_dir, "four_group_crash_counts_by_cell_imp.csv")
)

saveRDS(
  four_group_counts_by_cell_imp,
  file.path(out_dir, "four_group_crash_counts_by_cell_imp.rds")
)

write_csv(
  four_group_sparsity_summary,
  file.path(out_dir, "four_group_cell_sparsity_summary.csv")
)

# ------------------------------------------------------------------------------
# 7. A compact bundle for script 11
# ------------------------------------------------------------------------------

prep_bundle <- list(
  metadata = list(
    input_path = mi_driver_path,
    output_dir = out_dir,
    groups4 = groups4,
    two_car_types4 = two_car_types4,
    one_car_types4 = one_car_types4,
    cell_specs = cell_specs,
    rur_urb_keep_clear_only = rur_urb_keep_clear_only,
    created_at = as.character(Sys.time())
  ),
  cell_counts = four_group_counts_by_cell_imp,
  sparsity_summary = four_group_sparsity_summary,
  overall_counts_by_imp = overall_counts_by_imp,
  overall_counts_avg = overall_counts_avg,
  driver_bucket_summary_by_imp = driver_bucket_summary_by_imp
)

saveRDS(
  prep_bundle,
  file.path(out_dir, "four_group_mle_prep_bundle.rds")
)

# ------------------------------------------------------------------------------
# 8. Quick plot for sanity checking
# ------------------------------------------------------------------------------

plot_counts <- overall_counts_avg %>%
  filter(crash_type_4 != "missing_or_ambiguous") %>%
  mutate(
    crash_type_4 = factor(crash_type_4, levels = c(one_car_types4, two_car_types4))
  )

p_counts <- ggplot(plot_counts, aes(x = crash_type_4, y = n_crashes_avg)) +
  geom_col() +
  facet_wrap(~ crash_size, scales = "free_y") +
  labs(
    title = "Four-group crash counts prepared for MLE",
    subtitle = "Average across MI datasets; complete one- and two-car crash types only",
    x = "Four-group crash type",
    y = "Average number of crashes"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(
  file.path(fig_dir, "four_group_crash_counts_avg.png"),
  p_counts,
  width = 10,
  height = 5.5,
  dpi = 300
)

cat("\nFinished script 10 four-group MLE prep.\n")
cat("Main MLE input file:\n")
cat("  ", file.path(out_dir, "four_group_crash_counts_by_cell_imp.rds"), "\n")
