############################################################
# 06_prep_mi_data.R
#
# Purpose:
#   Combine 2015-2024 MIPER files with the main driver-level
#   Levitt-Porter analysis dataset.
#
# Main inputs:
#   data/processed/driver_analysis_lp_2015_2024.csv
#   data/raw/FARS{year}NationalCSV/MIPER.csv
#   data/raw/FARS{year}NationalCSV/miper.csv
#
# Main outputs:
#   data/processed/driver_miper_wide_2015_2024.csv
#   data/processed/driver_miper_long_2015_2024.csv
#   data/processed/mi_checks/...
#
# Notes:
#   - driver_miper_wide keeps one row per driver with P1-P10.
#   - driver_miper_long keeps one row per driver per imputation.
#   - For measured BAC, ALC_RES is in thousandths, so 80 = 0.080.
#   - MIPER P1-P10 are effectively in hundredths, so 8 = 0.08.
#   - For future MI analysis, analysis_bac_value uses:
#       actual measured BAC when valid;
#       imputed MIPER BAC only when measured BAC is missing/invalid.
############################################################

rm(list = ls())

library(tidyverse)

years <- 2015:2024

raw_dir <- "data/raw"
processed_dir <- "data/processed"
check_dir <- file.path(processed_dir, "mi_checks")

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(check_dir, recursive = TRUE, showWarnings = FALSE)

driver_path <- file.path(processed_dir, "driver_analysis_lp_2015_2024.csv")

p_cols <- paste0("P", 1:10)

############################################################
### 1. Helper functions
############################################################

find_miper_path <- function(year, raw_dir = "data/raw") {
  
  year_dir <- file.path(raw_dir, paste0("FARS", year, "NationalCSV"))
  
  if (!dir.exists(year_dir)) {
    stop("Could not find raw FARS directory for year ", year, ": ", year_dir)
  }
  
  possible_paths <- list.files(
    path = year_dir,
    pattern = "^miper\\.csv$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  if (length(possible_paths) == 0) {
    stop("Could not find MIPER/miper CSV file for year ", year, " in: ", year_dir)
  }
  
  if (length(possible_paths) > 1) {
    warning(
      "Found more than one MIPER-like file for year ", year,
      ". Using the first one: ", possible_paths[1]
    )
  }
  
  possible_paths[1]
}

read_miper_year <- function(year, raw_dir = "data/raw") {
  
  miper_path <- find_miper_path(year, raw_dir)
  
  message("Reading MIPER for ", year, ": ", miper_path)
  
  miper_year <- read.csv(miper_path) %>%
    as_tibble()
  
  # Standardize column names, just in case some files use lowercase names.
  names(miper_year) <- toupper(names(miper_year))
  
  required_cols <- c("ST_CASE", "VEH_NO", "PER_NO", p_cols)
  
  missing_cols <- setdiff(required_cols, names(miper_year))
  
  if (length(missing_cols) > 0) {
    stop(
      "MIPER file for year ", year, " is missing columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  miper_year %>%
    select(all_of(required_cols)) %>%
    mutate(
      YEAR_FILE = year,
      across(c(YEAR_FILE, ST_CASE, VEH_NO, PER_NO), as.integer),
      across(all_of(p_cols), as.integer)
    )
}

clean_measured_bac <- function(df, bac_col = "ALC_RES") {
  
  if (!(bac_col %in% names(df))) {
    stop("Could not find BAC column: ", bac_col)
  }
  
  df %>%
    mutate(
      measured_bac_raw = .data[[bac_col]],
      
      # Treat 0-940 as valid BAC values in thousandths.
      # 995-999 are missing/unknown codes.
      measured_bac_valid = case_when(
        is.na(measured_bac_raw) ~ FALSE,
        measured_bac_raw >= 0 & measured_bac_raw <= 940 ~ TRUE,
        TRUE ~ FALSE
      ),
      
      # Actual measured BAC on true BAC scale.
      # e.g., ALC_RES = 189 means BAC = 0.189.
      measured_bac_value = if_else(
        measured_bac_valid,
        measured_bac_raw / 1000,
        NA_real_
      ),
      
      # Convert measured BAC to MIPER scale.
      # Example: ALC_RES = 189 becomes 18.
      measured_bac_miper_scale = if_else(
        measured_bac_valid,
        as.integer(floor(measured_bac_raw / 10)),
        NA_integer_
      )
    )
}

############################################################
### 2. Read driver-level LP dataset
############################################################

driver_lp <- read.csv(driver_path) %>%
  as_tibble()

stopifnot(all(c("ST_CASE", "VEH_NO", "PER_NO") %in% names(driver_lp)))

# Make sure we have YEAR_FILE for joining.
# (ST_CASE numbers repeats across years, so YEAR_FILE is important)
if (!("YEAR_FILE" %in% names(driver_lp))) {
  if ("YEAR" %in% names(driver_lp)) {
    driver_lp <- driver_lp %>%
      mutate(YEAR_FILE = YEAR)
  } else {
    stop("Could not find YEAR_FILE or YEAR in driver_lp.")
  }
}

driver_lp <- driver_lp %>%
  mutate(
    across(c(YEAR_FILE, ST_CASE, VEH_NO, PER_NO), as.integer)
  ) %>%
  filter(YEAR_FILE %in% years)

############################################################
### 3. Read and combine MIPER files, 2015-2024
############################################################

miper_all <- map_dfr(years, read_miper_year, raw_dir = raw_dir)

############################################################
### 4. Check MIPER key uniqueness before joining
############################################################

# Should be empty:
miper_key_check <- miper_all %>%
  count(YEAR_FILE, ST_CASE, VEH_NO, PER_NO, name = "n_miper_rows") %>%
  filter(n_miper_rows > 1)

############################################################
### 5. Join MIPER to driver_lp
############################################################

driver_miper_wide <- driver_lp %>%
  left_join(
    miper_all,
    by = c("YEAR_FILE", "ST_CASE", "VEH_NO", "PER_NO")
  )

dim(driver_miper_wide)
dim(driver_lp) # should be 10 less columns

############################################################
### 6. Merge success checks
############################################################

merge_check_by_year <- driver_miper_wide %>%
  group_by(YEAR_FILE) %>%
  summarise(
    n_driver_rows = n(),
    n_matched_miper = sum(!is.na(P1)),
    n_unmatched_miper = sum(is.na(P1)),
    pct_matched_miper = mean(!is.na(P1)),
    .groups = "drop"
  )

print(merge_check_by_year, n = Inf)

write.csv(
  merge_check_by_year,
  file.path(check_dir, "miper_merge_check_by_year_2015_2024.csv"),
  row.names = FALSE
)

############################################################
### 7. Clean measured BAC and check measured-vs-MIPER match
############################################################

driver_miper_wide <- driver_miper_wide %>%
  clean_measured_bac(bac_col = "ALC_RES") %>%
  rowwise() %>%
  mutate(
    n_miper_values_present = sum(!is.na(c_across(all_of(p_cols)))),
    
    n_unique_miper_values = n_distinct(
      c_across(all_of(p_cols)),
      na.rm = TRUE
    ),
    
    # Do the 10 MIPER values vary?
    # For actual measured BAC cases, they should usually not vary.
    # For missing/unmeasured BAC cases, variation is expected.
    miper_values_vary = case_when(
      n_miper_values_present < 10 ~ NA,
      n_unique_miper_values > 1 ~ TRUE,
      TRUE ~ FALSE
    ),
    
    # If measured BAC is valid, all 10 MIPER values should equal
    # the measured BAC after converting ALC_RES to MIPER scale.
    all_10_equal_measured = case_when(
      !measured_bac_valid ~ NA,
      n_miper_values_present < 10 ~ NA,
      all(c_across(all_of(p_cols)) == measured_bac_miper_scale) ~ TRUE,
      TRUE ~ FALSE
    )
  ) %>%
  ungroup()

valid_measured_match_by_year <- driver_miper_wide %>%
  filter(measured_bac_valid) %>%
  group_by(YEAR_FILE) %>%
  summarise(
    n_valid_measured_bac = n(),
    n_matched_miper = sum(!is.na(P1)),
    n_all_10_equal_measured = sum(all_10_equal_measured == TRUE, na.rm = TRUE),
    n_not_all_equal_measured = sum(all_10_equal_measured == FALSE, na.rm = TRUE),
    pct_all_10_equal_measured = mean(all_10_equal_measured, na.rm = TRUE),
    .groups = "drop"
  )

print(valid_measured_match_by_year, n = Inf)

write.csv(
  valid_measured_match_by_year,
  file.path(check_dir, "valid_measured_bac_miper_match_by_year_2015_2024.csv"),
  row.names = FALSE
)


# Save mismatches for inspection.
# this file should have 0 rows.
measured_mismatches <- driver_miper_wide %>%
  filter(measured_bac_valid, all_10_equal_measured == FALSE) %>%
  select(
    YEAR_FILE, ST_CASE, VEH_NO, PER_NO,
    ALC_RES, measured_bac_value, measured_bac_miper_scale,
    all_of(p_cols),
    everything()
  )

# it has a few rows but fortunately none of the mismatches are severe.
write.csv(
  measured_mismatches,
  file.path(check_dir, "valid_measured_bac_miper_mismatches_2015_2024.csv"),
  row.names = FALSE
)

############################################################
### 8. Summary of measured vs imputed information in wide data
############################################################

mi_amount_by_year_wide <- driver_miper_wide %>%
  group_by(YEAR_FILE) %>%
  summarise(
    n_driver_rows = n(),
    n_matched_miper = sum(!is.na(P1)),
    
    n_valid_measured_bac = sum(measured_bac_valid),
    n_missing_or_invalid_measured_bac = sum(!measured_bac_valid),
    
    pct_valid_measured_bac = mean(measured_bac_valid),
    pct_missing_or_invalid_measured_bac = mean(!measured_bac_valid),
    
    n_miper_values_vary = sum(miper_values_vary == TRUE, na.rm = TRUE),
    n_miper_values_do_not_vary = sum(miper_values_vary == FALSE, na.rm = TRUE),
    pct_miper_values_vary = mean(miper_values_vary, na.rm = TRUE),
    
    .groups = "drop"
  )

mi_amount_total_wide <- driver_miper_wide %>%
  summarise(
    n_driver_rows = n(),
    n_matched_miper = sum(!is.na(P1)),
    
    n_valid_measured_bac = sum(measured_bac_valid),
    n_missing_or_invalid_measured_bac = sum(!measured_bac_valid),
    
    pct_valid_measured_bac = mean(measured_bac_valid),
    pct_missing_or_invalid_measured_bac = mean(!measured_bac_valid),
    
    n_miper_values_vary = sum(miper_values_vary == TRUE, na.rm = TRUE),
    n_miper_values_do_not_vary = sum(miper_values_vary == FALSE, na.rm = TRUE),
    pct_miper_values_vary = mean(miper_values_vary, na.rm = TRUE)
  )

print(mi_amount_by_year_wide, n = Inf, width = Inf)
print(mi_amount_total_wide, width = Inf)

write.csv(
  mi_amount_by_year_wide,
  file.path(check_dir, "mi_amount_by_year_wide_2015_2024.csv"),
  row.names = FALSE
)

write.csv(
  mi_amount_total_wide,
  file.path(check_dir, "mi_amount_total_wide_2015_2024.csv"),
  row.names = FALSE
)

############################################################
### 9. Create long MI analysis dataset
############################################################

driver_miper_long <- driver_miper_wide %>%
  pivot_longer(
    cols = all_of(p_cols),
    names_to = "imputation",
    values_to = "miper_bac_hundredths"
  ) %>%
  mutate(
    imputation = as.integer(str_remove(imputation, "^P")),
    
    # MIPER BAC on approximate BAC scale.
    # Example: miper_bac_hundredths = 8 means BAC = 0.08.
    miper_bac_value = miper_bac_hundredths / 100,
    
    # Main analysis BAC:
    # use exact measured BAC when available;
    # otherwise use imputed BAC from MIPER.
    analysis_bac_value = case_when(
      measured_bac_valid ~ measured_bac_value,
      !measured_bac_valid & !is.na(miper_bac_value) ~ miper_bac_value,
      TRUE ~ NA_real_
    ),
    
    analysis_bac_source = case_when(
      measured_bac_valid ~ "measured",
      !measured_bac_valid & !is.na(miper_bac_value) ~ "multiple_imputation",
      TRUE ~ "missing_no_miper"
    ),
    
    analysis_bac_source = factor(
      analysis_bac_source,
      levels = c("measured", "multiple_imputation", "missing_no_miper")
    ),
    
    # Threshold indicators to use in future.
    # These use analysis_bac_value, so they use exact measured BAC
    # when available and imputed BAC otherwise.
    mi_bac_ge_00 = case_when(
      is.na(analysis_bac_value) ~ NA_integer_,
      analysis_bac_value > 0 ~ 1L,
      TRUE ~ 0L
    ),
    
    mi_bac_ge_01 = case_when(
      is.na(analysis_bac_value) ~ NA_integer_,
      analysis_bac_value >= 0.01 ~ 1L,
      TRUE ~ 0L
    ),
    
    mi_bac_ge_05 = case_when(
      is.na(analysis_bac_value) ~ NA_integer_,
      analysis_bac_value >= 0.05 ~ 1L,
      TRUE ~ 0L
    ),
    
    mi_bac_ge_08 = case_when(
      is.na(analysis_bac_value) ~ NA_integer_,
      analysis_bac_value >= 0.08 ~ 1L,
      TRUE ~ 0L
    ),
    
    mi_bac_ge_10 = case_when(
      is.na(analysis_bac_value) ~ NA_integer_,
      analysis_bac_value >= 0.10 ~ 1L,
      TRUE ~ 0L
    )
  )

############################################################
### 10. Long-data summaries by imputation and year
############################################################

mi_source_summary_by_imp <- driver_miper_long %>%
  group_by(imputation, analysis_bac_source) %>%
  summarise(
    n = n(),
    .groups = "drop_last"
  ) %>%
  mutate(
    pct = n / sum(n)
  ) %>%
  ungroup()

mi_source_summary_by_year_imp <- driver_miper_long %>%
  group_by(YEAR_FILE, imputation, analysis_bac_source) %>%
  summarise(
    n = n(),
    .groups = "drop_last"
  ) %>%
  mutate(
    pct = n / sum(n)
  ) %>%
  ungroup()

mi_threshold_summary_by_imp <- driver_miper_long %>%
  group_by(imputation) %>%
  summarise(
    n_driver_rows = n(),
    
    n_measured = sum(analysis_bac_source == "measured"),
    n_multiple_imputation = sum(analysis_bac_source == "multiple_imputation"),
    n_missing_no_miper = sum(analysis_bac_source == "missing_no_miper"),
    
    pct_measured = mean(analysis_bac_source == "measured"),
    pct_multiple_imputation = mean(analysis_bac_source == "multiple_imputation"),
    pct_missing_no_miper = mean(analysis_bac_source == "missing_no_miper"),
    
    n_bac_gt_00 = sum(mi_bac_ge_00 == 1, na.rm = TRUE),
    n_bac_ge_01 = sum(mi_bac_ge_01 == 1, na.rm = TRUE),
    n_bac_ge_05 = sum(mi_bac_ge_05 == 1, na.rm = TRUE),
    n_bac_ge_08 = sum(mi_bac_ge_08 == 1, na.rm = TRUE),
    n_bac_ge_10 = sum(mi_bac_ge_10 == 1, na.rm = TRUE),
    
    pct_bac_gt_00 = mean(mi_bac_ge_00, na.rm = TRUE),
    pct_bac_ge_01 = mean(mi_bac_ge_01, na.rm = TRUE),
    pct_bac_ge_05 = mean(mi_bac_ge_05, na.rm = TRUE),
    pct_bac_ge_08 = mean(mi_bac_ge_08, na.rm = TRUE),
    pct_bac_ge_10 = mean(mi_bac_ge_10, na.rm = TRUE),
    
    .groups = "drop"
  )

mi_threshold_summary_by_year_imp <- driver_miper_long %>%
  group_by(YEAR_FILE, imputation) %>%
  summarise(
    n_driver_rows = n(),
    
    n_measured = sum(analysis_bac_source == "measured"),
    n_multiple_imputation = sum(analysis_bac_source == "multiple_imputation"),
    n_missing_no_miper = sum(analysis_bac_source == "missing_no_miper"),
    
    pct_measured = mean(analysis_bac_source == "measured"),
    pct_multiple_imputation = mean(analysis_bac_source == "multiple_imputation"),
    pct_missing_no_miper = mean(analysis_bac_source == "missing_no_miper"),
    
    n_bac_gt_00 = sum(mi_bac_ge_00 == 1, na.rm = TRUE),
    n_bac_ge_01 = sum(mi_bac_ge_01 == 1, na.rm = TRUE),
    n_bac_ge_05 = sum(mi_bac_ge_05 == 1, na.rm = TRUE),
    n_bac_ge_08 = sum(mi_bac_ge_08 == 1, na.rm = TRUE),
    n_bac_ge_10 = sum(mi_bac_ge_10 == 1, na.rm = TRUE),
    
    pct_bac_gt_00 = mean(mi_bac_ge_00, na.rm = TRUE),
    pct_bac_ge_01 = mean(mi_bac_ge_01, na.rm = TRUE),
    pct_bac_ge_05 = mean(mi_bac_ge_05, na.rm = TRUE),
    pct_bac_ge_08 = mean(mi_bac_ge_08, na.rm = TRUE),
    pct_bac_ge_10 = mean(mi_bac_ge_10, na.rm = TRUE),
    
    .groups = "drop"
  )

print(mi_source_summary_by_imp, n = Inf)
print(mi_threshold_summary_by_imp, n = Inf, width = Inf)

write.csv(
  mi_source_summary_by_imp,
  file.path(check_dir, "mi_source_summary_by_imp_2015_2024.csv"),
  row.names = FALSE
)

write.csv(
  mi_source_summary_by_year_imp,
  file.path(check_dir, "mi_source_summary_by_year_imp_2015_2024.csv"),
  row.names = FALSE
)

write.csv(
  mi_threshold_summary_by_imp,
  file.path(check_dir, "mi_threshold_summary_by_imp_2015_2024.csv"),
  row.names = FALSE
)

write.csv(
  mi_threshold_summary_by_year_imp,
  file.path(check_dir, "mi_threshold_summary_by_year_imp_2015_2024.csv"),
  row.names = FALSE
)

############################################################
### 11. Crash category summaries by imputation
############################################################
# These will be useful later if we want to compare how counts of
# DD, DS, SS, D, and S change across imputations under different
# BAC thresholds.
#
# This assumes driver_lp already has one_two_vehicle and
# VE_FORMS from previous scripts.

if (all(c("VE_FORMS", "ST_CASE") %in% names(driver_miper_long))) {
  
  lp_driver_counts_by_imp <- driver_miper_long %>%
    filter(VE_FORMS %in% c(1, 2)) %>%
    mutate(
      drunk_08 = mi_bac_ge_08,
      drunk_05 = mi_bac_ge_05,
      drunk_00 = mi_bac_ge_00
    ) %>%
    group_by(imputation, YEAR_FILE, ST_CASE, VE_FORMS) %>%
    summarise(
      n_drivers = n(),
      n_drunk_08 = sum(drunk_08 == 1, na.rm = TRUE),
      n_drunk_05 = sum(drunk_05 == 1, na.rm = TRUE),
      n_drunk_00 = sum(drunk_00 == 1, na.rm = TRUE),
      .groups = "drop"
    )
  
  lp_crash_type_summary_by_imp <- lp_driver_counts_by_imp %>%
    mutate(
      crash_type_08 = case_when(
        VE_FORMS == 1 & n_drunk_08 >= 1 ~ "D",
        VE_FORMS == 1 & n_drunk_08 == 0 ~ "S",
        VE_FORMS == 2 & n_drunk_08 == 2 ~ "DD",
        VE_FORMS == 2 & n_drunk_08 == 1 ~ "DS",
        VE_FORMS == 2 & n_drunk_08 == 0 ~ "SS",
        TRUE ~ NA_character_
      ),
      
      crash_type_05 = case_when(
        VE_FORMS == 1 & n_drunk_05 >= 1 ~ "D",
        VE_FORMS == 1 & n_drunk_05 == 0 ~ "S",
        VE_FORMS == 2 & n_drunk_05 == 2 ~ "DD",
        VE_FORMS == 2 & n_drunk_05 == 1 ~ "DS",
        VE_FORMS == 2 & n_drunk_05 == 0 ~ "SS",
        TRUE ~ NA_character_
      ),
      
      crash_type_00 = case_when(
        VE_FORMS == 1 & n_drunk_00 >= 1 ~ "D",
        VE_FORMS == 1 & n_drunk_00 == 0 ~ "S",
        VE_FORMS == 2 & n_drunk_00 == 2 ~ "DD",
        VE_FORMS == 2 & n_drunk_00 == 1 ~ "DS",
        VE_FORMS == 2 & n_drunk_00 == 0 ~ "SS",
        TRUE ~ NA_character_
      )
    ) %>%
    group_by(imputation) %>%
    summarise(
      DD_08 = sum(crash_type_08 == "DD", na.rm = TRUE),
      DS_08 = sum(crash_type_08 == "DS", na.rm = TRUE),
      SS_08 = sum(crash_type_08 == "SS", na.rm = TRUE),
      D_08  = sum(crash_type_08 == "D",  na.rm = TRUE),
      S_08  = sum(crash_type_08 == "S",  na.rm = TRUE),
      
      DD_05 = sum(crash_type_05 == "DD", na.rm = TRUE),
      DS_05 = sum(crash_type_05 == "DS", na.rm = TRUE),
      SS_05 = sum(crash_type_05 == "SS", na.rm = TRUE),
      D_05  = sum(crash_type_05 == "D",  na.rm = TRUE),
      S_05  = sum(crash_type_05 == "S",  na.rm = TRUE),
      
      DD_00 = sum(crash_type_00 == "DD", na.rm = TRUE),
      DS_00 = sum(crash_type_00 == "DS", na.rm = TRUE),
      SS_00 = sum(crash_type_00 == "SS", na.rm = TRUE),
      D_00  = sum(crash_type_00 == "D",  na.rm = TRUE),
      S_00  = sum(crash_type_00 == "S",  na.rm = TRUE),
      
      .groups = "drop"
    )
  
  print(lp_crash_type_summary_by_imp, n = Inf, width = Inf)
  
  write.csv(
    lp_driver_counts_by_imp,
    file.path(check_dir, "lp_driver_counts_by_imp_year_crash_2015_2024.csv"),
    row.names = FALSE
  )
  
  write.csv(
    lp_crash_type_summary_by_imp,
    file.path(check_dir, "lp_crash_type_summary_by_imp_2015_2024.csv"),
    row.names = FALSE
  )
}

############################################################
### 12. Save main datasets
############################################################

write.csv(
  driver_miper_wide,
  file.path(processed_dir, "driver_miper_wide_2015_2024.csv"),
  row.names = FALSE
)

write.csv(
  driver_miper_long,
  file.path(processed_dir, "driver_miper_long_2015_2024.csv"),
  row.names = FALSE
)

message("Saved wide MI driver dataset to: ",
        file.path(processed_dir, "driver_miper_wide_2015_2024.csv"))

message("Saved long MI driver dataset to: ",
        file.path(processed_dir, "driver_miper_long_2015_2024.csv"))

message("Saved checks and summaries to: ", check_dir)
