# Classify drivers in the ten multiply imputed FARS datasets into four BAC groups,
# then construct one- and two-vehicle fatal crash types. 
# Prepare the crash-level data for bootstrapping.
#
# Input:  data/processed/driver_miper_long_2015_2024.csv
# Outputs: data/processed/script03_four_group_prep/
#          four_group_crash_level_by_imp.csv
#          four_group_crash_counts_by_cell_imp.csv
#          four_group_crash_counts_by_cell_imp.rds

library(tidyverse)

input_path <- "data/processed/driver_miper_long_2015_2024.csv"
out_dir <- "data/processed/script03_four_group_prep"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

groups <- c("S", "L", "M", "H")

# Order each two-driver pair as S, L, M, H (e.g., LS becomes SL).
order_pair <- function(g1, g2) {
  i <- match(g1, groups)
  j <- match(g2, groups)
  paste0(groups[pmin(i, j)], groups[pmax(i, j)])
}

# Use measured BAC when available and imputed BAC otherwise
# Keep only rural/urban crashes
drivers <- read_csv(input_path, show_col_types = FALSE) %>%
  rename_with(toupper) %>%
  filter(RUR_URB %in% c(1, 2), VE_FORMS %in% c(1, 2)) %>%
  mutate(
    RUR_URB_GROUP_10 = if_else(RUR_URB == 1, "rural", "urban"),
    WEEKEND_10 = as.logical(WEEKEND_LP),
    CELL_RUR_URB_WEEKEND_10 = paste(RUR_URB_GROUP_10, WEEKEND_10, sep = "_"),
    CRASH_ID_10 = paste(YEAR_FILE, ST_CASE, sep = "_"),
    BAC_GROUP_4 = case_when(
      ANALYSIS_BAC_VALUE == 0 ~ "S",
      ANALYSIS_BAC_VALUE > 0 & ANALYSIS_BAC_VALUE < 0.05 ~ "L",
      ANALYSIS_BAC_VALUE >= 0.05 & ANALYSIS_BAC_VALUE < 0.08 ~ "M",
      ANALYSIS_BAC_VALUE >= 0.08 ~ "H",
      TRUE ~ NA_character_
    )
  )

# Collapse driver records to one row per crash and imputation. A one-vehicle
# crash needs one classified driver; a two-vehicle crash needs two.
crash_level_by_imp <- drivers %>%
  arrange(IMPUTATION, CRASH_ID_10, VEH_NO, PER_NO) %>%
  group_by(IMPUTATION, CRASH_ID_10, YEAR_FILE, ST_CASE, STATE,
           VE_FORMS, RUR_URB_GROUP_10, WEEKEND_10,
           CELL_RUR_URB_WEEKEND_10) %>%
  summarise(
    n_driver_rows = n(),
    n_bac_group_missing = sum(is.na(BAC_GROUP_4)),
    g1 = first(BAC_GROUP_4),
    g2 = nth(BAC_GROUP_4, 2, default = NA_character_),
    .groups = "drop"
  ) %>%
  mutate(
    crash_type_4 = case_when(
      VE_FORMS == 1 & n_driver_rows == 1 & n_bac_group_missing == 0 ~ g1,
      VE_FORMS == 2 & n_driver_rows == 2 & n_bac_group_missing == 0 ~ order_pair(g1, g2),
      TRUE ~ "missing_or_ambiguous"
    )
  )

write_csv(crash_level_by_imp,
          file.path(out_dir, "four_group_crash_level_by_imp.csv"))

# Count the four one-vehicle types and ten unordered two-vehicle types
# separately for each imputation and rural/urban-by-weekend cell.
cell_counts <- crash_level_by_imp %>%
  group_by(IMPUTATION, cell_id = CELL_RUR_URB_WEEKEND_10) %>%
  summarise(
    SS = sum(crash_type_4 == "SS"),
    SL = sum(crash_type_4 == "SL"),
    SM = sum(crash_type_4 == "SM"),
    SH = sum(crash_type_4 == "SH"),
    LL = sum(crash_type_4 == "LL"),
    LM = sum(crash_type_4 == "LM"),
    LH = sum(crash_type_4 == "LH"),
    MM = sum(crash_type_4 == "MM"),
    MH = sum(crash_type_4 == "MH"),
    HH = sum(crash_type_4 == "HH"),
    one_S = sum(crash_type_4 == "S"),
    one_L = sum(crash_type_4 == "L"),
    one_M = sum(crash_type_4 == "M"),
    one_H = sum(crash_type_4 == "H"),
    missing_two_vehicle = sum(VE_FORMS == 2 & crash_type_4 == "missing_or_ambiguous"),
    missing_one_vehicle = sum(VE_FORMS == 1 & crash_type_4 == "missing_or_ambiguous"),
    total_two_vehicle = sum(VE_FORMS == 2),
    total_one_vehicle = sum(VE_FORMS == 1),
    .groups = "drop"
  ) %>%
  mutate(
    n2 = SS + SL + SM + SH + LL + LM + LH + MM + MH + HH,
    n1 = one_S + one_L + one_M + one_H,
    total_complete_case = n1 + n2,
    cell_spec = "rur_urb_weekend",
    cell_var = "CELL_RUR_URB_WEEKEND_10",
    cell_spec_label = "rural/urban x weekend",
    .before = cell_id
  ) %>%
  arrange(IMPUTATION, cell_id)

# Save as CSV and RDS; the latter will be used for bootstrapping 
# in the four-group MLE
write_csv(cell_counts,
          file.path(out_dir, "four_group_crash_counts_by_cell_imp.csv"))
saveRDS(cell_counts,
        file.path(out_dir, "four_group_crash_counts_by_cell_imp.rds"))

message("Prepared ", n_distinct(crash_level_by_imp$CRASH_ID_10),
        " distinct crashes across ", n_distinct(cell_counts$IMPUTATION),
        " imputations.")
