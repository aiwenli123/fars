rm(list=ls())

library(tidyverse)

miacc <- read.csv("data/raw/FARS2024NationalCSV/MIACC.csv")
midrvacc <- read.csv("data/raw/FARS2024NationalCSV/MIDRVACC.csv")
miper <- read.csv("data/raw/FARS2024NationalCSV/MIPER.csv")

glimpse(miacc)
glimpse(midrvacc)
glimpse(miper)


#######################
### Notes
#######################
#
# MIACC: each A1, ..., A10 column is the imputed maximum BAC 
# among all relevant alcohol-imputed persons in the crash (so one row per crash)
#
# MIDRVACC: each A1, ..., A10 column is the imputed maximum 
# driver BAC in the crash (again, one row per crash)
#
# MIPER: one row per person involved in the crash (so each accident can have multiple rows). 
# Each P1, ..., P10 is the imputed BAC for that person.
#
# miacc:    36,297 rows
# midrvacc: 36,229 rows
# miper:    64,616 rows
#
# i think miper is best, bc it has ST_CASE, VEH_NO, and PER_NO, which are enough to 
# match with driver_lp (the main dataset from before).


#######################
### check if miper matches with driver_lp (just for 2024)
#######################

# Read in the main driver-level dataset from previous scripts
driver_lp <- read.csv("data/processed/driver_analysis_lp_2015_2024.csv")

glimpse(driver_lp)

# Check key columns exist
stopifnot(all(c("ST_CASE", "VEH_NO", "PER_NO") %in% names(driver_lp)))
stopifnot(all(c("ST_CASE", "VEH_NO", "PER_NO") %in% names(miper)))

# Define BAC column from driver_LP
bac_col <- "ALC_RES"

# Define imputed BAC columns from MIPER
p_cols <- paste0("P", 1:10)

# Since MIPER here is only from the 2024 raw FARS folder, restrict driver_lp to 2024
# before joining. Otherwise, older years can create incorrect or unmatched rows.
if ("YEAR" %in% names(driver_lp)) {
  driver_lp_2024 <- driver_lp %>%
    filter(YEAR == 2024)
} else if ("YEAR_FILE" %in% names(driver_lp)) {
  driver_lp_2024 <- driver_lp %>%
    filter(YEAR_FILE == 2024)
} else {
  stop("Could not find YEAR or YEAR_FILE in driver_lp. Please check names(driver_lp).")
}

# Merge using only 2024 driver rows
driver_miper_check <- driver_lp_2024 %>%
  left_join(
    miper,
    by = c("ST_CASE", "VEH_NO", "PER_NO")
  )

# Check merge success
merge_check <- driver_miper_check %>%
  summarise(
    n_driver_rows = n(),
    n_matched_miper = sum(!is.na(P1)),
    n_unmatched_miper = sum(is.na(P1)),
    pct_matched_miper = mean(!is.na(P1))
  )

print(merge_check)

# Next: create a clean measured BAC variable.
# NOTE: ALC_RES is in thousandths (e.g., 0.08 BAC is recorded as 80). 
# miper$P1-P10 is also in thousandths, but for some reason it's rounded down to the nearest hundredth.
# (so if actual BAC is 0.163 or 0.169, then miper will record 16.)
# So convert ALC_RES to MIPER scale using floor(ALC_RES / 10).
#
# Invalid/missing BAC codes are 995-999.

driver_miper_check <- driver_miper_check %>%
  mutate(
    measured_bac_raw = .data[[bac_col]],
    
    # Treat 0-940 as valid measured BAC values in thousandths.
    # Codes above this range are special/missing/unknown-type codes.
    measured_bac_valid = case_when(
      is.na(measured_bac_raw) ~ FALSE,
      measured_bac_raw >= 0 & measured_bac_raw <= 940 ~ TRUE,
      TRUE ~ FALSE
    ),
    
    measured_bac_miper_scale = if_else(
      measured_bac_valid,
      as.integer(floor(measured_bac_raw / 10)),
      NA_integer_
    )
  )

#######################
### row-level checks
#######################

driver_miper_check <- driver_miper_check %>%
  rowwise() %>%
  mutate(
    n_miper_values_present = sum(!is.na(c_across(all_of(p_cols)))),
    
    n_unique_miper_values = n_distinct(
      c_across(all_of(p_cols)),
      na.rm = TRUE
    ),
    
    # Do the 10 imputed values vary across imputations?
    miper_values_vary = case_when(
      n_miper_values_present < 10 ~ NA,
      n_unique_miper_values > 1 ~ TRUE,
      TRUE ~ FALSE
    ),
    
    # If measured BAC is valid, do all 10 MIPER values equal the measured BAC
    # after converting measured BAC to the MIPER scale?
    all_10_equal_measured = case_when(
      !measured_bac_valid ~ NA,
      n_miper_values_present < 10 ~ NA,
      all(c_across(all_of(p_cols)) == measured_bac_miper_scale) ~ TRUE,
      TRUE ~ FALSE
    )
  ) %>%
  ungroup()

#######################
### summary: valid measured BAC
#######################

valid_measured_summary <- driver_miper_check %>%
  filter(measured_bac_valid) %>%
  summarise(
    n_valid_measured_bac = n(),
    n_matched_miper = sum(!is.na(P1)),
    n_all_10_equal_measured = sum(all_10_equal_measured == TRUE, na.rm = TRUE),
    n_not_all_equal_measured = sum(all_10_equal_measured == FALSE, na.rm = TRUE),
    pct_all_10_equal_measured = mean(all_10_equal_measured, na.rm = TRUE)
  )

print(valid_measured_summary, width = Inf) 

# yay!
# note: gotta check for the other years