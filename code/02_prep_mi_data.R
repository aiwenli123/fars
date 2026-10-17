# Add NHTSA's ten imputed BAC values (MIPER) to the nighttime FARS driver
# dataset from Script 01. Keep measured BAC when valid, use an imputed value
# otherwise, and create one driver row per imputation. The resulting CSV is
# used by the two-group benchmarks and main four-group model.

library(tidyverse)

years <- 2015:2024
raw_dir <- "data/raw"
processed_dir <- "data/processed"
p_cols <- paste0("P", 1:10)
keys <- c("YEAR_FILE", "ST_CASE", "VEH_NO", "PER_NO")

# Read the driver-level sample prepared by Script 01.
drivers <- read_csv(
  file.path(processed_dir, "driver_analysis_lp_2015_2024.csv"),
  show_col_types = FALSE
)

# Locate and read the MIPER file, allowing capitalization to vary by year.
read_miper <- function(year) {
  year_dir <- file.path(raw_dir, paste0("FARS", year, "NationalCSV"))
  paths <- list.files(year_dir, pattern = "^miper\\.csv$",
                      ignore.case = TRUE, full.names = TRUE)
  if (length(paths) != 1) stop("Expected one MIPER file for ", year, ": ", year_dir)

  read_csv(
    paths[1],
    col_select = matches("^(ST_CASE|VEH_NO|PER_NO|P([1-9]|10))$", ignore.case = TRUE),
    col_types = cols(.default = col_double()),
    show_col_types = FALSE
  ) %>%
    rename_with(toupper) %>%
    select(ST_CASE, VEH_NO, PER_NO, all_of(p_cols)) %>%
    mutate(YEAR_FILE = year, .before = 1)
}

# Join and reshape one year at a time so the ten imputed datasets do not
# need to be stored in memory as separate full-year data frames.
prepare_year <- function(year) {
  year_drivers <- filter(drivers, YEAR_FILE == year)
  miper <- read_miper(year)

  # Each driver should have at most one MIPER record; otherwise the join
  # would duplicate drivers and change the crash composition.
  if (anyDuplicated(miper[keys])) stop("Duplicate MIPER driver keys in ", year)

  wide <- left_join(year_drivers, miper, by = keys)
  if (nrow(wide) != nrow(year_drivers)) stop("Driver count changed after MIPER join in ", year)

  wide %>%
    mutate(
      # ALC_RES records BAC in thousandths (80 means 0.080).
      # Retain the original script's measured-BAC validity range.
      measured_bac_valid = !is.na(ALC_RES) & between(ALC_RES, 0, 940)
    ) %>%
    pivot_longer(all_of(p_cols), names_to = "imputation",
                 values_to = "miper_bac_hundredths") %>%
    mutate(
      imputation = as.integer(sub("P", "", imputation)),

      # MIPER stores BAC in hundredths (8 means 0.08). Measured BAC
      # takes precedence and is unchanged across the ten imputations.
      analysis_bac_value = if_else(
        measured_bac_valid, ALC_RES / 1000, miper_bac_hundredths / 100
      ),
      analysis_bac_source = case_when(
        measured_bac_valid ~ "measured",
        !is.na(miper_bac_hundredths) ~ "multiple_imputation",
        TRUE ~ "missing_no_miper"
      ),

      # Indicators used by the two-group BAC benchmarks.
      mi_bac_ge_05 = if_else(is.na(analysis_bac_value), NA_integer_,
                             as.integer(analysis_bac_value >= 0.05)),
      mi_bac_ge_08 = if_else(is.na(analysis_bac_value), NA_integer_,
                             as.integer(analysis_bac_value >= 0.08))
    ) %>%
    select(-miper_bac_hundredths)
}

# Combine the completed datasets while retaining the original driver and
# crash variables needed for grouping by rural/urban status and weekend.
driver_miper_long <- map_dfr(years, prepare_year)

# Counts of drivers with valid BAC
measured_summary <- drivers %>%
  summarise(
    drivers = n(),
    valid_measured_bac = sum(!is.na(ALC_RES) & between(ALC_RES, 0, 940)),
    missing_or_invalid_bac = drivers - valid_measured_bac
  )
print(measured_summary)

# Save MI dataset
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
write_csv(driver_miper_long,
          file.path(processed_dir, "driver_miper_long_2015_2024.csv"))

message("Saved ", nrow(driver_miper_long), " driver-imputation rows.")
