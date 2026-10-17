# Prepare the 2015-2024 FARS nighttime fatal crash sample used in the paper.
# Keep crashes involving one or two vehicles in transport, with the number
# of driver records equal to the number of vehicles. 
# Note: Rural/urban exclusions occur in later scripts.

library(tidyverse)

years <- 2015:2024
raw_dir <- "data/raw"
processed_dir <- "data/processed"
results_dir <- "results"

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# Read only the coded fields needed for the analysis. Column capitalization
# varies across FARS years, so standardize names to uppercase.
read_fars <- function(year, file_type, variables) {
  path <- file.path(raw_dir, paste0("FARS", year, "NationalCSV"),
                    paste0(file_type, ".csv"))

  read_csv(
    path,
    col_select = matches(paste0("^(", paste(variables, collapse = "|"), ")$"),
                         ignore.case = TRUE),
    col_types = cols(.default = col_double()),
    show_col_types = FALSE
  ) %>%
    rename_with(toupper) %>%
    select(all_of(variables)) %>%
    mutate(YEAR_FILE = year, .before = 1)
}

# Build the eligible driver records one year at a time to limit memory use.
prepare_year <- function(year) {
  accidents <- read_fars(
    year, "accident",
    c("ST_CASE", "STATE", "DAY_WEEK", "HOUR", "VE_FORMS", "RUR_URB")
  ) %>%
    # Nighttime is 8:00 p.m. through 4:59 a.m.; exclude unknown hours.
    filter(HOUR %in% c(20:23, 0:4), VE_FORMS %in% c(1, 2)) %>%
    mutate(
      YEAR = year,
      # FARS codes Sunday as 1 and Saturday as 7.
      weekend_lp = DAY_WEEK %in% c(1, 7)
    )

  drivers <- read_fars(
    year, "person",
    c("ST_CASE", "VEH_NO", "PER_NO", "PER_TYP", "ALC_RES")
  ) %>%
    filter(PER_TYP == 1) %>%
    select(-PER_TYP)

  # Match the paper's sample rule: driver count equals vehicle count.
  # Preserve ALC_RES as reported; later scripts handle measured and imputed BAC.
  inner_join(drivers, accidents, by = c("YEAR_FILE", "ST_CASE")) %>%
    group_by(YEAR_FILE, ST_CASE) %>%
    filter(n() == first(VE_FORMS)) %>%
    ungroup()
}

driver_analysis_lp <- map_dfr(years, prepare_year)

# Table 1: count crashes and driver records before and after excluding
# unknown or other rural/urban classifications (only codes 1 and 2 are valid).
crashes <- driver_analysis_lp %>%
  distinct(YEAR_FILE, ST_CASE, VE_FORMS, RUR_URB)

primary_crashes <- crashes %>% filter(RUR_URB %in% c(1, 2))
excluded_crashes <- crashes %>% filter(!(RUR_URB %in% c(1, 2)))

sample_row <- function(label, data) {
  tibble(
    sample_component = label,
    crashes = nrow(data),
    driver_records = sum(data$VE_FORMS)
  )
}

sample_table <- bind_rows(
  sample_row("Initial dataset", crashes),
  sample_row("Initial: one-vehicle crashes", filter(crashes, VE_FORMS == 1)),
  sample_row("Initial: two-vehicle crashes", filter(crashes, VE_FORMS == 2)),
  sample_row("Excluded: other or unknown rural/urban status", excluded_crashes),
  sample_row("Primary dataset", primary_crashes),
  sample_row("Primary: one-vehicle crashes", filter(primary_crashes, VE_FORMS == 1)),
  sample_row("Primary: two-vehicle crashes", filter(primary_crashes, VE_FORMS == 2))
)

write_csv(driver_analysis_lp,
          file.path(processed_dir, "driver_analysis_lp_2015_2024.csv"))
write_csv(sample_table, file.path(results_dir, "sample_construction_table.csv"))

print(sample_table, n = Inf)
