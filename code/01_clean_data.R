# 01_clean_data.R

rm(list = ls())

library(tidyverse)

# folder with raw data
raw_dir <- "data/raw"

years <- 2015:2024

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

############################################################
### variables to keep
############################################################

### Crash-level variables from accident.csv
accident_vars <- c(
  # keys / identifiers
  "STATE", "STATENAME",
  "ST_CASE",
  "YEAR",
  
  # time
  "MONTH", "MONTHNAME",
  "DAY", "DAYNAME",
  "DAY_WEEK", "DAY_WEEKNAME",
  "HOUR", "HOURNAME",
  "MINUTE", "MINUTENAME",
  
  # location
  "COUNTY", "CITY",
  
  # vehicle/person counts
  "VE_TOTAL",      # total vehicles involved
  "VE_FORMS",      # motor vehicles in transport; key for 1- vs 2-vehicle crashes
  "PVH_INVL",      # parked/working vehicles involved
  "PERSONS",       # total persons
  "PERMVIT",       # persons in motor vehicles in transport
  "PEDS",          # number of pedestrians
  "PERNOTMVIT",    # persons not in motor vehicles in transport
  "FATALS",        # fatalities in crash
  
  # crash type / roadway / environment
  "HARM_EV", "HARM_EVNAME",
  "MAN_COLL", "MAN_COLLNAME",
  "RUR_URB", "RUR_URBNAME",
  "FUNC_SYS", "FUNC_SYSNAME",
  "LGT_COND", "LGT_CONDNAME",
  "WEATHER", "WEATHERNAME"
)

### Person-level variables from person.csv
person_vars <- c(
  # keys / identifiers
  "STATE", "STATENAME",
  "ST_CASE",
  "VEH_NO",
  "PER_NO",
  "VE_FORMS",
  
  # time cross-checks
  "MONTH", "MONTHNAME",
  "DAY", "DAYNAME",
  "HOUR", "HOURNAME",
  "MINUTE", "MINUTENAME",
  
  # driver demographics / role
  "AGE", "AGENAME",
  "SEX", "SEXNAME",
  "PER_TYP", "PER_TYPNAME",
  
  # injury / occupant info
  "INJ_SEV", "INJ_SEVNAME",
  "SEAT_POS", "SEAT_POSNAME",
  "REST_USE", "REST_USENAME",
  
  # alcohol variables
  "DRINKING", "DRINKINGNAME",
  "ALC_STATUS", "ALC_STATUSNAME",
  "ATST_TYP", "ATST_TYPNAME",
  "ALC_RES", "ALC_RESNAME",
  
  # drug variables, not central for LP but useful later
  "DRUGS", "DRUGSNAME",
  "DSTATUS", "DSTATUSNAME"
)

### Vehicle-level variables from vehicle.csv
vehicle_vars <- c(
  # keys / identifiers
  "STATE", "STATENAME",
  "ST_CASE",
  "VEH_NO",
  "VE_FORMS",
  
  # basic vehicle / unit info
  "NUMOCCS", "NUMOCCSNAME",
  "UNITTYPE", "UNITTYPENAME",
  "HIT_RUN", "HIT_RUNNAME",
  "MAKE", "MAKENAME",
  "MODEL",
  "BODY_TYP", "BODY_TYPNAME",
  "MOD_YEAR", "MOD_YEARNAME",
  
  # driver presence / vehicle-level drinking
  "DR_PRES", "DR_PRESNAME",
  "DR_DRINK", "DR_DRINKNAME",
  
  # prior driving history: key for LP bad driving record
  "PREV_ACC", "PREV_ACCNAME",
  "PREV_SPD", "PREV_SPDNAME",
  "PREV_OTH", "PREV_OTHNAME",
  "PREV_DWI", "PREV_DWINAME",
  "PREV_SUS1", "PREV_SUS1NAME", # missing in some years
  "PREV_SUS2", "PREV_SUS2NAME", # missing in some years
  "PREV_SUS3", "PREV_SUS3NAME", # missing in some years
  "FIRST_MO", "FIRST_MONAME",
  "FIRST_YR", "FIRST_YRNAME",
  "LAST_MO", "LAST_MONAME",
  "LAST_YR", "LAST_YRNAME"
)

############################################################
### helper functions
############################################################

make_fars_path <- function(year, file_type) {
  file.path(
    raw_dir,
    paste0("FARS", year, "NationalCSV"),
    paste0(file_type, ".csv")
  )
}

read_fars_file <- function(year, file_type, keep_vars) {
  
  path <- make_fars_path(year, file_type)
  
  if (!file.exists(path)) {
    stop("File does not exist: ", path)
  }
  
  df <- readr::read_csv(
    path,
    show_col_types = FALSE,
    guess_max = 100000
  ) %>%
    rename_with(toupper) # make all col names uppercase because FARS is inconsistent
  
  # Some person/vehicle files may not have YEAR, so add from folder name
  if (!"YEAR" %in% names(df)) {
    df <- df %>% mutate(YEAR = year)
  }
  
  # YEAR_FILE is useful for debugging in case YEAR is missing/wrong
  df <- df %>%
    mutate(YEAR_FILE = year, .before = 1)
  
  missing_vars <- setdiff(keep_vars, names(df))
  
  if (length(missing_vars) > 0) {
    message(
      "Missing in ", file_type, " ", year, ": ",
      paste(missing_vars, collapse = ", ")
    )
  }
  
  df %>%
    select(any_of(unique(c("YEAR_FILE", keep_vars))))
}

# Helper to clean prior-record count variables.
# Anything >= 98 seems to mean unknown for these variables.
clean_prior_count <- function(x) {
  x_num <- suppressWarnings(as.numeric(x))
  case_when(
    is.na(x_num) ~ NA_real_,
    x_num %in% c(98, 99, 998) ~ NA_real_,
    TRUE ~ x_num
  )
}

############################################################
### read data
############################################################

accident_10yr <- map_dfr(
  years,
  ~ read_fars_file(.x, "accident", accident_vars)
)

person_10yr <- map_dfr(
  years,
  ~ read_fars_file(.x, "person", person_vars)
)

vehicle_10yr <- map_dfr(
  years,
  ~ read_fars_file(.x, "vehicle", vehicle_vars)
)

############################################################
### quick structure checks
############################################################

glimpse(accident_10yr)
glimpse(person_10yr)
glimpse(vehicle_10yr)

accident_10yr %>% count(YEAR)
accident_10yr %>% count(YEAR_FILE, YEAR)
person_10yr %>% count(YEAR_FILE)
vehicle_10yr %>% count(YEAR_FILE)

# Missingness overview
colSums(is.na(accident_10yr))
colSums(is.na(person_10yr))
colSums(is.na(vehicle_10yr))

############################################################
### basic data checks
############################################################

# Accident file: number of fatal crashes by year
accident_10yr %>%
  count(YEAR, name = "n_crashes")

# Person file: person types (we want PER_TYP == 1 for drivers)
person_10yr %>%
  count(PER_TYP, PER_TYPNAME, sort = TRUE) %>%
  print(n = 50)

# Police drinking variable (0/1 if known, then 8/9 if unknown)
person_10yr %>%
  count(DRINKING, DRINKINGNAME, sort = TRUE) %>%
  print(n = 50)

# Alcohol test status (0/2 if not given/given, 1 if refused, 8/9 if unknown)
person_10yr %>%
  count(ALC_STATUS, ALC_STATUSNAME, sort = TRUE) %>%
  print(n = 50)

# Alcohol test type (95, 98, 99 unknown)
person_10yr %>%
  count(ATST_TYP, ATST_TYPNAME, sort = TRUE) %>%
  print(n = 50)

# BAC result examples (999, 997-995 are unknown. 998 = positive eading but no actual value)
person_10yr %>%
  count(ALC_RES, ALC_RESNAME, sort = TRUE) %>%
  arrange(desc(ALC_RES)) %>%
  print(n = 20)

# Driver presence
vehicle_10yr %>%
  count(DR_PRES, DR_PRESNAME, sort = TRUE) %>%
  print(n = 50)

# Vehicle-level drinking variable (only 0/1, no unknowns but will dive deeper into this later)
vehicle_10yr %>%
  count(DR_DRINK, DR_DRINKNAME, sort = TRUE) %>%
  print(n = 50)

# Prior-record variables: important for LP bad-record construction
prior_record_vars <- c(
  "PREV_ACC", "PREV_SPD", "PREV_OTH", 
  "PREV_DWI", "PREV_SUS1", "PREV_SUS2", "PREV_SUS3"
)

for (v in prior_record_vars) {
  message("\nCounts for ", v, ":")
  print(
    vehicle_10yr %>%
      count(.data[[v]], sort = TRUE) %>%
      print(n = 50)
  )
}

############################################################
### file for driver-level analysis
############################################################

unique(accident_10yr$HOUR)

# Accident-level core sample flags
accident_core <- accident_10yr %>%
  mutate(
    # 8pm through 4:59am, inclusive
    night_lp = (HOUR >= 20 | HOUR <= 4) & (HOUR != 99),
    
    # one- or two-motor-vehicle-in-transport crashes
    one_two_vehicle = VE_FORMS %in% c(1, 2),
    
    # weekend definition:
    # Monday = 2, Tuesday = 3, ..., Saturday = 7, Sunday = 1.
    weekend_lp = DAY_WEEK %in% c(1, 7),
    
    covid_period = case_when(
      YEAR %in% 2015:2019 ~ "2015-2019 pre-COVID",
      YEAR %in% 2020:2021 ~ "2020-2021 COVID",
      YEAR %in% 2022:2024 ~ "2022-2024 post-COVID",
      TRUE ~ NA_character_
    ),
    
    # cells useful later for LP Table 2-style specifications
    cell_hour = interaction(HOUR, drop = TRUE),
    cell_hour_year = interaction(HOUR, YEAR, drop = TRUE),
    cell_hour_year_weekend = interaction(HOUR, YEAR, weekend_lp, drop = TRUE),
    cell_hour_state_year = interaction(HOUR, STATE, YEAR, drop = TRUE),
    cell_hour_state_year_weekend = interaction(HOUR, STATE, YEAR, weekend_lp, drop = TRUE)
  )

unique(person_10yr$SEXNAME)
sort(unique(person_10yr$AGE), decreasing = TRUE)

# Driver rows from person file
drivers_person <- person_10yr %>%
  filter(PER_TYP == 1) %>%
  select(
    YEAR_FILE,
    ST_CASE,
    VEH_NO,
    PER_NO,
    AGE, AGENAME,
    SEX, SEXNAME,
    PER_TYP, PER_TYPNAME,
    INJ_SEV, INJ_SEVNAME,
    SEAT_POS, SEAT_POSNAME,
    REST_USE, REST_USENAME,
    DRINKING, DRINKINGNAME,
    ALC_STATUS, ALC_STATUSNAME,
    ATST_TYP, ATST_TYPNAME,
    ALC_RES, ALC_RESNAME,
    DRUGS, DRUGSNAME,
    DSTATUS, DSTATUSNAME
  ) %>%
  mutate(
    # LP age/gender variables
    # Treat FARS unknown age codes as missing.
    age_clean = case_when(
      is.na(AGE) ~ NA_real_,
      AGE >= 0 & AGE < 998 ~ as.numeric(AGE),
      TRUE ~ NA_real_
    ),
    young_driver_lp = case_when(
      !is.na(age_clean) & age_clean < 25 ~ 1,
      !is.na(age_clean) & age_clean >= 25 ~ 0,
      TRUE ~ NA_real_
    ),
    
    # 1 is Male, 2 is Female
    male_driver_lp = case_when(
      SEX == 1 ~ 1,
      SEX == 2 ~ 0,
      TRUE ~ NA_real_
    ),
    
    # Police drinking definition
    police_drinking_complete = case_when(
      DRINKING == 1 ~ 1,
      DRINKING == 0 ~ 0,
      TRUE ~ NA_real_
    ),
    
    # Alternative: unknown/not reported treated as not drinking
    police_drinking_unknown_as_no = case_when(
      DRINKING == 1 ~ 1,
      TRUE ~ 0
    ),
    
    # Valid observed BAC.
    # FARS ALC_RES appears as BAC * 1000, e.g., 80 = 0.080.
    # Special codes 995-999 are invalid / missing-style values.
    bac_valid = !is.na(ALC_RES) & ALC_RES >= 0 & ALC_RES < 995,
    bac_tested = bac_valid,
    bac_value = if_else(bac_valid, ALC_RES / 1000, NA_real_),
    
    bac_ge_08 = case_when(
      bac_valid & ALC_RES >= 80 ~ 1,
      bac_valid & ALC_RES < 80 ~ 0,
      TRUE ~ NA_real_
    ),
    
    bac_ge_10 = case_when(
      bac_valid & ALC_RES >= 100 ~ 1,
      bac_valid & ALC_RES < 100 ~ 0,
      TRUE ~ NA_real_
    )
  )

# Vehicle variables to join onto driver rows
vehicle_driver_vars <- vehicle_10yr %>%
  select(
    YEAR_FILE,
    ST_CASE,
    VEH_NO,
    
    # basic vehicle / unit info
    NUMOCCS, NUMOCCSNAME,
    UNITTYPE, UNITTYPENAME,
    HIT_RUN, HIT_RUNNAME,
    MAKE, MAKENAME,
    MODEL,
    BODY_TYP, BODY_TYPNAME,
    MOD_YEAR, MOD_YEARNAME,
    
    # driver presence / vehicle-level drinking
    DR_PRES, DR_PRESNAME,
    DR_DRINK, DR_DRINKNAME,

    # prior driving history variables
    PREV_ACC, PREV_ACCNAME,
    PREV_SPD, PREV_SPDNAME,
    PREV_OTH, PREV_OTHNAME,
    PREV_DWI, PREV_DWINAME,
    PREV_SUS1, PREV_SUS1NAME,
    PREV_SUS2, PREV_SUS2NAME,
    PREV_SUS3, PREV_SUS3NAME,
    FIRST_MO, FIRST_MONAME,
    FIRST_YR, FIRST_YRNAME,
    LAST_MO, LAST_MONAME,
    LAST_YR, LAST_YRNAME,
  ) %>%
  distinct(YEAR_FILE, ST_CASE, VEH_NO, .keep_all = TRUE) %>%
  mutate(
    # cleaned prior-record count variables
    prev_acc_n  = clean_prior_count(PREV_ACC),
    prev_spd_n  = clean_prior_count(PREV_SPD),
    prev_oth_n  = clean_prior_count(PREV_OTH),
    prev_dwi_n  = clean_prior_count(PREV_DWI),
    prev_sus1_n = clean_prior_count(PREV_SUS1),
    prev_sus2_n = clean_prior_count(PREV_SUS2),
    prev_sus3_n = clean_prior_count(PREV_SUS3),
    
    # LP-style minor blemishes:
    # previous crashes + previous speeding convictions + previous other moving violations
    minor_blemishes_lp = rowSums(
      across(c(prev_acc_n, prev_spd_n, prev_oth_n)),
      na.rm = TRUE
    ),
    
    # whether all minor variables were missing
    minor_all_missing = if_all(
      c(prev_acc_n, prev_spd_n, prev_oth_n),
      ~ is.na(.x)
    ),
    
    # LP-style major blemish: prior DWI conviction
    major_blemish_lp = case_when(
      if_any(
        c(prev_dwi_n, prev_sus1_n, prev_sus2_n, prev_sus3_n),
        ~ !is.na(.x) & .x > 0
      ) ~ 1,
      if_all(
        c(prev_dwi_n),
        ~ !is.na(.x)
      ) ~ 0,
      TRUE ~ NA_real_
    ),
    
    # Complete prior-record info flag
    prior_record_complete = if_all(
      c(prev_acc_n, prev_spd_n, prev_oth_n, 
        prev_dwi_n, prev_sus1_n, prev_sus2_n, prev_sus3_n),
      ~ !is.na(.x)
    ),
    
    # LP-style bad previous driving record:
    # 2+ minor blemishes OR 1+ major blemish
    bad_record_lp = case_when(
      minor_blemishes_lp >= 2 ~ 1,
      major_blemish_lp == 1 ~ 1,
      prior_record_complete ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  select(-minor_all_missing)

# Join crash-level and vehicle-level info onto driver rows
driver_analysis <- drivers_person %>%
  left_join(
    accident_core %>%
      select(
        YEAR_FILE, ST_CASE,
        YEAR, STATE, STATENAME,
        MONTH, MONTHNAME,
        DAY, DAYNAME,
        DAY_WEEK, DAY_WEEKNAME,
        HOUR, HOURNAME,
        MINUTE, MINUTENAME,
        COUNTY, CITY,
        VE_TOTAL, VE_FORMS,
        PVH_INVL, PERSONS, PERMVIT, PEDS, PERNOTMVIT,
        FATALS,
        night_lp, one_two_vehicle, weekend_lp, covid_period,
        cell_hour,
        cell_hour_year,
        cell_hour_year_weekend,
        cell_hour_state_year,
        cell_hour_state_year_weekend,
        HARM_EV, HARM_EVNAME,
        MAN_COLL, MAN_COLLNAME,
        RUR_URB, RUR_URBNAME,
        FUNC_SYS, FUNC_SYSNAME,
        LGT_COND, LGT_CONDNAME,
        WEATHER, WEATHERNAME,
      ),
    by = c("YEAR_FILE", "ST_CASE")
  ) %>%
  left_join(
    vehicle_driver_vars,
    by = c("YEAR_FILE", "ST_CASE", "VEH_NO")
  )

############################################################
### checks after joins
############################################################

glimpse(driver_analysis)

driver_analysis %>%
  count(YEAR, night_lp, one_two_vehicle)

driver_analysis %>%
  count(YEAR, VE_FORMS, sort = TRUE)

driver_analysis %>%
  count(YEAR, weekend_lp, sort = TRUE)

# Check merge coverage from vehicle file
driver_analysis %>%
  summarize(
    n_driver_rows = n(),
    n_missing_vehicle_match = sum(is.na(DR_PRES) & is.na(DR_DRINK)),
    pct_missing_vehicle_match = mean(is.na(DR_PRES) & is.na(DR_DRINK))
  )

# Check derived variables
driver_analysis %>%
  count(young_driver_lp, sort = TRUE)

driver_analysis %>%
  count(male_driver_lp, sort = TRUE)

driver_analysis %>%
  count(bad_record_lp, sort = TRUE)

driver_analysis %>%
  count(prior_record_complete, bad_record_lp, sort = TRUE)

driver_analysis %>%
  summarize(
    n = n(),
    pct_police_drinking_missing = mean(is.na(police_drinking_complete)),
    pct_bac_valid = mean(bac_valid, na.rm = TRUE),
    pct_age_missing = mean(is.na(age_clean)),
    pct_sex_missing = mean(is.na(male_driver_lp)),
    pct_bad_record_missing = mean(is.na(bad_record_lp))
  ) %>% 
  print(width = Inf)

############################################################
### driver counts by crash
############################################################

driver_counts_by_crash <- driver_analysis %>%
  group_by(YEAR, ST_CASE) %>%
  summarize(
    n_driver_rows = n(),
    VE_FORMS = first(VE_FORMS),
    night_lp = first(night_lp),
    one_two_vehicle = first(one_two_vehicle),
    weekend_lp = first(weekend_lp),
    covid_period = first(covid_period),
    .groups = "drop"
  )

driver_counts_by_crash %>%
  count(VE_FORMS, n_driver_rows, sort = TRUE)

# Clean Levitt-Porter style crash set
lp_clean_crashes <- driver_counts_by_crash %>%
  filter(
    night_lp,
    one_two_vehicle,
    n_driver_rows == VE_FORMS
  )

lp_clean_crashes %>%
  count(YEAR, VE_FORMS)

# Now make one row per driver, so 2 rows for 2-car crashes
driver_analysis_lp <- driver_analysis %>%
  semi_join(
    lp_clean_crashes %>% select(YEAR, ST_CASE),
    by = c("YEAR", "ST_CASE")
  )

driver_analysis_lp %>%
  count(YEAR, VE_FORMS)

driver_analysis_lp %>%
  count(YEAR, VE_FORMS, weekend_lp) %>% 
  print(n=40)

driver_analysis_lp %>%
  summarize(
    n_driver_rows = n(),
    n_crashes = n_distinct(paste(YEAR, ST_CASE)),
    pct_police_drinking_missing = mean(is.na(police_drinking_complete)),
    pct_bac_valid = mean(bac_valid, na.rm = TRUE),
    pct_bad_record_missing = mean(is.na(bad_record_lp)),
    pct_young_missing = mean(is.na(young_driver_lp)),
    pct_male_missing = mean(is.na(male_driver_lp))
  ) %>%
  print(width = Inf)

############################################################
### BAC completeness by state-year
############################################################

# This prepares for later BAC >= .08 / .10 analysis.
# LP used a restricted BAC sample based on high testing completeness among police-reported drinkers.
bac_completeness_state_year <- driver_analysis_lp %>%
  group_by(YEAR, STATE, STATENAME) %>%
  summarize(
    n_drivers = n(),
    n_police_drinking_yes = sum(police_drinking_complete == 1, na.rm = TRUE),
    n_police_drinking_yes_bac_valid = sum(police_drinking_complete == 1 & bac_valid, na.rm = TRUE),
    bac_complete_among_police_drinkers = case_when(
      n_police_drinking_yes > 0 ~ n_police_drinking_yes_bac_valid / n_police_drinking_yes,
      TRUE ~ NA_real_
    ),
    bac_complete_80 = bac_complete_among_police_drinkers >= 0.80,
    bac_complete_90 = bac_complete_among_police_drinkers >= 0.90,
    bac_complete_95 = bac_complete_among_police_drinkers >= 0.95,
    .groups = "drop"
  )

bac_completeness_state_year %>%
  arrange(YEAR, STATE) %>%
  print(n = 20, width = Inf)

# Add BAC-completeness flags back to driver-level LP data
driver_analysis_lp <- driver_analysis_lp %>%
  left_join(
    bac_completeness_state_year %>%
      select(
        YEAR, STATE,
        bac_complete_among_police_drinkers,
        bac_complete_80,
        bac_complete_90,
        bac_complete_95
      ),
    by = c("YEAR", "STATE")
  )

############################################################
### save files
############################################################

write_csv(accident_10yr, "data/processed/accident_2015_2024_core.csv")
write_csv(person_10yr, "data/processed/person_2015_2024_core.csv")
write_csv(vehicle_10yr, "data/processed/vehicle_2015_2024_core.csv")

write_csv(driver_analysis, "data/processed/driver_analysis_2015_2024.csv") # one row per driver
write_csv(driver_analysis_lp, "data/processed/driver_analysis_lp_2015_2024.csv") # still one row per driver, only keeps nighttime crashes, 1 or 2 cars. prolly the main file for later
write_csv(lp_clean_crashes, "data/processed/lp_clean_crashes_2015_2024.csv") # one row per crash
write_csv(bac_completeness_state_year, "data/processed/bac_completeness_state_year_2015_2024.csv")


############################################################################

#####################################
### check person$DRINKING vs. vehicle$DR_DRINK (this section is mainly just exploratory)
#####################################

# Keep driver rows from person file
drivers_person_check <- person_10yr %>%
  filter(PER_TYP == 1) %>%
  select(
    YEAR_FILE,
    ST_CASE,
    VEH_NO,
    PER_NO,
    PER_TYP,
    PER_TYPNAME,
    DRINKING,
    DRINKINGNAME,
    ALC_STATUS, ALC_STATUSNAME,
    ATST_TYP, ATST_TYPNAME,
    ALC_RES, ALC_RESNAME
  )

# Keep vehicle-level driver drinking variable if it exists
vehicle_drink_check <- vehicle_10yr %>%
  select(
    YEAR_FILE,
    ST_CASE,
    VEH_NO,
    any_of(c("DR_DRINK", "DR_DRINKNAME"))
  )

# Join person driver rows to vehicle rows
drink_compare <- drivers_person_check %>%
  left_join(
    vehicle_drink_check,
    by = c("YEAR_FILE", "ST_CASE", "VEH_NO")
  )

glimpse(drink_compare)

# See all combinations
drink_compare %>%
  count(DRINKING, DRINKINGNAME, DR_DRINK, DR_DRINKNAME, sort = TRUE) %>%
  print(n = 50)

### comparisons
drink_compare <- drink_compare %>%
  mutate(
    drink_match = case_when(
      is.na(DRINKING) & is.na(DR_DRINK) ~ TRUE,
      DRINKING == DR_DRINK ~ TRUE,
      TRUE ~ FALSE
    )
  )

drink_compare %>%
  count(drink_match)

# check by year
drink_compare %>%
  group_by(YEAR_FILE) %>%
  summarize(
    n_driver_rows = n(),
    n_missing_person_drinking = sum(is.na(DRINKING)),
    n_missing_vehicle_dr_drink = sum(is.na(DR_DRINK)),
    n_matches = sum(drink_match, na.rm = TRUE),
    n_mismatches = sum(!drink_match, na.rm = TRUE),
    match_rate = mean(drink_match, na.rm = TRUE),
    .groups = "drop"
  )

names(drink_compare)

drink_compare %>%
  filter(!drink_match) %>%
  select(
    YEAR_FILE, ST_CASE, VEH_NO, PER_NO,
    DRINKING, DRINKINGNAME,
    DR_DRINK, DR_DRINKNAME,
    ALC_STATUS, ALC_STATUSNAME,
    ATST_TYP, ATST_TYPNAME,
    ALC_RES, ALC_RESNAME
  ) %>%
  print(n = 100, width = Inf)

drink_compare %>%
  filter(!drink_match) %>%
  mutate(
    bac_valid = !is.na(ALC_RES) & ALC_RES >= 0 & ALC_RES < 995,
    bac_positive = bac_valid & ALC_RES > 0,
    bac_ge_08 = bac_valid & ALC_RES >= 80,
    bac_ge_10 = bac_valid & ALC_RES >= 100
  ) %>%
  count(
    DRINKING, DRINKINGNAME,
    DR_DRINK, DR_DRINKNAME,
    bac_valid, bac_positive, bac_ge_08,
    sort = TRUE
  ) %>%
  print(n = 100, length = Inf)


# check if join went ok
drivers_person_check %>%
  count(YEAR_FILE, ST_CASE, VEH_NO) %>%
  filter(n > 1)

vehicle_drink_check %>%
  count(YEAR_FILE, ST_CASE, VEH_NO) %>%
  filter(n > 1)

### save mismatch files for later reference if needed
drink_mismatch_summary <- drink_compare %>%
  filter(!drink_match) %>%
  mutate(
    bac_valid = !is.na(ALC_RES) & ALC_RES >= 0 & ALC_RES < 995,
    bac_positive = bac_valid & ALC_RES > 0,
    bac_ge_08 = bac_valid & ALC_RES >= 80,
    bac_ge_10 = bac_valid & ALC_RES >= 100
  ) %>%
  count(
    DRINKING, DRINKINGNAME,
    DR_DRINK, DR_DRINKNAME,
    bac_valid, bac_positive, bac_ge_08, bac_ge_10,
    sort = TRUE
  )

write_csv(
  drink_mismatch_summary,
  "data/processed/drinking_vs_dr_drink_mismatch_summary_2015_2024.csv"
)

drink_match_by_year <- drink_compare %>%
  mutate(
    drink_match = DRINKING == DR_DRINK
  ) %>%
  group_by(YEAR_FILE) %>%
  summarize(
    n_driver_rows = n(),
    n_matches = sum(drink_match, na.rm = TRUE),
    n_mismatches = sum(!drink_match, na.rm = TRUE),
    match_rate = mean(drink_match, na.rm = TRUE),
    person_drinking_yes_rate = mean(DRINKING == 1, na.rm = TRUE),
    dr_drink_yes_rate = mean(DR_DRINK == 1, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  drink_match_by_year,
  "data/processed/drinking_vs_dr_drink_match_by_year_2015_2024.csv"
)


# library(tidyverse)
# 
# driver_analysis_lp <- read_csv("data/processed/driver_analysis_lp_2015_2024.csv")
# names(driver_analysis_lp)
# 
# dim(driver_analysis_lp)
# table(driver_analysis_lp$ALC_RES)
# table(driver_analysis_lp$ALC_RESNAME)
