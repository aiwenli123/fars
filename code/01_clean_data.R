rm(list = ls())

library(tidyverse)

# folder with raw data
raw_dir <- "data/raw"

years <- 2015:2024

### Crash-level variables from accident.csv
accident_vars <- c(
  "STATE", "STATENAME",
  "ST_CASE",
  "YEAR",
  "MONTH", "DAY", "DAY_WEEK",
  "HOUR", "MINUTE",
  "VE_TOTAL",      # total vehicles involved
  "VE_FORMS",      # motor vehicles in transport; key for 1- vs 2-vehicle crashes
  "PVH_INVL",      # parked/working vehicles involved
  "PERSONS",       # total persons
  "PERMVIT",       # persons in motor vehicles in transport
  "PEDS",          # number of pedestrians
  "PERNOTMVIT",    # persons not in motor vehicles in transport
  "FATALS",        # fatalities in crash
  "HARM_EV", "HARM_EVNAME",
  "MAN_COLL", "MAN_COLLNAME",
  "RUR_URB", "RUR_URBNAME",
  "FUNC_SYS", "FUNC_SYSNAME",
  "LGT_COND", "LGT_CONDNAME",
  "WEATHER", "WEATHERNAME"
)

# Person-level variables from person.csv
person_vars <- c(
  "STATE", "STATENAME",
  "ST_CASE",
  "VEH_NO",
  "PER_NO",
  "VE_FORMS",
  "MONTH", "DAY",
  "HOUR", "MINUTE",
  "AGE", "AGENAME",
  "SEX", "SEXNAME",
  "PER_TYP", "PER_TYPNAME",
  "INJ_SEV", "INJ_SEVNAME",
  "SEAT_POS", "SEAT_POSNAME",
  "REST_USE", "REST_USENAME",
  "DRINKING", "DRINKINGNAME",
  "ALC_STATUS", "ALC_STATUSNAME",
  "ATST_TYP", "ATST_TYPNAME",
  "ALC_RES", "ALC_RESNAME",
  "DRUGS", "DRUGSNAME",
  "DSTATUS", "DSTATUSNAME"
)

# Vehicle-level variables from vehicle.csv
vehicle_vars <- c(
  "STATE", "STATENAME",
  "ST_CASE",
  "VEH_NO",
  "VE_FORMS",
  "NUMOCCS", "NUMOCCSNAME",
  "UNITTYPE", "UNITTYPENAME",
  "HIT_RUN", "HIT_RUNNAME",
  "MAKE", "MAKENAME",
  "MODEL",
  "BODY_TYP", "BODY_TYPNAME",
  "MOD_YEAR", "MOD_YEARNAME",
  "VPICMAKE", "VPICMAKENAME",
  "VPICMODEL", "VPICMODELNAME",
  "VPICBODYCLASS", "VPICBODYCLASSNAME",
  "GVWR_FROM", "GVWR_FROMNAME",
  "GVWR_TO", "GVWR_TONAME",
  "TRAV_SP", "TRAV_SPNAME",
  "ROLLOVER", "ROLLOVERNAME",
  "M_HARM", "M_HARMNAME",
  "DR_PRES", "DR_PRESNAME",
  "DR_DRINK", "DR_DRINKNAME"  # may or may not exist depending on year/file
)

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
    rename_with(toupper)
  
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

### read data
accident_10yr <- map_dfr(
  years,
  ~ read_fars_file(.x, "accident", accident_vars)
)
table(accident_10yr$YEAR)

person_10yr <- map_dfr(
  years,
  ~ read_fars_file(.x, "person", person_vars)
)

vehicle_10yr <- map_dfr(
  years,
  ~ read_fars_file(.x, "vehicle", vehicle_vars)
) # there are some missing columns in 2015-2019 but none of them are super important for replicating levitt-porter

glimpse(accident_10yr)
glimpse(person_10yr)
glimpse(vehicle_10yr)
vehicle_10yr %>%
  select(YEAR_FILE, ST_CASE, VEH_NO, VE_FORMS, UNITTYPE, BODY_TYP, DR_PRES, DR_DRINK) %>%
  glimpse()

accident_10yr %>% count(YEAR)
accident_10yr %>%
  count(YEAR_FILE, YEAR)
person_10yr %>% count(YEAR_FILE)
vehicle_10yr %>% count(YEAR_FILE)

# the ones that are missing in older vehicle datasets
vehicle_vars_optional <- c(
  "VPICMAKE", "VPICMAKENAME",
  "VPICMODEL", "VPICMODELNAME",
  "VPICBODYCLASS", "VPICBODYCLASSNAME",
  "GVWR_FROM", "GVWR_FROMNAME",
  "GVWR_TO", "GVWR_TONAME"
)

# remove these columns for now
vehicle_10yr <- vehicle_10yr %>%
  select(-any_of(vehicle_vars_optional))

colSums(is.na(accident_10yr))
colSums(is.na(person_10yr))
colSums(is.na(vehicle_10yr))

### checking data
# Accident file: number of fatal crashes by year
accident_10yr %>%
  count(YEAR, name = "n_crashes")

# Person file: person types
person_10yr %>%
  count(PER_TYP, PER_TYPNAME, sort = TRUE)

# Police drinking variable
person_10yr %>%
  count(DRINKING, DRINKINGNAME, sort = TRUE)

# Alcohol test status
person_10yr %>%
  count(ALC_STATUS, ALC_STATUSNAME, sort = TRUE)

# Alcohol test type
person_10yr %>%
  count(ATST_TYP, ATST_TYPNAME, sort = TRUE)

# BAC result examples
person_10yr %>%
  count(ALC_RES, ALC_RESNAME, sort = TRUE) %>%
  arrange((ALC_RES))

# Vehicle file: unit types
vehicle_10yr %>%
  count(UNITTYPE, UNITTYPENAME, sort = TRUE) # only one is 1: Motor Vehicle In-Transport (Inside or Outside the Trafficway)


#####################################
### file for driver-level analysis
#####################################


# Accident-level core sample flags
accident_core <- accident_10yr %>%
  mutate(
    night_lp = HOUR >= 20 | HOUR <= 4, # 8pm - 4:59am
    one_two_vehicle = VE_FORMS %in% c(1, 2),
    covid_period = case_when(
      YEAR %in% 2015:2019 ~ "2015-2019 pre-COVID",
      YEAR %in% 2020:2021 ~ "2020-2021 COVID",
      YEAR %in% 2022:2024 ~ "2022-2024 post-COVID",
      TRUE ~ NA_character_
    )
  )

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
    # Assuming ALC_RES is BAC * 1000 and special codes 995-999 are invalid.
    bac_valid = !is.na(ALC_RES) & ALC_RES >= 0 & ALC_RES < 995,
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

# Join crash-level info onto driver rows
driver_analysis <- drivers_person %>%
  left_join(
    accident_core %>%
      select(
        YEAR_FILE, ST_CASE,
        YEAR, STATE, STATENAME,
        MONTH, DAY, DAY_WEEK, HOUR, MINUTE,
        VE_TOTAL, VE_FORMS,
        PVH_INVL, PERSONS, PERMVIT, PEDS, PERNOTMVIT,
        FATALS,
        night_lp, one_two_vehicle, covid_period,
        HARM_EV, HARM_EVNAME,
        MAN_COLL, MAN_COLLNAME,
        RUR_URB, RUR_URBNAME,
        LGT_COND, LGT_CONDNAME,
        WEATHER, WEATHERNAME
      ),
    by = c("YEAR_FILE", "ST_CASE")
  )

glimpse(driver_analysis)

driver_analysis %>%
  count(YEAR, night_lp, one_two_vehicle)

driver_analysis %>%
  count(YEAR, VE_FORMS, sort = TRUE)

driver_analysis %>%
  count(YEAR, night_lp, one_two_vehicle)

### driver counts by crash
driver_counts_by_crash <- driver_analysis %>%
  group_by(YEAR, ST_CASE) %>%
  summarize(
    n_driver_rows = n(),
    VE_FORMS = first(VE_FORMS),
    night_lp = first(night_lp),
    one_two_vehicle = first(one_two_vehicle),
    covid_period = first(covid_period),
    .groups = "drop"
  )

driver_counts_by_crash %>%
  count(VE_FORMS, n_driver_rows, sort = TRUE)

# clean levitt-porter style crash set
lp_clean_crashes <- driver_counts_by_crash %>%
  filter(
    night_lp,
    one_two_vehicle,
    n_driver_rows == VE_FORMS
  )

lp_clean_crashes %>%
  count(YEAR, VE_FORMS)

# now make one row per driver, so 2x rows for 2-car crashes
driver_analysis_lp <- driver_analysis %>%
  semi_join(
    lp_clean_crashes %>% select(YEAR, ST_CASE),
    by = c("YEAR", "ST_CASE")
  )

driver_analysis_lp %>%
  count(YEAR, VE_FORMS)

### save files
write_csv(accident_10yr, "data/processed/accident_2015_2024_core.csv")
write_csv(person_10yr, "data/processed/person_2015_2024_core.csv")
write_csv(vehicle_10yr, "data/processed/vehicle_2015_2024_core.csv")
write_csv(driver_analysis_lp, "data/processed/driver_analysis_lp_2015_2024.csv")
