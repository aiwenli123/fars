# 00_read_data.R

rm(list=ls())

library(tidyverse)

person <- read.csv("data/raw/FARS2024NationalCSV/person.csv")
glimpse(person)

vehicle <- read.csv("data/raw/FARS2024NationalCSV/vehicle.csv")
glimpse(vehicle)

accident <- read.csv("data/raw/FARS2024NationalCSV/accident.csv")
glimpse(accident)

unique(accident$YEAR)

# manually created data dictionary
var_reference <- tribble(
  ~dataset, ~raw_column, ~meaning_for_replication, ~why_we_need_it,
  
  ############################################################
  # accident file: crash-level variables
  ############################################################
  
  "accident", "STATE", "State numeric code", "Used for state-year cells and possible state-level comparisons",
  "accident", "STATENAME", "State name", "Useful for readable summaries",
  "accident", "ST_CASE", "Unique crash identifier within year/state", "Primary key for joining accident, person, and vehicle files",
  "accident", "YEAR", "Crash year", "Needed for 2015-2024 replication, year-specific estimates, and COVID-period splits",
  
  # Time variables
  "accident", "MONTH", "Crash month", "Potential time control; useful for descriptive checks",
  "accident", "MONTHNAME", "Crash month label", "Readable version of MONTH if available",
  "accident", "DAY", "Day of month", "Potential time control; mostly for checking",
  "accident", "DAYNAME", "Day of month label", "Readable version of DAY if available",
  "accident", "DAY_WEEK", "Day of week", "Used to define weekend/night patterns for LP-style cells",
  "accident", "DAY_WEEKNAME", "Day of week label", "Readable version of DAY_WEEK; useful for verifying weekend definition",
  "accident", "HOUR", "Hour of crash", "Needed to restrict to Levitt-Porter nighttime window, roughly 8pm-5am",
  "accident", "HOURNAME", "Hour of crash label", "Readable version of HOUR; useful for checking unknown hour codes",
  "accident", "MINUTE", "Minute of crash", "Mostly for checking; not essential for first pass",
  "accident", "MINUTENAME", "Minute of crash label", "Readable version of MINUTE; useful for checking unknown minute codes",
  
  # Location variables
  "accident", "COUNTY", "County numeric code", "Useful for location checks and possible finer geographic robustness",
  "accident", "CITY", "City numeric code", "Useful for location checks and possible urban/local robustness",
  
  # Vehicle/person counts
  "accident", "VE_TOTAL", "Total vehicles involved", "Useful cross-check for one- and two-vehicle crash definitions",
  "accident", "VE_FORMS", "Number of motor vehicles in transport", "Primary variable for defining one-vehicle and two-vehicle crashes",
  "accident", "PVH_INVL", "Number of parked/working vehicles involved", "Useful for distinguishing motor vehicles in transport from parked/working vehicles",
  "accident", "PERSONS", "Number of persons involved", "Descriptive check",
  "accident", "PERMVIT", "Persons in motor vehicles in transport", "Descriptive check",
  "accident", "PEDS", "Number of pedestrians", "Useful for checking nonoccupant crashes",
  "accident", "PERNOTMVIT", "Persons not in motor vehicles in transport", "Useful for checking pedestrian/nonoccupant involvement",
  "accident", "FATALS", "Number of fatalities", "Crash severity; FARS crashes all have at least one fatality",
  
  # Crash type / roadway / environment
  "accident", "HARM_EV", "First harmful event code", "Useful for sample checks and possible exclusions",
  "accident", "HARM_EVNAME", "First harmful event label", "Readable version of HARM_EV",
  "accident", "MAN_COLL", "Manner of collision", "Useful for understanding two-vehicle crash type",
  "accident", "MAN_COLLNAME", "Manner of collision label", "Readable version of MAN_COLL",
  "accident", "RUR_URB", "Rural/urban indicator", "Optional heterogeneity check",
  "accident", "RUR_URBNAME", "Rural/urban label", "Readable version of RUR_URB",
  "accident", "FUNC_SYS", "Roadway functional system", "Optional roadway-control or descriptive variable",
  "accident", "FUNC_SYSNAME", "Roadway functional system label", "Readable version of FUNC_SYS",
  "accident", "LGT_COND", "Light condition", "Useful descriptive check for nighttime sample",
  "accident", "LGT_CONDNAME", "Light condition label", "Readable version of LGT_COND",
  "accident", "WEATHER", "Weather condition", "Optional descriptive/control variable",
  "accident", "WEATHERNAME", "Weather condition label", "Readable version of WEATHER",
  
  ############################################################
  # person file: person/driver-level variables
  ############################################################
  
  "person", "STATE", "State numeric code", "Join/check variable",
  "person", "STATENAME", "State name", "Readable summaries",
  "person", "ST_CASE", "Crash identifier", "Join key to accident file",
  "person", "VEH_NO", "Vehicle number within crash", "Join key to vehicle file; identifies driver's vehicle",
  "person", "PER_NO", "Person number within vehicle/crash", "Identifies person row",
  "person", "VE_FORMS", "Number of motor vehicles in transport", "Cross-check with accident file",
  
  # Time cross-checks
  "person", "MONTH", "Crash month", "Cross-check with accident file",
  "person", "MONTHNAME", "Crash month label", "Readable version of MONTH if available",
  "person", "DAY", "Crash day", "Cross-check with accident file",
  "person", "DAYNAME", "Crash day label", "Readable version of DAY if available",
  "person", "HOUR", "Crash hour", "Cross-check with accident file",
  "person", "HOURNAME", "Crash hour label", "Readable version of HOUR if available",
  "person", "MINUTE", "Crash minute", "Cross-check with accident file",
  "person", "MINUTENAME", "Crash minute label", "Readable version of MINUTE if available",
  
  # Driver demographics and role
  "person", "AGE", "Person age", "Needed for LP age extension, especially under-25 vs older drivers",
  "person", "AGENAME", "Person age label", "Readable version of AGE",
  "person", "SEX", "Person sex", "Needed for LP gender extension",
  "person", "SEXNAME", "Person sex label", "Readable version of SEX",
  "person", "PER_TYP", "Person type", "Needed to restrict to drivers of motor vehicles in transport",
  "person", "PER_TYPNAME", "Person type label", "Readable version of PER_TYP",
  
  # Injury / occupant information
  "person", "INJ_SEV", "Injury severity", "Useful for checking fatal/nonfatal driver status",
  "person", "INJ_SEVNAME", "Injury severity label", "Readable version of INJ_SEV",
  "person", "SEAT_POS", "Seat position", "Optional check; helps identify occupants vs nonoccupants",
  "person", "SEAT_POSNAME", "Seat position label", "Readable version of SEAT_POS",
  "person", "REST_USE", "Restraint use", "Potential covariate if using imputation-style predictors later",
  "person", "REST_USENAME", "Restraint use label", "Readable version of REST_USE",
  
  # Alcohol variables
  "person", "DRINKING", "Police-reported alcohol involvement", "Primary Levitt-Porter drinking definition for first replication",
  "person", "DRINKINGNAME", "Police-reported alcohol involvement label", "Readable version of DRINKING",
  "person", "ALC_STATUS", "Alcohol test status", "Needed to identify whether BAC test was given",
  "person", "ALC_STATUSNAME", "Alcohol test status label", "Readable version of ALC_STATUS",
  "person", "ATST_TYP", "Alcohol test type", "Blood/breath/urine/etc.; useful for BAC validity checks",
  "person", "ATST_TYPNAME", "Alcohol test type label", "Readable version of ATST_TYP",
  "person", "ALC_RES", "Alcohol test result / BAC result", "Needed for observed BAC >= .08 or >= .10 definitions",
  "person", "ALC_RESNAME", "Alcohol test result label", "Readable version of ALC_RES",
  
  # Drug variables, not core LP but useful later
  "person", "DRUGS", "Police-reported drug involvement", "Not central now, but useful for future robustness checks",
  "person", "DRUGSNAME", "Drug involvement label", "Readable version of DRUGS",
  "person", "DSTATUS", "Drug test status", "Not central now; useful later if considering drugs",
  "person", "DSTATUSNAME", "Drug test status label", "Readable version of DSTATUS",
  
  ############################################################
  # vehicle file: vehicle/driver-linked variables
  ############################################################
  
  "vehicle", "STATE", "State numeric code", "Join/check variable",
  "vehicle", "STATENAME", "State name", "Readable summaries",
  "vehicle", "ST_CASE", "Crash identifier", "Join key to accident/person files",
  "vehicle", "VEH_NO", "Vehicle number within crash", "Join key to person file",
  "vehicle", "VE_FORMS", "Number of motor vehicles in transport", "Cross-check with accident file",
  
  # Basic unit/vehicle descriptors
  "vehicle", "NUMOCCS", "Number of occupants", "Useful for descriptive checks",
  "vehicle", "NUMOCCSNAME", "Number of occupants label", "Readable version of NUMOCCS",
  "vehicle", "UNITTYPE", "Vehicle/unit type", "Helps confirm motor vehicle in transport",
  "vehicle", "UNITTYPENAME", "Vehicle/unit type label", "Readable version of UNITTYPE",
  "vehicle", "HIT_RUN", "Hit-and-run indicator", "Potential exclusion or robustness check",
  "vehicle", "HIT_RUNNAME", "Hit-and-run label", "Readable version of HIT_RUN",
  "vehicle", "MAKE", "Vehicle make code", "Optional vehicle descriptor",
  "vehicle", "MAKENAME", "Vehicle make label", "Readable vehicle descriptor",
  "vehicle", "MODEL", "Vehicle model code", "Optional vehicle descriptor",
  "vehicle", "BODY_TYP", "NCSA body type", "Useful for classifying cars/trucks/motorcycles",
  "vehicle", "BODY_TYPNAME", "NCSA body type label", "Readable version of BODY_TYP",
  "vehicle", "MOD_YEAR", "Vehicle model year", "Optional vehicle age/safety proxy",
  "vehicle", "MOD_YEARNAME", "Vehicle model year label", "Readable version of MOD_YEAR",

  # Driver presence and vehicle-level drinking variable
  "vehicle", "DR_PRES", "Driver presence", "Useful for checking whether the vehicle has a driver",
  "vehicle", "DR_PRESNAME", "Driver presence label", "Readable version of DR_PRES",
  "vehicle", "DR_DRINK", "Driver drinking, vehicle/driver-level version if present", "Useful check against person$DRINKING",
  "vehicle", "DR_DRINKNAME", "Driver drinking label if present", "Readable version of DR_DRINK",
 
  # Prior driving history variables: key for LP bad driving record extension
  "vehicle", "PREV_ACC", "Previous recorded crashes", "LP-style minor blemish; used to construct bad previous driving record",
  "vehicle", "PREV_ACCNAME", "Previous recorded crashes label", "Readable/check version of PREV_ACC",
  "vehicle", "PREV_SPD", "Previous speeding convictions", "LP-style minor blemish; used to construct bad previous driving record",
  "vehicle", "PREV_SPDNAME", "Previous speeding convictions label", "Readable/check version of PREV_SPD",
  "vehicle", "PREV_OTH", "Previous other moving violation convictions", "LP-style minor blemish; used to construct bad previous driving record",
  "vehicle", "PREV_OTHNAME", "Previous other moving violation convictions label", "Readable/check version of PREV_OTH",
  "vehicle", "PREV_DWI", "Previous DWI convictions", "LP-style major blemish; used to construct bad previous driving record",
  "vehicle", "PREV_DWINAME", "Previous DWI convictions label", "Readable/check version of PREV_DWI",
  "vehicle", "FIRST_MO", "Month of first crash/suspension/conviction", "Useful for checking lookback timing of prior-record information",
  "vehicle", "FIRST_MONAME", "Month of first prior event label", "Readable/check version of FIRST_MO",
  "vehicle", "FIRST_YR", "Year of first crash/suspension/conviction", "Useful for checking lookback timing of prior-record information",
  "vehicle", "FIRST_YRNAME", "Year of first prior event label", "Readable/check version of FIRST_YR",
  "vehicle", "LAST_MO", "Month of most recent crash/suspension/conviction", "Useful for checking lookback timing of prior-record information",
  "vehicle", "LAST_MONAME", "Month of most recent prior event label", "Readable/check version of LAST_MO",
  "vehicle", "LAST_YR", "Year of most recent crash/suspension/conviction", "Useful for checking lookback timing of prior-record information",
  "vehicle", "LAST_YRNAME", "Year of most recent prior event label", "Readable/check version of LAST_YR",

  ############################################################
  # derived variables: created in scripts, not raw FARS columns
  ############################################################
  
  "derived", "YEAR_FILE", "Year inferred from raw FARS folder/file", "Useful for debugging multi-year imports and joining files consistently",
  "derived", "night_lp", "Indicator for LP nighttime window, usually HOUR >= 20 or HOUR <= 4", "Restricts sample to 8pm-5am as in LP",
  "derived", "one_two_vehicle", "Indicator for VE_FORMS in 1 or 2", "Restricts sample to one- and two-vehicle fatal crashes",
  "derived", "weekend_lp", "Weekend indicator based on DAY_WEEK", "Needed for LP-style hour-state-year-weekend cells",
  "derived", "covid_period", "Pre-COVID, COVID, and post-COVID period label", "Useful for modern 2015-2024 comparison and sensitivity analysis",
  "derived", "police_drinking_complete", "Cleaned police drinking variable with yes/no and unknown as missing", "Main LP-style drinking definition using complete police-reported drinking status",
  "derived", "police_drinking_unknown_as_no", "Cleaned police drinking variable treating unknown/not reported as not drinking", "Sensitivity definition; useful because missing/unknown drinking classifications are common",
  "derived", "bac_valid", "Indicator for valid observed BAC result", "Needed before constructing BAC threshold definitions",
  "derived", "bac_value", "BAC value converted from FARS integer scale to decimal BAC", "Used for descriptive BAC checks and BAC threshold definitions",
  "derived", "bac_ge_08", "Indicator for observed BAC >= 0.08", "Modern legal-threshold drunk definition",
  "derived", "bac_ge_10", "Indicator for observed BAC >= 0.10", "Historical LP-comparable legal-threshold drunk definition",
  "derived", "bac_tested", "Indicator for whether a valid BAC test result is available", "Needed for BAC-completeness summaries by state-year",
  "derived", "bac_complete_state_year", "State-year indicator for sufficient BAC reporting completeness", "Needed to mimic LP's restricted BAC sample logic",
  "derived", "young_driver_lp", "Indicator for driver under age 25", "Needed for LP age extension and drinking-by-age interactions",
  "derived", "male_driver_lp", "Indicator for male driver", "Needed for LP gender extension and drinking-by-gender interactions",
  "derived", "minor_blemishes_lp", "Count of prior minor blemishes: prior crashes, speeding convictions, and other moving violations", "Intermediate variable for LP bad previous driving record definition",
  "derived", "major_blemish_lp", "Indicator for prior DWI", "Intermediate variable for LP bad previous driving record definition",
  "derived", "bad_record_lp", "Indicator for LP-style bad previous driving record", "Needed for LP bad-record extension and drinking-by-record interactions",
  "derived", "cell_hour", "Hour cell", "Used for Table 2-style increasingly disaggregated equal-mixing assumptions",
  "derived", "cell_hour_year", "Hour-by-year cell", "Used for Table 2-style increasingly disaggregated equal-mixing assumptions",
  "derived", "cell_hour_year_weekend", "Hour-by-year-by-weekend cell", "Used for Table 2-style increasingly disaggregated equal-mixing assumptions",
  "derived", "cell_hour_state_year", "Hour-by-state-by-year cell", "Used for Table 2-style increasingly disaggregated equal-mixing assumptions",
  "derived", "cell_hour_state_year_weekend", "Hour-by-state-by-year-by-weekend cell", "Closest to LP preferred Table 2 specification",
  "derived", "crash_drink_type", "Crash-level combination of driver drinking statuses, e.g., DD, DS, SS", "Core outcome/tally needed for LP two-car likelihood",
  "derived", "one_car_drink_type", "One-car crash drinking status, D or S", "Needed if estimating one-car relative risk lambda alongside two-car risk theta",
  
  ############################################################
  # external state-year policy variables: needed only for LP Table 8
  ############################################################
  
  "external_policy", "beer_tax_real", "State beer excise tax, inflation-adjusted", "Needed only if replicating LP policy regressions / Table 8",
  "external_policy", "first_offense_punishment", "Severity index for first DUI offense", "Needed only if replicating LP policy regressions / Table 8",
  "external_policy", "third_offense_punishment", "Severity index for third DUI offense", "Needed only if replicating LP policy regressions / Table 8",
  "external_policy", "police_per_capita", "Police officers per capita", "Needed only if replicating LP policy regressions / Table 8",
  "external_policy", "share_age_15_24", "State population share ages 15-24", "Control variable for LP policy regressions / Table 8",
  "external_policy", "share_age_25_54", "State population share ages 25-54", "Control variable for LP policy regressions / Table 8",
  "external_policy", "share_age_55_plus", "State population share ages 55+", "Control variable for LP policy regressions / Table 8",
  "external_policy", "pct_metro", "Percent of state population in metropolitan areas", "Control variable for LP policy regressions / Table 8",
  "external_policy", "pct_black", "Percent of state population Black", "Control variable for LP policy regressions / Table 8",
  "external_policy", "unemployment_rate", "State unemployment rate", "Control variable for LP policy regressions / Table 8"
)

var_reference

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

write_csv(var_reference, "data/processed/fars_variable_reference_for_replication.csv")

###################################
# most important variables for now:
# 
# accident: ST_CASE, YEAR, HOUR, DAY_WEEK, VE_FORMS
# person:   ST_CASE, VEH_NO, PER_TYP, DRINKING, ALC_STATUS, ATST_TYP, ALC_RES
# vehicle:  ST_CASE, VEH_NO, UNITTYPE, BODY_TYP
