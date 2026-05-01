rm(list=ls())

library(tidyverse)

person <- read.csv("data/raw/FARS2024NationalCSV/person.csv")
str(person)

vehicle <- read.csv("data/raw/FARS2024NationalCSV/vehicle.csv")
str(vehicle)

accident <- read.csv("data/raw/FARS2024NationalCSV/accident.csv")
str(accident)

unique(accident$YEAR)

# manually created data dictionary
var_reference <- tribble(
  ~dataset, ~raw_column, ~meaning_for_replication, ~why_we_need_it,
  
  "accident", "STATE", "State numeric code", "Used for state-year cells and possible state-level comparisons",
  "accident", "STATENAME", "State name", "Useful for readable summaries",
  "accident", "ST_CASE", "Unique crash identifier within year/state", "Primary key for joining accident, person, and vehicle files",
  "accident", "YEAR", "Crash year", "Needed for 2015-2024 and COVID-period splits",
  "accident", "MONTH", "Crash month", "Potential time control; useful for descriptive checks",
  "accident", "DAY", "Day of month", "Potential time control; mostly for checking",
  "accident", "DAY_WEEK", "Day of week", "Used to define weekend/night patterns if needed",
  "accident", "HOUR", "Hour of crash", "Needed to restrict to Levitt-Porter nighttime window, roughly 8pm-5am",
  "accident", "MINUTE", "Minute of crash", "Mostly for checking; not essential for first pass",
  "accident", "VE_TOTAL", "Total vehicles involved", "Useful cross-check for one- and two-vehicle crash definitions",
  "accident", "VE_FORMS", "Number of motor vehicles in transport", "Primary variable for defining one-vehicle and two-vehicle crashes",
  "accident", "PVH_INVL", "Number of parked/working vehicles involved", "Useful for distinguishing motor vehicles in transport from parked/working vehicles",
  "accident", "PERSONS", "Number of persons involved", "Descriptive check",
  "accident", "PERMVIT", "Persons in motor vehicles in transport", "Descriptive check",
  "accident", "PEDS", "Number of pedestrians", "Useful for checking nonoccupant crashes",
  "accident", "PERNOTMVIT", "Persons not in motor vehicles in transport", "Useful for checking pedestrian/nonoccupant involvement",
  "accident", "FATALS", "Number of fatalities", "Crash severity; FARS crashes all have at least one fatality",
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
  
  "person", "STATE", "State numeric code", "Join/check variable",
  "person", "STATENAME", "State name", "Readable summaries",
  "person", "ST_CASE", "Crash identifier", "Join key to accident file",
  "person", "VEH_NO", "Vehicle number within crash", "Join key to vehicle file; identifies driver’s vehicle",
  "person", "PER_NO", "Person number within vehicle/crash", "Identifies person row",
  "person", "VE_FORMS", "Number of motor vehicles in transport", "Cross-check with accident file",
  "person", "MONTH", "Crash month", "Cross-check with accident file",
  "person", "DAY", "Crash day", "Cross-check with accident file",
  "person", "HOUR", "Crash hour", "Cross-check with accident file",
  "person", "MINUTE", "Crash minute", "Cross-check with accident file",
  "person", "AGE", "Person age", "Useful for later subgroup analysis and replication extensions",
  "person", "AGENAME", "Person age label", "Readable version of AGE",
  "person", "SEX", "Person sex", "Useful for later subgroup analysis and replication extensions",
  "person", "SEXNAME", "Person sex label", "Readable version of SEX",
  "person", "PER_TYP", "Person type", "Needed to restrict to drivers of motor vehicles in transport",
  "person", "PER_TYPNAME", "Person type label", "Readable version of PER_TYP",
  "person", "INJ_SEV", "Injury severity", "Useful for checking fatal/nonfatal driver status",
  "person", "INJ_SEVNAME", "Injury severity label", "Readable version of INJ_SEV",
  "person", "SEAT_POS", "Seat position", "Optional check; helps identify occupants vs nonoccupants",
  "person", "SEAT_POSNAME", "Seat position label", "Readable version of SEAT_POS",
  "person", "REST_USE", "Restraint use", "Potential covariate if using imputation-style predictors later",
  "person", "REST_USENAME", "Restraint use label", "Readable version of REST_USE",
  "person", "DRINKING", "Police-reported alcohol involvement", "Primary Levitt-Porter drinking definition for first replication",
  "person", "DRINKINGNAME", "Police-reported alcohol involvement label", "Readable version of DRINKING",
  "person", "ALC_STATUS", "Alcohol test status", "Needed to identify whether BAC test was given",
  "person", "ALC_STATUSNAME", "Alcohol test status label", "Readable version of ALC_STATUS",
  "person", "ATST_TYP", "Alcohol test type", "Blood/breath/urine/etc.; useful for BAC validity checks",
  "person", "ATST_TYPNAME", "Alcohol test type label", "Readable version of ATST_TYP",
  "person", "ALC_RES", "Alcohol test result / BAC result", "Needed for observed BAC >= .08 or >= .10 definitions",
  "person", "ALC_RESNAME", "Alcohol test result label", "Readable version of ALC_RES",
  "person", "DRUGS", "Police-reported drug involvement", "Not central now, but useful for future robustness checks",
  "person", "DRUGSNAME", "Drug involvement label", "Readable version of DRUGS",
  "person", "DSTATUS", "Drug test status", "Not central now; useful later if considering drugs",
  "person", "DSTATUSNAME", "Drug test status label", "Readable version of DSTATUS",
  
  "vehicle", "STATE", "State numeric code", "Join/check variable",
  "vehicle", "STATENAME", "State name", "Readable summaries",
  "vehicle", "ST_CASE", "Crash identifier", "Join key to accident/person files",
  "vehicle", "VEH_NO", "Vehicle number within crash", "Join key to person file",
  "vehicle", "VE_FORMS", "Number of motor vehicles in transport", "Cross-check with accident file",
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
  "vehicle", "VPICMAKE", "vPIC make code", "Optional modern vehicle descriptor",
  "vehicle", "VPICMAKENAME", "vPIC make label", "Readable vPIC descriptor",
  "vehicle", "VPICMODEL", "vPIC model code", "Optional modern vehicle descriptor",
  "vehicle", "VPICMODELNAME", "vPIC model label", "Readable vPIC descriptor",
  "vehicle", "VPICBODYCLASS", "vPIC body class code", "Optional modern body class",
  "vehicle", "VPICBODYCLASSNAME", "vPIC body class label", "Readable vPIC body class",
  "vehicle", "GVWR_FROM", "Gross vehicle weight rating lower class", "Useful for broad vehicle class",
  "vehicle", "GVWR_FROMNAME", "GVWR lower class label", "Readable version of GVWR_FROM",
  "vehicle", "GVWR_TO", "Gross vehicle weight rating upper class", "Useful for broad vehicle class",
  "vehicle", "GVWR_TONAME", "GVWR upper class label", "Readable version of GVWR_TO",
  "vehicle", "TRAV_SP", "Travel speed", "Optional later robustness/control variable",
  "vehicle", "TRAV_SPNAME", "Travel speed label", "Readable version of TRAV_SP",
  "vehicle", "ROLLOVER", "Rollover indicator", "Optional crash-mechanism descriptor",
  "vehicle", "ROLLOVERNAME", "Rollover label", "Readable version of ROLLOVER",
  "vehicle", "M_HARM", "Most harmful event for vehicle", "Optional vehicle-level crash severity descriptor",
  "vehicle", "M_HARMNAME", "Most harmful event label", "Readable version of M_HARM",
  "vehicle", "DR_PRES", "Driver presence", "Useful for checking whether the vehicle has a driver",
  "vehicle", "DR_PRESNAME", "Driver presence label", "Readable version of DR_PRES",
  "vehicle", "DR_DRINK", "Driver drinking, vehicle/driver-level version if present", "Useful check against person$DRINKING",
  "vehicle", "DR_DRINKNAME", "Driver drinking label if present", "Readable version of DR_DRINK"
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
