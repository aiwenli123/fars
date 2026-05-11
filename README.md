# FARS Drunk Driving Replication Project

This repository contains code for preparing recent Fatality Analysis Reporting System (FARS) data for a Levitt–Porter-style drunk driving replication using 2015–2024 FARS National CSV files.

## Raw data

Raw FARS data are not included in this repository due to file size. To reproduce the analysis, download the FARS National CSV files from:

https://www.nhtsa.gov/file-downloads?p=nhtsa/downloads/FARS/

Place the raw yearly folders in:

```text
data/raw/
```

## Alcohol variable decision

For 2015–2024 FARS, `person$DRINKING` and `vehicle$DR_DRINK` do not match perfectly among driver rows. Mismatches are concentrated among cases where `person$DRINKING` is No, Unknown, or Not Reported.

In many of these cases, `vehicle$DR_DRINK` classifies the driver as drinking when observed BAC evidence is positive. This suggests that `DR_DRINK` may incorporate information beyond the raw police-reported `DRINKING` variable.

For this reason, we can treat `DRINKING`, `DR_DRINK`, and BAC-threshold definitions as separate alcohol classifications rather than assuming they are interchangeable.

The current plan is to compare results across several definitions of alcohol involvement:

- raw police-reported alcohol involvement using `person$DRINKING`;
- police-reported alcohol involvement with unknown or not reported values treated as no;
- vehicle-level drinking involvement using `vehicle$DR_DRINK`;
- observed BAC at or above 0.08;
- observed BAC at or above 0.10.

The BAC-threshold definitions are useful but should be interpreted carefully because observed BAC is missing for many drivers and BAC testing is unlikely to be random.

## Code scripts

### `00_read_data.R`

**Notes**

`00_read_data.R` is an initial exploration script. It reads in one year of raw FARS data, checks the structure of the accident, person, and vehicle files, and creates a manually curated variable reference table for the replication project.

This script is mainly useful for understanding the raw FARS files and documenting which variables are needed later.

**Inputs**

This script currently reads the 2024 FARS National CSV files:

```text
data/raw/FARS2024NationalCSV/person.csv
data/raw/FARS2024NationalCSV/vehicle.csv
data/raw/FARS2024NationalCSV/accident.csv
```

**Output**

The script saves the variable reference table to:

```text
data/processed/fars_variable_reference_for_replication.csv
```

### `01_clean_data.R`

**Notes**

`01_clean_data.R` is the main data cleaning script. It reads and combines the 2015–2024 FARS accident, person, and vehicle files, creates cleaned driver-level variables, constructs a Levitt-Porter-style analysis sample with alcohol variables and other covariates, and checks whether two FARS alcohol variables match.

It also checks and records BAC completeness by state-year in a separate dataset, to use for filtering out state-years with incomplete BAC testing (e.g., filtering out ones with < 95% complete BAC entries, as Levitt-Porter originally did).

**Inputs**

For each year from 2015 through 2024, the script expects the following files:

```text
data/raw/FARS{year}NationalCSV/accident.csv
data/raw/FARS{year}NationalCSV/person.csv
data/raw/FARS{year}NationalCSV/vehicle.csv
```

For example, for 2015:

```text
data/raw/FARS2015NationalCSV/accident.csv
data/raw/FARS2015NationalCSV/person.csv
data/raw/FARS2015NationalCSV/vehicle.csv
```

The full expected folder structure is:

```text
data/raw/
  FARS2015NationalCSV/
    accident.csv
    person.csv
    vehicle.csv
  FARS2016NationalCSV/
    accident.csv
    person.csv
    vehicle.csv
  ...
  FARS2024NationalCSV/
    accident.csv
    person.csv
    vehicle.csv
```

**Alcohol variables created**

The script creates several alcohol-involvement variables:

```text
police_drinking_complete
```

Equals 1 if `person$DRINKING == 1`, equals 0 if `person$DRINKING == 0`, and is missing otherwise.

```text
police_drinking_unknown_as_no
```

Equals 1 if `person$DRINKING == 1`, and equals 0 otherwise. This treats unknown or not reported alcohol involvement as no reported alcohol involvement.

```text
bac_valid
```

Indicates whether `person$ALC_RES` appears to contain a valid observed BAC result.

```text
bac_value
```

Observed BAC value, created as `ALC_RES / 1000` for valid BAC values.

```text
bac_ge_08
```

Equals 1 if observed BAC is at least 0.08 and equals 0 if observed BAC is valid but below 0.08.

```text
bac_ge_10
```

Equals 1 if observed BAC is at least 0.10 and equals 0 if observed BAC is valid but below 0.10.

**Driver characteristics variables created**

The script also creates variables for later Levitt-Porter extensions.

```text
young_driver_lp
```

Indicator for whether the driver is under age 25.

```text
male_driver_lp
```

Indicator for whether the driver is male.

```text
minor_blemishes_lp
```

Count of minor prior driving-record blemishes, based on previous crashes, previous speeding convictions, and previous other moving violations.

```text
major_blemishes_lp
```

Indicator for whether the driver has a major prior driving-record blemish. The current version considers prior DWI and license suspension/revocation variables (license suspension variables are missing for 2015-2017 though).

```text
bad_record_lp
```

Levitt-Porter-style bad previous driving record indicator. A driver is classified as having a bad prior record if they have either at least two minor blemishes, or at least one major blemish.

**Levitt-Porter-style sample**

The current cleaned analysis sample keeps:

- driver rows, using `PER_TYP == 1`;
- one- and two-vehicle crashes, using `VE_FORMS %in% c(1, 2)`;
- nighttime crashes, currently defined as crashes from 8:00pm through 4:59am;
- crashes where the number of driver rows equals the number of motor vehicles in transport.

The main cleaned analysis file has one row per driver. For two-vehicle crashes, there should be two driver rows per crash.

**Outputs**

The script saves cleaned core files:

```text
data/processed/accident_2015_2024_core.csv
data/processed/person_2015_2024_core.csv
data/processed/vehicle_2015_2024_core.csv
```

The script saves driver/crash analysis files:

```text
data/processed/driver_analysis_2015_2024.csv
data/processed/driver_analysis_lp_2015_2024.csv
data/processed/lp_clean_crashes_2015_2024.csv
data/processed/bac_completeness_state_year_2015_2024.csv
```

The main analysis file is:

```text
data/processed/driver_analysis_lp_2015_2024.csv
```

This file contains one row per driver in the cleaned one- and two-vehicle nighttime fatal crash sample.

The alcohol-variable diagnostic files are:

```text
data/processed/drinking_vs_dr_drink_mismatch_summary_2015_2024.csv
data/processed/drinking_vs_dr_drink_match_by_year_2015_2024.csv
```

These summarize how the person-level drinking variable `DRINKING` compares to the vehicle-level drinking variable `DR_DRINK`.

### `02_crash_type_counts.R`

Collapses the cleaned driver-level LP sample into crash-level drunk/sober crash types.

For two-vehicle crashes, the script classifies each crash as:

```text
DD = two drinking drivers
DS = one drinking driver and one sober driver
SS = two sober drivers
```

For one-vehicle crashes, the script classifies each crash as:

```text
D = one drinking driver
S = one sober driver
```

Crashes that cannot be classified because of missing or ambiguous driver drinking status are labeled:

```text
missing_or_ambiguous
```

The script then computes simple Levitt-Porter-style two-car relative risk estimates (theta) using the crash type counts.

**Main input**

```text
data/processed/driver_analysis_lp_2015_2024.csv,
```

created in '01_clean_data.R'.

The script compares five alcohol definitions:

```text
Police drinking, complete case
Police drinking, unknown as no
Vehicle DR_DRINK, complete case
Observed BAC >= 0.08
Observed BAC >= 0.10
```

over four different periods:

```text
2015-2024 all years
2015-2019 pre-COVID
2020-2021 COVID
2022-2024 post-COVID
```

This is a preliminary descriptive implementation of the Levitt-Porter two-car crash count logic. Later scripts may estimate richer versions with more disaggregated cells and additional assumptions.

**Outputs**

Crash-level files, counts, and 2-car crash theta estimates:

```text
data/processed/lp_crash_types_2015_2024.csv
data/processed/lp_crash_type_counts_2015_2024.csv
data/processed/lp_theta_simple_2015_2024.csv
data/processed/lp_one_vehicle_counts_2015_2024.csv
data/processed/lp_theta_simple_all_years_quick_table.csv
data/processed/lp_theta_simple_by_covid_period_quick_table.csv
```

Figures:

```text
figures/lp_crash_type_shares_by_definition_all_years.png
figures/lp_theta_by_definition_and_period.png
figures/lp_missing_by_definition_and_period.png
```

### '03_more_eda.R'

Creates cleaner descriptive tables and plots for the FARS / Levitt-Porter replication project. This script is mainly intended to make the outputs from scripts 01 and 02 easier to interpret. It does not redefine the main sample or estimate a new model. Instead, it organizes the current results into tables and figures that are useful for checking data quality, comparing alcohol definitions, and planning the next replication steps.

**Inputs**

From scripts 01 and 02:

```text
data/processed/driver_analysis_lp_2015_2024.csv
data/processed/lp_crash_types_2015_2024.csv
data/processed/lp_theta_simple_2015_2024.csv
data/processed/lp_one_vehicle_counts_2015_2024.csv
```

**Outputs**

The script saves 6 main CSV tables and four main plots in 

```text
data/processed/script03_tables/...
figures/script03_plots/...
```
