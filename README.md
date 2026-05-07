# FARS Drunk Driving Replication Project

This repository contains code for preparing recent Fatality Analysis Reporting System (FARS) data for a Levitt–Porter-style drunk driving replication using 2015–2024 FARS National CSV files.

The current goal is to build a cleaned driver-level analysis dataset for one- and two-vehicle fatal crashes, define several candidate alcohol-involvement measures, and prepare for replication analyses comparing fatal crash involvement among drunk and sober drivers.

## Raw data

Raw FARS data are not included in this repository due to file size. To reproduce the analysis, download the FARS National CSV files from:

https://www.nhtsa.gov/file-downloads?p=nhtsa/downloads/FARS/

Place the raw yearly folders in:

```text
data/raw/
```

## Code scripts

This repository currently contains two main R scripts:

```text
00_read_data.R
01_clean_data.R
```

The first script is for initial data inspection and variable documentation. The second script builds the cleaned 2015–2024 analysis datasets.

---

## `00_read_data.R`

### Purpose

`00_read_data.R` is an initial exploration script. It reads in one year of raw FARS data, checks the structure of the accident, person, and vehicle files, and creates a manually curated variable reference table for the replication project.

This script is mainly useful for understanding the raw FARS files and documenting which variables are needed later.

### Inputs

This script currently reads the 2024 FARS National CSV files:

```text
data/raw/FARS2024NationalCSV/person.csv
data/raw/FARS2024NationalCSV/vehicle.csv
data/raw/FARS2024NationalCSV/accident.csv
```

### Output

The script saves the variable reference table to:

```text
data/processed/fars_variable_reference_for_replication.csv
```

---

## `01_clean_data.R`

### Purpose

`01_clean_data.R` is the main data cleaning script. It reads and combines the 2015–2024 FARS accident, person, and vehicle files, creates driver-level alcohol variables, builds a cleaned Levitt-Porter-style analysis sample, and checks whether two FARS alcohol variables match.

### Inputs

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

### Alcohol variables created

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

### Levitt-Porter-style sample

The current cleaned analysis sample keeps:

- driver rows, using `PER_TYP == 1`;
- one- and two-vehicle crashes, using `VE_FORMS %in% c(1, 2)`;
- nighttime crashes, currently defined as crashes from 8:00pm through 4:59am;
- crashes where the number of driver rows equals the number of motor vehicles in transport.

The main cleaned analysis file has one row per driver. For two-vehicle crashes, there should be two driver rows per crash.

### Outputs

The script saves the following processed datasets:

```text
data/processed/accident_2015_2024_core.csv
data/processed/person_2015_2024_core.csv
data/processed/vehicle_2015_2024_core.csv
data/processed/driver_analysis_lp_2015_2024.csv
data/processed/drinking_vs_dr_drink_mismatch_summary_2015_2024.csv
data/processed/drinking_vs_dr_drink_match_by_year_2015_2024.csv
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

---

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

---

Next step: collapse two-vehicle crashes to crash-level drunk/sober pair types:

```text
DD = two drinking drivers
DS = one drinking driver and one sober driver
SS = two sober drivers
```

These crash-type counts can then be used to estimate the Levitt-Porter relative fatal-crash risk parameter under different alcohol definitions.
