# FARS Levitt-Porter Replication and BAC Extensions

This repository contains R code for preparing and analyzing 2015–2024 Fatality Analysis Reporting System (FARS) data. The project begins with a replication and adaptation of the Levitt and Porter (2001) fatal-crash model and then extends the analysis in two directions:
1. comparing estimates based on observed BAC with estimates based on the FARS multiple-imputation BAC files
2. replacing the binary sober/drinking classification with four BAC groups centered around the potential 0.05 BAC policy threshold.

The main quantities estimated are:
- `theta`: relative fatal two-car crash risk
- `lambda`: relative fatal one-car crash risk
- `N_hat`: the estimated drinking-to-sober road-exposure ratio in the binary model
- `r_L`, `r_M`, and `r_H`: estimated road-exposure ratios relative to the sober group in the four-group model.

The repository also contains data-quality diagnostics, regression analyses of the estimated exposure ratios, and comparisons between the binary and four-group models.

## Raw data

Raw FARS data are not included in this repository because of file size. Download the yearly FARS National CSV files from: https://www.nhtsa.gov/file-downloads?p=nhtsa/downloads/FARS/

Place the yearly folders under:

```text
data/raw/
```

The expected directory structure is:

```text
data/raw/
  FARS2015NationalCSV/
    accident.csv
    person.csv
    vehicle.csv
    MIPER.csv
  FARS2016NationalCSV/
    accident.csv
    person.csv
    vehicle.csv
    MIPER.csv
  ...
  FARS2024NationalCSV/
    accident.csv
    person.csv
    vehicle.csv
    MIPER.csv
```

`06_prep_mi_data.R` accepts either `MIPER.csv` or `miper.csv`, with capitalization ignored. The exploratory `05_mi_exploration.R` script additionally reads the 2024 files:

```text
data/raw/FARS2024NationalCSV/MIACC.csv
data/raw/FARS2024NationalCSV/MIDRVACC.csv
data/raw/FARS2024NationalCSV/MIPER.csv
```

## Core analysis definitions

### Levitt-Porter-style sample

The main cleaned analysis sample contains:

- driver rows, identified using `PER_TYP == 1`;
- fatal crashes involving one or two motor vehicles in transport, using `VE_FORMS %in% c(1, 2)`;
- crashes occurring from 8:00 p.m. through 4:59 a.m.;
- crashes in which the number of identified driver rows equals the number of motor vehicles in transport; and
- years 2015 through 2024.

The cleaned driver-level file has one row for a one-car crash and two rows for a two-car crash.

### Time definitions

The main period variable is:

```text
pre_covid  = 2015–2019
covid      = 2020–2021
post_covid = 2022–2024
```

The weekend indicator includes Saturday and Sunday:

```text
DAY_WEEK %in% c(7, 1)
```

### Rural and urban classifications

Based on FARS' documentation found at https://crashstats.nhtsa.dot.gov/#!/DocumentTypeList/23, the current scripts classify:

```text
RUR_URB == 1  as rural
RUR_URB == 2  as urban
RUR_URB in 6, 8, or 9 as other/unknown
```

Analyses using rural/urban cells generally keep only the clearly classified rural and urban observations.

## Alcohol and BAC definitions

Several alcohol classifications are used in the earlier scripts.

### Police-reported drinking

The complete-case police definition from the `person` file is:

```text
1 if DRINKING == 1
0 if DRINKING == 0
NA otherwise
```

An alternative definition treats unknown or unreported values as no:

```text
1 if DRINKING == 1
0 otherwise
```

### Vehicle-level drinking classification

The vehicle-level definition uses `DR_DRINK`:

```text
1 if DR_DRINK == 1
0 if DR_DRINK == 0
NA otherwise
```

`DRINKING` and `DR_DRINK` do not agree perfectly. The discrepancies are summarized by `01_clean_data.R` rather than assuming the two variables are interchangeable.

### Observed BAC thresholds

Observed BAC comes from `ALC_RES`, which records BAC in thousandths. For example:

```text
ALC_RES = 80 corresponds to BAC = 0.080
```

The binary BAC classifications are:

```text
BAC >= 0.05
BAC >= 0.08
BAC >= 0.10
```

Special missing or unknown codes are excluded from the observed-BAC classifications.

### FARS multiple-imputation BAC

The FARS `MIPER` files contain ten imputed BAC values, `P1` through `P10`, for each person. These values are recorded in hundredths. For example:

```text
P1 = 8 corresponds to BAC = 0.08
```

The main MI variable created in `06_prep_mi_data.R` is:

```text
analysis_bac_value
```

It uses:

```text
the exact observed BAC when a valid measured value is available;
the corresponding MIPER value when observed BAC is missing or invalid.
```

The ten resulting completed datasets are analyzed separately and then pooled.

### Four BAC groups

Scripts 09–12 use the following policy-relevant BAC groups:

```text
S: BAC = 0
L: 0 < BAC < 0.05
M: 0.05 <= BAC < 0.08
H: BAC >= 0.08
```

In the four-group model, `S` is the reference group.

## BAC-complete state-year restriction

The BAC-completeness variable follows the Levitt-Porter restriction. For each state-year, the script calculates:

```text
number of police-reported drinking drivers with valid observed BAC
------------------------------------------------------------------
total number of police-reported drinking drivers
```

A state-year is labeled `BAC_COMPLETE_95` when this proportion is at least 0.95. Note that this is not the same as 95% of _all_ drivers in the state-year having observed BAC.

## Model parameters

### Binary model

For a binary sober/drinking classification:

```text
theta
```

is the fatal two-car crash risk of the drinking group relative to the sober group.

```text
lambda
```

is the fatal one-car crash risk of the drinking group relative to the sober group.

```text
N_hat
```

is the estimated drinking-to-sober road-exposure ratio within a cell:

```text
N_hat = estimated drinking exposure / estimated sober exposure
```

The corresponding estimated drinking share is:

```text
N_hat / (1 + N_hat)
```

### Four-group model

The four-group model fixes the sober-group parameters at:

```text
theta_S = 1
lambda_S = 1
r_S = 1
```

It estimates:

```text
theta_L, theta_M, theta_H
lambda_L, lambda_M, lambda_H
r_Lc, r_Mc, r_Hc
```

where the `r` parameters are cell-specific exposure ratios relative to the sober group.

## Recommended script order

Run the scripts in the following order:

```text
00_read_data.R                         optional raw-data exploration
01_clean_data.R
02_crash_type_counts.R
03_more_eda.R                          optional descriptive outputs
04_lp_mle_period.R
04.5_validate_closed_form_theta.R
05_mi_exploration.R                    optional MI-file exploration
06_prep_mi_data.R
07_mi_lp_mle.R
08_regression_analysis.R
09_mi_diagnostics.R
10_four_group_prep.R
11_four_group_mle.R
12_four_group_results.R
```

`07_mi_lp_mle.R` analyzes one cell-specification family per run. Its current default is:

```r
location_spec <- "LOCATION_WEEKEND"
```

For the complete cross-specification analyses in Scripts 08 and 09, run Script 07 three times, changing the option to:

```r
location_spec <- "STATE"
location_spec <- "RUR_URB"
location_spec <- "LOCATION_WEEKEND"
```

Each run is saved to a separate output folder, so the results do not overwrite one another.

# Code scripts

## `00_read_data.R`

`00_read_data.R` is an initial exploration script. It reads in one year of raw FARS data, checks the structure of the accident, person, and vehicle files, and creates a manually curated variable reference table for the replication project.

This script is mainly useful for understanding the raw FARS files and documenting which variables are needed later.

## `01_clean_data.R`

`01_clean_data.R` is the main data cleaning script. It reads and combines the 2015–2024 FARS accident, person, and vehicle files, creates cleaned driver-level variables, constructs a Levitt-Porter-style analysis sample with alcohol variables and other covariates, and checks whether two FARS alcohol variables match.

It also checks and records BAC completeness by state-year in a separate dataset, to use for filtering out state-years with incomplete BAC testing (e.g., filtering out ones with < 95% complete BAC entries, as Levitt-Porter originally did).

**Driver characteristics variables created**

The script also creates other driver variables for potential extensions based on Levitt & Porter, 2001:

- `young_driver_lp`: Indicator for whether the driver is under age 25.
- `male_driver_lp`: Indicator for whether the driver is male.
- `minor_blemishes_lp`: Count of minor prior driving-record blemishes, based on previous crashes, previous speeding convictions, and previous other moving violations.
- `major_blemishes_lp`: Indicator for whether the driver has a major prior driving-record blemish. The current version considers prior DWI and license suspension/revocation variables (note that license suspension variables are missing for 2015-2017).
- `bad_record_lp`: Levitt-Porter-style bad previous driving record indicator. A driver is classified as having a bad prior record if they have either at least two minor blemishes, or at least one major blemish.

## `02_crash_type_counts.R`

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

Crashes that cannot be classified because of missing or ambiguous driver drinking status are labeled `missing_or_ambiguous`.

The script then computes simple Levitt-Porter-style two-car relative risk estimates (theta) using the crash type counts:

```text
R = DS^2 / (DD * SS)
```

and:

```text
theta = ((R - 2) + sqrt(R^2 - 4R)) / 2, if R >= 4
theta = 1,                              if R < 4
```

The script compares six alcohol definitions:

```text
Police drinking, complete case
Police drinking, unknown as no
Vehicle DR_DRINK, complete case
Observed BAC >= 0.05
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

## `03_more_eda.R`

Creates cleaner descriptive tables and plots for the LP replication project. This script is mainly intended to make the outputs from scripts 01 and 02 easier to interpret. It does not redefine the main sample or estimate a new model. Instead, it organizes the current results into tables and figures that are useful for checking data quality, comparing alcohol definitions, and planning the next replication steps.

## `04_lp_mle_period.R`

Estimates a Levitt-Porter-style maximum likelihood model using the cleaned 2015-2024 FARS sample.

This script moves beyond the simple crash-count estimates from `02_crash_type_counts.R` by estimating Levitt-Porter-style risk parameters using maximum likelihood (Levitt-Porter also used MLE in their 2001 paper). The script is designed to reproduce the spirit of LP's main Table 2 results, while adapting the setup to recent FARS data and the 2015-2024 period.

The script estimates two main parameters:

```text
theta
lambda
```

The script first estimates `theta` using only two-car crashes. It then uses the estimated cell-level drinking-to-sober driver prevalence, `N_hat`, from the two-car model to estimate `lambda` from the one-car crash counts.

Unlike the original Levitt-Porter paper, this script uses period-based cells rather than year-based cells. The periods are:

```text
pre_covid  = 2015-2019
covid      = 2020-2021
post_covid = 2022-2024
```

This is partly because the recent sample has fewer years than the original Levitt-Porter setting, so period-based cells help avoid making the cell structure too sparse.

The script runs each specification on two samples:

```text
All state-years
BAC-complete 95% state-years
```

The BAC-complete sample keeps only state-years where BAC reporting among police-reported drinking drivers is at least 95% complete; same as what Levitt-Porter did.

**Cell-size checks**

The script runs each model twice:

```text
min_n2 = 1
min_n2 = 5
```

where `n2` is the number of complete two-car crashes in a cell. The `min_n2 = 1` version keeps all cells with at least one complete two-car crash. The `min_n2 = 5` version drops very sparse cells and is useful as a sensitivity check, since highly detailed cell specifications (e.g., `hour x state x period x weekend`) can create many cells with very small DD/DS/SS counts.

## `04.5_validate_closed_form_theta.R`

Validates the binary `theta` estimation code in the simplest possible setting:

```text
one all-data cell
two-car crashes only
```

For every alcohol definition and both script 04 sample variants, it compares:

1. the Levitt-Porter closed-form estimate;
2. an independently rerun profile-likelihood optimizer; and
3. the corresponding all-data estimate saved by script 04.

The script also verifies that the `min_n2 = 1` and `min_n2 = 5` estimates are identical when there is only one pooled cell. Cases with `DD = 0` or `SS = 0` are labeled `not_comparable`, because the finite closed-form ratio is undefined.

## `05_mi_exploration.R`

Explores the FARS multiple imputation BAC files for 2024 just to explore their structure. This script reads the three alcohol imputation files downloaded from FARS:

```text
MIACC
MIDRVACC
MIPER
```

We conclude that `MIPER` is the most useful file for this project because it has driver/person-level identifiers, which allow the imputed BAC values to be joined directly to the cleaned driver-level Levitt-Porter sample:

```text
ST_CASE
VEH_NO
PER_NO
```

## `06_prep_mi_data.R`

Prepares driver-level multiple imputation BAC datasets for 2015-2024. This script reads and combines the yearly `MIPER` files from 2015 through 2024, joins them to the cleaned driver-level Levitt-Porter sample, checks merge success, and creates analysis-ready datasets for future multiple imputation analyses.

It also creates a combined BAC variable that uses the exact observed BAC value when valid observed BAC is available, and uses the corresponding MIPER imputed BAC value only when observed BAC is missing or invalid. This is useful because observed BAC is more precise (three significant digits), while MIPER is only recorded to two significant digits.

## `07_mi_lp_mle.R`

Estimates the LP-style MLE using the multiple imputation BAC data prepared in `06_prep_mi_data.R`.

This script uses the same MLE method as in `04_lp_mle_period.R`, and extends on script 04 by comparing raw measured BAC estimates with estimates based on FARS MIPER imputed BAC.

The script runs the analysis for three definitions of `drunk`:

```text
BAC >= 0.05
BAC >= 0.08
BAC >= 0.10
```

**Analysis samples**

The script compares four analysis samples:

```text
A. Raw BAC, all state-years
B. Raw BAC, 95% complete state-years
C. MI BAC, same 95% complete state-years
D. MI BAC, all state-years
```

Comparing the four samples helps show whether differences in `theta` and `lambda` come from the 95% completeness restriction, from imputing missing BAC, or from adding back all state-years after imputation.

The cell specifications are selected using:

```r
location_spec <- "STATE", "RUR_URB", or "LOCATION_WEEKEND"
```

#### `STATE`

```text
all data
hour
hour x period
hour x period x weekend
hour x state x period
hour x state x period x weekend
```

#### `RUR_URB`

```text
all data
hour
hour x weekend
hour x rural/urban
hour x rural/urban x weekend
hour x rural/urban x period
hour x rural/urban x period x weekend
```

#### `LOCATION_WEEKEND`

```text
all data
rural/urban
state
rural/urban x state
rural/urban x weekend
state x weekend
rural/urban x state x weekend
```

## `08_regression_analysis.R` 

Runs exploratory regressions using the cell-level exposure ratios estimated in `07_mi_lp_mle.R`. The main goal is to understand which cell characteristics are associated with higher implied drunk-driver prevalence, and whether MI changes these exposure estimates systematically. The main outcome is: 

```text
log(N_hat)
```

which is the LP-implied drunk-to-sober driver exposure ratio within a cell. The script also analyzes:

```text
Delta log(N_hat) = log(N_hat_MI) - log(N_hat_raw)
```

for matched raw vs. MI cells. The script uses raw BAC, 95% complete state-years vs. MI BAC, all state-years.

The script runs both:

```text
unweighted OLS
OLS weighted by the number of complete two-car crashes, n2
```

MI regressions are estimated separately by imputation and pooled using Rubin-style combining rules.

The regressions are descriptive and are not interpreted causally.

## `09_mi_diagnostics.R` 

Runs diagnostics/EDA to understand why MI-based LP estimates may differ from raw-BAC estimates. The script compares drivers with observed raw BAC values to drivers with missing raw BAC values filled in using FARS MIPER multiple imputation. It checks whether MI-filled missing-BAC drivers appear systematically less alcohol-involved than drivers with observed BAC. It also looks at four buckets of BAC levels:

```text
S: BAC = 0
L: 0 < BAC < 0.05
M: 0.05 <= BAC < 0.08
H: BAC >= 0.08
```

## `10_four_group_prep.R`

Prepares crash-count data for the four-BAC-group LP-style MLE. The script uses `ANALYSIS_BAC_VALUE`, so measured BAC is retained when available and MIPER BAC is used otherwise.

The main sample uses MI BAC for all state-years. With the current option:

```r
rur_urb_keep_clear_only <- TRUE
```

the prepared analysis data retain only observations classified clearly as rural or urban, including for the broad `all_data` and `weekend` specifications.

### Cell specifications

To limit sparsity, the current script prepares only:

```text
all data
weekend
rural/urban
rural/urban x weekend
```

It reports sparsity under:

```text
min_n2 = 1, 5, 10, and 25
```

but does not fit the MLE yet.

## `11_four_group_mle.R`

Estimates the four-BAC-group LP-style model separately for each imputation and cell specification.

### Two-car stage

The script estimates:

```text
theta_L
theta_M
theta_H
```

and cell-specific:

```text
r_Lc
r_Mc
r_Hc
```

with:

```text
theta_S = 1
r_S = 1
```

### One-car stage

Using the exposure ratios estimated in the two-car stage, the script estimates:

```text
lambda_L
lambda_M
lambda_H
```

with:

```text
lambda_S = 1
```

The current two-car cell threshold is:

```r
min_n2 <- 1
```

After fitting all ten imputations, global `theta` and `lambda` parameters are pooled on the log scale using Rubin-style within- and between-imputation variance calculations.

## `12_four_group_results.R`

Creates figures and comparison tables from the four-group MLE.

The script:

- plots pooled four-group estimates across all cell specifications;
- displays imputation-level variability;
- highlights the preferred rural/urban-by-weekend specification;
- compares the four-group model with the binary models from Script 07; and
- summarizes optimizer and Hessian diagnostics.
