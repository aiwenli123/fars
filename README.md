# FARS Levitt–Porter Replication and BAC Extensions

This repository contains R code for estimating relative fatal crash risk across blood alcohol concentration (BAC) groups using the 2015–2024 Fatality Analysis Reporting System (FARS). We extend the Levitt–Porter (2001) model to four BAC groups, focusing on drivers with BAC from 0.05 to below 0.08, without requiring a separate sample of non-crash drivers.

## Repository structure

```text
.
├── code/
│   ├── 01_clean_data.R
│   ├── 02_prep_mi_data.R
│   ├── 03_four_group_prep.R
│   ├── 04_two_group_bootstrap.R
│   ├── 05_four_group_bootstrap_parallel.R
│   ├── 06_assumption2_bootstrap_parallel.R
│   ├── 07_assumption3_sensitivity.R
│   ├── 08_results_visualization.R
│   ├── 09a_simulation_full.R
│   ├── 09b_combine_simulation_batches.R
│   └── hpc/
│       └── fars_bootstrap.sbatch
├── data/
│   ├── raw/
│   └── processed/
├── results/
└── figures/
```

- `code/`: Data preparation, estimation, bootstrap inference, sensitivity analyses, simulation validation, and visualization.
- `data/raw/`: Original FARS files from [NHTSA](https://www.nhtsa.gov/file-downloads?p=nhtsa/downloads/FARS/), not included in this repository.
- `data/processed/`: Cleaned driver-level datasets and crash-type classifications.
- `results/`: Model estimates, confidence intervals, sensitivity results, and simulation summaries.
- `figures/`: Figures generated for the paper.

## Analysis workflow

Scripts are numbered in their intended order:

| Scripts | Purpose |
|---|---|
| `01`–`03` | Prepare nighttime FARS data, incorporate ten BAC imputations, and classify crashes into four BAC groups. |
| `04` | Estimate two-group benchmarks at BAC cutoffs of 0.05 and 0.08. |
| `05` | Estimate the main four-group model with crash-level bootstrap inference. |
| `06`–`07` | Evaluate alternative assumptions about driver mixing, crash seriousness, and avoidance ability. |
| `08` | Generate the crash-count summary and main results figure. |
| `09a`–`09b` | Simulate fatal crash datasets, evaluate bootstrap confidence interval coverage, and combine simulation batches. |

## Reproducibility

The analysis uses R and the `tidyverse` package. Scripts should be run from the repository root, in numerical order.

1. Download the 2015–2024 FARS accident, person, and MIPER files from NHTSA. Extract each year's files into `data/raw/FARS{year}NationalCSV/`.
2. Run Scripts `01`–`08` to reproduce the main estimates, confidence intervals, sensitivity analyses, and figures. Scripts `04`–`06` involve computationally intensive bootstrapping.
3. Run Script `09a` as a 100-task Slurm job array using `code/hpc/fars_bootstrap.sbatch`, after configuring the cluster account and project directory. Each task simulates 10 datasets with 200 bootstrap samples per dataset.
4. Run Script `09b` after all Slurm tasks finish to combine the 1,000 simulated datasets and calculate empirical confidence interval coverage.

Scripts `04`–`06` can reuse existing bootstrap outputs to avoid unnecessary recomputation. The simulation scripts were run on the University of Pennsylvania's PARCC Betty cluster.

All generated paths for scripts `01`–`08` are relative to the repository root.
