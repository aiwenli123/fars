# Apply the Assumption 3 sensitivity scenarios to the four-group
# relative-risk estimates and their confidence intervals.
# The scenario parameters are fixed, so no models or bootstraps are refitted.
#
# Input:  results/four_group_bootstrap/four_group_bootstrap_combined_results.csv
# Output: results/assumption3_sensitivity/assumption3_scenario_results.csv

# rm(list=ls())
library(tidyverse)

input_path <- "results/four_group_bootstrap/four_group_bootstrap_combined_results.csv"
out_dir <- "results/assumption3_sensitivity"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Scenario-specific relative inability to avoid crashes (mu) and error
# seriousness (delta). 
scenarios <- tribble(
  ~scenario, ~scenario_label, ~BAC_group, ~mu, ~delta,
  1, "Baseline",                                       "L", 1.00, 1.00,
  1, "Baseline",                                       "M", 1.00, 1.00,
  1, "Baseline",                                       "H", 1.00, 1.00,
  2, "25% greater inability to avoid",                 "L", 1.25, 1.00,
  2, "25% greater inability to avoid",                 "M", 1.25, 1.00,
  2, "25% greater inability to avoid",                 "H", 1.25, 1.00,
  3, "Avoidance inability rises with BAC",             "L", 1.05, 1.00,
  3, "Avoidance inability rises with BAC",             "M", 1.15, 1.00,
  3, "Avoidance inability rises with BAC",             "H", 1.25, 1.00,
  4, "Avoidance inability and seriousness rise equally", "L", 1.05, 1.05,
  4, "Avoidance inability and seriousness rise equally", "M", 1.15, 1.15,
  4, "Avoidance inability and seriousness rise equally", "H", 1.25, 1.25,
  5, "Seriousness rises slightly faster",              "L", 1.05, 1.10,
  5, "Seriousness rises slightly faster",              "M", 1.15, 1.20,
  5, "Seriousness rises slightly faster",              "H", 1.25, 1.30
)

# Put the two- and one-vehicle estimates (and confidence limits) side by side.
primary_results <- read_csv(input_path, show_col_types = FALSE) %>%
  filter(cell_spec == "rur_urb_weekend", stage %in% c("two_car", "one_car")) %>%
  transmute(BAC_group = str_extract(parameter, "[LMH]$"),
            stage, estimate, ci_low, ci_high) %>%
  pivot_wider(names_from = stage, values_from = c(estimate, ci_low, ci_high))

# From Appendix D: theta_unadj = (mu / delta) * theta_adj and
# lambda_unadj = mu * lambda_adj. Apply the same factors to confidence limits.
scenario_results <- scenarios %>%
  left_join(primary_results, by = "BAC_group") %>%
  mutate(
    two_vehicle_risk = (mu / delta) * estimate_two_car,
    two_vehicle_ci_low = (mu / delta) * ci_low_two_car,
    two_vehicle_ci_high = (mu / delta) * ci_high_two_car,
    one_vehicle_risk = mu * estimate_one_car,
    one_vehicle_ci_low = mu * ci_low_one_car,
    one_vehicle_ci_high = mu * ci_high_one_car,
    BAC_group_label = recode(BAC_group,
                             "L" = "0 < BAC < 0.05",
                             "M" = "0.05 <= BAC < 0.08",
                             "H" = "BAC >= 0.08")
  ) %>%
  select(scenario, scenario_label, BAC_group, BAC_group_label, mu, delta,
         two_vehicle_risk, two_vehicle_ci_low, two_vehicle_ci_high,
         one_vehicle_risk, one_vehicle_ci_low, one_vehicle_ci_high) %>%
  arrange(scenario, factor(BAC_group, levels = c("L", "M", "H")))

# Save the complete estimates and intervals
write_csv(scenario_results, file.path(out_dir, "assumption3_scenario_results.csv"))
print(scenario_results, n = Inf, width = Inf)
