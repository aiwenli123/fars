################################################################################
# Calculate the Assumption 3 seriousness and crash-avoidance scenarios from the
# four-group bootstrap results.
################################################################################

rm(list = ls())

library(tidyverse)

input_path <- "results/four_group_bootstrap/four_group_bootstrap_combined_results.csv"
out_dir <- "results/assumption3_sensitivity"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Fixed values used in Table 6 of the paper.
scenarios <- tribble(
  ~scenario, ~scenario_label, ~BAC_group, ~mu, ~delta,
  1, "Baseline", "L", 1.00, 1.00,
  1, "Baseline", "M", 1.00, 1.00,
  1, "Baseline", "H", 1.00, 1.00,
  2, "25% greater inability to avoid", "L", 1.25, 1.00,
  2, "25% greater inability to avoid", "M", 1.25, 1.00,
  2, "25% greater inability to avoid", "H", 1.25, 1.00,
  3, "Avoidance inability rises with BAC", "L", 1.05, 1.00,
  3, "Avoidance inability rises with BAC", "M", 1.15, 1.00,
  3, "Avoidance inability rises with BAC", "H", 1.25, 1.00,
  4, "Avoidance inability and seriousness rise equally", "L", 1.05, 1.05,
  4, "Avoidance inability and seriousness rise equally", "M", 1.15, 1.15,
  4, "Avoidance inability and seriousness rise equally", "H", 1.25, 1.25,
  5, "Seriousness rises slightly faster", "L", 1.05, 1.10,
  5, "Seriousness rises slightly faster", "M", 1.15, 1.20,
  5, "Seriousness rises slightly faster", "H", 1.25, 1.30
)

# Put the two- and one-vehicle estimates from the four-group model on the same row.
primary_results <- read_csv(input_path, show_col_types = FALSE) %>%
  filter(cell_spec == "rur_urb_weekend", stage %in% c("two_car", "one_car")) %>%
  transmute(BAC_group = str_extract(parameter, "[LMH]$"), stage, estimate, ci_low, ci_high) %>%
  pivot_wider(names_from = stage, values_from = c(estimate, ci_low, ci_high))

# Fixed mu and delta values to rescale the estimates and confidence limits
# model doesn't need to be refitted
scenario_results <- scenarios %>%
  left_join(primary_results, by = "BAC_group") %>%
  mutate(
    two_vehicle_risk = (mu / delta) * estimate_two_car,
    two_vehicle_ci_low = (mu / delta) * ci_low_two_car,
    two_vehicle_ci_high = (mu / delta) * ci_high_two_car,
    one_vehicle_risk = mu * estimate_one_car,
    one_vehicle_ci_low = mu * ci_low_one_car,
    one_vehicle_ci_high = mu * ci_high_one_car,
    BAC_group_label = recode(BAC_group, "L" = "0 < BAC < 0.05", "M" = "0.05 <= BAC < 0.08", "H" = "BAC >= 0.08")
  ) %>%
  select(scenario, scenario_label, BAC_group, BAC_group_label, mu, delta, two_vehicle_risk, two_vehicle_ci_low, two_vehicle_ci_high, one_vehicle_risk, one_vehicle_ci_low, one_vehicle_ci_high) %>%
  arrange(scenario, factor(BAC_group, levels = c("L", "M", "H")))

# More compact table, cleaner to read
formatted_table <- scenario_results %>%
  transmute(
    Scenario = scenario,
    `BAC group` = BAC_group_label,
    mu = sprintf("%.2f", mu),
    delta = sprintf("%.2f", delta),
    `Two-vehicle risk` = sprintf("%.2f", two_vehicle_risk),
    `One-vehicle risk` = sprintf("%.2f", one_vehicle_risk)
  )

write_csv(scenario_results, file.path(out_dir, "assumption3_scenario_results.csv"))
write_csv(formatted_table, file.path(out_dir, "assumption3_table_formatted.csv"))
