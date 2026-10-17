# Record crash-type counts and make a figure for the four-group model results.
# The table created here summarizes four-group crash counts over ten BAC imputations;
# The figure here plots the pooled four-group risks and 95% confidence intervals.
#
# Inputs:  data/processed/script03_four_group_prep/four_group_crash_counts_by_cell_imp.csv
#          results/four_group_bootstrap/four_group_bootstrap_combined_results.csv
# Outputs: results/results_visualization/mean_crash_type_counts_by_cell.csv
#          figures/results_visualization/four_group_bootstrap_relative_risk.{png,pdf}

# rm(list=ls())
library(tidyverse)

counts_path <- "data/processed/script03_four_group_prep/four_group_crash_counts_by_cell_imp.csv"
risk_path <- "results/four_group_bootstrap/four_group_bootstrap_combined_results.csv"
table_dir <- "results/results_visualization"
figure_dir <- "figures/results_visualization"

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# Average counts in each of the four rural/urban-by-weekend cells.
cell_labels <- c(
  rural_FALSE = "Rural weekday",
  rural_TRUE  = "Rural weekend",
  urban_FALSE = "Urban weekday",
  urban_TRUE  = "Urban weekend"
)
one_types <- paste0("one_", c("S", "L", "M", "H"))
two_types <- c("SS", "SL", "SM", "SH", "LL", "LM", "LH", "MM", "MH", "HH")

mean_counts <- read_csv(counts_path, show_col_types = FALSE) %>%
  filter(cell_spec == "rur_urb_weekend") %>%
  select(IMPUTATION, cell_id, all_of(c(one_types, two_types))) %>%
  pivot_longer(all_of(c(one_types, two_types)),
               names_to = "count_variable", values_to = "count") %>%
  mutate(
    crash_size = if_else(str_starts(count_variable, "one_"),
                         "One vehicle", "Two vehicles"),
    crash_type = str_remove(count_variable, "^one_"),
    cell = recode(cell_id, !!!cell_labels),
    cell = factor(cell, levels = unname(cell_labels)),
    crash_size = factor(crash_size, levels = c("One vehicle", "Two vehicles")),
    crash_type = factor(crash_type, levels = c("S", "L", "M", "H", two_types))
  ) %>%
  group_by(cell, crash_size, crash_type) %>%
  summarise(mean_count = mean(count), .groups = "drop") %>%
  select(crash_size, crash_type, cell, mean_count) %>%
  pivot_wider(names_from = cell, values_from = mean_count) %>%
  arrange(crash_size, crash_type)

write_csv(mean_counts, file.path(table_dir, "mean_crash_type_counts_by_cell.csv"))

# Plot of four-group estimates with 95% Rubin intervals
risk_colors <- c("Two-vehicle crashes" = "#0072B2",
                 "One-vehicle crashes" = "#D55E00")

plot_data <- read_csv(risk_path, show_col_types = FALSE) %>%
  filter(cell_spec == "rur_urb_weekend") %>%
  mutate(
    bac_group = case_when(
      str_ends(parameter, "_L") ~ "0 < BAC < 0.05",
      str_ends(parameter, "_M") ~ "0.05 ≤ BAC < 0.08",
      str_ends(parameter, "_H") ~ "BAC ≥ 0.08"
    ),
    crash_size = recode(stage,
                        two_car = "Two-vehicle crashes",
                        one_car = "One-vehicle crashes"),
    bac_group = factor(bac_group,
                       levels = c("0 < BAC < 0.05", "0.05 ≤ BAC < 0.08", "BAC ≥ 0.08")),
    crash_size = factor(crash_size, levels = names(risk_colors))
  )

risk_figure <- ggplot(plot_data,
                      aes(x = estimate, y = bac_group, color = crash_size)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.5, color = "grey45") +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0,
                linewidth = 0.8, position = position_dodge(width = 0.42)) +
  geom_point(size = 2.8, position = position_dodge(width = 0.42)) +
  scale_color_manual(values = risk_colors) +
  scale_x_continuous(expand = expansion(mult = c(0.04, 0.08))) +
  labs(
    title = "Four-group relative fatal-crash risk",
    subtitle = "Pooled estimates across 10 completed datasets with 95% confidence intervals",
    x = "Relative fatal-crash risk",
    y = "BAC group",
    caption = "Drivers with BAC = 0 are the reference group; the dashed line marks equal risk."
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "grey15"),
    axis.title.x = element_text(face = "bold"),
    legend.position = "top",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "grey30"),
    plot.caption = element_text(color = "grey35", hjust = 0),
    plot.margin = margin(10, 16, 10, 10)
  )

ggsave(file.path(figure_dir, "four_group_bootstrap_relative_risk.png"),
       risk_figure, width = 8, height = 5.2, dpi = 400, bg = "white")
ggsave(file.path(figure_dir, "four_group_bootstrap_relative_risk.pdf"),
       risk_figure, width = 8, height = 5.2, device = "pdf")
