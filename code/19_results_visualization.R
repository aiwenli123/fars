################################################################################
# Summarize counts of crash types by analysis cell and plot the two- and four-group
# bootstrap estimates used in the paper.
################################################################################

rm(list = ls())

library(tidyverse)

cell_counts_path <- "data/processed/script10_four_group_prep/four_group_crash_counts_by_cell_imp.csv"
two_group_results_path <- "results/two_group_bootstrap/two_group_bootstrap_combined_results.csv"
four_group_results_path <- "results/four_group_bootstrap/four_group_bootstrap_combined_results.csv"

table_dir <- "results/results_visualization"
figure_dir <- "figures/results_visualization"
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# Summarize the crash-type counts separately within the four cells.
# Means are probably the main summaries, but medians are also computed here.
cell_labels <- c(
  "rural_FALSE" = "Rural weekday",
  "rural_TRUE" = "Rural weekend",
  "urban_FALSE" = "Urban weekday",
  "urban_TRUE" = "Urban weekend"
)

one_vehicle_types <- c("one_S", "one_L", "one_M", "one_H")
two_vehicle_types <- c("SS", "SL", "SM", "SH", "LL", "LM", "LH", "MM", "MH", "HH")

cell_count_summary <- read_csv(cell_counts_path, show_col_types = FALSE) %>%
  filter(cell_spec == "rur_urb_weekend") %>%
  select(IMPUTATION, cell_id, all_of(c(one_vehicle_types, two_vehicle_types))) %>%
  pivot_longer(all_of(c(one_vehicle_types, two_vehicle_types)), names_to = "count_variable", values_to = "count") %>%
  mutate(
    crash_size = if_else(str_starts(count_variable, "one_"), "One vehicle", "Two vehicles"),
    crash_type = str_remove(count_variable, "^one_"),
    cell = recode(cell_id, !!!cell_labels),
    cell = factor(cell, levels = unname(cell_labels)),
    crash_size = factor(crash_size, levels = c("One vehicle", "Two vehicles")),
    crash_type = factor(crash_type, levels = c("S", "L", "M", "H", two_vehicle_types))
  ) %>%
  group_by(cell, crash_size, crash_type) %>%
  summarize(mean_count = mean(count), median_count = median(count), .groups = "drop") %>%
  arrange(crash_size, crash_type, cell)

mean_count_table <- cell_count_summary %>%
  select(crash_size, crash_type, cell, mean_count) %>%
  pivot_wider(names_from = cell, values_from = mean_count) %>%
  arrange(crash_size, crash_type)

median_count_table <- cell_count_summary %>%
  select(crash_size, crash_type, cell, median_count) %>%
  pivot_wider(names_from = cell, values_from = median_count) %>%
  arrange(crash_size, crash_type)

write_csv(cell_count_summary, file.path(table_dir, "crash_type_counts_by_cell_summary.csv"))
write_csv(mean_count_table, file.path(table_dir, "mean_crash_type_counts_by_cell.csv"))
write_csv(median_count_table, file.path(table_dir, "median_crash_type_counts_by_cell.csv"))

# Shared visual style for the two relative-risk figures.
risk_colors <- c("Two-vehicle crashes" = "#0072B2", "One-vehicle crashes" = "#D55E00")

risk_figure_theme <- theme_minimal(base_size = 12) +
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

# Two-group benchmark figure.
two_group_plot_data <- read_csv(two_group_results_path, show_col_types = FALSE) %>%
  mutate(
    crash_size = recode(parameter, "theta" = "Two-vehicle crashes", "lambda" = "One-vehicle crashes"),
    threshold_label = recode(threshold_label, "BAC >= 0.05" = "BAC ≥ 0.05", "BAC >= 0.08" = "BAC ≥ 0.08"),
    threshold_label = factor(threshold_label, levels = c("BAC ≥ 0.05", "BAC ≥ 0.08")),
    crash_size = factor(crash_size, levels = names(risk_colors))
  )

two_group_figure <- ggplot(two_group_plot_data, aes(x = estimate, y = threshold_label, color = crash_size)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.5, color = "grey45") +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0, linewidth = 0.8, position = position_dodge(width = 0.42)) +
  geom_point(size = 2.8, position = position_dodge(width = 0.42)) +
  scale_color_manual(values = risk_colors) +
  scale_x_continuous(expand = expansion(mult = c(0.04, 0.08))) +
  labs(
    title = "Two-group relative fatal-crash risk",
    subtitle = "Pooled estimates across 10 completed datasets with 95% confidence intervals",
    x = "Relative fatal-crash risk",
    y = "Drinking definition",
    caption = "The reference group consists of drivers below the indicated BAC cutoff; the dashed line marks equal risk."
  ) +
  risk_figure_theme

ggsave(file.path(figure_dir, "two_group_bootstrap_relative_risk.png"), two_group_figure, width = 8, height = 4.8, dpi = 400, bg = "white")
ggsave(file.path(figure_dir, "two_group_bootstrap_relative_risk.pdf"), two_group_figure, width = 8, height = 4.8, device = "pdf")

# Four-group primary-model figure.
four_group_plot_data <- read_csv(four_group_results_path, show_col_types = FALSE) %>%
  mutate(
    bac_group = case_when(
      str_ends(parameter, "_L") ~ "0 < BAC < 0.05",
      str_ends(parameter, "_M") ~ "0.05 ≤ BAC < 0.08",
      str_ends(parameter, "_H") ~ "BAC ≥ 0.08"
    ),
    crash_size = recode(stage, "two_car" = "Two-vehicle crashes", "one_car" = "One-vehicle crashes"),
    bac_group = factor(bac_group, levels = c("0 < BAC < 0.05", "0.05 ≤ BAC < 0.08", "BAC ≥ 0.08")),
    crash_size = factor(crash_size, levels = names(risk_colors))
  )

four_group_figure <- ggplot(four_group_plot_data, aes(x = estimate, y = bac_group, color = crash_size)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.5, color = "grey45") +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0, linewidth = 0.8, position = position_dodge(width = 0.42)) +
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
  risk_figure_theme

ggsave(file.path(figure_dir, "four_group_bootstrap_relative_risk.png"), four_group_figure, width = 8, height = 5.2, dpi = 400, bg = "white")
ggsave(file.path(figure_dir, "four_group_bootstrap_relative_risk.pdf"), four_group_figure, width = 8, height = 5.2, device = "pdf")

