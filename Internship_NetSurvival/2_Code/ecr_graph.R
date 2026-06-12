# Load required libraries
library(tidyverse)
library(scales)

# Read the dataset
results_df <- read_csv("outputs/tables/final_results_complete.csv")

# Convert time_t to a categorical factor for stratification
results_df <- results_df %>%
  mutate(time_t_factor = factor(time_t, labels = paste("Year", unique(sort(time_t)))))

# Generate the plot
ggplot(results_df, aes(x = pct_cancer, y = ecr, color = time_t_factor, group = time_t_factor)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "darkred", alpha = 0.7) +
  geom_line(linewidth = 1.2, alpha = 0.8) +
  geom_point(size = 2.5) +
  scale_color_viridis_d(option = "plasma", end = 0.8) + # Professional color palette
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  scale_y_continuous(limits = c(0.5, 1.05), breaks = seq(0, 1, by = 0.2)) +
  labs(
    title = "Estimator Reliability Under Competing Risks",
    subtitle = "ECR vs. Cohort Cancer Incidence (pct_cancer)",
    x = "Proportion of Cohort Experiencing Cancer Death",
    y = "Empirical Coverage Rate (ECR)",
    color = "Follow-up Time"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )
