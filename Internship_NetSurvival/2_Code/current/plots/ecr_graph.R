library(tidyverse)
library(scales)

results_df <- read_csv("current/outputs/tables/final_results_complete.csv")

# Convert time_t to a categorical factor for stratification
results_df <- results_df %>%
  mutate(time_t_factor = factor(time_t, labels = paste("Year", unique(sort(time_t)))))

# Generate the plot
ggplot(results_df, aes(x = pct_cancer, y = ecr_unconditional)) +
  geom_hline(aes(yintercept = 0.95, linetype = "Nominal Target (95%)"), color = "darkred", alpha = 0.7) +
  geom_line(aes(color = time_t_factor, group = time_t_factor), linewidth = 1.2, alpha = 0.8) +
  geom_point(aes(color = time_t_factor), size = 2.5) +
  
  # Scales
  scale_color_viridis_d(option = "plasma", end = 0.8) + 
  scale_linetype_manual(values = c("Nominal Target (95%)" = "dashed")) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  scale_y_continuous(limits = c(0.5, 1.0), breaks = seq(0.6, 1.0, by = 0.2)) +
  
  # Labels
  labs(
    title = "Estimator Reliability Under Competing Risks",
    x = "Proportion of Cohort Experiencing Cancer Death",
    y = "Empirical Coverage Rate (ECR)",
    color = "Follow-up Time",
    linetype = NULL 
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.3),
    
    plot.subtitle = element_text(size = 12, hjust = 0.4),
    axis.title.x  = element_text(size = 12),
    axis.title.y  = element_text(size = 12),
    axis.text.x   = element_text(size = 12),
    axis.text.y   = element_text(size = 12),
    legend.title  = element_text(size = 12),
    legend.text   = element_text(size = 12),
    
    # Legend layout alignment
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.box.just = "center",
    legend.spacing.x = unit(0.2, "cm"),
    
    panel.grid.minor = element_blank()
  )