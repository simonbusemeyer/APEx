library(tidyverse)
library(scales)

# Convert time_t to a categorical factor safely
results_df <- results_df %>%
  mutate(time_t_factor = factor(time_t, labels = paste("Year", sort(unique(time_t)))))

ggplot(results_df, aes(x = pct_cancer, y = ecr_unconditional)) +
  geom_hline(aes(yintercept = 0.95, linetype = "Nominal Target (95%)"), color = "darkred", alpha = 0.8) +
  # Map both color and shape to the follow-up year factor, and remove geom_line()
  geom_point(aes(color = time_t_factor, shape = time_t_factor), size = 3) +
  
  # Scales
  scale_color_viridis_d(option = "plasma", end = 0.8) + 
  scale_shape_discrete() + # Automatically assigns distinct professional shapes
  scale_linetype_manual(values = c("Nominal Target (95%)" = "dashed")) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  scale_y_continuous(breaks = seq(0.6, 1.0, by = 0.05)) +
  coord_cartesian(ylim = c(0.6, 1.0)) + 
  
  # Labels
  labs(
    title = "Estimator Reliability Under Competing Risks",
    x = "Proportion of Cohort Experiencing Cancer Death vs Non-Cancer Death",
    y = "Empirical Coverage Rate",
    color = "Follow-up",
    shape = "Follow-up", # Combined legend by using identical labels
    linetype = NULL 
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    
    # Maintain uniform font sizes at the bottom
    axis.title.x = element_text(size = 14),
    axis.text.x  = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 14),
    
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )