library(ggplot2)
library(dplyr)

# Focusing on the vulnerable 85+ age group as an example
ggplot(metrics, aes(x = time_t)) +
  # The "Box" (IQR) - bounded by Median +/- (IQR/2) assuming symmetry, 
  # or you can plot actual 25th/75th percentiles if you add them to compute_metrics.R
  geom_ribbon(aes(ymin = pmax(0, med_n_at_risk_o85 - (iqr_n_at_risk_o85/2)), 
                  ymax = med_n_at_risk_o85 + (iqr_n_at_risk_o85/2)), 
              fill = "#2980b9", alpha = 0.4) +
  # The "Median Line"
  geom_line(aes(y = med_n_at_risk_o85), color = "#2c3e50", linewidth = 1) +
  facet_wrap(~ lambda, scales = "free_y") +
  theme_minimal() +
  labs(title = "Non-Parametric Attrition Track (85+ Age Group)",
       subtitle = "Solid line = Median; Shaded band = Approximated IQR (Middle 50% of simulations)",
       x = "Time (Years)", y = "Median Patients at Risk")