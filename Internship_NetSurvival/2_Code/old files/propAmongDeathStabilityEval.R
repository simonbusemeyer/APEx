# ---------------------------------------------------------
# 9. STABILITY OF CAUSE-OF-DEATH PROPORTIONS
# ---------------------------------------------------------
# Ensure your variables are defined (from your previous code)
total_deaths <- sum(simulated_data$status == 1)
cancer_deaths <- sum(simulated_data$cause == 1)
other_deaths <- sum(simulated_data$status == 1 & simulated_data$cause == 0)

# 1. Calculate the proportions (Estimate)
p_cancer <- cancer_deaths / total_deaths
p_other  <- other_deaths / total_deaths

# 2. Calculate the Absolute Standard Error (SE) using Binomial formula
se_cancer <- sqrt((p_cancer * (1 - p_cancer)) / total_deaths)
se_other  <- sqrt((p_other * (1 - p_other)) / total_deaths)

# 3. Calculate 95% Confidence Intervals (Normal Approximation / Wald)
ci_lower_cancer <- pmax(0, p_cancer - (1.96 * se_cancer))
ci_upper_cancer <- pmin(1, p_cancer + (1.96 * se_cancer))

ci_lower_other <- pmax(0, p_other - (1.96 * se_other))
ci_upper_other <- pmin(1, p_other + (1.96 * se_other))

# 4. Compile into a comprehensive evaluation table
prop_eval_table <- data.frame(
  Metric = c("Prop_Cancer_among_dead", "Prop_Other_among_dead"),
  Count  = c(cancer_deaths, other_deaths),
  Total_Deaths = c(total_deaths, total_deaths),
  Estimate = c(p_cancer, p_other),
  Absolute_SE = c(se_cancer, se_other),
  CI_Width = c(ci_upper_cancer - ci_lower_cancer, ci_upper_other - ci_lower_other),
  RSE_Percent = c((se_cancer / p_cancer) * 100, (se_other / p_other) * 100)
)

# Round the numerical columns for clean viewing
prop_eval_table[, 4:7] <- round(prop_eval_table[, 4:7], 4)
print(prop_eval_table)