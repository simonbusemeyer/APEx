# ---------------------------------------------------------
# 10. EMPIRICAL STABILITY ACROSS POOLED SIMULATIONS
# ---------------------------------------------------------
# 1. Create a data frame to store the proportions for each simulation
N_sim <- max(all_simulated_data$sim_id)
sim_stats <- data.frame(
  sim_id       = 1:N_sim,
  total_deaths = numeric(N_sim),
  p_cancer     = numeric(N_sim),
  p_other      = numeric(N_sim)
)

# 2. Loop through each simulation to calculate its specific proportions
for(i in 1:N_sim){
  dat <- subset(all_simulated_data, sim_id == i)
  
  td <- sum(dat$status == 1)
  cd <- sum(dat$status == 1 & dat$cause == 1)
  od <- sum(dat$status == 1 & dat$cause == 0)
  
  sim_stats$total_deaths[i] <- td
  # Avoid division by zero if a simulation randomly had 0 deaths
  sim_stats$p_cancer[i] <- ifelse(td > 0, cd / td, NA)
  sim_stats$p_other[i]    <- ifelse(td > 0, od / td, NA)
}

# 3. Calculate empirical metrics across the N simulations
mean_p_cancer <- mean(sim_stats$p_cancer, na.rm = TRUE)
sd_p_cancer   <- sd(sim_stats$p_cancer, na.rm = TRUE)

# Calculate Empirical 95% Confidence Intervals using the 2.5th and 97.5th percentiles
ci_lower <- quantile(sim_stats$p_cancer, 0.025, na.rm = TRUE)
ci_upper <- quantile(sim_stats$p_cancer, 0.975, na.rm = TRUE)

# Calculate Empirical RSE (Coefficient of Variation)
empirical_rse <- (sd_p_cancer / mean_p_cancer) * 100

# 4. Compile the evaluation table
pooled_eval <- data.frame(
  Metric                = "Empirical_Prop_Cancer",
  Mean_Estimate         = round(mean_p_cancer, 4),
  Empirical_SD          = round(sd_p_cancer, 4),
  Empirical_CI_Lower    = round(as.numeric(ci_lower), 4),
  Empirical_CI_Upper    = round(as.numeric(ci_upper), 4),
  Empirical_CI_Width    = round(as.numeric(ci_upper - ci_lower), 4),
  Empirical_RSE_Percent = round(empirical_rse, 2)
)

print(pooled_eval)

# Optional: Visualize the distribution of the proportions
hist(
  sim_stats$p_cancer, 
  breaks = 15, 
  col = "lightcoral", 
  main = "Distribution of Cancer Death Proportions Across Simulations",
  xlab = "Proportion of Cancer Deaths"
)
abline(v = mean_p_cancer, col = "red", lwd = 2, lty = 2)