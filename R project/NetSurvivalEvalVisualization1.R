# ---------------------------------------------------------
# 6. OVERLAY PLOT: INDIVIDUAL PP, POOLED PP, & THEORETICAL
# ---------------------------------------------------------

# Set up a blank plot canvas
plot(
  0, type = "n",
  xlim = c(0, params$max_time * 365.241),
  ylim = c(0.7, 1), # Adjust depending on your mortality rate
  xlab = "Time since diagnosis (Days)",
  ylab = "Net Survival Probability",
  main = paste("Pohar-Perme vs Theoretical (", N_sim, " Simulations)")
)
grid()

# ---------------------------------------------------------
# 1. Plot Individual Pohar-Perme Curves (Loop)
# ---------------------------------------------------------
for (i in 1:N_sim) {
  data_sim <- subset(all_simulated_data, sim_id == i)
  
  # Calculate PP for this specific iteration
  pp_sim <- rs.surv(
    Surv(observed_time_days, status) ~ 1,
    data = data_sim,
    ratetable = survexp.usr,
    rmap = list(age = age_days, sex = sex, race = race, year = year_diagnosis),
    method = "pohar-perme"
  )
  
  # Plot as a step function ("s") with semi-transparency
  lines(
    pp_sim$time, pp_sim$surv, 
    type = "s", 
    col = rgb(0.2, 0.5, 0.8, alpha = 0.3), 
    lwd = 1
  )
}

# ---------------------------------------------------------
# 2. Plot Pooled Pohar-Perme Curve
# ---------------------------------------------------------
pp_pooled <- rs.surv(
  Surv(observed_time_days, status) ~ 1,
  data = all_simulated_data,
  ratetable = survexp.usr,
  rmap = list(age = age_days, sex = sex, race = race, year = year_diagnosis),
  method = "pohar-perme"
)

# Plot pooled curve in bold red
lines(
  pp_pooled$time, pp_pooled$surv, 
  type = "s", 
  col = "red", 
  lwd = 3
)

# ---------------------------------------------------------
# 3. Plot Marginal Theoretical Curve
# ---------------------------------------------------------
# Define a time sequence in years (for the math) and days (for the plot x-axis)
t_seq_years <- seq(0, params$max_time, by = 0.05)
t_seq_days <- t_seq_years * 365.241

# Calculate the theoretical survival curve for EVERY patient individually
ind_theo_curves <- sapply(1:nrow(all_simulated_data), function(i) {
  exp(-params$lambda * t_seq_years * exp(
    params$beta_sex * all_simulated_data$sex_num[i] +
      params$beta_age * all_simulated_data$ageStand[i] +
      params$beta_X   * all_simulated_data$race_num[i]
  ))
})

# Average the individual curves to get the true marginal population survival
surv_theo_mean <- rowMeans(ind_theo_curves)

# Plot theoretical curve as a dashed black line
lines(
  t_seq_days, surv_theo_mean, 
  col = "black", 
  lwd = 3, 
  lty = 2
)

# ---------------------------------------------------------
# Add Legend
# ---------------------------------------------------------
legend(
  "bottomleft", 
  legend = c("Individual PP Runs", "Pooled PP Estimate", "Marginal Theoretical S(t)"),
  col = c(rgb(0.2, 0.5, 0.8, alpha = 0.5), "red", "black"), 
  lwd = c(2, 3, 3), 
  lty = c(1, 1, 2),
  bty = "n"
)