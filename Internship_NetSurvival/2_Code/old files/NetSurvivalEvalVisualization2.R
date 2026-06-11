# ---------------------------------------------------------
# 5. VISUALIZING POHAR-PERME VS KM VS THEORETICAL
# ---------------------------------------------------------
message("Calculating estimates for all simulations. This may take a moment...")

# 1. Define the theoretical survival function (from UtilsExponential)
survtheo <- function(t, lambda, beta_sex, beta_age, sex, ageStand) {
  exp(-lambda * t * exp(beta_sex * sex + beta_age * ageStand))
}

# 2. Calculate the average theoretical net survival curve across the whole cohort
t_seq_years <- seq(0, params$max_time, by = 0.1)

surv_theo_matrix <- sapply(1:nrow(all_simulated_data), function(i) {
  survtheo(
    t         = t_seq_years,
    lambda    = params$lambda,
    beta_sex  = params$beta_sex, 
    beta_age  = params$beta_age,
 #   beta_X    = params$beta_X,
    sex       = all_simulated_data$sex_num[i],
    ageStand  = all_simulated_data$ageStand[i]
 #   X         = all_simulated_data$race_num[i]
  )
})

surv_theo_mean <- rowMeans(surv_theo_matrix)

# 3. Calculate the Pooled Pohar-Perme Estimate (Real-world estimate)
pp_pooled <- rs.surv(
  Surv(observed_time_days, status) ~ 1,
  data = all_simulated_data,
  ratetable = survexp.us,
  rmap = list(age = age_days, sex = sex, year = year_diagnosis),
  method = "pohar-perme"
)

# 4. Calculate the Pooled Kaplan-Meier Estimate (Hypothetical truth)
#km_pooled <- survfit(Surv(hypo_time_days, hypothetical_status) ~ 1, data = all_simulated_data)

# 5. Set up the blank plot
plot(
  0, type = "n",
  xlim = c(0, params$max_time),
  ylim = c(0.8, 1), # Adjust this lower bound depending on your lambda
  xlab = "Time since diagnosis (Years)",
  ylab = "Net Survival Probability",
main = paste("Net Survival: PP vs Theoretical (", N_sim, " runs) | Age ~ U(80-89) ", sep="")
)
grid()

# 6. Loop through and plot each individual simulation's Pohar-Perme curve
for (i in 1:N_sim) {
  data_sim <- subset(all_simulated_data, sim_id == i)
  
  pp_sim <- rs.surv(
    Surv(observed_time_days, status) ~ 1,
    data = data_sim,
    ratetable = survexp.us,
    rmap = list(age = age_days, sex = sex, year = year_diagnosis),
    method = "pohar-perme"
  )
  
  # Individual PP curves in faint blue
  lines(pp_sim$time / 365.241, pp_sim$surv, col = rgb(0.2, 0.5, 0.8, alpha = 0.2), lwd = 1, type = "s")
}

# 7. Overlay the Pooled Pohar-Perme curve (in bold Red)
lines(pp_pooled$time / 365.241, pp_pooled$surv, col = "red", lwd = 3, type = "s")

# 8. Overlay the Pooled Kaplan-Meier curve (in bold Green)
#lines(km_pooled$time / 365.241, km_pooled$surv, col = "green3", lwd = 3, type = "s")

# 9. Overlay the Theoretical Excess Hazard curve (in bold dashed Black)
lines(t_seq_years, surv_theo_mean, col = "black", lwd = 3, lty = 2)

# 10. Add a comprehensive legend
legend(
  "bottomleft", 
  legend = c(
    "Individual PP Estimates", 
    "Pooled PP Estimate", 
    "Theoretical Net Survival Curve S(t) = exp(-lambda * t * exp(beta * Z)) "
  ),
  col = c(rgb(0.2, 0.5, 0.8, alpha = 0.5), "red", "black"), 
  lwd = c(2, 3, 3, 3), 
  lty = c(1, 1, 1, 2), 
  bty = "n",
  cex = 0.9 # Slightly scale down text to fit cleanly
)