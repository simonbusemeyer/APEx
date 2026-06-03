# Calculate the crude probability of death
crude_mortality <- cmp.rel(
  Surv(observed_time_days, status) ~ 1, 
  data = simulated_data, 
  ratetable = survexp.usr, 
  rmap = list(age = age_days, sex = sex, race = race, year = year_diagnosis)
)

# Extract the estimates and their variances at exactly 5 and 10 years
summary(crude_mortality, times = c(5 * 365.241, 10 * 365.241))



# Step 2: Plot the curves
plot(
  crude_mortality,
  col = c("lightcoral", "lightblue"), # Colors for cancer vs. other causes
  conf.int = c(1, 2),                 # 1 = Cancer CI, 2 = Other causes CI
  xscale = 365.241,                   # Convert the x-axis from days to years
  xlab = "Time since diagnosis (Years)",
  ylab = "Crude Probability of Death",
  main = "Crude Mortality (Cancer vs. Other Causes)",
  lwd = 2                             # Line width
)

# Step 3: Add a clear legend
legend(
  "topleft",
  legend = c("Death from Cancer", "Death from Other Causes"),
  col = c("lightcoral", "lightblue"),
  lwd = 2,
  lty = 1,
  bty = "n"
)