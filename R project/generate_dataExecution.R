library(survival)
library(relsurv)

set.seed(12345)

params <- list(
  lambda         = 0.03,
  # Baseline excess hazard
  age_option     = "A",
  # Age distribution option
  n              = 10000,
  # Number of patients to simulate
  max_time       = 5,
  # Administrative censoring at 5 years
  prop_female    = 0.5,
  # 50% women in the dataset
  prop_x0        = 0.5,
  # Proportion of X=0 in the dataset
  year.start_min = 2000,
  # Minimum diagnosis year
  year.start_max = 2005,
  # Maximum diagnosis year
  beta_sex       = -1.5,
  # Effect of sex on excess mortality
  beta_age       = 0.02,
  # Effect of age on excess mortality
  beta_X         = 0.01,
  # Effect of binary covariate X on excess mortality
  borne_a        = 6       # Maximum random censoring time (uniform[0, 6])
)

simulated_data <- generate_data(
  lambda         = params$lambda,
  age_option     = params$age_option,
  n              = params$n,
  max_time       = params$max_time,
  prop_female    = params$prop_female,
  prop_x0        = params$prop_x0,
  year.start_min = params$year.start_min,
  year.start_max = params$year.start_max,
  beta_sex       = params$beta_sex,
  beta_age       = params$beta_age,
  beta_X         = params$beta_X,
  borne_a        = params$borne_a
)

head(simulated_data)

#proportion of events vs. censored observations
table(simulated_data$status)

#distributions of age, sex, and X
summary(simulated_data$age)
hist(simulated_data$age)
prop.table(table(simulated_data$sex))
prop.table(table(simulated_data$race_num))

#time conversion
simulated_data$observed_time_days <- simulated_data$observed_time * 365.241
simulated_data$hypo_time_days     <- simulated_data$hypothetical_time * 365.241
simulated_data$age_days           <- simulated_data$age * 365.241

# model estimates
KM_estimate <- survfit(Surv(hypo_time_days, hypothetical_status) ~ 1, data = simulated_data)

netSurv_estimate <- rs.surv(
  Surv(observed_time_days, status) ~ 1,
  data = simulated_data,
  ratetable = survexp.usr,
  rmap = list(
    age = age_days,
    sex = sex,
    race = race,
    year = year_diagnosis
  ),
  method = "pohar-perme"
)

# survival curves
plot(
  netSurv_estimate,
  conf.int = TRUE,
  col = "blue",
  lwd = 2,
  xscale = 365.241,
  xlab = "Time since diagnosis (Years)",
  ylab = "Net Survival",
  main = "Estimated vs Theoretical Net Survival",
  ylim = c(0, 1)
)

lines(
  KM_estimate,
  conf.int = TRUE,
  col = "red",
  lwd = 2,
  lty = 2,
  xscale = 365.241
)

grid()
legend(
  "bottomright",
  legend = c("Pohar-Perme Net Survival", "Theoretical Net Survival (KM)"),
  col = c("blue", "red"),
  lwd = c(2, 2),
  lty = c(1, 2),
  bty = "n"
)

#nessie
#Create age groups
breaks <- pretty(simulated_data$age, n = 5)
simulated_data$agegr <- cut(simulated_data$age,
                            breaks = breaks,
                            include.lowest = TRUE)

nessie_output <- nessie(
  Surv(observed_time_days, status) ~ sex + agegr,
  data = simulated_data,
  ratetable = survexp.usr,
  times = seq(0, params$max_time, 1),
  rmap = list(age = age_days, sex = sex, year = year_diagnosis)
)

print(nessie_output)

# additional calculations
total_deaths <- sum(simulated_data$status == 1)
cancer_deaths <- sum(simulated_data$cause == 1)
other_deaths <- sum(simulated_data$status == 1 &
                      simulated_data$cause == 0)

#proportions|dead
prop_cancer_among_dead <- cancer_deaths / total_deaths
prop_other_among_dead  <- other_deaths / total_deaths
