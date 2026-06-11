library(survival)
library(relsurv)
library(ggplot2)
library(dplyr)
library(tidyr)

# 1. Wrapper function to extract 5-year Net Survival
get_5yr_ns <- function(mod_params) {
  # Generate data with modified parameters
  sim_data <- generate_data(
    lambda = mod_params$lambda, age_option = mod_params$age_option, n = mod_params$n,
    max_time = mod_params$max_time, prop_female = mod_params$prop_female, prop_x0 = mod_params$prop_x0,
    year.start_min = mod_params$year.start_min, year.start_max = mod_params$year.start_max,
    beta_sex = mod_params$beta_sex, beta_age = mod_params$beta_age,
    beta_X = mod_params$beta_X, borne_a = mod_params$borne_a
  )
  
  sim_data$observed_time_days <- sim_data$observed_time * 365.241
  sim_data$age_days <- sim_data$age * 365.241
  
  # Fit net survival
  ns_est <- rs.surv(
    Surv(observed_time_days, status) ~ 1,
    data = sim_data, ratetable = survexp.usr,
    rmap = list(age = age_days, sex = sex, year = year_diagnosis),
    method = "pohar-perme"
  )
  
  # Extract survival probability closest to 5 years (1826 days)
  target_time <- 5 * 365.241
  idx <- which.min(abs(ns_est$time - target_time))
  return(ns_est$surv[idx])
}

# 2. Define Base Case and Bounds (+/- 20% for continuous variables)
base_params <- params
base_ns <- get_5yr_ns(base_params)

sa_vars <- list(
  lambda   = c(0.005, 0.055),
  beta_sex = c(-3.0, 0),
  beta_age = c(0.01, 0.03),
  beta_X   = c(0.008, 0.012)
)

# 3. Run One-Way Sensitivity Analysis
results <- data.frame(Parameter = character(), Low_Val = numeric(), High_Val = numeric(), stringsAsFactors = FALSE)

for (var in names(sa_vars)) {
  # Test Lower Bound
  params_low <- base_params
  params_low[[var]] <- sa_vars[[var]][1]
  ns_low <- get_5yr_ns(params_low)
  
  # Test Upper Bound
  params_high <- base_params
  params_high[[var]] <- sa_vars[[var]][2]
  ns_high <- get_5yr_ns(params_high)
  
  results <- rbind(results, data.frame(Parameter = var, Low_Val = ns_low, High_Val = ns_high))
}

# 4. Prepare data for Tornado Diagram
results <- results %>%
  mutate(
    Spread = abs(High_Val - Low_Val),
    Parameter = reorder(Parameter, Spread) # Order by impact magnitude
  ) %>%
  pivot_longer(cols = c(Low_Val, High_Val), names_to = "Bound", values_to = "NetSurvival")

# 5. Plot using ggplot2
ggplot(results, aes(x = Parameter, y = NetSurvival, fill = Bound)) +
  geom_bar(stat = "identity", position = "identity", width = 0.5, alpha = 0.8) +
  geom_hline(yintercept = base_ns, linetype = "dashed", color = "black", size = 1) +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Tornado Diagram: Deterministic Sensitivity Analysis",
    subtitle = paste("Impact of parameter variation on 5-Year Net Survival (Base Case =", round(base_ns, 3), ")"),
    x = "Parameters",
    y = "5-Year Net Survival Probability"
  ) +
  scale_fill_manual(values = c("Low_Val" = "steelblue", "High_Val" = "darkred"))