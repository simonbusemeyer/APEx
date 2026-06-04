library(survival)
library(relsurv)

# ---------------------------------------------------------
# 1. SIMULATION PARAMETERS
# ---------------------------------------------------------
params <- list(
  lambda         = 0.1,
  age_option     = "C",    
  n              = 1000,   
  max_time       = 10,      
  prop_female    = 0.5,    
 prop_x0        = 0.5,    
  year.start_min = 2008,   
  year.start_max = 2010,   
  beta_sex       = -1.5,   
  beta_age       = 0.0,   
  beta_X         = 0.0,   
  borne_a        = 15       
)
# Number of simulations to run
N_sim <- 20

# Base seed for reproducibility
base_seed <- 12345

# Pre-allocate a list to store the generated datasets
sim_list <- vector("list", N_sim)

# ---------------------------------------------------------
# 2. DATA GENERATION LOOP
# ---------------------------------------------------------
for (i in 1:N_sim) {
  # Set a unique seed for each iteration
  set.seed(base_seed + i)
  
  # Generate the data for iteration i
  temp_data <- generate_data(
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
  
  # Add a simulation ID to differentiate patients 
  # from each dataset once they are aggregated
  temp_data$sim_id <- i
  
  # Time conversions (to days) useful for the estimators
  temp_data$observed_time_days <- temp_data$observed_time * 365.241
  temp_data$hypo_time_days     <- temp_data$hypothetical_time * 365.241
  temp_data$age_days           <- temp_data$age * 365.241
  
  # Store the current dataset in the list
  sim_list[[i]] <- temp_data
}

# ---------------------------------------------------------
# 3. DATA AGGREGATION
# ---------------------------------------------------------
# Combine all datasets into one large DataFrame
all_simulated_data <- do.call(rbind, sim_list)

# Quick check of the aggregated dataset
message("Total rows: ", nrow(all_simulated_data))
head(all_simulated_data)

# ---------------------------------------------------------
# 4. EXAMPLE OF GLOBAL OR PER-SIMULATION ANALYSIS
# ---------------------------------------------------------
# Global proportion of deaths vs censored
table(all_simulated_data$status)

# If you want to run the Kaplan-Meier or Pohar-Perme estimator 
# on the 1st simulation only:
data_sim1 <- subset(all_simulated_data, sim_id == 1)

# Estimators (executed on sim_id == 1)
KM_estimate <- survfit(Surv(hypo_time_days, hypothetical_status) ~ 1, data = data_sim1)

netSurv_estimate <- rs.surv(
  Surv(observed_time_days, status) ~ 1,
  data = data_sim1,
  ratetable = survexp.usr,
  rmap = list(
    age = age_days,
    sex = sex,
    race = race,
    year = year_diagnosis
  ),
  method = "pohar-perme"
)

# Plot for the first simulation
plot(
  netSurv_estimate,
  conf.int = TRUE,
  col = "blue",
  lwd = 2,
  xscale = 365.241,
  xlab = "Time since diagnosis (Years)",
  ylab = "Net Survival",
  main = "Estimated vs Theoretical Net Survival (Sim 1)",
  ylim = c(0.7, 1)
)

lines(
  KM_estimate,
  conf.int = TRUE,
  col = "red",
  lwd = 2,
  lty = 2,
  xscale = 365.241
)

legend(
  "bottomright",
  legend = c("Pohar-Perme Net Survival", "Theoretical Net Survival (KM)"),
  col = c("blue", "red"),
  lwd = c(2, 2),
  lty = c(1, 2),
  bty = "n"
)