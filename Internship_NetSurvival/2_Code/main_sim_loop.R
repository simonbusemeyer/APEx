#Simulation Metrics Loop (run last)

library(survival)
library(relsurv)

source("generate_dataModified_ng.R")
source("analyze_one_ng.R")
source("compute_metrics.R")

# --- 0. Parameters --- ####
n_patients <- 2000
age_option <- "A"
lambda <- 0.01 
beta_age <- 0.02
beta_sex <- 0
max_time <- 4
max_time_days <- max_time * 365.241
year.start_min <- 2008
year.start_max <- 2010
prop_female <- 0
N_files <- 20

#lambda/borne_a pairs
scenarios <- data.frame(
  lambda = c(0.05, 0.10, 0.20, 0.30, 0.50, 1.00, 1.50),
  borne_a = c(Inf, Inf, 30, 12, 5, 3, 2)
)



all_metrics <- list()

#reconciles the output directory
if(!dir.exists("outputs/tables")) dir.create("outputs/tables", recursive = TRUE)

start <- proc.time()

#cycles through each scenario
for (k in 1:nrow(scenarios)) {
  lambda_scenario <- scenarios$lambda[k]
  borne_a_scenario <- scenarios$borne_a[k]
  cat(sprintf("scenario %d/%d: lambda = %.2f, borne_a = %.0f\n", k, nrow(scenarios), lambda_scenario, borne_a_scenario))
  
  results_scenarios <- vector("list", N_files)
  
  # cycles an inner loop through N_files datasets
  for (j in 1:N_files) {
    df <- generate_data(
      lambda = lambda_scenario,     
      age_option = age_option,   
      n = n_patients,        
      max_time = max_time, 
      prop_female = prop_female,     
      year.start_min = year.start_min, 
      year.start_max = year.start_max,
      beta_sex = beta_sex,         
      beta_age = beta_age, 
      borne_a = borne_a_scenario         
    )
    
    results_scenarios[[j]] <- analyze_one(df, lambda = lambda_scenario, beta_age = beta_age, times_years = c(1, 2, 3))
  }
  
  # calculates aggregate Bias, RMSE, and ECR of all datasets
  all_metrics[[k]] <- compute_metrics(results_list = results_scenarios, lambda_val = lambda_scenario)
}

elapsed <- proc.time() - start
elapsed

# stacks each scenario in one df
final_results <- do.call(rbind, all_metrics)
