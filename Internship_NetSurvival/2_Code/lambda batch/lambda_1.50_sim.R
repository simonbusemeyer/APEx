# lambda_1.50_sim.R

library(future.apply)
library(data.table)

# --- Scenario Parameters ---
lambda_scenario <- 1.50
borne_a_scenario <- Inf

cat(sprintf("Running scenario: lambda = %.2f, borne_a = %.0f\n", lambda_scenario, borne_a_scenario))
results_scenarios <- vector("list", N_files)
df <- vector("list", N_files)

# Setup parallel backend (leaves 1 core free for the OS)
plan(multisession, workers = availableCores() - 1)

library(future)

# 1. Verify how many cores the OS has made available to R
cat("Total Cores Available:", availableCores(), "\n")

# 2. Verify how many workers the future plan has actually registered
cat("Active Parallel Workers:", nbrOfWorkers(), "\n")

# Pre-allocate list for datasets
df <- vector("list", N_files)

set.seed(12345) #use the same N_files dfs for each lambda

for (j in 1:N_files) {
  df[[j]] <- generate_data(
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
  
  # --- CHANGED: Assign simulation ID for downstream visualization subsetting ---
  df[[j]]$sim_id <- j
}

# --- 2. Analyze Data in Parallel ---
# future_lapply automatically distributes the data and execution across CPU cores
results_scenarios <- future_lapply(1:N_files, function(j) {
  analyze_one(df[[j]], lambda = lambda_scenario, beta_age = beta_age, times_years = c(1, 2, 3))
}, future.seed = TRUE)

# --- 3. Aggregate Data Efficiently ---
if(!dir.exists("outputs/data")) dir.create("outputs/data", recursive = TRUE)

# CHANGED: Use data.table::rbindlist instead of do.call(rbind, ...) for exponential speedup
all_scenario_data <- rbindlist(df)
saveRDS(all_scenario_data, file = sprintf("outputs/data/simulated_cohort_lambda_%.2f.rds", lambda_scenario))

# --- 4. Calculate and Save Metrics ---
metrics <- compute_metrics(results_list = results_scenarios, lambda_val = lambda_scenario, borne_a_val = borne_a_scenario)
saveRDS(metrics, file = sprintf("outputs/tables/metrics_lambda_%.2f.rds", lambda_scenario))

# Close parallel backend to free up resources
plan(sequential)
