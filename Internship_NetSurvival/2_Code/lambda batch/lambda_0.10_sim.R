# lambda_0.10_sim.R

# --- Scenario Parameters ---
lambda_scenario <- 0.10
borne_a_scenario <- Inf

cat(sprintf("Running scenario: lambda = %.2f, borne_a = %.0f\n", lambda_scenario, borne_a_scenario))

results_scenarios <- vector("list", N_files)
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
  
  results_scenarios[[j]] <- analyze_one(df[[j]], lambda = lambda_scenario, beta_age = beta_age, times_years = c(1, 2, 3))
}

# --- CHANGED: Save the aggregated raw cohort data to disk ---
if(!dir.exists("outputs/data")) dir.create("outputs/data", recursive = TRUE)
all_scenario_data <- do.call(rbind, df)
saveRDS(all_scenario_data, file = sprintf("outputs/data/simulated_cohort_lambda_%.2f.rds", lambda_scenario))
# -----------------------------------------------------------

# Calculate metrics and save independently
metrics <- compute_metrics(results_list = results_scenarios, lambda_val = lambda_scenario, borne_a_val = borne_a_scenario)
saveRDS(metrics, file = sprintf("outputs/tables/metrics_lambda_%.2f.rds", lambda_scenario))