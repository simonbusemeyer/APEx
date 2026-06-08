library(survival)
source("generate_dataModified.R")
source("analyze_one.R")

# --- 0. Parameters ---
n_patients <- 2000 #
lambda <- 0.05 
beta_age <- 0.2
max_time <- 4
max_time_days <- max_time * 365.241
borne_a <- max_time_days

# --- 2. Timing check (run once on 1 dataset) ---
start <- proc.time()

df <- generate_data(
  lambda, age_option = "A", n = n_patients, max_time = max_time_days,
  prop_female = 0,
  year.start_min = 2008, year.start_max = 2010,
  beta_sex = 0, beta_age = beta_age, borne_a = borne_a
)

res <- analyze_one(df, lambda, beta_age, times_years = c(1, 2, 3))

elapsed <- proc.time() - start
elapsed

start <- proc.time()

# --- 3. Censoring calibration ---

lambdas <- c(0.5) 

# Dataframe to store calibration results and pilot metrics
calibration_results <- data.frame(
  lambda = numeric(),
  calibrated_borne_a = numeric(),
  observed_censoring_rate = numeric(),
  pct_cancer_pilot = numeric()
)

set.seed(12345)

cat("Starting censoring calibration...\n")

for (lam in lambdas) {
  current_borne_a <- max_time_days
  calibrated <- FALSE
  iteration <- 1
  
  while (!calibrated) {
    censoring_rates <- numeric(1)
    pct_cancer_rates <- numeric(1)
    
    # 1. Generate 10 pilot datasets
    for (i in 1:1) {
      df <- generate_data(
        lambda = lam, 
        age_option = "A", 
        n = n_patients, 
        max_time = max_time_days,
        prop_female = 0,
        year.start_min = 2008, 
        year.start_max = 2010,
        beta_sex = 0, 
        beta_age = beta_age, 
        borne_a = current_borne_a
      )
      
      # Extract metrics per dataset
      censoring_rates[i] <- mean(df$event_type == "censored")
      pct_cancer_rates[i]  <- mean(df$event_type == "cancer")
    }
    
    # 2. Compute mean rates across the 10 datasets
    mean_censoring <- mean(censoring_rates)
    mean_cancer    <- mean(pct_cancer_rates)
    
    # 3. Check boundaries and adjust borne_a
    if (mean_censoring > 0.32) {
      current_borne_a <- current_borne_a * 1.1 
      iteration <- iteration + 1
    } else if (mean_censoring < 0.28) {
      current_borne_a <- current_borne_a * 0.9
      iteration <- iteration + 1
    } else {
      calibrated <- TRUE
      
      # Save successful calibration parameters
      calibration_results <- rbind(
        calibration_results, 
        data.frame(
          lambda = lam, 
          calibrated_borne_a = current_borne_a, 
          observed_censoring_rate = mean_censoring,
          pct_cancer_pilot = mean_cancer
        )
      )
      
      # Print output for the required Task 4 table
      cat(sprintf("Success: Lambda = %.2f | borne_a = %7.1f | Censoring = %4.1f%% | Cancer = %4.1f%% | (Iterations: %d)\n", 
                  lam, current_borne_a, mean_censoring * 100, mean_cancer * 100, iteration))
    }
  }
}

elapsed <- proc.time() - start
elapsed

# View the final calibrated table
print(calibration_results)