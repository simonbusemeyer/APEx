
# censoring_calibration.R
source("generate_dataModified_ng.R")

calibrate_censoring_grid <- function(lambdas, n_patients, max_time_days, beta_age, n_pilots = 1, max_iter = 50) {
  # Dataframe to store calibration results and pilot metrics
  calibration_results <- data.frame(
    lambda = numeric(),
    calibrated_borne_a = numeric(),
    observed_censoring_rate = numeric(),
    pct_cancer_pilot = numeric()
  )
  
  cat("Starting censoring calibration...\n")
  
  for (lam in lambdas) {
    current_borne_a <- max_time_days
    calibrated <- FALSE
    iteration <- 1
    step_count <- 1 # Internal counter for the fail-safe
    
    while (!calibrated) {
      censoring_rates <- numeric(n_pilots)
      pct_cancer_rates <- numeric(n_pilots)
      
      # 1. Generate pilot datasets
      for (i in 1:n_pilots) {
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
        pct_cancer <- mean(df$event_type[df$status == 1] == "cancer")
      }
      
      # 2. Compute mean rates across the datasets
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
        
        # Save calibration parameters
        calibration_results <- rbind(
          calibration_results, 
          data.frame(
            lambda = lam, 
            calibrated_borne_a = current_borne_a, 
            observed_censoring_rate = mean_censoring,
            pct_cancer_pilot = mean_cancer
          )
        )
        
        #Task 4 table
        cat(sprintf("Success: Lambda = %.2f | borne_a = %7.1f | Censoring = %4.1f%% | Cancer = %4.1f%% | (Iterations: %d)\n", 
                    lam, current_borne_a, mean_censoring * 100, mean_cancer * 100, iteration))
      }
      
      # 4. FAIL-SAFE: Break if max iterations reached without calibrating
      if (!calibrated && step_count >= max_iter) {
        cat(sprintf("Failed to converge: Lambda = %.2f | Iteration limit reached. Returning borne_a = Inf\n", lam))
        
        calibration_results <- rbind(
          calibration_results, 
          data.frame(
            lambda = lam, 
            calibrated_borne_a = Inf, 
            observed_censoring_rate = mean_censoring, # Logs the closest it got
            pct_cancer_pilot = mean_cancer
          )
        )
        break # Exits the while loop, moves to next lambda
      }
      
      step_count <- step_count + 1
    }
  }
  
  return(calibration_results)
}