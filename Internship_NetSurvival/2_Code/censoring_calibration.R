
# --- 3. Censoring calibration ---

calibrate_censoring <- function(lambda_val, max_time_days, n_patients, beta_age) {
  current_borne_a <- max_time_days
  target_reached <- FALSE
  
  cat(sprintf("Calibrating borne_a for lambda = %.3f...\n", lambda_val))
  
  while (!target_reached) {
    censoring_rates <- numeric(100)
    
    # Generate 100 pilot datasets to test current_borne_a
    for (i in 1:100) {
      df_pilot <- generate_data(
        lambda = lambda_val, 
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
      
      # Calculate proportion of censored events
      censoring_rates[i] <- mean(df_pilot$event_type == "censored")
    }
    
    mean_cens_rate <- mean(censoring_rates)
    
    # Adjust borne_a based on the observed mean censoring rate
    if (mean_cens_rate > 0.32) {
      current_borne_a <- current_borne_a * 1.1  # Too much censoring -> increase borne_a
    } else if (mean_cens_rate < 0.28) {
      current_borne_a <- current_borne_a * 0.9  # Too little censoring -> decrease borne_a
    } else {
      target_reached <- TRUE
      cat(sprintf("Success! Calibrated borne_a: %.2f (Censoring Rate: %.1f%%)\n", 
                  current_borne_a, mean_cens_rate * 100))
    }
  }
  
  return(list(borne_a = current_borne_a, obs_censoring = mean_cens_rate))
}


# ---------------------------------------------------------
# Example implementation for Task 4 (Lambda Grid Integration)
# ---------------------------------------------------------
start <- proc.time()

lambda_grid <- c(0.3)
calibration_results <- list()
#
for (k in seq_along(lambda_grid)) {
res <- calibrate_censoring(
     lambda_val = lambda_grid[k],
     max_time_days = max_time_days,
     n_patients = n_patients,
     beta_age = beta_age
   )
   calibration_results[[k]] <- data.frame(
     lambda = lambda_grid[k],
     calibrated_borne_a = res$borne_a,
     observed_censoring = res$obs_censoring
   )
 }
elapsed <- proc.time() - start

 final_calibration_table <- do.call(rbind, calibration_results)
 print(final_calibration_table)