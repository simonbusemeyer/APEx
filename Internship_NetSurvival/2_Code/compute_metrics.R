# =============================================================================
# compute_metrics.R
# =============================================================================

compute_metrics <- function(results_list, lambda_val, borne_a_val){
  # Combine J dataframes
  all_res <- do.call(rbind, results_list)
  
  # Check for time points  
  times_eval <- sort(unique(all_res$time))
  n_times <- length(times_eval)
  
  # Create the metrics df
  metrics <- data.frame(
    lambda = lambda_val,
    borne_a = borne_a_val,
    censoring_rate = mean(all_res$cens_rate, na.rm = TRUE),
    pct_cancer = mean(all_res$pct_cancer, na.rm = TRUE),
    time_t = times_eval,
    
    # New columns for absolute counts
    n_total = numeric(n_times),
    n_success = numeric(n_times),
    
    # New columns for reliability
    convergence_rate = numeric(n_times),
    failure_rate = numeric(n_times),
    
    # New columns for bifurcated metrics
    bias_conditional = numeric(n_times),
    rmse_conditional = numeric(n_times),
    ecr_conditional = numeric(n_times),
    ecr_unconditional = numeric(n_times)
  )
  
  # Calculate metrics by time t
  for (i in 1:n_times) {
    res_t <- all_res[all_res$time == times_eval[i],]
    
    metrics$n_total[i] <- nrow(res_t)
    metrics$n_success[i] <- sum(!is.na(res_t$diff))
    
    metrics$convergence_rate[i] <- mean(!is.na(res_t$diff))
    metrics$failure_rate[i] <- mean(is.na(res_t$diff))
    
    #Calculate metrics conditional on successful convergence
    metrics$bias_conditional[i] <- mean(res_t$diff, na.rm = TRUE)
    metrics$rmse_conditional[i] <- sqrt(mean(res_t$diff^2, na.rm = TRUE))
    metrics$ecr_conditional[i] <- mean(res_t$covered, na.rm = TRUE)
    
    #Calculate Unconditional Coverage where NA is a failure to cover
    metrics$ecr_unconditional[i] <- sum(res_t$covered %in% TRUE) / nrow(res_t)
  }
  
  return(metrics)
}