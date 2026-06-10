# =============================================================================
# compute_metrics.R
# =============================================================================

#combine J dfs
compute_metrics <- function(results_list, lambda_val){
  all_res <- do.call(rbind, results_list)
  
#check for time points  
  times_eval <- sort(unique(all_res$time))
  n_times <- length(times_eval)
  
  #create the metrics df
  metrics <- data.frame(
    lambda = lambda_val,
    pct_cancer = mean(all_res$pct_cancer, na.rm = TRUE),
    time_t = times_eval,
    bias = numeric(n_times),
    rmse = numeric(n_times),
    ecr = numeric(n_times)
  )
  #metrics by time t
  for (i in 1:n_times) {
    res_t <- all_res[all_res$time == times_eval[i],]
    metrics$bias[i] <- mean(res_t$diff, na.rm = TRUE)
    metrics$rmse[i] <- sqrt(mean(res_t$diff^2, na.rm = TRUE))
    metrics$ecr[i] <- mean(res_t$covered, na.rm = TRUE)
  }
  return(metrics)
}