# =============================================================================
# compute_metrics.R
# =============================================================================

compute_metrics <- function(results_list, lambda_val, borne_a_val) {
  # 1. Faster binding than do.call(rbind, ...)
  dt <- data.table::rbindlist(results_list)
  
  # 2. Native C-level grouped aggregation by time
  metrics <- dt[, .(
    absolute_bias = abs(mean(diff, na.rm = TRUE)),
    bias          = mean(diff, na.rm = TRUE),
    rmse          = sqrt(mean(diff^2, na.rm = TRUE)),
    estimate_sd       = sd(net_surv_pp, na.rm = TRUE),
    estimation_error_sd = sd(diff, na.rm = TRUE),
    mean_se       = mean(se, na.rm = TRUE),
    ecr           = mean(covered, na.rm = TRUE)
  ), by = .(time_t = time)]
  
  metrics[, se_calibration_ratio :=
            mean_se / estimation_error_sd]
  
  # 3. Assign global scenario averages in place
  metrics[, `:=`(
    lambda          = lambda_val,
    borne_a         = borne_a_val,
    censoring_rate  = mean(dt$cens_rate, na.rm = TRUE),
    pct_cancer      = mean(dt$pct_cancer, na.rm = TRUE),
    n_deaths_cancer = mean(dt$n_deaths_cancer, na.rm = TRUE),
    n_deaths_other  = mean(dt$n_deaths_other, na.rm = TRUE)
  )]
  
  # Reorder
  data.table::setcolorder(metrics, c("lambda", "borne_a", "censoring_rate", "pct_cancer", 
                                     "n_deaths_cancer", "n_deaths_other", "time_t"))
  
  return(as.data.frame(metrics))
}