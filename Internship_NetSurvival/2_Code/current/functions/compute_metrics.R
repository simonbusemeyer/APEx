# =============================================================================
# compute_metrics.R
# =============================================================================

compute_metrics <- function(results_list, lambda_val, borne_a_val) {
  # 1. Faster binding than do.call(rbind, ...)
  dt <- data.table::rbindlist(results_list)
  
  metrics <- dt[, .(
    absolute_bias       = abs(mean(diff, na.rm = TRUE)),
    bias                = mean(diff, na.rm = TRUE),
    rmse                = sqrt(mean(diff^2, na.rm = TRUE)),
    estimate_sd         = sd(net_surv_pp, na.rm = TRUE),
    estimation_error_sd = sd(diff, na.rm = TRUE),
    mean_se             = mean(se, na.rm = TRUE),
    ecr                 = mean(covered, na.rm = TRUE),
    
    # --- RISK-SET DIAGNOSTICS ---
    # Overall At Risk
    mean_n_at_risk      = mean(n_at_risk, na.rm = TRUE),
    sd_n_at_risk        = sd(n_at_risk, na.rm = TRUE),
    med_n_at_risk       = median(n_at_risk, na.rm = TRUE),
    iqr_n_at_risk       = IQR(n_at_risk, na.rm = TRUE),
    
    # Age < 65 At Risk
    mean_n_at_risk_u65  = mean(n_at_risk_u65, na.rm = TRUE),
    sd_n_at_risk_u65    = sd(n_at_risk_u65, na.rm = TRUE),
    med_n_at_risk_u65   = median(n_at_risk_u65, na.rm = TRUE),
    iqr_n_at_risk_u65   = IQR(n_at_risk_u65, na.rm = TRUE),
    
    # Age 65-85 At Risk
    mean_n_at_risk_65_85= mean(n_at_risk_65_85, na.rm = TRUE),
    sd_n_at_risk_65_85  = sd(n_at_risk_65_85, na.rm = TRUE),
    med_n_at_risk_65_85 = median(n_at_risk_65_85, na.rm = TRUE),
    iqr_n_at_risk_65_85 = IQR(n_at_risk_65_85, na.rm = TRUE),
    
    # Age 85+ At Risk
    mean_n_at_risk_o85  = mean(n_at_risk_o85, na.rm = TRUE),
    sd_n_at_risk_o85    = sd(n_at_risk_o85, na.rm = TRUE),
    med_n_at_risk_o85   = median(n_at_risk_o85, na.rm = TRUE),
    iqr_n_at_risk_o85   = IQR(n_at_risk_o85, na.rm = TRUE),
    
    # Cumulative Cancer Events
    mean_cum_cancer     = mean(cum_cancer, na.rm = TRUE),
    sd_cum_cancer       = sd(cum_cancer, na.rm = TRUE),
    med_cum_cancer      = median(cum_cancer, na.rm = TRUE),
    iqr_cum_cancer      = IQR(cum_cancer, na.rm = TRUE),
    
    # Cumulative Censored
    mean_cum_censored   = mean(cum_censored, na.rm = TRUE),
    sd_cum_censored     = sd(cum_censored, na.rm = TRUE),
    med_cum_censored    = median(cum_censored, na.rm = TRUE),
    iqr_cum_censored    = IQR(cum_censored, na.rm = TRUE),
    
    # Cumulative Other Deaths
    mean_cum_other      = mean(cum_other, na.rm = TRUE),
    sd_cum_other        = sd(cum_other, na.rm = TRUE),
    med_cum_other       = median(cum_other, na.rm = TRUE),
    iqr_cum_other       = IQR(cum_other, na.rm = TRUE)
    
  ), by = .(
    pp_age_class,
    time_t = time)]
  
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