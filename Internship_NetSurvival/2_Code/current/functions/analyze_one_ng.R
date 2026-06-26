# =============================================================================
# analyze_one.R
# =============================================================================

analyze_one <- function(df, lambda, beta_age, max_time) {times_years <- sort(unique(c(
  1, 
  4, 
  ceiling((max_time + 1)/2), 
  # pmax guarantees sequence elements are bounded above 1
  #pmax(1, seq(max_time - 6, max_time, by = 0.5))
  max_time
)))
  # 1. Convert analysis times from years to days
  times_days <- times_years * 365.241
  
  # 2. Fit the Pohar-Perme estimator
  pp_fit <- rs.surv(
    Surv(observed_time*365.241, status) ~ 1, 
    data = df, 
    ratetable = survexp.us, 
    rmap = list(age = age * 365.241, sex = sex, year = year_diagnosis), 
    method = "pohar-perme", 
    add.times = times_days
  )
  
  #plot(pp_fit)
  
  pp_summary <- summary(pp_fit, times = times_days, extend = TRUE) 
  net_surv_pp    <- pp_summary$surv
  net_surv_lower <- pp_summary$lower
  net_surv_upper <- pp_summary$upper
  
  net_surv_se <- pp_summary$std.err
  
  net_surv_theo <- function(t){
    mean(exp(-lambda * t * exp(beta_age * df$ageCentre)))
  }
  net_surv_theo_3points <- sapply(times_years, net_surv_theo) # YEARS INSTEAD OF DAYS HERE!
  
  diff       <- net_surv_pp - net_surv_theo_3points
  covered    <- (net_surv_lower <= net_surv_theo_3points) & (net_surv_theo_3points <= net_surv_upper)
  
  # --- RISK-SET DIAGNOSTICS CALCULATIONS ---
  # Define age classes matching the Luo/clinical partitions
  df$age_class <- cut(df$age, breaks = c(0, 65, 85, Inf), right = FALSE, labels = c("<65", "65-85", "85+"))
  
  n_at_risk_overall <- numeric(length(times_years))
  n_at_risk_u65     <- numeric(length(times_years))
  n_at_risk_65_85   <- numeric(length(times_years))
  n_at_risk_o85     <- numeric(length(times_years))
  cum_cancer        <- numeric(length(times_years))
  cum_censored      <- numeric(length(times_years))
  cum_other         <- numeric(length(times_years))
  
  for (i in seq_along(times_years)) {
    t <- times_years[i]
    
    # Still at risk (survival time >= horizon)
    n_at_risk_overall[i] <- sum(df$observed_time >= t)
    n_at_risk_u65[i]     <- sum(df$observed_time >= t & df$age_class == "<65")
    n_at_risk_65_85[i]   <- sum(df$observed_time >= t & df$age_class == "65-85")
    n_at_risk_o85[i]     <- sum(df$observed_time >= t & df$age_class == "85+")
    
    # Cumulative events BEFORE the horizon
    cum_cancer[i]   <- sum(df$observed_time < t & df$event_type == "cancer")
    cum_censored[i] <- sum(df$observed_time < t & df$event_type == "censored")
    cum_other[i]    <- sum(df$observed_time < t & df$event_type == "other")
  }
  
  # Global metrics
  pct_cancer      <- mean(df$event_type[df$status == 1] == "cancer")
  n_deaths_cancer <- sum(df$event_type[df$status == 1] == "cancer")
  n_deaths_other  <- sum(df$event_type[df$status == 1] == "other")
  cens_rate       <- mean(df$status == 0)
  
  res <- data.frame(
    time              = times_years,
    net_surv_pp       = net_surv_pp,
    se                = net_surv_se,
    net_surv_lower    = net_surv_lower,
    net_surv_upper    = net_surv_upper,
    net_surv_theo     = net_surv_theo_3points,
    diff              = diff,
    covered           = covered,
    n_at_risk         = n_at_risk_overall,
    n_at_risk_u65     = n_at_risk_u65,
    n_at_risk_65_85   = n_at_risk_65_85,
    n_at_risk_o85     = n_at_risk_o85,
    cum_cancer        = cum_cancer,
    cum_censored      = cum_censored,
    cum_other         = cum_other,
    n_deaths_cancer   = n_deaths_cancer,
    n_deaths_other    = n_deaths_other,
    pct_cancer        = pct_cancer,
    cens_rate         = cens_rate
  )
  
  return(res)
}