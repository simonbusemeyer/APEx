# =============================================================================
# analyze_one.R
# =============================================================================

analyze_one <- function(df, lambda, beta_age, max_time) {times_years <- c(1, ceiling((max_time + 1)/2), max_time)  
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
  
  #if conf.type="log" (default), std.err is on the log scale.
  # Multiplying by net_surv_pp converts it to the standard error of S(t).
  raw_se <- pp_summary$std.err
  net_surv_se <- raw_se * net_surv_pp
  
  net_surv_theo <- function(t){
    mean(exp(-lambda * t * exp(beta_age * df$ageCentre)))
  }
  net_surv_theo_3points <- sapply(times_years, net_surv_theo) # YEARS INSTEAD OF DAYS HERE!
  
  diff       <- net_surv_pp - net_surv_theo_3points
  covered    <- (net_surv_lower <= net_surv_theo_3points) & (net_surv_theo_3points <= net_surv_upper)
  pct_cancer <- mean(df$event_type[df$status == 1] == "cancer")
  n_deaths_cancer <- sum((df$event_type[df$status == 1] == "cancer"))
  n_deaths_other <- sum((df$event_type[df$status == 1] == "other"))
  cens_rate  <- mean(df$status == 0)
  
  res <- data.frame(
    time       = times_years,
    net_surv_pp       = net_surv_pp,
    se                = net_surv_se,
    net_surv_lower    = net_surv_lower,
    net_surv_upper    = net_surv_upper,
    net_surv_theo     = net_surv_theo_3points,
    diff              = diff,
    covered           = covered,
    n_deaths_cancer   = n_deaths_cancer,
    n_deaths_other    = n_deaths_other,
    pct_cancer        = pct_cancer,
    cens_rate         = cens_rate
  )
  
  return(res)
}