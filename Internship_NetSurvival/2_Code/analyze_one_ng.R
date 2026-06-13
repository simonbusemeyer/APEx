# =============================================================================
# analyze_one.R
# =============================================================================

analyze_one <- function(df, lambda, beta_age, times_years = c(1, 3, 4)) {
  
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
  
  net_surv_theo <- function(t){
    mean(exp(-lambda * t * exp(beta_age * df$ageCentre)))
  }
  net_surv_theo_3points <- sapply(times_years, net_surv_theo) # YEARS INSTEAD OF DAYS HERE!
  
  diff       <- net_surv_pp - net_surv_theo_3points
  covered    <- (net_surv_lower <= net_surv_theo_3points) & (net_surv_theo_3points <= net_surv_upper)
  pct_cancer <- mean(df$event_type[df$status == 1] == "cancer")
  cens_rate  <- mean(df$status == 0)
  
  res <- data.frame(
    time       = times_years,
    net_surv_pp       = net_surv_pp,
    net_surv_lower    = net_surv_lower,
    net_surv_upper    = net_surv_upper,
    net_surv_theo     = net_surv_theo_3points,
    diff       = diff,
    covered    = covered,
    pct_cancer = pct_cancer,
    cens_rate  = cens_rate
  )
  
  return(res)
}