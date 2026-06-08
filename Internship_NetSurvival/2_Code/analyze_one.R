# =============================================================================
# analyze_one.R
# =============================================================================

library(survival)
library(relsurv)

analyze_one <- function(df, lambda, beta_age, times_years = c(1, 2, 3)) {
  
  # 1. Convert analysis times from years to days
  times_days <- times_years * 365.241
  
  # 2. Fit the Pohar-Perme estimator
  pp_fit <- rs.surv(
    Surv(observed_time, status) ~ 1, 
    data = df, 
    ratetable = survexp.us, 
    rmap = list(age = age * 365.241, sex = sex, year = year_diagnosis), 
    method = "pohar-perme"
  )
  
  pp_summary <- summary(pp_fit, times = times_days) #
  s_pp    <- pp_summary$surv
  s_lower <- pp_summary$lower
  s_upper <- pp_summary$upper
  
  s_theo <- sapply(times_days, function(t) {
    mean(exp(-lambda * t * exp(beta_age * df$ageStand)))
  })
  
  diff       <- s_pp - s_theo
  covered    <- (s_lower <= s_theo) & (s_theo <= s_upper)
  pct_cancer <- mean(df$event_type == "cancer")
  
  res <- data.frame(
    time       = times_years,
    s_pp       = s_pp,
    s_lower    = s_lower,
    s_upper    = s_upper,
    s_theo     = s_theo,
    diff       = diff,
    covered    = covered,
    pct_cancer = pct_cancer
  )
  
  return(res)
}