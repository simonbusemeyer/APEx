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
  # Convert analysis times from years to days
  times_days <- times_years * 365.241
  
 # define age classes
 # PP estimation
  df$pp_age_class <- cut(
    df$age,
    breaks = c(0, 65, Inf),
    right = FALSE,
    labels = c("<65", ">=65")
  )
  
  # diagnostics
  df$diag_age_class <- cut(
    df$age,
    breaks = c(0, 65, 85, Inf),
    right = FALSE,
    labels = c("<65", "65-85", "85+")
  )
  
  df$pp_age_class   <- droplevels(df$pp_age_class)
  df$diag_age_class <- droplevels(df$diag_age_class)
  
  
  # Fit the Pohar-Perme estimator stratified by age class
  pp_fit <- rs.surv(
    Surv(observed_time*365.241, status) ~ pp_age_class, #age class ^ <65; >=65
    data = df, 
    ratetable = survexp.us, 
    rmap = list(age = age * 365.241, sex = sex, year = year_diagnosis), 
    method = "pohar-perme", 
    add.times = times_days
  )
  
  #plot(pp_fit)
  
  pp_summary <- summary(pp_fit, times = times_days, extend = TRUE) 
  
  #Put PP summary into a data frame
  pp_df <- data.frame(
    time_days = pp_summary$time,
    time = pp_summary$time / 365.241,
    pp_age_class = 
      as.character(pp_summary$strata),
    net_surv_pp    = pp_summary$surv,
    se             = pp_summary$std.err,
    net_surv_lower = pp_summary$lower,
    net_surv_upper = pp_summary$upper
  )
  
  pp_df$pp_age_class <- sub("^pp_age_class=", "", pp_df$pp_age_class)
  
  pp_df$time <- times_years[
    match(
      round(pp_df$time, 8),
      round(times_years, 8)
    )
  ]
  
  #age specific theoretical net survival 
  theo_list <- lapply(split(df, df$pp_age_class), function(dfg){
    
    data.frame(
      pp_age_class = unique(dfg$pp_age_class),
      time = times_years,
      net_surv_theo = sapply(times_years, function(t){
        mean(exp(-lambda * t * exp(beta_age * dfg$ageCentre)))
      })
    )
  })
  
  theo_df <- do.call(rbind, theo_list)
  
 # Risk-set and event diagnostics by <65, 65-85, 85+
  diag_grid <- expand.grid(
    diag_age_class = levels(df$diag_age_class),
    time = times_years,
    stringsAsFactors = FALSE
  )
  
  diag_df <- do.call(rbind, lapply(seq_len(nrow(diag_grid)), function(i) {
    
    age_group <- diag_grid$diag_age_class[i]
    t <- diag_grid$time[i]
    
    dfg <- df[df$diag_age_class == age_group, ]
    
    n_deaths_cancer <- sum(dfg$status == 1 & dfg$event_type == "cancer")
    n_deaths_other  <- sum(dfg$status == 1 & dfg$event_type == "other")
    n_deaths_total  <- n_deaths_cancer + n_deaths_other
    
    pct_cancer <- ifelse(
      n_deaths_total > 0,
      n_deaths_cancer / n_deaths_total,
      NA_real_
    )
    
    data.frame(
      diag_age_class  = age_group,
      time            = t,
      n_patients      = nrow(dfg),
      n_at_risk       = sum(dfg$observed_time >= t),
      cum_cancer      = sum(dfg$observed_time < t & dfg$event_type == "cancer"),
      cum_censored    = sum(dfg$observed_time < t & dfg$event_type == "censored"),
      cum_other       = sum(dfg$observed_time < t & dfg$event_type == "other"),
      n_deaths_cancer = n_deaths_cancer,
      n_deaths_other  = n_deaths_other,
      pct_cancer      = pct_cancer,
      cens_rate       = mean(dfg$status == 0)
    )
  }))
  
  # Convert diagnostics wide
  diag_wide <- reshape(
    diag_df,
    idvar = "time",
    timevar = "diag_age_class",
    direction = "wide"
  )
  
  names(diag_wide) <- gsub("\\.", "_", names(diag_wide))
  names(diag_wide) <- gsub("<65", "u65", names(diag_wide), fixed = TRUE)
  names(diag_wide) <- gsub("65-85", "65_85", names(diag_wide), fixed = TRUE)
  names(diag_wide) <- gsub("85+", "o85", names(diag_wide), fixed = TRUE)
  
  # Reconstruct the global/unstratified totals
  diag_wide$n_at_risk <- diag_wide$n_at_risk_u65 + diag_wide$n_at_risk_65_85 + diag_wide$n_at_risk_o85
  diag_wide$cum_cancer <- diag_wide$cum_cancer_u65 + diag_wide$cum_cancer_65_85 + diag_wide$cum_cancer_o85
  diag_wide$cum_censored <- diag_wide$cum_censored_u65 + diag_wide$cum_censored_65_85 + diag_wide$cum_censored_o85
  diag_wide$cum_other    <- diag_wide$cum_other_u65 + diag_wide$cum_other_65_85 + diag_wide$cum_other_o85
  
  # merge estimates and diagnostics
  res <- merge(
    pp_df,
    theo_df,
    by = c("pp_age_class", "time"),
    all.x = TRUE
  )
  
  res <- merge(
    res,
    diag_wide,
    by = "time",
    all.x = TRUE
  )
  
  # compute diff and coverage
  res$diff <- res$net_surv_pp - res$net_surv_theo
  
  res$covered <- with(
    res,
    net_surv_lower <= net_surv_theo & net_surv_theo <= net_surv_upper
  )
  
  res <- res[order(res$pp_age_class, res$time), ]
  
  rownames(res) <- NULL
  
  
  return(res)
}