# =============================================================================
# generate_data.R
# Fonctions de génération de données de survie
# =============================================================================

generate_data <- function(lambda,
                          age_option,
                          n,
                          max_time,
                          prop_female,
                          year.start_min,
                          year.start_max,
                          beta_sex,
                          beta_age,
                          borne_a) {
  
  # Covariables generation: Age
  if (age_option == "A") {
    # Option A: age ~ Uniform[80, 90)
    age <- runif(n, min = 80, max = 90)
    
  } else if (age_option == "C") {
    # Option C: age ~ Uniform[15, 40)
    age <- runif(n, min = 15, max = 40)
    
  } else if (age_option == "D") {
    # Option D: age ~ Uniform[50, 75)
    age <- trunc(runif(n, min = 50, max = 75))
    
  } else if (age_option == "E") {
    # Option E: age ~ Uniform[50, 60)
    age <- runif(n, min = 50, max = 60)
    
  } else if (age_option == "F") {
    # Option F: age ~ Discrete-Uniform[80, 90)
    age <- trunc(runif(n, min = 80, max = 90))
    
  } else if (age_option == "Luo") {
    # Option Luo: Empirical distribution from Luo et al. 2023 (Prostate Cancer)
    # <65 (46.2%), 65-85 (51.9%), >85 (2.0%)
    u <- runif(n)
    age <- numeric(n)
    
    idx1 <- u <= 0.462
    idx2 <- u > 0.462 & u <= 0.981
    idx3 <- u > 0.981
    
    # clinical bounds: [50-65), [65-85), [85-95)
    age[idx1] <- runif(sum(idx1), min = 50, max = 65)
    age[idx2] <- runif(sum(idx2), min = 65, max = 85)
    age[idx3] <- runif(sum(idx3), min = 85, max = 95)
    
  } else if (age_option == "LuoTrunc") {
    # Option LuoTrunc: Truncated Empirical distribution from Luo et al. 2023 (Prostate Cancer)
    # trunc(<65 (46.2%), 65-85 (51.9%), >85 (2.0%))
    u <- runif(n)
    age <- numeric(n)
    
    idx1 <- u <= 0.462
    idx2 <- u > 0.462 & u <= 0.981
    idx3 <- u > 0.981
    
    # clinical bounds: [50-65), [65-85), [85-95)
    age[idx1] <- trunc(runif(sum(idx1), min = 50, max = 65))
    age[idx2] <- trunc(runif(sum(idx2), min = 65, max = 85))
    age[idx3] <- trunc(runif(sum(idx3), min = 85, max = 95))
    
  } else {
    # Catch any invalid inputs
    stop("age_option must be 'A', 'C', 'D', 'E', 'F', 'Luo'or 'LuoTrunc'")
  }
  
  # Dynamic Calculation of Theoretical Means
  # Define standard bounds for the uniform distributions: c(min, max)
  bounds <- switch(age_option,
                   "A"        = c(80, 90),
                   "C"        = c(15, 40),
                   "D"        = c(50, 75),
                   "E"        = c(50, 60),
                   "F"        = c(80, 90),
                   "Luo"      = NULL, # Handled separately below
                   "LuoTrunc" = NULL)
  
  if (!is.null(bounds)) {
    if (age_option %in% c("D", "F")) {
      # For discrete uniform distributions via trunc(runif(min, max)):
      # The values are integers from min to (max - 1).
      # The mean of integers from a to b is (a + b) / 2.
      theoretical_mean <- (bounds[1] + (bounds[2] - 1)) / 2
    } else {
      # For continuous uniform distributions: (min + max) / 2
      theoretical_mean <- mean(bounds)
    }
  } else {
    # For Luo and LuoTrunc empirical mixture distributions:
    # Expected value E[X] = sum( P(Interval) * E[X | Interval] )
    weights <- c(0.462, 0.519, 0.020) # the paper has a rounding error of +0.001
    
    if (age_option == "Luo") {
      # Continuous sub-intervals
      means <- c(mean(c(50, 65)), mean(c(65, 85)), mean(c(85, 95)))
    } else if (age_option == "LuoTrunc") {
      # Truncated (discrete) sub-intervals: trunc(runif(min, max))
      # Means are calculated using the true integer ranges: [50,64], [65,84], [85,94]
      means <- c((50 + 64)/2, (65 + 84)/2, (85 + 94)/2)
    }
    
    theoretical_mean <- sum(weights * means)
  }
  
  ageCentre <- age - theoretical_mean
  
  uSex <- runif(n)
  pFem <- prop_female
  sex  <- rbinom(n, size = 1, prob = prop_female)   # hommes=0, femmes=1
  
  if (year.start_min != year.start_max) {
    year.start <- as.Date(paste0(
      sample(year.start_min:year.start_max, n, replace = TRUE),
      "-01-01"
    )) +
      sample(0:364, n, replace = TRUE)
  } else {
    year.start <- as.Date(paste0(year.start_min, "-01-01")) +
      sample(0:364, n, replace = TRUE)
  }
  
  # T_E generation
  tempuS <- runif(n)
  exp.betaz <- exp(beta_sex * sex + beta_age * ageCentre)
  
  ui <- runif(n)
  tpsSpe <- -log(ui) / (lambda * exp.betaz)
  
  # T_P generation
  tpsGene <- rep(0, n)
  TauxAtt <- rep(NA, n)
  sexNom <- ifelse(sex == 0, "male", "female")
  
  f1 <- function(i) {
    U <- runif(1)
    H_target <- -log(U) 
    H_accum <- 0        
    tpsGene_val <- 0    
    
    # Calculate how much time remains in the starting age block
    current_age_floor <- trunc(age[i])
    fraction_first_year <- 1 - (age[i] - current_age_floor) 
    
    i.age <- which(attr(survexp.us, which = "dimnames")[[1]] == current_age_floor)
    i.sex <- which(attr(survexp.us, which = "dimnames")[[2]] == sexNom[i])
    i.year <- which(attr(survexp.us, which = "dimnames")[[3]] == format(year.start[i], "%Y"))
    
    if (length(i.age) == 0 || length(i.sex) == 0 || length(i.year) == 0) {
      return(NA)  
    }
    
    max.i.age <- length(attributes(survexp.us)$dimnames[[1]])
    max.i.year <- length(attributes(survexp.us)$dimnames[[3]])
    
    is_first_year <- TRUE # Flag to track if we are in the starting fractional year
    
    while (TRUE) {
      h_daily <- survexp.us[i.age, i.sex, i.year]
      
      # Year block size depends on whether it's the first partial year or a full year
      year_fraction <- ifelse(is_first_year, fraction_first_year, 1)
      H_step <- h_daily * 365.241 * year_fraction 
      
      if (H_accum + H_step >= H_target) {
        # Calculate the exact fraction of the year lived before the event
        H_remaining <- H_target - H_accum
        
        # Time lived in this interval is proportional to remaining hazard vs daily rate
        tpsGene_val <- tpsGene_val + (H_remaining / (h_daily * 365.241))
        break
      } else {
        # Patient survives this interval; accumulate hazard and time
        H_accum <- H_accum + H_step
        tpsGene_val <- tpsGene_val + year_fraction
        
        # Advance ratetable indices (capped at max available limits)
        i.age <- min(i.age + 1, max.i.age)
        i.year <- min(i.year + 1, max.i.year)
        is_first_year <- FALSE # All subsequent loops are full 1-year blocks
      }
    }
    return(tpsGene_val)
  }
  
  tpsGene <- sapply(1:n, f1)
  
  # tpsGene <- tpsGene * 365.241
  # tpsSpe  <- tpsSpe * 365.241
  tpsSurv <- pmin(tpsGene, tpsSpe)
  
  # CENSORING
  borne_a = borne_a
  
  # censoring time
  # tpsCens <- runif(n, min = 0, max = borne_a)
  if (is.finite(borne_a)) {
    tpsCens <- runif(n, min = 0, max = borne_a)
  } else {
    tpsCens <- rep(Inf, n)
  }
  
  temps <- pmin(tpsCens, tpsSurv)
  temps2 <- pmin(tpsCens, tpsSpe)
  
  statut     <- ifelse(temps == tpsCens, 0, 1)
  
  # cause of death
  cause <- ifelse(tpsSurv == tpsSpe & temps != tpsCens, 1, 0)
  
  # hypothetical world
  cause2     <- ifelse(temps2 == tpsCens, 0, 1)
  
  # administrative censoring
  
  statut[temps > max_time]   <-  0
  cause[temps > max_time] <- 0
  # hypothetical world
  cause2[temps2 > max_time]  <-  0
  
  temps[temps > max_time]    <- max_time
  # hypothetical world
  temps2[temps2 > max_time]  <- max_time
  
  # Added event_type 
  event_type <- ifelse(temps == tpsCens | temps == max_time, "censored",
                       ifelse(temps == tpsSpe, "cancer", "other"))
  
  result <- data.frame(
    patient_id = 1:n,
    age = age,
    ageCentre = ageCentre,
    sex_num = sex,
    sex = factor(sexNom, levels = c("male", "female")),
    tpsCens = tpsCens,
    tpsGene = tpsGene,
    tpsSpe = tpsSpe,
    year_diagnosis = year.start,
    observed_time = temps,
    status = statut,
    cause = cause,
    event_type = event_type,
    hypothetical_time = temps2,
    hypothetical_status = cause2,
    admin_cens = max_time
  )
  
  return(result)
}