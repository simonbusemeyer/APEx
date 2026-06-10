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
    # Option A: age ~ Uniform[80, 89]
    age <- runif(n, min = 80, max = 89)
    
  } else if (age_option == "C") {
    # Option C: age ~ Uniform[15, 39]
    age <- runif(n, min = 15, max = 39)
    
  } else if (age_option == "D") {
    # Option D: age ~ Uniform[50, 74]
    age <- runif(n, min = 50, max = 74)
    
  } else {
    # Catch any invalid inputs
    stop("age_option must be 'A', 'C', or 'D'")
  }
  
  ageMoyen <- mean(age)
  # ageStand <- (age - ageMoyen) / sd(age)
  ageCentre <- (age - ageMoyen)
  
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
    uAtt <- runif(1)
    
    i.age <- which(attr(survexp.us, which = "dimnames")[[1]] == trunc(age[i]))
    i.sex <- which(attr(survexp.us, which = "dimnames")[[2]] == sexNom[i])
    i.year <- which(attr(survexp.us, which = "dimnames")[[3]] == format(year.start[i], "%Y"))
    
    if (length(i.age) == 0 ||
        length(i.sex) == 0 || length(i.year) == 0) {
      return(NA)  # Retourner NA si indices invalides
    }
    
    max.i.age <- length(attributes(survexp.us)$dimnames[[1]])
    max.i.year <- length(attributes(survexp.us)$dimnames[[3]])
    TauxAtt[i] <- 1 - exp(-365.241 * survexp.us[i.age, i.sex, i.year])
    
    if (uAtt <= TauxAtt[i]) {
      tpsG <- runif(1)
      tpsGene[i] <- tpsG
    } else{
      while (uAtt > TauxAtt[i]) {
        tpsGene[i] <- tpsGene[i] + 1
        i.age <- min(i.age + 1, max.i.age) # happy birthday !
        i.year <- min(i.year + 1, max.i.year) # and one calendar year more...
        
        uAtt <- runif(1)
        TauxAtt[i] <- 1 - exp(-365.241 * survexp.us[i.age, i.sex, i.year])
      }
      tpsG <- runif(1)
      tpsGene[i] <- tpsGene[i] + tpsG
    }
    return(tpsGene[i])
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
    age = round(age, 1),
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