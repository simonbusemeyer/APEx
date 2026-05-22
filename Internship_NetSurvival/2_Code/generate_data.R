# =============================================================================
# generate_data.R
# Fonctions de génération de données de survie
# =============================================================================

generate_data <- function(k, pho, alpha, n, max_time, prop_female, prop_x0, year.start_min, year.start_max, beta_sex, beta_age, beta_X, borne_a) {
  # Covariables generation
  age_1 <- runif(n * .25, min = 30, max = 65)
  age_2 <- runif(n * .35, min = 65, max = 75)
  age_3 <- runif(n * .40, min = 75, max = 80)
  age <- c(age_1, age_2, age_3)
  ageMoyen <- mean(age)
  ageStand <- (age - ageMoyen) / sd(age)
  
  uSex <- runif(n)
  pFem <- prop_female
  sex  <- rbinom(n, size = 1, prob = prop_female)   # hommes=0, femmes=1
  
  uX   <- runif(n)
  X    <- rbinom(n, size = 1, prob = 1 - prop_x0)   # X = 0 ou 1
  
  
  if (year.start_min != year.start_max) {
    year.start <- as.Date(paste0(sample(year.start_min:year.start_max, n, replace = TRUE), "-01-01")) +
      sample(0:364, n, replace = TRUE)
  } else {
    year.start <- as.Date(paste0(year.start_min, "-01-01")) +
      sample(0:364, n, replace = TRUE)
  }
  
  # T_E generation
  tempuS <- runif(n)
  exp.betaz <- exp(
    beta_sex * sex + beta_age * ageStand + beta_X * X
  )
  
  ui <- runif(n)
  tpsSpe <- 1 / pho * (alpha * ((1 / (1 - ui))^(1 / (alpha * exp.betaz)) - 1))^(1 / k) 
  
  # T_P generation
  tpsGene <- rep(0, n)
  TauxAtt <- rep(NA, n)
  sexNom <- ifelse(sex == 0, "male", "female")
  Xrace <- ifelse(X == 0, "white", "black")
  
  
  f1 <- function(i){
    uAtt <- runif(1)
    
    i.age <- which(attr(survexp.usr, which = "dimnames")[[1]] == trunc(age[i]))
    i.sex <- which(attr(survexp.usr, which = "dimnames")[[2]] == sexNom[i])
    i.race <- which(attr(survexp.usr, which = "dimnames")[[3]] == Xrace[i])
    i.year <- which(attr(survexp.usr, which = "dimnames")[[4]] == format(year.start[i], "%Y"))
    
    if(length(i.age) == 0 || length(i.sex) == 0 || length(i.year) == 0) {
      return(NA)  # Retourner NA si indices invalides
    }
    
    max.i.age <- length(attributes(survexp.usr)$dimnames[[1]])
    max.i.year <- length(attributes(survexp.usr)$dimnames[[4]])
    TauxAtt[i] <- 1-exp(-365.24*survexp.usr[i.age,i.sex, i.race, i.year])
    
    if(uAtt <= TauxAtt[i]){
      tpsG <- runif(1)
      tpsGene[i] <- tpsG
    }else{
      while(uAtt > TauxAtt[i]){
        tpsGene[i] <- tpsGene[i] + 1
        i.age <- min(i.age + 1, max.i.age) # happy birthday !
        i.year <- min(i.year + 1, max.i.year) # and one calendar year more...
        
        uAtt <- runif(1)
        TauxAtt[i] <- 1-exp(-365.24*survexp.usr[i.age,i.sex, i.race, i.year])
      }
      tpsG <- runif(1)
      tpsGene[i] <- tpsGene[i] + tpsG
    }
    return(tpsGene[i])
  }
  tpsGene <- sapply(1:n,f1)
  
  tpsSurv <- pmin(tpsGene, tpsSpe)
  
  # CENSORING
  borne_a = borne_a
  
  # censoring time
  tpsCens <- runif(n, min=0, max=borne_a)
  
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
  
  result <- data.frame(
    patient_id = 1:n,
    age = round(age, 1),
    ageStand = ageStand,
    sex_num = sex,
    sex = factor(sexNom,  levels = c("male", "female")),
    race = factor(Xrace, levels = c("white", "black")),
    race_num = X,
    tpsGene = tpsGene,
    tpsSpe = tpsSpe,
    year_diagnosis = year.start,
    observed_time = temps,
    status = statut, 
    hypothetical_time = temps2,
    hypothetical_status = cause2,
    admin_cens = max_time
  )
  
  return(result)
}