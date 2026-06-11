# =============================================================================
# generate_data.R
# Fonctions de génération de données de survie
# =============================================================================
#is this weibull?

k <- 1
pho <- 1
alpha <- 1
n <- 10000
max_time <- 50
prop_female <- 0.5
prop_x0 <- 0.3
year.start_min <- 2000
year.start_max <- 2010
beta_sex <- 0.3
beta_age <- 0.2
beta_X <- 0.4
borne_a <- 1
  
  # --- Phase 1: Covariables generation ---
  
  # Generate ages for 25% of the cohort using a uniform distribution between 30 and 65 years
  age_1 <- runif(n * .25, min = 30, max = 65)
  # Generate ages for 35% of the cohort using a uniform distribution between 65 and 75 years
  age_2 <- runif(n * .35, min = 65, max = 75)
  # Generate ages for 40% of the cohort using a uniform distribution between 75 and 80 years
  age_3 <- runif(n * .40, min = 75, max = 80)
  # Combine the three age groups into a single vector of length n
  age <- c(age_1, age_2, age_3)
  # Calculate the mean age of the entire generated cohort
  ageMoyen <- mean(age)
  # Standardize the age (Z-score calculation) for numerical stability in regression models
  ageStand <- (age - ageMoyen) / sd(age)
  
  par(mfrow = c(1, 2), mar = c(5, 4, 4, 2) + 0.1)
  
  # Histogram 1: Raw Age Distribution
  hist(age, 
       breaks = 40, 
       col = "darkslategray3", 
       border = "white",
       main = "Raw Age Distribution\n(Mixed Uniform Cohorts)", 
       xlab = "Age (Years)", 
       ylab = "Frequency")
  abline(v = ageMoyen, col = "firebrick", lwd = 2, lty = 2) # Mark the sample mean
  text(ageMoyen - 2, n*0.02, paste("Mean:", round(ageMoyen, 1)), col = "firebrick", pos = 2)
  
  # Histogram 2: Standardized Age Distribution
  hist(ageStand, 
       breaks = 40, 
       col = "coral1", 
       border = "white",
       main = "Standardized Age Distribution\n(Z-score Scaled)", 
       xlab = "Z-score Value", 
       ylab = "Frequency")
  abline(v = 0, col = "darkblue", lwd = 2, lty = 2) # Mark Z-score mean (0)
  
  
  
  # Generate n random uniform values (0 to 1) for sex assignment (though unused directly below)
  uSex <- runif(n)
  # Store the target proportion of females
  pFem <- prop_female
  # Generate binary sex variable using binomial distribution: 0 = male, 1 = female
  sex  <- rbinom(n, size = 1, prob = prop_female)   
  
  # Generate n random uniform values (0 to 1) for X assignment (unused directly below)
  uX   <- runif(n)
  # Generate binary variable X using binomial distribution: 1 - prop_x0 gives the probability of X=1
  X    <- rbinom(n, size = 1, prob = 1 - prop_x0)   
  
  # Check if there is a range of possible starting years
  if (year.start_min != year.start_max) {
    # Randomly select a year within the range, set to Jan 1st, then add a random number of days (0-364) to get a random date
    year.start <- as.Date(paste0(sample(year.start_min:year.start_max, n, replace = TRUE), "-01-01")) +
      sample(0:364, n, replace = TRUE)
  } else {
    # If min and max year are the same, just use that single year, set to Jan 1st, and add a random number of days
    year.start <- as.Date(paste0(year.start_min, "-01-01")) +
      sample(0:364, n, replace = TRUE)
  }
  
  # --- Phase 2: T_E generation (Time to death from cancer) ---
  
  # Generate n random uniform values (0 to 1) for survival time generation
  tempuS <- runif(n)
  # Calculate the linear predictor (exponential of the sum of covariate effects multiplied by their betas)
  exp.betaz <- exp(
    beta_sex * sex + beta_age * ageStand + beta_X * X
  )
  
  # Generate n random uniform probabilities for the inverse probability transform method
  ui <- runif(n)
  # Use the inverse cumulative hazard function of the specific complex distribution to calculate exact time to cancer death
  tpsSpe <- 1 / pho * (alpha * ((1 / (1 - ui))^(1 / (alpha * exp.betaz)) - 1))^(1 / k)  
  
  
  # --- Phase 3: T_P generation (Time to death from background/general causes) ---
  
  # Initialize a vector of 0s to store the background survival times for each patient
  tpsGene <- rep(0, n)
  # Initialize a vector of NAs to temporarily store the expected annual hazard rate for each patient
  TauxAtt <- rep(NA, n)
  # Convert numeric sex (0/1) to character strings ("male"/"female") to match life table formatting
  sexNom <- ifelse(sex == 0, "male", "female")
  # Convert numeric X (0/1) to race strings ("white"/"black") to match the US life table formatting
  Xrace <- ifelse(X == 0, "white", "black")
  
  # Define a custom function to calculate background death time for a single patient 'i'
  f1 <- function(i){
    # Roll a random probability to test against the patient's annual risk of dying
    uAtt <- runif(1)
    
    # Find the row index in the life table matching the patient's truncated (integer) age
    i.age <- which(attr(survexp.usr, which = "dimnames")[[1]] == trunc(age[i]))
    # Find the column index matching the patient's sex
    i.sex <- which(attr(survexp.usr, which = "dimnames")[[2]] == sexNom[i])
    # Find the depth index matching the patient's race
    i.race <- which(attr(survexp.usr, which = "dimnames")[[3]] == Xrace[i])
    # Find the dimension index matching the calendar year of diagnosis
    i.year <- which(attr(survexp.usr, which = "dimnames")[[4]] == format(year.start[i], "%Y"))
    
    # Safety check: if any index is missing (e.g., age not in table), return NA for this patient
    if(length(i.age) == 0 || length(i.sex) == 0 || length(i.year) == 0) {
      return(NA)  
    }
    
    # Store the maximum age index available in the life table to prevent out-of-bounds errors
    max.i.age <- length(attributes(survexp.usr)$dimnames[[1]])
    # Store the maximum year index available in the life table
    max.i.year <- length(attributes(survexp.usr)$dimnames[[4]])
    # Convert the daily hazard rate from the life table into an annual probability of dying
    TauxAtt[i] <- 1-exp(-365.24*survexp.usr[i.age,i.sex, i.race, i.year])
    
    # Test if the patient died this year (random probability is less than their risk)
    if(uAtt <= TauxAtt[i]){
      # If they died, assign a random fraction of the current year as their exact time of death
      tpsG <- runif(1)
      tpsGene[i] <- tpsG
    } else {
      # If they survived the year, enter a while loop to keep simulating year-by-year until they die
      while(uAtt > TauxAtt[i]){
        # Add 1 full year to their survival time
        tpsGene[i] <- tpsGene[i] + 1
        # Increase their age index by 1 (happy birthday!)
        i.age <- min(i.age + 1, max.i.age) 
        # Increase their calendar year index by 1
        i.year <- min(i.year + 1, max.i.year) 
        
        # Roll a new random probability for the new year
        uAtt <- runif(1)
        # Look up their new, updated annual risk of dying using their older age and new calendar year
        TauxAtt[i] <- 1-exp(-365.24*survexp.usr[i.age,i.sex, i.race, i.year])
      }
      # Once they finally fail the survival test, add a random fraction of their final year to the total time
      tpsG <- runif(1)
      tpsGene[i] <- tpsGene[i] + tpsG
    }
    # Return the total accumulated time to background death
    return(tpsGene[i])
  }
  # Apply the function f1 to all n patients, returning a vector of background survival times
  tpsGene <- sapply(1:n,f1)
  
  # --- Phase 4: Competing Risks & Censoring Resolution ---
  
  # The biological survival time is whichever happens first: background death or cancer death
  tpsSurv <- pmin(tpsGene, tpsSpe)
  
  # Reassign the borne_a (max dropout time) variable to itself
  borne_a = borne_a
  
  # Generate a random dropout/censoring time for every patient between time 0 and borne_a
  tpsCens <- runif(n, min=0, max=borne_a)
  
  # The actual observed real-world time is the minimum of their death time and their dropout time
  temps <- pmin(tpsCens, tpsSurv)
  # The hypothetical time (for net survival truth) compares only cancer death and dropout (ignoring background)
  temps2 <- pmin(tpsCens, tpsSpe)
  
  # Status is 0 if they dropped out (censored), or 1 if they were observed dying
  statut     <- ifelse(temps == tpsCens, 0, 1)
  
  # Cause is 1 (cancer) if their biological death was from cancer AND they didn't drop out first; otherwise 0
  cause <- ifelse(tpsSurv == tpsSpe & temps != tpsCens, 1, 0)
  
  # In the hypothetical world, status is 0 if they dropped out, or 1 if they died of cancer
  cause2     <- ifelse(temps2 == tpsCens, 0, 1)
  
  # --- Phase 5: Administrative Censoring ---
  
  # If real-world observed time exceeds the max study time, mark them as censored (0)
  statut[temps > max_time]   <-  0
  # If real-world observed time exceeds max study time, their cause is also marked as 0
  cause[temps > max_time] <- 0
  # If hypothetical time exceeds max study time, mark as censored (0) in the hypothetical world
  cause2[temps2 > max_time]  <-  0
  
  # Hard-cap all real-world observed times that exceed the max study time to exactly max_time
  temps[temps > max_time]    <- max_time
  # Hard-cap all hypothetical times that exceed the max study time to exactly max_time
  temps2[temps2 > max_time]  <- max_time
  
  # --- Phase 6: Output Assembly ---
  
  # Bind all generated and calculated vectors into a single data.frame
  result <- data.frame(
    patient_id = 1:n,                                      # Unique ID for each patient
    age = round(age, 1),                                   # Raw age rounded to 1 decimal
    ageStand = ageStand,                                   # Standardized age used in the model
    sex_num = sex,                                         # Binary sex (0=male, 1=female)
    sex = factor(sexNom,  levels = c("male", "female")),   # Factored string representation of sex
    race = factor(Xrace, levels = c("white", "black")),    # Factored string representation of race/X
    race_num = X,                                          # Binary representation of X
    tpsGene = tpsGene,                                     # True time to background death
    tpsSpe = tpsSpe,                                       # True time to cancer death
    year_diagnosis = year.start,                           # Date of diagnosis
    observed_time = temps,                                 # Final observed survival time
    status = statut,                                       # Final observed status (1=dead, 0=censored)
    hypothetical_time = temps2,                            # Pure time to cancer death with censoring
    hypothetical_status = cause2,                          # Status for hypothetical world (1=cancer death, 0=censored)
    admin_cens = max_time                                  # The max study follow-up time applied
  )
  
  # Output the final assembled dataset
  