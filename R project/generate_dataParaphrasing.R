k <- 1.5
pho <- 0.05
alpha <- 100
n <- 10000
max_time <- 50
prop_female <- 0.5
prop_x0 <- 0.3
year.start_min <- 2000
year.start_max <- 2010
beta_sex <- 0.3
beta_age <- 0.2
beta_X <- 0.4
borne_a <- 10

# --- Phase 1: Covariables generation ---

age_1 <- runif(n * .25, min = 30, max = 65)
age_2 <- runif(n * .35, min = 65, max = 75)
age_3 <- runif(n * .40, min = 75, max = 80)
age <- c(age_1, age_2, age_3)
ageMoyen <- mean(age)
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

uSex <- runif(n)
pFem <- prop_female
sex  <- rbinom(n, size = 1, prob = prop_female)   

uX   <- runif(n)
X    <- rbinom(n, size = 1, prob = 1 - prop_x0)   

# continuous enrollement?
if (year.start_min != year.start_max) {
  year.start <- as.Date(paste0(sample(year.start_min:year.start_max, n, replace = TRUE), "-01-01")) +
    sample(0:364, n, replace = TRUE)
} else {
  year.start <- as.Date(paste0(year.start_min, "-01-01")) +
    sample(0:364, n, replace = TRUE)
}

# --- Phase 2: T_E generation (Time to death from cancer) ---
# why is this not T_O generation?
# calculating the risk multiplier before TE generation?
tempuS <- runif(n)
exp.betaz <- exp(
  beta_sex * sex + beta_age * ageStand + beta_X * X
)

# generalized weibull extension
ui <- runif(n)
tpsSpe <- 1 / pho * (alpha * ((1 / (1 - ui))^(1 / (alpha * exp.betaz)) - 1))^(1 / k)  
hist(tpsSpe)


# --- Phase 3: T_P generation (Time to death from background/general causes) ---

tpsGene <- rep(0, n)
TauxAtt <- rep(NA, n)
sexNom <- ifelse(sex == 0, "male", "female")
Xrace <- ifelse(X == 0, "white", "black")

f1 <- function(i){
  # the random chance of competing risk
  uAtt <- runif(1)
  
  # lookup the coresponding demographics
  i.age <- which(attr(survexp.usr, which = "dimnames")[[1]] == trunc(age[i]))
  i.sex <- which(attr(survexp.usr, which = "dimnames")[[2]] == sexNom[i])
  i.race <- which(attr(survexp.usr, which = "dimnames")[[3]] == Xrace[i])
  i.year <- which(attr(survexp.usr, which = "dimnames")[[4]] == format(year.start[i], "%Y"))
  
#return NA if the index is missing
  if(length(i.age) == 0 || length(i.sex) == 0 || length(i.year) == 0) {
    return(NA)  
  }
  
  # retrieve the position of the table maximums
  max.i.age <- length(attributes(survexp.usr)$dimnames[[1]])
  max.i.year <- length(attributes(survexp.usr)$dimnames[[4]])
  # Convert the daily hazard rate from the life table into an annual probability of dying
  TauxAtt[i] <- 1-exp(-365.24*survexp.usr[i.age,i.sex, i.race, i.year])
  
  #  if random probability is less than risk of dying then they are deemed dead
  if(uAtt <= TauxAtt[i]){
    # randomize time of death during the year
    tpsG <- runif(1)
    tpsGene[i] <- tpsG
  } else {
    # begin a while loop for those who keep living
    while(uAtt > TauxAtt[i]){
      # count up indexes
      tpsGene[i] <- tpsGene[i] + 1
      i.age <- min(i.age + 1, max.i.age) 
      i.year <- min(i.year + 1, max.i.year) 
      
      # iterate for the next year
      uAtt <- runif(1)
      TauxAtt[i] <- 1-exp(-365.24*survexp.usr[i.age,i.sex, i.race, i.year])
    }
    tpsG <- runif(1)
    tpsGene[i] <- tpsGene[i] + tpsG
  }
  # Return the total accumulated time to background death
  return(tpsGene[i])
}
# run the function for all n patients
tpsGene <- sapply(1:n,f1)

hist(tpsGene)

# --- Phase 4: Competing Risks & Censoring Resolution ---

# take the minimum of the patient's population and cancer survival time 
tpsSurv <- pmin(tpsGene, tpsSpe)
hist(tpsSurv)

# simulate random loss to followup time
borne_a = borne_a
tpsCens <- runif(n, min=0, max=borne_a)

# observed patient time is their survival time or the point they dropped out before.
temps <- pmin(tpsCens, tpsSurv)
hist(temps)
# The hypothetical time observed for net survival where patients can only die of cancer type or drop out first.
temps2 <- pmin(tpsCens, tpsSpe)
hist(temps2)

# classic status indicator equals 0 if the patient was censored or 1 if they were observed to die of any cause.
statut     <- ifelse(temps == tpsCens, 0, 1)

#indicate if patient died of cancer by making sure survival time is equal to excess survival time and that total time is not equal to censored time 
cause <- ifelse(tpsSurv == tpsSpe & temps != tpsCens, 1, 0)

# hypothetical indicator if patient could only die of cancer or be censored
cause2     <- ifelse(temps2 == tpsCens, 0, 1)

# --- Phase 5: Administrative Censoring ---

# if patient time is greater than max time then mark censored
statut[temps > max_time]   <-  0
# if patient time exceeds max study time, then they are deemed not dead of cancer
cause[temps > max_time] <- 0
cause2[temps2 > max_time]  <-  0

# any patient time greater than the max study time is set to shortened to the max time
temps[temps > max_time]    <- max_time
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