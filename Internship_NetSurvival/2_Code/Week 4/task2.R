library(survival)
library(relsurv)

simulate_censoring_rate <- function(lambda, beta_age, borne_a, max_time_years, n = 2000) {
  
  age <- runif(n, min = 80, max = 89)
  ageStand <- (age - mean(age)) / sd(age)
  sexNom <- rep("male", n)
  year.start <- as.Date("2008-01-01") + sample(0:364, n, replace = TRUE)
  
  exp.betaz <- exp(beta_age * ageStand)
  tpsSpe <- -log(runif(n)) / (lambda * exp.betaz)
  
  tpsGene <- numeric(n)
  for (i in 1:n) {
    uAtt <- runif(1)
    
    i.age <- match(trunc(age[i]), attr(survexp.us, "dimnames")[[1]])
    i.sex <- match(sexNom[i], attr(survexp.us, "dimnames")[[2]])
    i.year <- match(format(year.start[i], "%Y"), attr(survexp.us, "dimnames")[[3]])
    
    max.i.age <- length(attr(survexp.us, "dimnames")[[1]])
    max.i.year <- length(attr(survexp.us, "dimnames")[[3]])
    
    TauxAtt <- 1 - exp(-365.241 * survexp.us[i.age, i.sex, i.year])
    
    if (uAtt <= TauxAtt) {
      tpsGene[i] <- runif(1)
    } else {
      while (uAtt > TauxAtt) {
        tpsGene[i] <- tpsGene[i] + 1
        i.age <- min(i.age + 1, max.i.age)
        i.year <- min(i.year + 1, max.i.year)
        uAtt <- runif(1)
        TauxAtt <- 1 - exp(-365.24 * survexp.us[i.age, i.sex, i.year])
      }
      tpsGene[i] <- tpsGene[i] + runif(1)
    }
  }
  
  #Censoring Logic & Administrative Boundary
  tpsSurv <- pmin(tpsGene, tpsSpe)
  tpsCens <- runif(n, min = 0, max = borne_a)
  
  temps <- pmin(tpsCens, tpsSurv)
  temps[temps > max_time_years] <- max_time_years
  
  #Calculate Censored Proportion
  # An event is censored if observation ended due to uniform C OR administrative boundary
  is_censored <- (temps == tpsCens) | (temps == max_time_years)
  
  return(mean(is_censored))
}


#Calibration Loop
calibrate_scenario <- function(lambda_val, beta_age_val, max_time_years) {
  borne_a <- max_time_years 
  
  cat(sprintf("--- Calibrating scenario: lambda = %.4f ---\n", lambda_val))
  
  repeat {
    # Generate 100 datasets and compute mean censoring rate
    rates <- replicate(100, simulate_censoring_rate(lambda = lambda_val, 
                                                    beta_age = beta_age_val, 
                                                    borne_a = borne_a, 
                                                    max_time_years = max_time_years))
    
    mean_cens <- mean(rates)
    cat(sprintf("  borne_a: %7.2f | Censoring Rate: %5.2f%%\n", borne_a, mean_cens * 100))
    
    # Adjust boundary parameters
    if (mean_cens > 0.32) {
      borne_a <- borne_a * 1.1
    } else if (mean_cens < 0.28) {
      borne_a <- borne_a * 0.9
    } else {
      cat(sprintf("  SUCCESS: Target reached. Saved borne_a = %.2f\n\n", borne_a))
      break
    }
  }
  
  return(borne_a)
}

# Define your specific lambda scenarios (update as needed)
lambda_scenarios <- c(0.01, 0.03, 0.05, 0.10) 
beta_age_val <- 0.05
max_fup_years <- 4 # Extracted from Task 1

# Initialize an empty dataframe
results_table <- data.frame(
  Lambda = numeric(),
  Borne_A_Years = numeric(),
  Censoring_Rate_Pct = numeric()
)

set.seed(12345)

for (l in lambda_scenarios) {
  # 1. Calibrate the boundary for this lambda
  final_borne_a <- calibrate_scenario(lambda_val = l, 
                                      beta_age_val = beta_age_val, 
                                      max_time_years = max_fup_years)
  
  # 2. Run a final verification batch to get the exact observed rate to report
  final_rates <- replicate(100, simulate_censoring_rate(lambda = l, 
                                                        beta_age = beta_age_val, 
                                                        borne_a = final_borne_a, 
                                                        max_time_years = max_fup_years))
  
  observed_rate <- mean(final_rates) * 100
  
  # 3. Store the results
  results_table <- rbind(results_table, data.frame(
    Lambda = l,
    Borne_A_Years = round(final_borne_a, 2),
    Censoring_Rate_Pct = sprintf("%.2f%%", observed_rate)
  ))
}

# Print the final Markdown table
cat("\n### Final Calibration Table\n\n")
print(knitr::kable(results_table, format = "markdown", 
                   col.names = c("λ (Excess Hazard)", "Calibrated borne_a (Years)", "Observed Censoring Rate (%)"),
                   align = "c"))
