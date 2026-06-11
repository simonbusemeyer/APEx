# ---------------------------------------------------------
# 5. VISUALIZING POHAR-PERME VS THEORETICAL
# Standalone execution adapted for 1 value of lambda
# ---------------------------------------------------------
library(survival)
library(relsurv)

# Source data generation 
source("generate_dataModified_ng.R")

message("Setting up parameters and generating data for all simulations. This may take a moment...")

# --- 0. Global Parameters ---
n_patients <- 2000
age_option <- "A"
beta_age <- 0.02
beta_sex <- 0
max_time <- 4
max_time_days <- max_time * 365.241
year.start_min <- 2008
year.start_max <- 2010
prop_female <- 0
N_files <- 100
N_plot <- 25 # Number of curves to randomly sample for the plot

# Select 1 value of lambda
lambda_scenario <- 0.05
borne_a_scenario <- Inf

# --- 1. Generate & Pool Data for all N_files ---
set.seed(12345)
df_list <- vector("list", N_files)

for (j in 1:N_files) {
  df_list[[j]] <- generate_data(
    lambda = lambda_scenario,     
    age_option = age_option,   
    n = n_patients,        
    max_time = max_time, 
    prop_female = prop_female,     
    year.start_min = year.start_min, 
    year.start_max = year.start_max,
    beta_sex = beta_sex,         
    beta_age = beta_age, 
    borne_a = borne_a_scenario         
  )
  # Add tracking ID for random sampling later
  df_list[[j]]$sim_id <- j 
}
all_simulated_data <- do.call(rbind, df_list)

# --- 1.5 Convert Time Units to Days for relsurv ---
all_simulated_data$age_days <- all_simulated_data$age * 365.241
all_simulated_data$observed_time_days <- all_simulated_data$observed_time * 365.241

# --- 2. Calculate Theoretical Curve ---
survtheo <- function(t, lambda, beta_sex, beta_age, sex, ageCentre) {
  exp(-lambda * t * exp(beta_sex * sex + beta_age * ageCentre))
}

t_seq_years <- seq(0, max_time, by = 0.1)

# Averages over the full pooled simulated cohort
surv_theo_matrix <- sapply(1:nrow(all_simulated_data), function(i) {
  survtheo(
    t         = t_seq_years,
    lambda    = lambda_scenario,
    beta_sex  = beta_sex, 
    beta_age  = beta_age,
    sex       = all_simulated_data$sex_num[i],
    ageCentre = all_simulated_data$ageCentre[i] 
  )
})
surv_theo_mean <- rowMeans(surv_theo_matrix)

# --- 3. Calculate Pooled Pohar-Perme (All N_files) ---
pp_pooled <- rs.surv(
  Surv(observed_time_days, status) ~ 1,
  data = all_simulated_data,
  ratetable = survexp.us,
  rmap = list(age = age_days, sex = sex, year = year_diagnosis),
  method = "pohar-perme"
)

# --- 4. Plot Setup ---
plot(
  0, type = "n",
  xlim = c(0, max_time),
  ylim = c(0.7, 1), # Adjusted lower bound to reasonably fit lambda=0.05
  xlab = "Time since diagnosis (Years)",
  ylab = "Net Survival Probability",
  main = paste0("Net Survival: PP vs Theoretical (\u03bb = ", lambda_scenario, ")\n",
                N_files, " pooled runs | ", N_plot, " random samples plotted")
)
grid()

# --- 5. Plot a Random Sample of Individual Pohar-Perme Curves ---
set.seed(42) # Reproducible sample mapping
sampled_sims <- sample(1:N_files, N_plot)

for (i in sampled_sims) {
  data_sim <- subset(all_simulated_data, sim_id == i)
  
  # Safety check to ensure the subset isn't empty before fitting
  if(nrow(data_sim) > 0) {
    pp_sim <- rs.surv(
      Surv(observed_time_days, status) ~ 1,
      data = data_sim,
      ratetable = survexp.us,
      rmap = list(age = age_days, sex = sex, year = year_diagnosis),
      method = "pohar-perme"
    )
    
    # Individual PP curves in faint blue
    lines(pp_sim$time / 365.241, pp_sim$surv, col = rgb(0.2, 0.5, 0.8, alpha = 0.3), lwd = 1, type = "s")
  }
}

# --- 6. Overlay Pooled Pohar-Perme & Theoretical Curves ---
lines(pp_pooled$time / 365.241, pp_pooled$surv, col = "red", lwd = 3, type = "s")
lines(t_seq_years, surv_theo_mean, col = "black", lwd = 3, lty = 2)

# --- 7. Comprehensive Legend ---
legend(
  "bottomleft", 
  legend = c(
    paste0("Individual PP Estimates (n=", N_plot, ")"), 
    paste0("Pooled PP Estimate (N=", N_files, ")"), 
    "Theoretical Net Survival Curve S(t)"
  ),
  col = c(rgb(0.2, 0.5, 0.8, alpha = 0.5), "red", "black"), 
  lwd = c(2, 3, 3), 
  lty = c(1, 1, 2), 
  bty = "n",
  cex = 0.85
)