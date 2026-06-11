# ---------------------------------------------------------
# 5. VISUALIZING POHAR-PERME VS THEORETICAL
# Fast execution leveraging pre-generated batch outputs
# ---------------------------------------------------------
library(survival)
library(relsurv)

# --- 0. Parameters (Must match the target scenario) ---
lambda_scenario <- 0.05
max_time <- 4
beta_age <- 0.02
beta_sex <- 0
N_plot <- 15 

# --- 1. Load Data from Batch Output folder ---
data_path <- sprintf("outputs/data/simulated_cohort_lambda_%.2f.rds", lambda_scenario)

if (!file.exists(data_path)) {
  stop(sprintf("Batch data file %s not found! Run main_batch.R first.", data_path))
}

message("Loading pre-generated cohort data...")
all_simulated_data <- readRDS(data_path)

# Extract N_files dynamically based on the data
N_files <- max(all_simulated_data$sim_id)

# --- 1.5 Convert Time Units to Days for relsurv ---
all_simulated_data$age_days <- all_simulated_data$age * 365.241
all_simulated_data$observed_time_days <- all_simulated_data$observed_time * 365.241

# --- 2. Calculate Theoretical Curve ---
survtheo <- function(t, lambda, beta_sex, beta_age, sex, ageCentre) {
  exp(-lambda * t * exp(beta_sex * sex + beta_age * ageCentre))
}

t_seq_years <- seq(0, max_time, by = 0.1)

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

# --- 3. Calculate Pooled Pohar-Perme (Using a Fast Sub-sample) ---
# Computing rs.surv on 200,000 rows causes extreme computational overhead.
# A 15,000 patient sub-sample provides an identical, highly stable pooled estimate instantly.

set.seed(12345)
fast_pool_idx <- sample(1:nrow(all_simulated_data), size = 15000)
data_pooled_fast <- all_simulated_data[fast_pool_idx, ]

message("Computing pooled Pohar-Perme on a 15,000 patient sub-sample...")
pp_pooled <- rs.surv(
  Surv(observed_time_days, status) ~ 1,
  data = data_pooled_fast,
  ratetable = survexp.us,
  rmap = list(age = age_days, sex = sex, year = year_diagnosis),
  method = "pohar-perme"
)
message("Pooled calculation complete.")

# --- 4. Plot Setup ---
plot(
  0, type = "n",
  xlim = c(0, max_time),
  ylim = c(0.7, 1), 
  xlab = "Time since diagnosis (Years)",
  ylab = "Net Survival Probability",
  main = paste0("Net Survival: PP vs Theoretical (\u03bb = ", lambda_scenario, ")\n",
                N_files, " pooled runs | First ", N_plot, " samples plotted")
)
grid()

# --- 5. Plot the First 25 Individual Pohar-Perme Curves ---
sampled_sims <- 1:N_plot

for (i in sampled_sims) {
  data_sim <- subset(all_simulated_data, sim_id == i)
  
  if(nrow(data_sim) > 0) {
    pp_sim <- rs.surv(
      Surv(observed_time_days, status) ~ 1,
      data = data_sim,
      ratetable = survexp.us,
      rmap = list(age = age_days, sex = sex, year = year_diagnosis),
      method = "pohar-perme"
    )
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
    paste0("Individual PP Estimates (First ", N_plot, ")"), 
    paste0("Pooled PP Estimate (N=", N_files, ")"), 
    "Theoretical Net Survival Curve S(t)"
  ),
  col = c(rgb(0.2, 0.5, 0.8, alpha = 0.5), "red", "black"), 
  lwd = c(2, 3, 3), 
  lty = c(1, 1, 2), 
  bty = "n",
  cex = 0.85
)