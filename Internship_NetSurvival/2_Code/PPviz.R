# ---------------------------------------------------------
# VISUALIZING POHAR-PERME VS THEORETICAL
# ---------------------------------------------------------
library(survival)
library(relsurv)

# Parameters (Must match the target scenario)
lambda_scenario <- 0.05
max_time <- 4
beta_age <- 0.02
beta_sex <- 0
N_plot <- 15 

#Load Data from Batch Output folder
data_path <- sprintf("outputs/data/simulated_cohort_lambda_%.2f.rds", lambda_scenario)
metrics_path <- sprintf("outputs/tables/metrics_lambda_%.2f.rds", lambda_scenario)

if (!file.exists(data_path)) {
  stop(sprintf("Batch data file %s not found! Run main_batch.R first.", data_path))
}

all_simulated_data <- readRDS(data_path)
metrics_data       <- readRDS(metrics_path)

# Extract N_files dynamically based on the data
N_files <- max(all_simulated_data$sim_id)

#extract and round pct_cancer
pct_cancer_viz <- round(metrics_data$pct_cancer[1] * 100, 1)

# Convert Time Units to Days for relsurv
all_simulated_data$age_days <- all_simulated_data$age * 365.241
all_simulated_data$observed_time_days <- all_simulated_data$observed_time * 365.241

#Calculate Theoretical Curve
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

# Calculate Pooled Pohar-Perme (Using a Fast Sub-sample)

set.seed(12345)
fast_pool_idx <- sample(1:nrow(all_simulated_data), size = 15000)
data_pooled_fast <- all_simulated_data[fast_pool_idx, ]

pp_pooled <- rs.surv(
  Surv(observed_time_days, status) ~ 1,
  data = data_pooled_fast,
  ratetable = survexp.us,
  rmap = list(age = age_days, sex = sex, year = year_diagnosis),
  method = "pohar-perme"
)
message("Pooled calculation complete.")

#Plot Setup
plot(
  0, type = "n",
  xlim = c(0, max_time),
  ylim = c(0.3, 1.1), 
  xlab = "Time since diagnosis (Years)",
  ylab = "Net Survival Probability",
  main = paste0("Net Survival: PP vs Theoretical \n Proportion of deaths due to cancer: ", pct_cancer_viz, "%")
)
grid()

# randomly select 15 individual Pohar-Perme curves
set.seed(54321)
unique_sims <- unique(all_simulated_data$sim_id)
sampled_sims <- sample(unique_sims, size = min(N_plot, length(unique_sims)), replace = FALSE)

for (i in sampled_sims) {
  data_sim <- subset(all_simulated_data, sim_id == i)
  
  if(nrow(data_sim) > 0) {
    pp_sim <- rs.surv(
      Surv(observed_time_days, status) ~ 1,
      data = data_sim,
      ratetable = survexp.us,
      rmap = list(age = age_days, sex = sex, year = year_diagnosis),
      method = "pohar-perme", add.times = c(1,3,4) * 365.241
    )
    lines(pp_sim$time / 365.241, pp_sim$surv, col = rgb(0.2, 0.5, 0.8, alpha = 0.3), lwd = 1, type = "s")
  }
}

# Overlay Pooled Pohar-Perme & Theoretical Curves
lines(pp_pooled$time / 365.241, pp_pooled$surv, col = "red", lwd = 3, type = "s")
lines(t_seq_years, surv_theo_mean, col = "black", lwd = 3, lty = 2)

#Legend ---
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