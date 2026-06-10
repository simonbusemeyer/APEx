library(survival)
source("generate_dataModified.R")
source("analyze_one.R")
source("censoring_calibration.R")

# --- 0. Parameters ---
n_patients <- 2000 #
lambda <- 0.5 
beta_age <- 0.2
max_time <- 4
max_time_days <- max_time * 365.241
borne_a <- max_time_days

# --- 2. Timing check (run once on 1 dataset) ---
start <- proc.time()

df <- generate_data(
  lambda, age_option = "A", n = n_patients, max_time = max_time_days,
  prop_female = 0,
  year.start_min = 2008, year.start_max = 2010,
  beta_sex = 0, beta_age = beta_age, borne_a = borne_a
)

res <- analyze_one(df, lambda, beta_age, times_years = c(1, 2, 3))

elapsed <- proc.time() - start
elapsed

start <- proc.time()

# --- 3. Censoring calibration ---

start <- proc.time()

lambdas <- c(0.1) 
set.seed(12345)

# Call the external function
calibration_results <- calibrate_censoring_grid(
  lambdas = lambdas,
  n_patients = n_patients,
  max_time_days = max_time_days,
  beta_age = beta_age,
  n_pilots = 10
)

elapsed <- proc.time() - start
print(elapsed)

# View the final calibrated table
print(calibration_results)