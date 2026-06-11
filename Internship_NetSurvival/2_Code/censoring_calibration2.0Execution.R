# 1. Load required libraries and source the necessary files
library(survival)
source("generate_dataModified_ng.R")  # Needs to be sourced so generate_data() is available
source("censoring_calibration2.0.R")     # The file containing your updated function

# 2. Define the base parameters (from main_ng.R)
n_patients <- 2000 
beta_age <- 0.02
max_time_years <- 4
max_time_days <- max_time_years * 365.241

# 3. Define the lambdas you want to calibrate
lambdas_to_test <- c(0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0, 1.5)

# 4. Execute the calibration
# n_pilots = 100 matches your original N_files = 100 for stable averaging
# max_iter = 50 activates the new fail-safe to prevent infinite loops
calibration_summary <- calibrate_censoring_grid(
  lambdas = lambdas_to_test,
  n_patients = n_patients,
  max_time_days = max_time_days,
  beta_age = beta_age,
  n_pilots = 25, 
  max_iter = 40
)

# 5. View the final results table
print("=== FINAL CALIBRATION RESULTS ===")
print(calibration_summary)

# Optional: Save the results to a CSV so you don't lose them
# write.csv(calibration_summary, "calibration_results.csv", row.names = FALSE)