# Run this once outside of Quarto
source("plots/PPviz_parallel_function.R")
dir.create("plots_data", showWarnings = FALSE)

lambda_files <- list.files("data", pattern = "^simulated_cohort_lambda_.*\\.rds$", full.names = FALSE)
lambda_values <- sort(as.numeric(stringr::str_remove(stringr::str_remove(lambda_files, "^simulated_cohort_lambda_"), "\\.rds$")))

for (lambda_val in lambda_values) {
  message("Pre-computing plot data for Lambda = ", lambda_val)
  
  # Capture the invisible list returned by your function
  plot_data <- plot_ppviz_parallel(
    lambda_scenario = lambda_val,
    data_dir = "data",
    metrics_dir = "tables",
    use_parallel = TRUE
  )
  
  # Save purely the lightweight coordinates
  saveRDS(plot_data, sprintf("plots_data/ppviz_data_lambda_%.4f.rds", lambda_val))
  
  # Force garbage collection between iterations
  rm(plot_data)
  gc()
}