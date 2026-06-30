library(here)

source_path <- here("current", "plots", "age_stratified_estimate_plots", "PPviz_parallel_function_age_strata.R")
if (!file.exists(source_path)) stop("CRITICAL: Function script not found at ", source_path)
source(source_path)

data_dir    <- here("current", "outputs", "data")
metrics_dir <- here("current", "outputs", "tables")
save_dir    <- here("current", "outputs", "plots_data_age_strata")

dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

lambda_files <- list.files(data_dir, pattern = "^simulated_cohort_lambda_.*\\.parquet$", full.names = FALSE)
lambda_values <- sort(as.numeric(stringr::str_remove(stringr::str_remove(lambda_files, "^simulated_cohort_lambda_"), "\\.parquet$")))

for (lambda_val in lambda_values) {
  message("Pre-computing age strata plot for Lambda = ", lambda_val)
  
  plot_data <- plot_ppviz_parallel(
    lambda_scenario = lambda_val,
    data_dir = data_dir,
    metrics_dir = metrics_dir,
    use_parallel = TRUE
  )
  
  save_path <- file.path(save_dir, sprintf("ppviz_age_strata_lambda_%.4f.rds", lambda_val))
  saveRDS(plot_data, save_path)
  
  rm(plot_data)
  gc()
}

message("SUCCESS: Saved all pre-computed age strata plot objects directly into '", save_dir, "'")