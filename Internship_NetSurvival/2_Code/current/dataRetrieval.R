library(data.table)

# 1. Identify all cohort .rds files
dataset_files <- list.files(
  path = "current/outputs/data",
  pattern = "simulated_cohort_lambda_.*\\.rds$", 
  full.names = TRUE
)

# 2. Read files and add the lambda value
all_cohorts_list <- lapply(dataset_files, function(file_path) {
  # Load the data
  dt <- readRDS(file_path)
  
  extracted_lambda <- as.numeric(gsub(".*lambda_([0-9.]+)\\.rds", "\\1", basename(file_path)))
  
  dt$lambda <- extracted_lambda
  
  return(dt)
})

# 3. Bind into one master dataset
master_cohort_data <- rbindlist(all_cohorts_list, fill = TRUE)

# 4. Verify again
cat("Total dimensions of the master dataset:\n")
print(dim(master_cohort_data))

cat("\nVerification of lambda scenarios:\n")
print(table(master_cohort_data$lambda))