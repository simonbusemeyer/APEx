library(arrow)
library(dplyr)

# Maps the directory without loading it into RAM
cohort_dataset <- arrow::open_dataset(
  sources = "current/outputs/data", 
  format = "parquet"
)

# Performs the aggregation safely at the C++ level
lambda_verification <- cohort_dataset %>%
  group_by(lambda) %>%
  summarize(total_rows = n()) %>%
  collect()

print(lambda_verification)