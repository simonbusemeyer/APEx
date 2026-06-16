# main_batch.R
rm(list=ls())

library(survival)
library(relsurv)

source("current/functions/generate_dataModified_ng.R")
source("current/functions/analyze_one_ng.R")
source("current/functions/compute_metrics.R")

# Global Parameters
n_patients <- 2000
age_option <- "E"
beta_age <- 0.02
beta_sex <- 0
max_time <- 10
max_time_days <- max_time * 365.241
year.start_min <- 2008
year.start_max <- 2010
prop_female <- 0
N_files <- 100

#set.seed(12345) #placed inside each scenario instead for better comparison

# Reconcile output directory
if(!dir.exists("current/outputs/tables")) dir.create("current/outputs/tables", recursive = TRUE)

start <- proc.time()

# Execute Scenarios
# Each sourced file will run its loop and save an .rds file to outputs/tables/
source("current/lambda_batch/lambda_0.01_sim.R")#optional
source("current/lambda_batch/lambda_0.02_sim.R")#optional
source("current/lambda_batch/lambda_0.03_sim.R")#optional
source("current/lambda_batch/lambda_0.04_sim.R")#optional
source("current/lambda_batch/lambda_0.05_sim.R")
source("current/lambda_batch/lambda_0.07_sim.R") #optional
source("current/lambda_batch/lambda_0.10_sim.R")
source("current/lambda_batch/lambda_0.20_sim.R")
source("current/lambda_batch/lambda_0.30_sim.R")
source("current/lambda_batch/lambda_0.50_sim.R")
#source("lambda_batch/lambda_1.00_sim.R") #invalid results at time 3,4 due to censoring calibration
#source("lambda_batch/lambda_1.50_sim.R") #invalid results at time 2,3,4 due to censoring calibration

elapsed <- proc.time() - start
print(elapsed)

# Aggregate Results
# Read all individual scenario outputs and bind them into one dataframe
rds_files <- list.files("current/outputs/tables", pattern = "metrics_lambda_.*\\.rds$", full.names = TRUE)
all_metrics_list <- lapply(rds_files, readRDS)

final_results <- do.call(rbind, all_metrics_list)

#Save the finalized bound dataframe
saveRDS(final_results, "current/outputs/tables/final_results_complete.rds")

write.csv(final_results, 
          file = "current/outputs/tables/final_results_complete.csv", 
          row.names = FALSE)