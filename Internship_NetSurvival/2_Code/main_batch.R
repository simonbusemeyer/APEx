# main_batch.R
rm(list=ls())

library(survival)
library(relsurv)

source("generate_dataModified_ng.R")
source("analyze_one_ng.R")
source("compute_metrics.R")

# --- 0. Global Parameters --- ####
n_patients <- 2000
age_option <- "A"
beta_age <- 0.02
beta_sex <- 0
max_time <- 4
max_time_days <- max_time * 365.241
year.start_min <- 2008
year.start_max <- 2010
prop_female <- 0
N_files <- 10

#set.seed(12345)

# Reconcile output directory
if(!dir.exists("outputs/tables")) dir.create("outputs/tables", recursive = TRUE)

start <- proc.time()

# --- 1. Execute Scenarios --- ####
# Each sourced file will run its loop and save an .rds file to outputs/tables/
source("lambda batch/lambda_0.05_sim.R")
#source("lambda batch/lambda_0.10_sim.R")
#source("lambda batch/lambda_0.20_sim.R")
#source("lambda batch/lambda_0.30_sim.R")
#source("lambda batch/lambda_0.50_sim.R")
#source("lambda batch/lambda_0.1.00_sim.R")
#source("lambda batch/lambda_1.50_sim.R")

elapsed <- proc.time() - start
print(elapsed)

# --- 2. Aggregate Results --- ####
# Read all individual scenario outputs and bind them into one dataframe
rds_files <- list.files("outputs/tables", pattern = "metrics_lambda_.*\\.rds$", full.names = TRUE)
all_metrics_list <- lapply(rds_files, readRDS)

final_results <- do.call(rbind, all_metrics_list)

# Optional: Save the finalized bound dataframe
saveRDS(final_results, "outputs/tables/final_results_complete.rds")