library(survival)
library(relsurv)

source("generate_dataModified_ng.R")
source("analyze_one_ng.R")

# --- 0. Parameters --- ####
n_patients <- 2000 
lambda <- 0.01 
beta_age <- 0.02
max_time <- 4
max_time_days <- max_time * 365.241
borne_a <- Inf
N_files <- 100


# --- 2. Generating 1 dataset --- ####
df <- generate_data(
  lambda, age_option = "A", n = n_patients, max_time = max_time,#max_time_days,
  prop_female = 0,
  year.start_min = 2008, year.start_max = 2010,
  beta_sex = 0, beta_age = beta_age, borne_a = borne_a
)


# --- 2. (Task 1) nessie --- ####
#Create age groups
breaks <- pretty(df$age, n = 5)
df$agegr <- cut(df$age, breaks = breaks, include.lowest = TRUE)

nessie_output <- nessie(
  Surv(observed_time*365.241, status) ~ sex + agegr,
  data = df,
  ratetable = survexp.us,
  times = seq(0, max_time, 1),
  rmap = list(age = age*365.241, sex = sex, year = year_diagnosis)
)

# 0     1     2     3     4 c.exp.surv
# sexfemale=0,agegr(82,84] 448 410.7 373.1 335.6 298.4        6.7
# sexfemale=0,agegr(84,86] 441 395.5 351.1 307.8 266.0        5.8
# sexfemale=0,agegr(86,88] 453 395.4 340.7 289.0 241.0        5.0
# sexfemale=0,agegr(88,90] 194 165.6 139.2 114.9  93.1        4.5
# sexfemale=0,agegr[80,82] 464 432.5 400.7 368.3 335.6        7.7

# chosen follow-up = 4 years

# --- 3. (Task 2) censoring calibration and proportion of cacner deaths --- ####
datasets <- list()
censoring_rate <- numeric(length = N_files)

borne_a <- 30
for(j in 1:N_files){
  datasets[[j]] <- generate_data(
    lambda, age_option = "A", n = n_patients, max_time = max_time,#max_time_days,
    prop_female = 0,
    year.start_min = 2008, year.start_max = 2010,
    beta_sex = 0, beta_age = beta_age, borne_a = borne_a
  )
  
  censoring_rate[j] <- sum(datasets[[j]]$status==0) / nrow(datasets[[j]])
}

(mean(censoring_rate))

# borne_a = 10 ---- censoring_rate = 0.606575
# borne_a = 20 ---- censoring_rate = 0.56168
# borne_a = 100 ---- censoring_rate = 0.5289

# => PB !

# quick checks
nb_total_deaths <- sum(datasets[[1]]$status==1) 
(cancer_deaths_prop <- sum(datasets[[1]]$status==1 & datasets[[1]]$observed_time==datasets[[1]]$tpsSpe)/nb_total_deaths)
# 0.3213909

lambdas <- c(0.05, 0.1, 0.2, 0.3, 0.5, 1.0, 1.5)
prop_cancer <- lapply(lambdas, function(lam) {
  d <- generate_data(lambda = lam, age_option = "A", n = 2000,
                     max_time = 4, prop_female = 0,
                     year.start_min = 2008, year.start_max = 2010,
                     beta_sex = 0, beta_age = 0.02, borne_a = Inf)
  data.frame(
    lambda      = lam,
    prop_cancer = mean(d$cause == 1),
    cens_rate   = mean(d$status == 0)
  )
})
do.call(rbind, prop_cancer)
#lambda prop_cancer cens_rate
#1   0.05      0.1440    0.5400
#2   0.10      0.2730    0.4310
#3   0.20      0.4510    0.2915
#4   0.30      0.5915    0.2045
#5   0.50      0.7410    0.1015
#6   1.00      0.8920    0.0100
#7   1.50      0.9355    0.0035


# Only administrative censoring for lambda =0.05 and 0.10
# For other values of lambda

# lambda = 0.2 
datasets <- list()
censoring_rate <- numeric(length = N_files)

borne_a <- 30
for(j in 1:N_files){
  datasets[[j]] <- generate_data(
    lambda=0.2, age_option = "A", n = n_patients, max_time = max_time,
    prop_female = 0,
    year.start_min = 2008, year.start_max = 2010,
    beta_sex = 0, beta_age = beta_age, borne_a = borne_a
  )
  
  censoring_rate[j] <- sum(datasets[[j]]$status==0) / nrow(datasets[[j]])
}

(mean(censoring_rate)) # 0.3252

# lambda = 0.3 
datasets <- list()
censoring_rate <- numeric(length = N_files)

borne_a <- 12
for(j in 1:N_files){
  datasets[[j]] <- generate_data(
    lambda=0.3, age_option = "A", n = n_patients, max_time = max_time,
    prop_female = 0,
    year.start_min = 2008, year.start_max = 2010,
    beta_sex = 0, beta_age = beta_age, borne_a = borne_a
  )
  
  censoring_rate[j] <- sum(datasets[[j]]$status==0) / nrow(datasets[[j]])
}

(mean(censoring_rate)) # [1] 0.29357

# borne_a = 20 -- cens_rate = 0.25
# borne_a = 15 -- cens_rate = 0.27

# lambda = 0.5 
datasets <- list()
censoring_rate <- numeric(length = N_files)

borne_a <- 5
for(j in 1:N_files){
  datasets[[j]] <- generate_data(
    lambda=0.5, age_option = "A", n = n_patients, max_time = max_time,
    prop_female = 0,
    year.start_min = 2008, year.start_max = 2010,
    beta_sex = 0, beta_age = beta_age, borne_a = borne_a
  )
  
  censoring_rate[j] <- sum(datasets[[j]]$status==0) / nrow(datasets[[j]])
}

(mean(censoring_rate)) # 0.319595

# borne_a = 10 -- cens_rate = 0.20

# lambda = 1 
datasets <- list()
censoring_rate <- numeric(length = N_files)

borne_a <- 3
for(j in 1:N_files){
  datasets[[j]] <- generate_data(
    lambda=1, age_option = "A", n = n_patients, max_time = max_time,
    prop_female = 0,
    year.start_min = 2008, year.start_max = 2010,
    beta_sex = 0, beta_age = beta_age, borne_a = borne_a
  )
  
  censoring_rate[j] <- sum(datasets[[j]]$status==0) / nrow(datasets[[j]])
}

(mean(censoring_rate)) # 0.29164

# lambda = 1.50 
datasets <- list()
censoring_rate <- numeric(length = N_files)

borne_a <- 2
for(j in 1:N_files){
  datasets[[j]] <- generate_data(
    lambda=1.50, age_option = "A", n = n_patients, max_time = max_time,
    prop_female = 0,
    year.start_min = 2008, year.start_max = 2010,
    beta_sex = 0, beta_age = beta_age, borne_a = borne_a
  )
  
  censoring_rate[j] <- sum(datasets[[j]]$status==0) / nrow(datasets[[j]])
}

(mean(censoring_rate)) # 0.300955

# values for borne_a and censoring rate for each lambda in table
#see grid in amubox
# correct the md file
# another table with these values (possible incl in manuscript)

# --- 4. Timing check (run once on 1 dataset) --- ####
start <- proc.time()

res <- analyze_one(df, lambda, beta_age, times_years = c(1, 2, 3))

elapsed <- proc.time() - start
elapsed
# utilisateur     système      écoulé 
#     0.77        0.05        0.81

