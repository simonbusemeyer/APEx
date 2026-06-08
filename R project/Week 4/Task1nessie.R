library(survival)
library(relsurv)

set.seed(12345)
n <- 2000
age_days <- runif(n, 80, 89) * 365.241
sex <- rep("male", n)
diagnosis_year <- as.Date("2008-01-01")

time <- rep(100, n) 
status <- rep(0, n)

pilot_df <- data.frame(
  age = age_days, 
  sex = sex, 
  year = diagnosis_year, 
  time = time, 
  status = status
)

pilot_df$age_strata <- cut(pilot_df$age / 365.241, breaks = c(80, 81, 82, 83, 84, 85, 86, 87, 88, 89))

nessie_out <- nessie(
  Surv(time, status) ~ age_strata,
  data = pilot_df,
  ratetable = survexp.us,
  rmap = list(age = age, sex = sex, year = year)
)

