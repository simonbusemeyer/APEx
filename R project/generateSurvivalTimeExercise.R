# Task:
#   Using the approach from Bender et al. (2005), generate survival data with the following specifications:
#   
#   •	N = 20000 patients
#   •	One binary covariate X ~ Bernoulli(0.5)
#   •	Exponential baseline hazard: λ₀ = 0.02
#   •	Effect of X: HR = 2 (so β = log(2) ≈ 0.693)
#   •	Administrative censoring at t = 5 years
  
library(survival)

# 1. parameter configuration
N <- 20000
lambda <- 0.02
beta_X <- log(2)


# 2. covariate generation
X <- rbinom(N, size = 1, prob = 0.5)


# 3. Bender 2005 survival time generation (exponential baseline)
U <- runif(N)
T_val <- -log(U) / (lambda * exp(beta_X * X))
hist(T_val, breaks = 25)


# 4. administrative censoring
cens_time <- 5
observed_time <- pmin(T_val, cens_time)
status <- as.integer(T_val <= cens_time)


# 5. output dataset
sim_df <- data.frame(X, observed_time, status)
head(sim_df)


#   Verification steps:
#     
#     •	Plot Kaplan-Meier curves for X=0 vs X=1 groups
#   •	Fit a Cox model and report the estimated HR (hint: it should be close to 2)
#   •	Report the proportion of events (deaths) vs censored observations


# 1. KM plot by variable X

km_fit <- survfit(Surv(observed_time, status) ~ X, data = sim_df)
plot(km_fit, xlab = "Time (Years)", col = c("blue", "red"), ylab = "Survival Probability")
legend("bottomleft", legend = c("X = 0", "X = 1"), 
       col = c("blue", "red"), lwd = 2, bty = "n")


# 2. cox validation
cox_model <- coxph(Surv(observed_time, status) ~ X, data=sim_df)
summary(cox_model)


# 3. proportion of events (1) vs censored obs (0)
print(round(prop.table(table(sim_df$status)), 2))


