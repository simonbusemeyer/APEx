# ---------------------------------------------------------
# 8. ADVANCED METRICS FOR COMPETENT STUDY DESIGN (CORRECTED)
# ---------------------------------------------------------
# Define the time points (Years 1 to max_time in days)
years_to_check <- 1:params$max_time
yearly_times_days <- years_to_check * 365.241

# 1. ISOLATE CAUSE 1
# cmp.rel returns a list of curves. [[1]] is Disease, [[2]] is Population.
cause1_curve <- crude_mortality[[1]]

# 2. Find the exact indices in the time vector for these years
# Using the specific $time vector inside the cause 1 curve
indices <- findInterval(yearly_times_days, cause1_curve$time)

# Handle potential edge cases where a time is before the very first event
indices[indices == 0] <- 1

# 3. Extract Estimate and Variance directly from the Cause 1 list
# Note: These are 1D vectors, not matrices, so we don't need the "[, 1]" syntax
est_cancer <- cause1_curve$est[indices]
var_cancer <- cause1_curve$var[indices]

# 4. Calculate Absolute Standard Error (SE)
se_cancer <- sqrt(var_cancer)

# 5. Calculate 95% Confidence Intervals (Linear approximation)
lower_ci <- est_cancer - (1.96 * se_cancer)
upper_ci <- est_cancer + (1.96 * se_cancer)

# Bound the CIs logically between 0 and 1
lower_ci <- pmax(0, lower_ci)
upper_ci <- pmin(1, upper_ci)

# Calculate CI Width
ci_width <- upper_ci - lower_ci

# 6. Build the comprehensive evaluation table
eval_table <- data.frame(
  Year        = years_to_check,
  Estimate    = est_cancer,
  Variance    = var_cancer,
  Absolute_SE = se_cancer,
  CI_Width    = ci_width
)

# 7. Calculate Logit RSE
p    <- eval_table$Estimate
se_p <- eval_table$Absolute_SE

# Avoid division by zero if p is exactly 0 or 1
p[p <= 0] <- 1e-6
p[p >= 1] <- 1 - 1e-6

logit_p    <- log(p / (1 - p))
se_logit_p <- se_p / (p * (1 - p))

eval_table$Logit_RSE_Percent <- abs(se_logit_p / logit_p) * 100

# 8. Format the table for readability (rounding)
eval_table_rounded <- data.frame(
  Year          = eval_table$Year,
  Estimate      = round(eval_table$Estimate, 4),
  Absolute_SE   = round(eval_table$Absolute_SE, 4),
  CI_Width      = round(eval_table$CI_Width, 4),
  Logit_RSE_pct = round(eval_table$Logit_RSE_Percent, 2)
)

print(eval_table_rounded)