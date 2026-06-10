# Code Structure

---

## File organization

```
simulation_study/
├── main.R
├── functions/
│   ├── generate_data.R
│   ├── analyze_one.R
│   └── compute_metrics.R
└── outputs/
    ├── plots/
    └── tables/
```

---

## Functions — what goes where

### `generate_data.R`
→ **Your existing generation code**, wrapped in a function.
Inputs: n, lambda, beta_age, borne_a, max_time_days
Output: 1 dataframe 

**Do not forget:** `event_type` that is needed to compute % cancer deaths

```r
event_type <- ifelse(obs_time == C,   "censored",
              ifelse(obs_time == T_E, "cancer", "other"))
```

---

### `analyze_one.R`
→ **Your existing PP call**, wrapped in a function.

Inputs: df (1 dataset), lambda, beta_age, times_years

Output: dataframe with columns:
  time, s_pp, s_lower, s_upper, s_theo, diff, covered, pct_cancer

New things to put here vs your existing code:

1. Theoretical net survival:


2. Extract CI from rs.surv:

```r
pp_summary <- summary(pp_fit, times = c(1,2,3) * 365.241)
# pp_summary$surv / pp_summary$lower / pp_summary$upper
```

3. Return all columns including:

```r
diff       = s_pp - s_theo
covered    = (s_lower <= s_theo) & (s_theo <= s_upper)
pct_cancer = mean(df$event_type == "cancer")
```

Note: `diff` and `covered` are computed **per simulation**.

`compute_metrics` then aggregates them over 1000 simulations:
- mean(diff)         → Bias
- sqrt(mean(diff²))  → RMSE
- mean(covered)      → ECR

---

### `compute_metrics.R`
→ **New function** — takes results from 1000 simulations, returns metrics table.
Input: list of J=1000 dataframes (each from analyze_one)
Output: dataframe with bias, RMSE, ECR per time point

```r
compute_metrics <- function(results_list, lambda) {
  all_res <- do.call(rbind, results_list)
  # group by time, compute mean(diff), sqrt(mean(diff^2)), mean(covered)
  # add lambda and pct_cancer columns
}
```

---

## `main.R` — skeleton only

```r
source("functions/generate_data.R")
source("functions/analyze_one.R")
source("functions/compute_metrics.R")

# --- 0. Parameters ---
# n, n_sim, beta_age, times_years, etc.

# --- 1. nessie (run once, sets max_time_years) ---

# --- 2. Timing check (run once on 1 dataset) ---
# Wrap your existing generate + analyze calls like this:
#
# start   <- proc.time()
# df      <- generate_data(...)   # ← your call
# res     <- analyze_one(df, ...) # ← your call
# elapsed <- proc.time() - start
# cat("Time for 1 dataset:", round(elapsed["elapsed"], 2), "seconds\n")
#
# If > 10 seconds: stop and send me an email pls

# --- 3. Censoring calibration (run once per lambda) ---
# For each lambda value:
#   Start with borne_a = max_follow_up_days
#   Repeat:
#     Generate 100 datasets with current borne_a
#     Compute censoring_rate = mean % censored across 100 datasets
#     If censoring_rate > 32% : borne_a is too small → multiply by 1.1
#     If censoring_rate < 28% : borne_a is too large → multiply by 0.9
#     If between 28% and 32%  : stop, save this borne_a
#   Use the saved borne_a for all 1000 simulations of this scenario
# Reminder: large borne_a → C can be large → less censoring

# --- 4. Verification plots (run before full simulation) ---
# For each scenario: 1 dataset, plot PP vs S_theo
# Check: % cancer, % other, % censored sum to 100%
# Check: sex = male only, age in [80,89]
# Save to outputs/plots/

# --- 5. Full simulation loop ---
for (k in scenarios) {
  for (j in 1:n_sim) {
    df            <- generate_data(...)
    results[[j]]  <- analyze_one(df, ...)
  }
  metrics[[k]] <- compute_metrics(results, lambda[k])
}

# --- 6. Save results ---
# write.csv → outputs/tables/metrics_week4.csv
```

---

## Note on verification plots (section 4)

Before launching 1000 simulations, for each scenario plot on 1 dataset:
- PP curve (blue)
- S_theo curve (red dashed)

If PP curve goes UP → follow-up too long → go back to nessie output.
If PP and S_theo far apart → check event_type proportions first.