# Week 4 — Instructions

---

## Fixed parameters

| Parameter         | Value                         |
|-------------------|-------------------------------|
| n per dataset     | 2000                          |
| n simulations     | 1000                          |
| Age               | Uniform(80, 89)               |
| Sex               | **Males only**                    |
| Ratetable         | `survexp.us` — males only       |
| Excess hazard     | Exponential(λ)                |
| β\_age             | [VALUE]                       |
| Follow-up         | determined by nessie (Task 1) |
| Censoring rate    | 30%                           |
| Times for metrics | 1, 2, 3 years                 |

---

## Task 1 — Follow-up duration via `nessie`

Run `nessie` once on a pilot dataset (covariates only — no event times needed).

**Rule:** set max follow-up = floor( min(`c.exp.surv`) across age strata )

**Why:** Pohar-Perme weights = 1/S*(t). When t > c.exp.surv for the oldest
patients, S*(t) → 0, weights → infinity → the net survival curve goes UP.
You have seen this in your plots. This is why follow-up must not exceed
c.exp.surv for any stratum.

Report: full nessie output + chosen follow-up with justification.

---

## Task 2 — Censoring calibration

C ~ Uniform(0, borne_a). borne_a must be calibrated **separately for each λ**.

```
For each scenario (each λ):

  Initialize: borne_a = max_follow_up_days

  Repeat:
    1. Generate 100 datasets using current borne_a
    2. censoring_rate = mean % censored across 100 datasets

    3. If censoring_rate > 32%  →  borne_a too small  →  borne_a = borne_a × 1.1
       If censoring_rate < 28%  →  borne_a too large  →  borne_a = borne_a × 0.9
       If 28% ≤ censoring_rate ≤ 32%  →  STOP

  Save this borne_a → use it for all 1000 simulations of this scenario
```

Note: large borne_a = C can be large = less censoring. Small borne_a = more censoring.

Report: table with λ, calibrated borne_a, observed censoring rate.

---

## Task 3 — Time your code

```r
start   <- proc.time()
# one full generate + analyze call on 1 dataset
elapsed <- proc.time() - start
cat("Time for 1 dataset:", round(elapsed["elapsed"], 2), "seconds\n")
```

If > 10 seconds: stop and send me an email before launching 1000 simulations.

---

## Task 4 — Lambda grid

Propose 6 values of λ. For each, run 100 pilot datasets and report
observed % cancer deaths before launching the full 1000 simulations.

Target range: ~10% to ~80% cancer deaths across the 6 scenarios.

| Scenario | λ | % cancer (pilot, 100 datasets) | borne_a |
|----------|---|--------------------------------|---------|
| S1       |   |                                |         |
| S2       |   |                                |         |
| S3       |   |                                |         |
| S4       |   |                                |         |
| S5       |   |                                |         |
| S6       |   |                                |         |

---

## Performance metrics

For J = 1000 simulations, at each time t:

- S_theo(t) for simulation $j = (1/n) Σ exp(−λ · t · exp(β\_{age} · age_i))$ computed on the n=2000 patients of that simulation
- [L_j(t), U_j(t)] = 95% CI from rs.surv ` (pp_summary$lower, pp_summary$upper)`

Bias(t)  = (1/J) Σ\_j [ S_PP_j(t) − S_theo_j(t) ]

RMSE(t)  = sqrt( (1/J) Σ\_j [ S_PP_j(t) − S_theo_j(t) ]² )

ECR(t)   = (1/J) Σ\_j 1[ L_j(t) ≤ S_theo_j(t) ≤ U_j(t) ]
→ expected value under correct specification: 0.95

---

## Expected results table

| lambda | pct_cancer | time | bias | rmse | ecr |
|--------|------------|------|------|------|-----|
| ...    | ...        | 1    | ...  | ...  | ... |
| ...    | ...        | 2    | ...  | ...  | ... |
| ...    | ...        | 3    | ...  | ...  | ... |

Save as .csv in outputs/tables/

---

## End of day

Upload to shared AmuBox by 5:30pm:

- Reorganized script
- nessie output + chosen follow-up
- Timing result
- Calibrated borne_a per scenario

---

## Week schedule

| Day       | Task                                                                                       |
|-----------|--------------------------------------------------------------------------------------------|
| Monday    | Tasks 1–3 only (nessie + censoring + timing). Pls share code by 5:30pm                     |
| Tuesday   | Task 4 (lambda grid + pilots) after my feedback. Hoping full simulations if code validated |
| Wednesday | Results table + interpret metrics + start slides                                           |
| Friday am | Slides review                                                                              |

## Tuesday instructions

The calibration is done. You will find on the shared cloud:

- the parameters table (lambda, borne_a, prop_cancer, cens_rate)
- corrected generate_data, main and analyze_one scripts

### Your tasks today (before and after the meeting)

1. **Timing check** — measure on YOUR computer:

```r
start <- proc.time()

res <- analyze_one(df, lambda, beta_age, times_years = c(1, 2, 3))

elapsed <- proc.time() - start
```

Report the result before going further.

1. **Implement compute_metrics**: using the formulas in this document or in your written notes for 

   \- bias

   \- RMSE  
   - ECR   
     
   at t = 1, 2, 3 years across J simulations
2. **Write the simulation loop** for one scenario (i.e., one value of lambda):

```r
results <- vector("list", N_files)
for (j in 1:N_files) {
  df         <- generate_data(lambda = ..., borne_a = ...)
  results[[j]] <- analyze_one(df, ...)
}
metrics <- compute_metrics(results, lambda = ...)
```

1. **Run all 7 scenarios** sequentially using the parameters table. Save results to outputs/tables/metrics.csv

Upload your script to the shared cloud by tomorrow. We'll try to look at them at 2 pm.