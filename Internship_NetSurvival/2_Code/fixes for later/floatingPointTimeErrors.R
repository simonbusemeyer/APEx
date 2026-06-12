# ---------------------------------------------------------
# 1. ESTABLISH BIOLOGICAL TRUTH (The Competing Risks)
# ---------------------------------------------------------
# Determine what happened to the patient in a world without censoring
true_time <- pmin(tpsGene, tpsSpe)
true_cause <- ifelse(tpsSpe < tpsGene, "cancer", "other")

# ---------------------------------------------------------
# 2. APPLY OBSERVATIONAL LIMITS (The Censoring Mechanisms)
# ---------------------------------------------------------
# The patient is observed until they die, are randomly censored, or hit administrative max_time
observed_time <- pmin(true_time, tpsCens, max_time)

# ---------------------------------------------------------
# 3. CLASSIFY FINAL OBSERVED STATUS
# ---------------------------------------------------------
# Status is 1 (dead) ONLY if their true death occurred before or exactly when observation ended
status <- as.integer(true_time <= tpsCens & true_time <= max_time)

# Map the final event type based on the strict status
event_type <- ifelse(status == 0, "censored", true_cause)

# Binary cause indicator for the Pohar-Perme estimator (1 = cancer, 0 = other/censored)
cause <- ifelse(event_type == "cancer", 1, 0)