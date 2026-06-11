# =============================================================================
# utils.R
# Fonctions utilitaires partagées
# =============================================================================

library(survival)
library(relsurv)

# -----------------------------------------------------------------------------
# Check : résumé d'un dataset
# -----------------------------------------------------------------------------
check_data <- function(data) {
  message("Dimensions : ", nrow(data), " lignes x ", ncol(data), " colonnes")
  message("Taux de censure : ",
          round(100 * mean(data$status == 0), 1), "%")
  message("Suivi médian : ", round(median(data$observed_time), 2))
  print(summary(data))
  invisible(data)
}

# -----------------------------------------------------------------------------
# Courbe de survie théorique (Exponentielle)
# -----------------------------------------------------------------------------
survtheo <- function(t, lambda, beta_sex, beta_age, beta_X,
                     sex, ageStand, X) {
  # S(t) = exp(-lambda * t * exp(beta * Z))
  exp(-lambda * t * exp(beta_sex * sex + beta_age * ageStand + beta_X * X))
}

# -----------------------------------------------------------------------------
# Plot 1 : KM par dataset + courbe théorique moyennée
# -----------------------------------------------------------------------------
plot_km_theo <- function(datasets, scenario, max_time) {
  
  t_seq <- seq(0, max_time, by = 0.1)
  
  surv_theo_mean <- rowMeans(
    sapply(datasets, function(df) {
      ind_curves <- sapply(1:nrow(df), function(i) {
        survtheo(
          t         = t_seq,
          lambda    = scenario$lambda,
          beta_sex  = scenario$beta_sex,
          beta_age  = scenario$beta_age,
          beta_X    = scenario$beta_X,
          sex       = df$sex_num[i],
          ageStand  = df$ageStand[i],
          X         = df$race_num[i]
        )
      })
      rowMeans(ind_curves)
    })
  )
  
  plot(0, type = "n", ylab = "Survival probability", xlab = "Time (years)",
       xlim = c(0, max_time), ylim = c(0.7, 1),
       main = "KM_hyp_world par dataset + courbe SN théorique moyenne")
  grid()
  
  for (data in datasets) {
    km <- survfit(Surv(hypothetical_time, hypothetical_status) ~ 1, data = data)
    lines(km$time, km$surv, col = rgb(0.7, 0.7, 0.7, 0.4))
  }
  
  lines(t_seq, surv_theo_mean, col = "blue", lwd = 2)
  legend("topright", legend = c("KM_hyp_world", "SN_théo moyenne"),
         col = c(rgb(0.7, 0.7, 0.7, 0.4), "blue"), lwd = c(1, 2))
}

# -----------------------------------------------------------------------------
# Check : analyse des événements sur les N datasets
# -----------------------------------------------------------------------------
check_events <- function(datasets) {
  
  # -----------------------------------------------------------------------------
  # Construire la cause de décès pour chaque dataset
  # -----------------------------------------------------------------------------
  get_cause <- function(data) {
    data$cause <- ifelse(data$status == 0, "censure",
                         ifelse(abs(data$observed_time - data$tpsSpe) < 1e-10, "cancer",
                                "autres"))
    return(data)
  }
  
  # -----------------------------------------------------------------------------
  # Proportions par dataset
  # -----------------------------------------------------------------------------
  props <- sapply(datasets, function(data) {
    data  <- get_cause(data)
    n_dc  <- sum(data$status == 1)
    c(
      prop_cancer  = sum(data$cause == "cancer")  / n_dc,
      prop_autres  = sum(data$cause == "autres")  / n_dc,
      n_cancer     = sum(data$cause == "cancer"),
      n_autres     = sum(data$cause == "autres"),
      n_censure    = sum(data$cause == "censure"),
      n_total_dc   = n_dc
    )
  })
  
  # -----------------------------------------------------------------------------
  # Boxplots proportions
  # -----------------------------------------------------------------------------
  par(mfrow = c(1, 2))
  
  boxplot(props["prop_cancer", ],
          main = "Proportion DC cancer / total DC",
          ylab = "Proportion",
          col  = "lightcoral")
  grid()
  
  boxplot(props["prop_autres", ],
          main = "Proportion DC autres causes / total DC",
          ylab = "Proportion",
          col  = "lightblue")
  grid()
  
  # -----------------------------------------------------------------------------
  # Histogrammes sur données poolées
  # -----------------------------------------------------------------------------
  all_data <- do.call(rbind, lapply(datasets, get_cause))
  
  par(mfrow = c(1, 2))
  
  hist(all_data$tpsSpe[all_data$cause == "cancer"],
       main  = "Temps DC cancer (poolé)",
       xlab  = "Temps (années)",
       col   = "lightcoral",
       border = "white")
  grid()
  
  hist(all_data$tpsGene[all_data$cause == "autres"],
       main  = "Temps DC autres causes (poolé)",
       xlab  = "Temps (années)",
       col   = "lightblue",
       border = "white")
  grid()
  
  par(mfrow = c(1, 1))
  
  # -----------------------------------------------------------------------------
  # Résumé numérique
  # -----------------------------------------------------------------------------
  message("--- Résumé événements ---")
  message("DC cancer    — moyenne : ", round(mean(props["n_cancer", ]), 1),
          " (", round(mean(props["prop_cancer", ]) * 100, 1), "%)")
  message("DC autres    — moyenne : ", round(mean(props["n_autres", ]), 1),
          " (", round(mean(props["prop_autres", ]) * 100, 1), "%)")
  message("Censures     — moyenne : ", round(mean(props["n_censure", ]), 1))
  message("Total DC     — moyenne : ", round(mean(props["n_total_dc", ]), 1))
  
  invisible(list(props = props, data_pooled = all_data))
}