# =============================================================================
# batch.R
# Fonction batch : appelle generate_data
# selon différents scénarios
# =============================================================================

library(here)
source(here::here("R/generate_data.R"))

# -----------------------------------------------------------------------------
# Scénarios : liste de paramètres
# -----------------------------------------------------------------------------
scenarios <- list(
  scenario_1 = list(
    k = 2,
    pho = 0.5,
    alpha = 0.2,
    n = 2000,
    max_time = 5,
    prop_female = 0.5, prop_x0 = 0.5,
    year.start_min = 2000, year.start_max = 2002,
    beta_sex = 1.2, beta_age = 0.3, beta_X = -0.2,
    borne_a = 18
  ),
  scenario_2 = list(
    k = 2,
    pho = 0.1,
    alpha = 0.2,
    n = 2000,
    max_time = 5,
    prop_female = 0.5, prop_x0 = 0.5,
    year.start_min = 2000, year.start_max = 2002,
    beta_sex = 1.2, beta_age = 0.3, beta_X = -0.2,
    borne_a = 50
  )
)

# -----------------------------------------------------------------------------
# Fonction batch
# -----------------------------------------------------------------------------

run_batch <- function(scenario, generate_fct = generate_data, N = 100, seed = 123) {
  
  set.seed(seed)
  
  datasets <- lapply(1:N, function(sc_name) {
    data <- do.call(generate_fct, scenario)
    data$sim_id <- 1
    return(data)
  })
  
  names(datasets) <- paste0("sim_", 1:N)
  return(datasets)
}

save_batch <- function(datasets, scenario_name = "scenario_1", 
                       path = "data/processed/") {
  filepath <- file.path(path, paste0(scenario_name, ".rds"))
  saveRDS(datasets, file = filepath)
  message("Sauvegardé : ", filepath)
}