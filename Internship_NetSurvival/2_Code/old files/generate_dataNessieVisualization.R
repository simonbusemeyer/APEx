library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)

plot_nessie <- function(nessie_obj) {
  # 1. Extract the matrix from the nessie object
  mat <- nessie_obj$mata
  
  # 2. Convert to a dataframe and pull the row names into a 'Strata' column
  df <- as.data.frame(mat) %>%
    rownames_to_column(var = "Strata")
  
  # 3. Drop the 'c.exp.surv' column (it messes up the time series line)
  df_plot <- df %>% select(-c.exp.surv)
  
  # 4. Pivot from Wide to Long format
  df_long <- df_plot %>%
    pivot_longer(
      cols = -Strata, 
      names_to = "Time", 
      values_to = "Expected_N"
    ) %>%
    mutate(Time = as.numeric(Time)) # Ensure X-axis is treated as continuous numbers
  
  # 5. Build a dynamic ggplot
  p <- ggplot(df_long, aes(x = Time, y = Expected_N, color = Strata, group = Strata)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2.5) +
    theme_minimal(base_size = 14) +
    labs(
      title = "Expected Net Sample Size Over Time",
      subtitle = "Projected remaining patients if only population mortality applies",
      x = "Follow-up Time (Years)",
      y = "Expected Number of Patients",
      color = "Patient Subgroup"
    ) +
    theme(
      legend.position = "bottom",
      legend.direction = "vertical",
      panel.grid.minor = element_blank() # Cleans up the background
    )
  
  return(p)
}

# Execute the function on your saved output
plot_nessie(nessie_output)