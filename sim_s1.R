library(pracma)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(parallel)
library(fda)
begin = Sys.time()
source("~/FDA_with_accelerated_longitudinal_designs/ldm.R")
source("~/FDA_with_accelerated_longitudinal_designs/kerFctn.R")
source("~/FDA_with_accelerated_longitudinal_designs/Rfunctions.R")

sim_S1 <- function(method, rep_cond, iter) {
  set.seed(iter)
  results <- list()
  n_subj <- 250
  tfine <- seq(0, 19, length.out = 240)
  
  # Simulate true data and sparse sampling
  sim_data <- simulate_exponential_decay(n_subjects = n_subj, x = tfine, noise = 0.1)
  sim_al_data <- simulate_accelerated_longitudinal(Y_full = sim_data$trajectories, n_time_points = 3)
  
  age <- 0:19
  
  # Generate initial values depending on method
  if (method == "normal") {
    z02 <- generate_initial_values(
      Lt_list = sim_al_data$Lt,
      Ly_list = sim_al_data$Ly,
      age_grid = age,
      m = rep_cond,
      method = "normal"
    )
  } else if (method == "replicate") {
    z02 <- generate_initial_values(
      Lt_list = sim_data$Lt_true,
      Ly_list = sim_data$Ly_true,
      age_grid = age,
      m = rep_cond,
      method = "replicate"
    )
  } else {
    stop("Invalid method. Use 'normal' or 'replicate'.")
  }
  
  # Fit models
  acd_fda_2 <- ldm(sim_al_data$Ly, sim_al_data$Lt, z02, age,
                   optns = list(
                     M = nrow(z02),
                     regular = FALSE, cores = 24
                   )
  )
  
  # Evaluate final model
  results <- evaluate_simulations(tfine = tfine, sim_al_data = sim_al_data, sim_data = sim_data, acd_fda_2 = acd_fda_2)
  rm(sim_data, sim_al_data, acd_fda_2)
  
  # Return results
  return(list(
    ie_opt = results$ie_opt,
    ie_all = results$ie_all,
    ise_opt = results$ise_opt,
    ise_all = results$ise_all,
    mse_opt = results$mse_opt,
    mse_all = results$mse_all,
    bias_opt = results$bias_opt,
    bias_all = results$bias_all
  ))
}

methods <- c("normal", "replicate")
rep_conds <- c(4, 10, 100)
iters <- 1:50

# Initialize result tibble
result_tbl <- tibble(iter = iters)

# Loop through method and repetition conditions
for (method in methods) {
  for (m in rep_conds) {
    
    # Create column name for results
    col_name <- paste0(method, "_m", m)
    
    # Print start of this condition
    message("Starting simulations for method = ", method, ", m = ", m)
    
    # Run sim_S1 for all iterations with progress message
    result_tbl[[col_name]] <- map(iters, function(i) {
      message("  Iteration ", i, " for ", method, "_m", m)
      
      # Wrap sim_S1 in tryCatch to handle errors
      tryCatch(
        sim_S1(method = method, rep_cond = m, iter = i),
        error = function(e) {
          message("    Error at iter ", i, ": ", conditionMessage(e))
          return(NA)
        }
      )
    })
  }
}

saveRDS(result_tbl, file = "~/FDA_with_accelerated_longitudinal_designs/result_tbl2")
result_tbl2 <- readRDS(file = "~/FDA_with_accelerated_longitudinal_designs/result_tbl2")

end = Sys.time()
end - begin
