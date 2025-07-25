library(pracma)
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

run_fda_lcs_simulation <- function(i = 1, n_subj = 100, n_time_points = 2,
                                   data_gen = "homogeneous",
                                   age = 0:19, tfine = seq(0, 19, length.out = 240), method = "normal", rep_cond = 4,
                                   cores = 24) {
  set.seed(i)
  # Simulate data
  if (data_gen == "homogeneous") {
    sim_data <- simulate_exponential_decay(n_subjects = n_subj, x = tfine, noise = 0.1)
  } else if (data_gen == "heterogeneous") {
    sim_data <- simulate_exponential_decay(n_subjects = n_subj, x = tfine, noise = 0.1, sd_L = 5, sd_a = 5)
  } else if (data_gen == "functional_heterogeneous") {
    sim_data <- simulate_exponential_decay(n_subjects = n_subj, x = tfine, noise = 0.1, subgrp = 0.5)
  } else {
    stop("Unknown data_gen type: must be 'homogeneous', 'heterogeneous', or 'functional_heterogeneous'")
  }
  
  sim_al_data <- simulate_accelerated_longitudinal(Y_full = sim_data$trajectories, n_time_points = n_time_points)
  
  # FDA model
  z02 <- generate_initial_values(
    Lt_list = sim_al_data$Lt,
    Ly_list = sim_al_data$Ly,
    age_grid = age,
    m = rep_cond,
    method = method
  )
  
  acd_fda_2 <- ldm(sim_al_data$Ly, sim_al_data$Lt, z02, age,
                   optns = list(M = nrow(z02), cores = cores, regular = FALSE))
  
  results <- evaluate_simulations(tfine = tfine, sim_al_data = sim_al_data, sim_data = sim_data, acd_fda_2 = acd_fda_2, plot = FALSE)
  
  # LCS-CT-SSM
  data <- convert_to_full_subject_list(sim_al_data)
  LCS_SSM <- run_full_kalman_simulation(data, sim_data)
  
  results2 <- suppressWarnings(evaluate_lcsssm(LCS_SSM, sim_data, sim_al_data, plot = FALSE))
  
  # Comparison
  metrics_df <- suppressWarnings(compare_fda_lcs(results, results2, LCS_SSM, plot = FALSE))
  
  list(
    Metrics = metrics_df$metrics,
    DB_FDA = metrics_df$DB_FDA,
    DB_LCS = metrics_df$DB_LCS
  )
}

n_subj_vec <- c(100, 250, 500)
n_time_points_vec <- c(2, 3)
data_gen_vec <- c("homogeneous", "heterogeneous", "functional_heterogeneous")
iters <- 1:10

# Create a results tibble
result_tbl2 <- tibble(iter = iters)

# Loop over all parameter combinations
for (subject in n_subj_vec) {
  for (time_point in n_time_points_vec) {
    for(data_method in data_gen_vec){
      
      # Create column name for results
      col_name <- paste0(data_method, "_n", subject, "_tp", time_point)
      
      # Print start of this condition
      message("Starting simulations for data method = ", data_method, ", n = ", subject, ", tp = ", time_point)
      
      # Run run_fda_lcs_simulation for all iterations with progress message
      result_tbl2[[col_name]] <- map(iters, function(i) {
        message("  Iteration ", i, " for ", data_method, "_n", subject, "_tp", time_point)
        
        # Wrap run_fda_lcs_simulation in tryCatch to handle errors
        tryCatch(
          run_fda_lcs_simulation(i = i, n_subj = subject, n_time_points = time_point, data_gen = data_method, cores = 24),
          error = function(e) {
            message("    Error at iter ", i, ": ", conditionMessage(e))
            return(NA)
          }
        )
      })
    }
  }
}
   
# Save results
saveRDS(result_tbl2, file = "~/FDA_with_accelerated_longitudinal_designs/result_tbl2.rds")

end <- Sys.time()
end - begin
