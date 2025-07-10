simulate_exponential_decay <- function(n_subjects = 500, 
                                       x = seq(0, 19, length.out = 1000),
                                       mu_L = 50, sd_L = 1,
                                       mu_a = 50, sd_a = 1,
                                       mu_k = 0.5, sd_k = 0.1,
                                       noise = 0.1) {
  # Bounds for asymptote check
  L_lower <- 0.5 * mu_L
  L_upper <- 1.5 * mu_L
  
  # Initialize storage
  Y <- matrix(NA, nrow = length(x), ncol = n_subjects)
  Y_no_noise <- matrix(NA, nrow = length(x), ncol = n_subjects)
  a_vals <- numeric(n_subjects)
  L_vals <- numeric(n_subjects)
  k_vals <- numeric(n_subjects)
  
  for (i in seq_len(n_subjects)) {
    repeat {
      # Sample parameters
      a_i <- rnorm(1, mean = mu_a, sd = sd_a)
      L_i <- rnorm(1, mean = mu_L, sd = sd_L)
      k_i <- rnorm(1, mean = mu_k, sd = sd_k)
      
      # Simulate trajectory
      traj <- L_i - a_i * exp(-abs(k_i) * x)
      asymptote_val <- traj[length(x)]
      
      # Check the asymptote
      if (asymptote_val >= L_lower && asymptote_val <= L_upper) {
        Y_no_noise[, i] <- traj
        Y[, i] <- traj + rnorm(length(x), sd = noise)
        a_vals[i] <- a_i
        L_vals[i] <- L_i
        k_vals[i] <- k_i
        break
      }
    }
  }
  
  # Return as list
  Ly_true <- lapply(seq_len(n_subjects), function(i) Y_no_noise[, i])
  Lt_true <- lapply(seq_len(n_subjects), function(i) x)
  
  return(list(
    trajectories = Y,
    x = x,
    a = a_vals,
    k = k_vals,
    L = L_vals,
    Y_no_noise = Y_no_noise,
    Lt_true = Lt_true,
    Ly_true = Ly_true
  ))
}


simulate_accelerated_longitudinal <- function(Y_full, 
                                              time_range = c(0, 19),
                                              n_time_points = 3,
                                              time_grid_length = 240,
                                              x = 0.05) {
  n_subjects <- ncol(Y_full)
  
  # Step 1: Subsample to monthly measurements
  original_grid <- seq(time_range[1], time_range[2], length.out = nrow(Y_full))
  monthly_grid <- seq(time_range[1], time_range[2], length.out = time_grid_length)
  
  # Find closest points on original grid for monthly grid
  subsample_idx <- sapply(monthly_grid, function(t) which.min(abs(original_grid - t)))
  Y_full_subsampled <- Y_full[subsample_idx, ]
  x_full <- monthly_grid
  
  # Step 2: Precompute subject IDs for fixed start
  n_fixed <- round(x * n_subjects)
  fixed_start_subjects <- sample(seq_len(n_subjects), n_fixed)
  
  Lt <- vector("list", n_subjects)
  Ly <- vector("list", n_subjects)
  
  for (i in seq_len(n_subjects)) {
    sampled_times <- numeric(n_time_points)
    
    if (i %in% fixed_start_subjects) {
      # Force first time point to be time_range[1]
      sampled_times[1] <- time_range[1]
    } else {
      # Random first time point
      max_start <- time_range[2] - (n_time_points - 1) * 0.92
      sampled_times[1] <- runif(1, min = time_range[1], max = max_start)
    }
    
    for (j in 2:n_time_points) {
      min_time <- sampled_times[j - 1] + 0.9
      max_time <- sampled_times[j - 1] + 2
      if (min_time >= time_range[2]) {
        sampled_times <- sampled_times[1:(j - 1)]
        break
      }
      sampled_times[j] <- runif(1, min = min_time, max = min(max_time, time_range[2]))
    }
    
    # Match to closest x_full grid points
    closest_idx <- sapply(sampled_times, function(t) which.min(abs(x_full - t)))
    Lt[[i]] <- x_full[closest_idx]
    Ly[[i]] <- Y_full_subsampled[closest_idx, i]
  }
  
  return(list(Lt = Lt, Ly = Ly, x_full = x_full, Y_full = Y_full_subsampled))
}

generate_initial_values <- function(Lt_list, Ly_list, age_grid, m = 1, method = c("normal", "replicate"), eps = 1e-6) {
  method <- match.arg(method)
  n_subj <- length(Ly_list)
  
  # Unlist Lt and Ly to extract values at time 0
  Lt_unlisted <- unlist(Lt_list)
  Ly_unlisted <- unlist(Ly_list)
  
  # Values where time is (almost) zero
  Ly_at_t0 <- Ly_unlisted[abs(Lt_unlisted - 0) < eps]
  
  # If no exact or near-zero match, find closest one
  if (length(Ly_at_t0) == 0) {
    closest_index <- which.min(abs(Lt_unlisted - 0))
    Ly_at_t0 <- Ly_unlisted[closest_index]
  }
  
  if (method == "normal") {
    mu <- mean(Ly_at_t0)
    sigma2 <- var(Ly_at_t0)
    z_vals <- rnorm(n_subj * m, mean = mu, sd = sqrt(sigma2))
  } else if (method == "replicate") {
    # Repeat the vector or single value to get n_subj * m values
    z_vals <- rep(Ly_at_t0, length.out = n_subj * m)
  }
  
  # Final z0 matrix: columns = (initial value, age)
  z0 <- cbind(z_vals, rep(age_grid[1], n_subj * m))
  return(z0)
}

compute_incremental_avg_mse <- function(sim_vals, obs_vals) {
  sim_vals <- sim_vals[, colSums(is.na(sim_vals)) == 0, drop = FALSE]
  # Compute individual MSEs
  mse <- colMeans((sim_vals - obs_vals)^2)

  # Sort simulations by MSE
  sorted_idx <- order(mse)
  sorted_sim_vals <- sim_vals[, sorted_idx]
  
  n <- ncol(sorted_sim_vals)
  cumulative_avg_mse <- numeric(n)
  
  # Iteratively average top-N and compute MSE to obs_vals
  cumulative_sum <- rep(0, length(obs_vals))
  for (i in 1:n) {
    cumulative_sum <- cumulative_sum + sorted_sim_vals[, i]
    avg_curve <- cumulative_sum / i
    cumulative_avg_mse[i] <- mean((avg_curve - obs_vals)^2)
  }
  
  # Define x and y
  x <- 1:n
  y <- cumulative_avg_mse
  
  # Line connecting first and last point
  point1 <- c(x[1], y[1])
  point2 <- c(x[length(x)], y[length(y)])
  
  # Function to compute perpendicular distance from a point to a line
  perpendicular_dist <- function(x0, y0, x1, y1, x2, y2) {
    num <- abs((y2 - y1)*x0 - (x2 - x1)*y0 + x2*y1 - y2*x1)
    denom <- sqrt((y2 - y1)^2 + (x2 - x1)^2)
    return(num / denom)
  }
  
  # Compute distances of all points to the line
  distances <- mapply(function(x0, y0) {
    perpendicular_dist(x0, y0, point1[1], point1[2], point2[1], point2[2])
  }, x, y)
  
  elbow_index <- which.max(distances)
  
  return(list(
    cumulative_avg_mse = cumulative_avg_mse,
    elbow_index = elbow_index
  ))
}

plot_metric <- function(metric_type, df = metrics_long) {
  df_sub <- df %>% filter(Metric_Type == metric_type)
  
  # Calculate outliers for each group
  df_outliers <- df_sub %>%
    group_by(Version) %>%
    mutate(
      Q1 = quantile(Value, 0.25, na.rm = TRUE),
      Q3 = quantile(Value, 0.75, na.rm = TRUE),
      IQR_val = Q3 - Q1,
      upper_whisker = Q3 + 1.5 * IQR_val,
      lower_whisker = Q1 - 1.5 * IQR_val,
      is_outlier = Value > upper_whisker | Value < lower_whisker
    )
  
  # Plot
  ggplot(df_outliers, aes(x = Version, y = Value)) +
    geom_boxplot(outlier.shape = NA, fill = "lightgray") +
    geom_jitter(aes(color = is_outlier), width = 0.1, size = 2, alpha = 0.6) +
    geom_text(data = filter(df_outliers, is_outlier),
              aes(label = subject),
              hjust = -0.3, color = "red", size = 3.5) +
    scale_color_manual(values = c("black", "red"), guide = "none") +
    labs(title = paste(metric_type, ": Opt vs All Simulated Trajectories"),
         y = metric_type, x = "") +
    theme_minimal()
}

evaluate_simulations <- function(tfine, sim_al_data, sim_data, acd_fda_2, plot = FALSE, plot_subjects = NULL) {

  # --- Step 1: Setup grid and basis ---
  n_obs <- tfine[length(tfine)]
  Y_estimated <- acd_fda_2$path   # dim: subjects × time points
  time_est <- seq(0, n_obs, length.out = ncol(Y_estimated))  # 1 × T
  
  # --- Step 2: Smooth each subject using smooth.spline ---
  n_subjects <- nrow(Y_estimated)
  fd_eval_mat <- matrix(NA, nrow = length(tfine), ncol = n_subjects)
  
  error_count <- 0  # Counter for errors
  
  for (i in 1:n_subjects) {
    y <- Y_estimated[i, ]
    
    tryCatch({
      smooth_fit <- smooth.spline(x = time_est, y = y, cv = FALSE)
      fd_eval_mat[, i] <- predict(smooth_fit, x = tfine)$y
    }, error = function(e) {
      error_count <<- error_count + 1  # Increment error counter
    })
  }
  
  # Final message
  if (error_count > 0) {
    message("  Skipped ", error_count, " simulated trajectories due to smooth.spline errors.")
  }
  
  # Calculate mean of the simulated trajectories
  mean_fd <- rowMeans(fd_eval_mat, na.rm = TRUE)
  
  # Step 3: Match and compute Ly_sim and Ly_sim_avg
  age_grid <- tfine
  n_sim <- ncol(fd_eval_mat)
  n_sub <- length(sim_al_data$Ly)
  
  Ly_sim <- vector("list", n_sub)
  Ly_sim_avg <- vector("list", n_sub)
  mse_vec <- numeric(n_sub)
  bias_vec <- numeric(n_sub)
  mse_all <- numeric(n_sub)
  bias_all <- numeric(n_sub)
  
  for (i in seq_len(n_sub)) {
    obs_time <- sim_al_data$Lt[[i]]
    obs_vals <- sim_al_data$Ly[[i]]
    pos <- match(round(obs_time, 3), round(age_grid, 3))
    
    sim_vals <- fd_eval_mat[pos, ]
    sim_mean_vals <- mean_fd[pos]
    
    # Compute subject-wise MSEs across all sim trajectories
    mse <- colMeans((sim_vals - obs_vals)^2)
    sorted_idx <- order(mse)
    
    # Get incremental MSE and elbow index
    inc_mse <- compute_incremental_avg_mse(sim_vals, obs_vals)
    elbow_index <- inc_mse$elbow_index
    top_idx <- sorted_idx[1:elbow_index]
    
    # Store simulations and average
    Ly_sim[[i]] <- fd_eval_mat[, top_idx]
    Ly_sim_avg[[i]] <- rowMeans(Ly_sim[[i]])
    
    # Store optimized MSE and bias
    mse_vec[i] <- inc_mse$cumulative_avg_mse[elbow_index]
    bias_vec[i] <- mean(obs_vals - rowMeans(sim_vals[, top_idx, drop = FALSE]))
    
    # Store MSE and bias from overall mean_fd
    mse_all[i] <- mean((sim_mean_vals - obs_vals)^2)
    bias_all[i] <- mean(obs_vals - sim_mean_vals)
  }
  
  # --- Step 4a: ISE computation ---
  calc_ise <- function(true_phi, est_phi) {
    if (is.vector(true_phi)) true_phi <- matrix(true_phi, ncol = 1)
    if (is.vector(est_phi)) est_phi <- matrix(est_phi, ncol = 1)
    sapply(1:ncol(true_phi), function(j) {
      diff_sq <- (est_phi[, j] - true_phi[, j])^2
      trapz(tfine, diff_sq)
    })
  }
  
  ise_vec <- sapply(seq_along(Ly_sim_avg), function(i) {
    est_phi <- Ly_sim_avg[[i]]
    true_phi <- sim_data$Y_no_noise[, i]
    calc_ise(true_phi, est_phi)
  })
  
  mean_fd_mat <- matrix(mean_fd, nrow = length(tfine), ncol = ncol(sim_data$Y_no_noise))
  ise_all <- calc_ise(true_phi = sim_data$Y_no_noise, est_phi = mean_fd_mat)

  # --- Step 4b: IE computation ---
  calc_ie <- function(true_phi, est_phi) {
    if (is.vector(true_phi)) true_phi <- matrix(true_phi, ncol = 1)
    if (is.vector(est_phi)) est_phi <- matrix(est_phi, ncol = 1)
    sapply(1:ncol(true_phi), function(j) {
      diff <- true_phi[, j] - est_phi[, j]
      trapz(tfine, diff)
    })
  }
  
  ie_vec <- sapply(seq_along(Ly_sim_avg), function(i) {
    est_phi <- Ly_sim_avg[[i]]
    true_phi <- sim_data$Y_no_noise[, i]
    calc_ie(true_phi, est_phi)
  })

  ie_all <- calc_ie(true_phi = sim_data$Y_no_noise, est_phi = mean_fd_mat)

  # --- Step 5: Plots ---
  if (plot) {
    
    # --- Plot smoothed trajectories ---
    matplot(tfine, fd_eval_mat, type = "l", lty = 1, 
            xlab = "Age", ylab = "Value", main = "Smoothed Estimated Trajectories")
    
    # --- Subject-level trajectory plots ---
    if (is.null(plot_subjects)) {
      plot_subjects <- sample(seq_along(Ly_sim), 5)
    }
    
    for (i in plot_subjects) {
      sim_trajs <- Ly_sim[[i]]
      true_traj <- sim_data$Y_no_noise[, i]
      sparse_time <- sim_al_data$Lt[[i]]
      sparse_vals <- sim_al_data$Ly[[i]]
      avg_traj <- Ly_sim_avg[[i]]
      
      sim_df <- as.data.frame(sim_trajs)
      sim_df$Time <- age_grid
      sim_long <- pivot_longer(sim_df, -Time, names_to = "Trajectory", values_to = "Value")
      true_df <- data.frame(Time = age_grid, Value = true_traj)
      obs_df <- data.frame(Time = sparse_time, Value = sparse_vals)
      avg_df <- data.frame(Time = age_grid, Value = avg_traj)
      
      p <- ggplot() +
        geom_line(data = sim_long, aes(x = Time, y = Value, group = Trajectory),
                  color = "steelblue", alpha = 0.4) +
        geom_line(data = true_df, aes(x = Time, y = Value, color = "True"), size = 1) +
        geom_line(data = avg_df, aes(x = Time, y = Value, color = "Avg Simulated"), size = 1) +
        geom_point(data = obs_df, aes(x = Time, y = Value, color = "Observed"), size = 2) +
        labs(title = paste("Subject", i, ": Simulated Trajectories, True, Observed, and Average"),
             x = "Age", y = "Value", color = "Legend") +
        theme_minimal() +
        scale_color_manual(values = c("True" = "red", "Avg Simulated" = "blue", "Observed" = "black")) +
        theme(legend.position = "bottom")
      
      print(p)
    }
    
    # --- Step 6: SE summary boxplot ---
    metrics_df <- data.frame(
      subject = seq_along(ie_vec),
      IE_Opt = ie_vec,
      IE_All = ie_all,
      ISE_Opt = ise_vec,
      ISE_All = ise_all,
      MSE_Opt = mse_vec,
      MSE_All = mse_all,
      Bias_Opt = bias_vec,
      Bias_All = bias_all
    )
    
    # Convert to long format for grouped plotting
    metrics_long <- metrics_df %>%
      pivot_longer(-subject, names_to = "Metric", values_to = "Value") %>%
      separate(Metric, into = c("Metric_Type", "Version"), sep = "_")
    
    print(plot_metric("ISE", df = metrics_long))
    print(plot_metric("IE", df = metrics_long))
    print(plot_metric("MSE", df = metrics_long))
    print(plot_metric("Bias", df = metrics_long))
    
  }
  
  
  # --- Return results ---
  return(list(
    Ly_sim = Ly_sim,
    Ly_sim_avg = Ly_sim_avg,
    ie_opt = ie_vec,
    ie_all = ie_all,
    ise_opt = ise_vec,
    ise_all = ise_all,
    mse_opt = mse_vec,
    mse_all = mse_all,
    bias_opt = bias_vec,
    bias_all = bias_all
  ))
}

