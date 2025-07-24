simulate_exponential_decay <- function(n_subjects = 500, 
                                       x = seq(0, 19, length.out = 1000),
                                       mu_L = 50, sd_L = 1,
                                       mu_a = 50, sd_a = 1,
                                       mu_k = 0.5, sd_k = 0.1,
                                       noise = 0.1,
                                       subgrp = 0) {
  
  # Validate subgrp input
  if (subgrp < 0 || subgrp > 1) stop("`subgrp` must be between 0 and 1.")
  
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
  
  # --- Apply slow decay to subgrp% of subjects ---
  # Generalized slow decay for arbitrary base
  apply_slow_decay <- function(base, x, A = 0.06, k = 1, center = 10) {
    S <- 1 / (1 + exp(-k * (x - center)))
    base - A * (x - center)^2 * S
  }

  # Determine which subjects to modify
  n_modify <- floor(subgrp * n_subjects)
  modify_idx <- sample(n_subjects, n_modify)

  for (i in modify_idx) {
    Y_no_noise[, i] <- apply_slow_decay(Y_no_noise[, i], x)
    Y[, i] <- Y_no_noise[, i] + rnorm(length(x), sd = noise)
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
    Ly_true = Ly_true,
    slow_decay_subjects = modify_idx
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
  
  subsample_idx <- sapply(monthly_grid, function(t) which.min(abs(original_grid - t)))
  Y_full_subsampled <- Y_full[subsample_idx, ]
  x_full <- monthly_grid
  
  # Step 2: Precompute subject IDs for fixed start and fixed end
  n_fixed <- round(x * n_subjects)
  fixed_start_subjects <- sample(seq_len(n_subjects), n_fixed)
  remaining_subjects <- setdiff(seq_len(n_subjects), fixed_start_subjects)
  fixed_end_subjects <- sample(remaining_subjects, n_fixed)
  
  Lt <- vector("list", n_subjects)
  Ly <- vector("list", n_subjects)
  
  for (i in seq_len(n_subjects)) {
    sampled_times <- numeric(n_time_points)
    
    if (i %in% fixed_start_subjects) {
      # Forward sampling with fixed start
      sampled_times[1] <- time_range[1]
      for (j in 2:n_time_points) {
        min_time <- sampled_times[j - 1] + 0.9
        max_time <- sampled_times[j - 1] + 2
        sampled_times[j] <- runif(1, min = min_time, max = min(max_time, time_range[2]))
      }
    } else if (i %in% fixed_end_subjects) {
      # Backward sampling with fixed end
      sampled_times[n_time_points] <- time_range[2]
      for (j in (n_time_points - 1):1) {
        max_time <- sampled_times[j + 1] - 0.9
        min_time <- sampled_times[j + 1] - 2
        sampled_times[j] <- runif(1, max(min_time, time_range[1]), max_time)
      }
    } else {
      # Random forward sampling (default)
      max_start <- time_range[2] - (n_time_points - 1) * 2
      sampled_times[1] <- runif(1, min = time_range[1], max = max_start)
      for (j in 2:n_time_points) {
        min_time <- sampled_times[j - 1] + 0.9
        max_time <- sampled_times[j - 1] + 2
        if (min_time >= time_range[2]) {
          sampled_times <- sampled_times[1:(j - 1)]
          break
        }
        sampled_times[j] <- runif(1, min = min_time, max = min(max_time, time_range[2]))
      }
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
  Ly_at_t0 <- Ly_unlisted[abs(Lt_unlisted - age_grid[1]) < eps]
  
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
  
  # Convert to wide format for paired t-test
  df_wide <- df_sub %>%
    dplyr::select(subject, Version, Value) %>%
    pivot_wider(names_from = Version, values_from = Value)
  
  # Paired t-test
  ttest_res <- t.test(df_wide$Opt, df_wide$All, paired = TRUE)
  t_stat <- round(ttest_res$statistic, 2)
  p_val <- signif(ttest_res$p.value, 3)
  test_label <- paste0("t = ", t_stat, ", p = ", p_val)
  
  # Calculate outliers
  df_outliers <- df_sub %>%
    group_by(Version) %>%
    mutate(
      Q1 = quantile(Value, 0.25, na.rm = TRUE),
      Q3 = quantile(Value, 0.75, na.rm = TRUE),
      IQR_val = Q3 - Q1,
      upper_whisker = Q3 + 1.5 * IQR_val,
      lower_whisker = Q1 - 1.5 * IQR_val,
      is_outlier = Value > upper_whisker | Value < lower_whisker
    ) %>%
    ungroup()
  
  # Plot
  ggplot(df_outliers, aes(x = Version, y = Value)) +
    geom_boxplot(outlier.shape = NA, fill = "lightgray") +
    geom_jitter(aes(color = is_outlier), width = 0.1, size = 2, alpha = 0.6) +
    geom_text(data = filter(df_outliers, is_outlier),
              aes(label = subject),
              hjust = -0.3, color = "red", size = 3.5) +
    scale_color_manual(values = c("black", "red"), guide = "none") +
    labs(title = paste(metric_type, ": Opt vs All Simulated Trajectories"),
         subtitle = test_label,
         y = metric_type, x = "") +
    theme_minimal()
}

evaluate_simulations <- function(tfine, sim_al_data, sim_data, acd_fda_2, plot = FALSE, plot_subjects = NULL) {
  
  # --- Step 1: Setup grid and basis ---
  n_obs <- tfine[length(tfine)]
  Y_estimated <- acd_fda_2$path
  time_est <- seq(tfine[1], n_obs, length.out = ncol(Y_estimated))
  
  # --- Step 2: Smooth each subject ---
  n_sim_traj <- nrow(Y_estimated)
  fd_eval_mat <- matrix(NA, nrow = length(tfine), ncol = n_sim_traj)
  error_count <- 0
  
  for (i in 1:n_sim_traj) {
    y <- Y_estimated[i, ]
    tryCatch({
      smooth_fit <- smooth.spline(x = time_est, y = y, cv = FALSE)
      fd_eval_mat[, i] <- predict(smooth_fit, x = tfine)$y
    }, error = function(e) {
      error_count <<- error_count + 1
    })
  }
  if (error_count > 0) {
    message("  Skipped ", error_count, " simulated trajectories due to smooth.spline errors.")
  }
  
  mean_fd <- rowMeans(fd_eval_mat, na.rm = TRUE)
  
  # Step 3: Matching and simulated fits
  age_grid <- tfine
  n_sim <- ncol(fd_eval_mat)
  n_sub <- length(sim_al_data$Ly)
  
  Ly_sim <- vector("list", n_sub)
  Ly_sim_avg <- vector("list", n_sub)
  mse_vec <- numeric(n_sub)
  bias_vec <- numeric(n_sub)
  mse_all <- numeric(n_sub)
  bias_all <- numeric(n_sub)
  bound_vec <- matrix(NA, nrow = length(tfine), ncol = n_sub)
  
  for (i in seq_len(n_sub)) {
    obs_time <- sim_al_data$Lt[[i]]
    obs_vals <- sim_al_data$Ly[[i]]
    if (length(obs_vals) == 0 && is.numeric(obs_vals)) {
      # Assign NA and skip if obs_vals is numeric(0)
      Ly_sim[[i]] <- NA
      Ly_sim_avg[[i]] <- NA
      mse_vec[i] <- NA
      bias_vec[i] <- NA
      mse_all[i] <- NA
      bias_all[i] <- NA
      bound_vec[, id] <- rep(NA, length(tfine))
      next
    }
    pos <- sapply(obs_time, function(t) {
      rounded_grid <- round(age_grid, 3)
      idx <- which(rounded_grid == round(t, 3))
      
      if (length(idx) > 0) {
        return(idx[1])
      } else {
        # Find the closest two neighbors
        diffs <- rounded_grid - round(t, 3)
        before <- max(which(diffs < 0))
        after <- min(which(diffs > 0))
        
        if (!is.na(before) && !is.na(after)) {
          # Randomly select between before and after
          return(sample(c(before, after), 1))
        } else if (!is.na(before)) {
          return(before)
        } else if (!is.na(after)) {
          return(after)
        } else {
          return(NA_integer_)
        }
      }
    })
    
    sim_vals <- fd_eval_mat[pos, ]
    sim_mean_vals <- mean_fd[pos]
    
    if (length(obs_vals) == 1) {
      # Not enough data to compute meaningful metrics
      Ly_sim[[i]] <- NA
      Ly_sim_avg[[i]] <- NA
      mse_vec[i] <- NA
      bias_vec[i] <- NA
      mse_all[i] <- NA
      bias_all[i] <- NA
      bound_vec[, id] <- rep(NA, length(tfine))
      next
    } 
    
    mse <- colMeans((sim_vals - obs_vals)^2)
    sorted_idx <- order(mse)
    inc_mse <- compute_incremental_avg_mse(sim_vals, obs_vals)
    elbow_index <- inc_mse$elbow_index
    top_idx <- sorted_idx[1:elbow_index]
    
    Ly_sim[[i]] <- fd_eval_mat[, top_idx]
    Ly_sim_avg[[i]] <- rowMeans(Ly_sim[[i]])
    mse_vec[i] <- inc_mse$cumulative_avg_mse[elbow_index]
    bias_vec[i] <- mean(obs_vals - rowMeans(sim_vals[, top_idx, drop = FALSE]))
    mse_all[i] <- mean((sim_mean_vals - obs_vals)^2)
    bias_all[i] <- mean(obs_vals - sim_mean_vals)
    bound_vec[, i] <- apply(Ly_sim[[i]], 1, function(x) {
      if (all(is.finite(x))) {
        quantile(x, 0.95, na.rm = TRUE) - quantile(x, 0.05, na.rm = TRUE)
      } else {
        NA  
      }
    })
  }
  
  # Initialize ISE/IE variables
  ise_vec <- NULL
  ise_all <- NULL
  ie_vec <- NULL
  ie_all <- NULL
  
  if (!is.null(sim_data)) {
    # Step 4a: ISE
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
    ise_all <- calc_ise(sim_data$Y_no_noise, mean_fd_mat)
    
    # Step 4b: IE
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
    ie_all <- calc_ie(sim_data$Y_no_noise, mean_fd_mat)
  }
  
  # --- Step 5: Plots ---
  if (plot) {

    est_df <- as.data.frame(fd_eval_mat) %>%
      mutate(Time = tfine) %>%
      pivot_longer(-Time, names_to = "Subject", values_to = "y") %>%
      mutate(group = "Estimated")

    obs_df <- map2_dfr(sim_al_data$Ly, sim_al_data$Lt,
                       ~ tibble(Time = .y, y = .x),
                       .id = "Subject") %>%
      mutate(group = "Observed")

    df_all <- bind_rows(est_df, obs_df)

    trajectory_plot <- ggplot(df_all, aes(x = Time, y = y, group = Subject, color = group)) +
      geom_line(data = filter(df_all, group == "Estimated"), 
                alpha = 0.4, linewidth = 0.4) +
      geom_line(data = filter(df_all, group == "Observed"), 
                linewidth = 0.4) +
      scale_color_manual(values = c("Estimated" = "steelblue",
                                    "Observed" = "red")) +
      labs(x = "Time", y = "Functional Value", color = "Trajectory Type",
           title = "Accelerated Longitudinal Sampling:\nObserved vs Estimated Trajectories") +
      theme_minimal(base_size = 14)
    
    print(trajectory_plot)
    
    if (is.null(plot_subjects)) {
      valid_indices <- which(sapply(Ly_sim, function(x) length(x) > 1))
      plot_subjects <- sample(valid_indices, 5)
    }
    
    for (i in plot_subjects) {
      sim_trajs <- Ly_sim[[i]]
      avg_traj <- Ly_sim_avg[[i]]
      sparse_time <- sim_al_data$Lt[[i]]
      sparse_vals <- sim_al_data$Ly[[i]]
      obs_df <- data.frame(Time = sparse_time, Value = sparse_vals)
      avg_df <- data.frame(Time = age_grid, Value = avg_traj)
      
      sim_df <- as.data.frame(sim_trajs)
      sim_df$Time <- age_grid
      sim_long <- pivot_longer(sim_df, -Time, names_to = "Trajectory", values_to = "Value")
      
      p <- ggplot() +
        geom_line(data = sim_long, aes(x = Time, y = Value, group = Trajectory), color = "steelblue", alpha = 0.4) +
        geom_point(data = obs_df, aes(x = Time, y = Value, color = "Observed"), size = 2) +
        geom_line(data = avg_df, aes(x = Time, y = Value, color = "Avg Simulated"), size = 1)
      
      if (!is.null(sim_data)) {
        true_df <- data.frame(Time = age_grid, Value = sim_data$Y_no_noise[, i])
        p <- p + geom_line(data = true_df, aes(x = Time, y = Value, color = "True"), size = 1)
      }
      
      p <- p +
        labs(title = paste("Subject", i, ": Simulated Trajectories"),
             x = "Age", y = "Value", color = "Legend") +
        scale_color_manual(values = c("True" = "red", "Avg Simulated" = "blue", "Observed" = "black")) +
        theme_minimal() +
        theme(legend.position = "bottom")
      
      print(p)
    }
    
    # --- Step 6: SE summary boxplot ---
    metrics_df <- data.frame(
      subject = seq_len(n_sub),
      MSE_Opt = mse_vec,
      MSE_All = mse_all,
      Bias_Opt = bias_vec,
      Bias_All = bias_all
    )
    
    if (!is.null(sim_data)) {
      metrics_df$IE_Opt <- ie_vec
      metrics_df$IE_All <- ie_all
      metrics_df$ISE_Opt <- ise_vec
      metrics_df$ISE_All <- ise_all
    }
    
    metrics_long <- metrics_df %>%
      pivot_longer(-subject, names_to = "Metric", values_to = "Value") %>%
      separate(Metric, into = c("Metric_Type", "Version"), sep = "_")
    
    if (!is.null(sim_data)) {
      print(plot_metric("ISE", df = metrics_long))
      print(plot_metric("IE", df = metrics_long))
    }
    
    print(plot_metric("MSE", df = metrics_long))
    print(plot_metric("Bias", df = metrics_long))
  }
  
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
    bias_all = bias_all,
    d_bound = bound_vec
  ))
}

library(OpenMx)
library(dplyr)
library(parallel)

# Function to build and fit multigroup model
fit_multisubj_ct_model <- function(data) {
  opmxL <- list()
  
  # Matrix A. Latent dynamics
  opmxL$amat_ct <- mxMatrix(name = "A", "Full", 2, 2, free = c(TRUE,FALSE,FALSE,FALSE),
                            values = c(-.2,0,1,0),
                            dimnames = list( c("y0", "yA"), c("y0", "yA") ),
                            labels = c("b_y", NA,NA,NA),
                            lbound = c(-1, NA,NA,NA),
                            ubound = c(0, NA,NA,NA))
  
  # Matrix B. Effects of time-varying covariates on latent variables (null in our model)
  opmxL$bmat <- mxMatrix(name = "B", "Zero", 2, 1)
  
  # Matrix C. Measurement model linking latent and observed variables
  opmxL$cmat <- mxMatrix(name = "C", "Full", 1, 2, free = FALSE,
                         values = c(1,0), 
                         dimnames = list( c("y"), c("y0", "yA") ),
                         labels = c(NA, NA)  )
  
  # Matrix D. Effects of time-varying covariates on observed variables (null in our model)
  opmxL$dmat <- mxMatrix("Zero", 1, 1, name = "D")
  
  # Matrix Q. Covariance matrix for the dynamic error (i.e., latent innovations)
  opmxL$qmat <- mxMatrix("Zero", 2, 2, name = "Q")
  
  # Matrix R. Covariance matrix for the measurement error
  opmxL$rmat <- mxMatrix("Diag", 1, 1, TRUE, 2,
                         name = "R", labels = "MerY")
  
  # Matrix x0. Means of the initial latent levels
  opmxL$xmat <- mxMatrix(name = "x0", "Full", 2, 1, free = TRUE,
                         values = c(12, 7),
                         labels = c("y0mn", "yAmn"))
  
  # Matrix P0. Covariance of the initial latent levels
  opmxL$pmat <- mxMatrix(name = "P0", "Symm", 2, 2, TRUE,
                         values = c(25, 3, .4),
                         labels = c("y0v", "y0Acv", "yAv"),
                         lbound = c(0, NA, 0))
  
  # Matrix u. Values of the time-varying covariates (null in our model)
  opmxL$umat <- mxMatrix("Zero", 1, 1, name = "u")
  
  # Matrix t. Variable in the data set with the times at which measurement happened.
  opmxL$tmat <- mxMatrix('Full', 1, 1, name='time', labels='data.age')
  
  opmxL$modL_ct <- with(opmxL, list(amat_ct, bmat, cmat, dmat,
                                    qmat, rmat, xmat, pmat,
                                    umat, tmat))
  
  opmxL$expODE <- mxExpectationStateSpaceContinuousTime(A = "A", B = "B",
                                                        C = "C", D = "D",
                                                        Q = "Q", R = "R",
                                                        x0 = "x0", P0 = "P0",
                                                        u = "u", t = "time")
  
  genMxIndModels_ct <- function(x, dwork, modNames_ct) {
    DataSetForSubjectK <- dwork[[x]]
    ctModel <- opmxL$modL_ct
    indivmodels <- mxModel(name = modNames_ct[x],
                           ctModel,
                           opmxL$expODE,
                           mxFitFunctionML(),
                           mxData(DataSetForSubjectK, type ='raw')  )  }
  
  # Create individual models
  cases <- seq(1,length(data))
  modNames_ct <- paste0("i", cases, "ODE") 
  indivmodels_ct <- lapply(cases, genMxIndModels_ct, data, modNames_ct)
  # Create multiple-subject model
  multiSubjODE <- mxModel(name = "MultiODE", indivmodels_ct,
                          mxFitFunctionMultigroup(modNames_ct))
  # Fit the model
  multiSubjODERun <- mxRun(multiSubjODE, silent = TRUE)
  multiSubjODERun
}

# Function to compute Kalman-smoothed trajectories for a subject
simulate_one_kalman_path <- function(id) {
  # Get Kalman scores
  ks_scores <- suppressMessages(
    OpenMx::mxKalmanScores(model = multiSubjODERun@submodels[[id]], data = data[[id]])
  )
  y_KS <- ks_scores$xSmoothed[-1, 1]  # Drop first row
  
  # Compute confidence intervals
  errors <- sqrt(ks_scores$PSmoothed[1, 1, 2:(length(y_KS) + 1)])
  alpha <- 0.05
  qZ <- qnorm(1 - alpha / 2)
  upper <- y_KS + qZ * errors
  lower <- y_KS - qZ * errors
  
  # Construct output dataframe
  df <- data.frame(
    id = id,
    age = data[[id]]$age,
    y_observed = data[[id]]$y,
    y_KS = y_KS,
    CI_upper = upper,
    CI_lower = lower,
    y_true = sim_data$Y_no_noise[, id]
  )
  
  return(df)
}


# Main execution
run_full_kalman_simulation <- function(data, parallel = TRUE) {
  multiSubjODERun <<- fit_multisubj_ct_model(data)
  n_subj <- length(data)
  
  if (parallel) {
    num_cores <- detectCores() - 1
    cl <- makeCluster(num_cores)
    clusterExport(cl, varlist = c("multiSubjODERun", "data", "simulate_one_kalman_path", "sim_data"), envir = environment())
    LCS_SSM <- parLapply(cl, 1:n_subj, simulate_one_kalman_path)
    stopCluster(cl)
  } else {
    LCS_SSM <- lapply(1:n_subj, simulate_one_kalman_path)
  }
  
  LCS_SSM
}

convert_to_full_subject_list <- function(sim_al_data) {
  n_subjects <- length(sim_al_data$Lt)
  x_full <- sim_al_data$x_full
  
  full_data_list <- vector("list", n_subjects)
  
  for (i in seq_len(n_subjects)) {
    # Initialize vectors for all time points as NA
    y_full <- rep(NA_real_, length(x_full))
    
    # Get indices where the subject has data
    match_idx <- match(sim_al_data$Lt[[i]], x_full)
    y_full[match_idx] <- sim_al_data$Ly[[i]]
    
    # Construct data.frame for the subject
    subject_df <- tibble::tibble(
      id = i,
      y = y_full,
      age = x_full
    )
    
    full_data_list[[i]] <- subject_df
  }
  
  return(full_data_list)
}

evaluate_lcsssm <- function(LCS_SSM, sim_data, sim_al_data, plot = FALSE, plot_subjects = NULL) {
  library(dplyr)
  library(ggplot2)
  library(tibble)
  library(pracma)
  library(tidyr)
  
  n_subj <- length(LCS_SSM)
  tfine <- sim_al_data$x_full
  
  # Helper functions
  calc_ise <- function(true_phi, est_phi, tfine) {
    if (is.vector(true_phi)) true_phi <- matrix(true_phi, ncol = 1)
    if (is.vector(est_phi)) est_phi <- matrix(est_phi, ncol = 1)
    sapply(1:ncol(true_phi), function(j) {
      diff_sq <- (est_phi[, j] - true_phi[, j])^2
      trapz(tfine, diff_sq)
    })
  }
  
  calc_ie <- function(true_phi, est_phi, tfine) {
    if (is.vector(true_phi)) true_phi <- matrix(true_phi, ncol = 1)
    if (is.vector(est_phi)) est_phi <- matrix(est_phi, ncol = 1)
    sapply(1:ncol(true_phi), function(j) {
      diff <- true_phi[, j] - est_phi[, j]
      trapz(tfine, diff)
    })
  }
  
  mse_vec <- numeric(n_subj)
  bias_vec <- numeric(n_subj)
  ise_vec <- numeric(n_subj)
  ie_vec <- numeric(n_subj)
  bound_vec <- matrix(NA, nrow = length(tfine), ncol = n_subj)
  
  for (id in seq_len(n_subj)) {
    ks <- LCS_SSM[[id]]$y_KS
    true_traj <- sim_data$Y_no_noise[, id]
    
    ise_vec[id] <- calc_ise(true_traj, ks, tfine)
    ie_vec[id] <- calc_ie(true_traj, ks, tfine)
    
    obs_ages <- LCS_SSM[[id]]$age[!is.na(LCS_SSM[[id]]$y_observed)]
    obs_vals <- LCS_SSM[[id]]$y_observed[!is.na(LCS_SSM[[id]]$y_observed)]
    if (length(obs_vals) == 0) next
    
    true_vals_at_obs <- approx(x = tfine, y = true_traj, xout = obs_ages, rule = 2)$y
    ks_vals_at_obs   <- approx(x = tfine, y = ks, xout = obs_ages, rule = 2)$y
    
    mse_vec[id] <- mean((true_vals_at_obs - obs_vals)^2)
    bias_vec[id] <- mean(ks_vals_at_obs - obs_vals)
    bound_vec[, id] <- with(LCS_SSM[[id]], ifelse(is.finite(CI_upper) & is.finite(CI_lower), CI_upper - CI_lower, NA))
  }
  
  if (plot) {
    # 1. Subject-level plot for selected subjects
    if (is.null(plot_subjects)) {
      plot_subjects <- sample(seq_len(n_subj), min(5, n_subj))
    }
    
    for (id in plot_subjects) {
      df_subj <- LCS_SSM[[id]]

      df_plot <- df_subj |>
        dplyr::mutate(
          `Kalman Upper` = CI_upper,
          `Kalman Lower` = CI_lower,
          `Kalman Smoothed` = y_KS,
          `True` = y_true,
          `Observed` = y_observed
        )
      
      p <- ggplot(data = df_plot, aes(x = age)) +
        
        # Confidence ribbon
        geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper, fill = "95% CI"),
                    alpha = 0.15, show.legend = TRUE) +
        
        # Kalman CI bounds
        geom_line(aes(y = CI_upper, color = "95% CI Bounds"), size = 0.2) +
        geom_line(aes(y = CI_lower, color = "95% CI Bounds"), size = 0.2) +
        
        # Kalman-smoothed line
        geom_line(aes(y = y_KS, color = "Kalman Smoothed"),
                  size = 0.8, linetype = "dashed") +
        
        # Observed points and line
        geom_point(aes(y = y_observed, color = "Observed"), size = 2) +
        geom_line(data = df_plot[!is.na(df_plot$y_observed), ],
                  aes(y = y_observed, color = "Observed"), size = 0.5) +
        
        # True latent line
        geom_line(aes(y = y_true, color = "True"), size = 0.8) +
        
        # Manual color and fill scales
        scale_color_manual(name = "Legend",
                           values = c("Observed" = "black",
                                      "Kalman Smoothed" = "blue",
                                      "True" = "red",
                                      "95% CI Bounds" = "blue")) +
        scale_fill_manual(name = "Legend",
                          values = c("95% CI" = "blue")) +
        
        scale_x_continuous(breaks = seq(0, 20, 1)) +
        ylab("y scores") + 
        xlab("Age") +
        ggtitle(paste("Subject", id, "- Kalman Smoothed vs Observed and True")) +
        theme_minimal() +
        theme(legend.position = "bottom",
              legend.title = element_blank())
      
      print(p)
    }
    
    # 2. Error metric boxplots
    results_df <- data.frame(
      Subject = seq_along(ise_vec),
      ISE = ise_vec,
      IE = ie_vec,
      MSE = mse_vec,
      Bias = bias_vec
    )
    
    df_long <- results_df %>%
      pivot_longer(cols = -Subject, names_to = "Metric", values_to = "Value")
    
    df_outliers <- df_long %>%
      group_by(Metric) %>%
      mutate(
        Q1 = quantile(Value, 0.25, na.rm = TRUE),
        Q3 = quantile(Value, 0.75, na.rm = TRUE),
        IQR_val = Q3 - Q1,
        upper_whisker = Q3 + 1.5 * IQR_val,
        lower_whisker = Q1 - 1.5 * IQR_val,
        is_outlier = Value > upper_whisker | Value < lower_whisker
      ) %>%
      ungroup()
    
    metric_levels <- c("MSE", "Bias", "ISE", "IE")
    df_outliers$Metric <- factor(df_outliers$Metric, levels = metric_levels)
    
    boxplot_panel <- ggplot(df_outliers, aes(x = "", y = Value)) +
      geom_boxplot(outlier.shape = NA, fill = "lightgray") +
      geom_jitter(aes(color = is_outlier), width = 0.1, size = 2, alpha = 0.6) +
      geom_text(
        data = filter(df_outliers, is_outlier),
        aes(label = Subject),
        hjust = -0.3, color = "red", size = 3.5
      ) +
      scale_color_manual(values = c("black", "red"), guide = "none") +
      facet_wrap(~ Metric, scales = "free_y", ncol = 2) +
      labs(
        title = "Subject-Level Error Metrics",
        y = "Value", x = NULL
      ) +
      theme_minimal() +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    
    print(boxplot_panel)
  }
  
  return(list(
    mse = mse_vec,
    bias = bias_vec,
    ise = ise_vec,
    ie = ie_vec,
    d_bound = bound_vec
  ))
}

plot_metric_s2 <- function(metric_type, df = metrics_long) {
  df_sub <- df %>% filter(Metric_Type == metric_type)
  
  # Convert to wide format for paired t-test
  df_wide <- df_sub %>%
    dplyr::select(subject, Version, Value) %>%
    tidyr::pivot_wider(names_from = Version, values_from = Value)
  
  # Paired t-test
  ttest_res <- t.test(df_wide$FDA, df_wide$LCS, paired = TRUE)
  t_stat <- round(ttest_res$statistic, 2)
  p_val <- signif(ttest_res$p.value, 3)
  test_label <- paste0("t = ", t_stat, ", p = ", p_val)
  
  # Calculate outliers
  df_outliers <- df_sub %>%
    group_by(Version) %>%
    mutate(
      Q1 = quantile(Value, 0.25, na.rm = TRUE),
      Q3 = quantile(Value, 0.75, na.rm = TRUE),
      IQR_val = Q3 - Q1,
      upper_whisker = Q3 + 1.5 * IQR_val,
      lower_whisker = Q1 - 1.5 * IQR_val,
      is_outlier = Value > upper_whisker | Value < lower_whisker
    ) %>%
    ungroup()
  
  # Plot
  ggplot(df_outliers, aes(x = Version, y = Value)) +
    geom_boxplot(outlier.shape = NA, fill = "gray90", color = "black") +
    geom_jitter(aes(color = is_outlier), width = 0.1, size = 2, alpha = 0.6) +
    geom_text(data = filter(df_outliers, is_outlier),
              aes(label = subject),
              hjust = -0.3, color = "red", size = 3.5) +
    scale_color_manual(values = c("black", "red"), guide = "none") +
    labs(title = paste(metric_type, ": FDA vs LCS"),
         subtitle = test_label,
         y = metric_type, x = "") +
    theme_minimal(base_size = 14)
}

plot_d_bound_comparison <- function(d_bound_fda, d_bound_lcs, n_time) {
  
  compute_summary_df <- function(d_bound_matrix, method_name, n_time) {
    data <- data.frame(
      Time = n_time,
      Mean = NA_real_,
      CI_lower = NA_real_,
      CI_upper = NA_real_,
      Method = method_name
    )
    
    for (t in 1:length(n_time)) {
      vals <- d_bound_matrix[t, ]
      vals <- vals[is.finite(vals)]
      n_vals <- length(vals)
      
      if (n_vals > 0) {
        mean_val <- mean(vals)
        se <- if (n_vals > 1) sd(vals) / sqrt(n_vals) else 0
        t_crit <- if (n_vals > 1) qt(0.975, df = n_vals - 1) else NA
        
        data$Mean[t] <- mean_val
        data$CI_lower[t] <- if (!is.na(t_crit)) mean_val - t_crit * se else NA
        data$CI_upper[t] <- if (!is.na(t_crit)) mean_val + t_crit * se else NA
      }
    }
    
    return(data)
  }
  
  df_fda <- compute_summary_df(d_bound_fda, "FDA", n_time)
  df_lcs <- compute_summary_df(d_bound_lcs, "LCS", n_time)
  plot_df <- rbind(df_fda, df_lcs)
  
  ggplot(plot_df, aes(x = Time, y = Mean, color = Method, fill = Method)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2, linetype = 0) +
    labs(title = "95% CI of Error Bound Width: FDA vs LCS",
         x = "Age", y = "Width of Error Bound") +
    scale_color_manual(values = c("FDA" = "blue", "LCS" = "darkred")) +
    scale_fill_manual(values = c("FDA" = "blue", "LCS" = "darkred")) +
    theme_minimal(base_size = 14)
}


compare_fda_lcs <- function(results, results2, LCS_SSM, 
                            plot = TRUE, plot_subjects = NULL) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  n_subj <- length(LCS_SSM)
  tfine <- LCS_SSM[[1]]$age
  
  metrics_df <- data.frame(
    subject = seq_len(n_subj),
    MSE_FDA = results$mse_opt,
    MSE_LCS = results2$mse,
    Bias_FDA = results$bias_opt,
    Bias_LCS = results2$bias,
    IE_FDA = results$ie_opt,
    IE_LCS = results2$ie,
    ISE_FDA = results$ise_opt,
    ISE_LCS = results2$ise
  )
  
  # 1. Plot Subject Trajectories
  if (plot) {
    if (is.null(plot_subjects)) {
      set.seed(123)
      plot_subjects <- sample(seq_len(n_subj), 5)
    }
    
    for (id in plot_subjects) {
      fda_id_vec <- results$Ly_sim_avg[[id]]
      
      df_subj <- LCS_SSM[[id]] %>%
        mutate(Functional_Avg = fda_id_vec)
      
      p <- ggplot(df_subj, aes(x = age)) +
        geom_point(aes(y = y_observed, color = "Observed"), size = 2) +
        geom_line(aes(y = y_observed, color = "Observed"), size = 0.5) +
        geom_line(aes(y = y_KS, color = "LCS"), 
                  size = 0.8, linetype = "dashed") +
        geom_line(aes(y = Functional_Avg, color = "FDA SDE"), 
                  size = 0.8, linetype = "longdash") +
        geom_line(aes(y = y_true, color = "True"), 
                  size = 0.8, linetype = "solid") +
        scale_color_manual(
          name = "Legend",
          values = c("Observed" = "black",
                     "LCS" = "blue",
                     "FDA SDE" = "darkgreen",
                     "True" = "red")
        ) +
        scale_x_continuous(breaks = seq(0, 20, 1)) +
        labs(x = "Age", y = "y scores", 
             title = paste("Subject", id, "- Observed vs LCS vs FDA SDE vs True")) +
        theme_minimal(base_size = 12) +
        theme(legend.position = "bottom")
      
      print(p)
    }  
    
    
    metrics_long <- metrics_df %>%
      pivot_longer(-subject, names_to = "Metric", values_to = "Value") %>%
      separate(Metric, into = c("Metric_Type", "Version"), sep = "_")
    
    print(plot_metric_s2("MSE", metrics_long))
    print(plot_metric_s2("Bias", metrics_long))
    print(plot_metric_s2("IE", metrics_long))
    print(plot_metric_s2("ISE", metrics_long))
    
    error_bound_plot <- plot_d_bound_comparison(results$d_bound, results2$d_bound, tfine)
    print(error_bound_plot)
  }
  
  return(list(metrics = metrics_df,
              DB_FDA = results$d_bound,
              DB_LCS = results2$d_bound))
}
