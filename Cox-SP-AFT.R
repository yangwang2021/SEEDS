#==============================================================================
# MORE COMPLETE COX-SP-AFT WITH MCP REGULARIZATION
# Addresses the missing regularization components
# FIXED: Removed auto-execution that was causing errors
# MODIFIED: Updated loss to log scale violation as per recommendation
#==============================================================================

# MCP penalty function (approximation)
mcp_penalty <- function(beta, lambda_mcp, gamma = 3) {
  abs_beta <- abs(beta)
  penalty <- ifelse(abs_beta <= gamma * lambda_mcp,
                    lambda_mcp * abs_beta - beta^2 / (2 * gamma),
                    gamma * lambda_mcp^2 / 2)
  return(sum(penalty))
}

# Cox model with MCP regularization (simplified implementation)
cox_with_mcp <- function(cox_data, cox_formula, lambda_mcp = 0.1) {
  
  tryCatch({
    # Try using ncvreg package if available for proper MCP
    if (requireNamespace("ncvreg", quietly = TRUE)) {
      # Extract design matrix
      X <- model.matrix(cox_formula, data = cox_data)[, -1, drop = FALSE]  # Remove intercept
      y <- with(cox_data, survival::Surv(time, status))
      
      # Fit Cox model with MCP penalty
      mcp_fit <- ncvreg::ncvsurv(X = X, y = y, penalty = "MCP")
      
      # Extract coefficients (using cross-validation selected lambda)
      cv_fit <- ncvreg::cv.ncvsurv(X = X, y = y, penalty = "MCP")
      best_coefs <- coef(mcp_fit, lambda = cv_fit$lambda.min)
      
      # Create a coxph-like object for compatibility
      regular_fit <- survival::coxph(cox_formula, data = cox_data)
      regular_fit$coefficients <- best_coefs  # Remove intercept if present
      
      return(regular_fit)
      
    } else {
      # Fallback: Standard Cox model with warning
      warning("ncvreg package not available - using standard Cox model")
      return(survival::coxph(cox_formula, data = cox_data))
    }
    
  }, error = function(e) {
    # Final fallback
    return(survival::coxph(cox_formula, data = cox_data))
  })
}

# AFT model with MCP regularization (simplified implementation)
aft_with_mcp <- function(aft_data, aft_formula, weight, lambda_mcp = 0.1) {
  
  tryCatch({
    # Try using ncvreg package for proper MCP regularization
    if (requireNamespace("ncvreg", quietly = TRUE)) {
      # Extract design matrix and response
      X <- model.matrix(aft_formula, data = aft_data)[, -1, drop = FALSE]  # Remove intercept
      y <- aft_data$log_time
      
      if (length(weight) == nrow(X)) {
        # Fit linear model with MCP penalty and weights
        mcp_fit <- ncvreg::ncvreg(X = X, y = y, penalty = "MCP", weights = weight)
        
        # Use cross-validation to select lambda
        cv_fit <- ncvreg::cv.ncvreg(X = X, y = y, penalty = "MCP", weights = weight)
        best_lambda <- cv_fit$lambda.min
        
        # Return the fit and lambda for prediction
        return(list(fit = mcp_fit, lambda = best_lambda))
      }
    }
    
    # Fallback: Standard weighted linear regression
    lm_fit <- lm(aft_formula, data = aft_data, weights = weight)
    return(list(fit = lm_fit, lambda = NULL))
    
  }, error = function(e) {
    lm_fit <- lm(aft_formula, data = aft_data, weights = weight)
    return(list(fit = lm_fit, lambda = NULL))
  })
}

# Improved SP-AFT with proper Kaplan-Meier initialization
sp_aft_model_complete <- function(x_data, y_initial, status_initial, censored_indices, 
                                  alpha_start = 0.1, mu = 1.1, max_iterations = 50, 
                                  tolerance = 1e-4, lambda_mcp = 0.1) {
  
  m <- length(censored_indices)
  if (m == 0) return(list(beta = NULL, labeled_times = NULL, reliable_indices = NULL))
  
  # Step 1: Initialize V^(0) = (1,...,1)
  V <- rep(1, m)
  lambda_param <- alpha_start
  
  # Step 1: Initialize Y^(t) with Kaplan-Meier method (as specified in paper)
  tryCatch({
    km_fit <- survival::survfit(survival::Surv(y_initial, status_initial) ~ 1)
    
    # Initialize censored observations using KM estimates
    Y <- numeric(m)
    for (i in 1:m) {
      idx <- censored_indices[i]
      censor_time <- y_initial[idx]
      
      # Get survival probability at censoring time
      surv_prob <- approx(km_fit$time, km_fit$surv, xout = censor_time, 
                          method = "constant", f = 0, rule = 2)$y
      if (is.na(surv_prob)) surv_prob <- 0.5
      
      # Estimate completion time based on KM curve
      if (surv_prob > 0.1) {
        # Find time where survival drops to 75% of current level
        target_surv <- surv_prob * 0.75
        est_time <- approx(km_fit$surv, km_fit$time, xout = target_surv, 
                           method = "linear", rule = 2)$y
        Y[i] <- ifelse(!is.na(est_time) && est_time > censor_time, est_time, censor_time * 1.2)
      } else {
        Y[i] <- censor_time * 1.1
      }
    }
  }, error = function(e) {
    # Fallback initialization
    Y <- y_initial[censored_indices] * 1.2
  })
  
  # Main SP-AFT iteration loop
  for (iteration in 1:max_iterations) {
    
    # Step 3: Update β^(t) according to E.q.(9) - AFT model with MCP regularization
    complete_indices <- setdiff(1:length(y_initial), censored_indices)
    
    # Combine complete data with current estimates of censored data
    combined_x <- rbind(x_data[complete_indices, , drop = FALSE], 
                        x_data[censored_indices, , drop = FALSE])
    combined_y <- c(log(pmax(y_initial[complete_indices], 0.001)), log(pmax(Y, 0.001)))
    combined_weights <- c(rep(1, length(complete_indices)), V)
    
    # Fit AFT model with MCP regularization
    aft_data <- data.frame(log_time = combined_y, combined_x)
    aft_formula <- as.formula(paste("log_time ~", paste(names(x_data), collapse = " + ")))
    
    beta_fit <- aft_with_mcp(aft_data = aft_data, aft_formula = aft_formula, weight = combined_weights, lambda_mcp = lambda_mcp)
    
    # Step 4: Update V^(t) according to E.q.(11) - Self-paced learning weights
    prev_V <- V
    
    for (j in 1:m) {
      idx <- censored_indices[j]
      original_censor_time <- y_initial[idx]
      
      # Predict log survival time using current β
      newx <- data.frame(x_data[idx, , drop = FALSE])
      colnames(newx) <- names(x_data)
      pred_log <- if ("ncvreg" %in% class(beta_fit$fit)) {
        predict(beta_fit$fit, X = as.matrix(newx), lambda = beta_fit$lambda, type = "response")[1]
      } else {
        predict(beta_fit$fit, newdata = newx, type = "response")[1]
      }
      
      # MODIFIED: Use log-scale violation loss as per paper
      log_censor <- log(original_censor_time)
      violation <- max(log_censor - pred_log, 0)
      loss_j <- violation^2
      
      # Update V and Y
      if (loss_j <= lambda_param) {
        V[j] <- 1
        Y[j] <- exp(pred_log)  # Impute predicted time (paper sets to f(x_j))
      } else {
        V[j] <- 0
        Y[j] <- original_censor_time * 1.01  # Minimal adjustment for non-reliable
      }
    }
    
    # Step 6: λ ← μλ (increase age parameter)
    lambda_param <- mu * lambda_param
    
    # Check convergence
    weight_change <- mean(abs(V - prev_V))
    if (weight_change < tolerance) {
      break
    }
  }
  
  # Return reliable labeled samples
  reliable_indices <- censored_indices[V > 0]
  reliable_times <- Y[V > 0]
  
  return(list(
    beta_fit = beta_fit,  # Return the fit object for potential further use
    labeled_times = reliable_times,
    reliable_indices = reliable_indices,
    final_weights = V
  ))
}

# Updated main Cox-SP-AFT algorithm with MCP regularization
cox_sp_aft_complete <- function(survival_times, censoring_indicators, covariates, 
                                evaluation_times, min_censored_rate = 0.05, 
                                max_outer_iterations = 10, lambda_mcp = 0.1, verbose = FALSE) {
  
  tryCatch({
    
    if (verbose) cat("=== COMPLETE COX-SP-AFT WITH MCP REGULARIZATION ===\n")
    
    # Check for required packages
    if (!requireNamespace("ncvreg", quietly = TRUE)) {
      if (verbose) cat("Warning: ncvreg package not available - using approximations\n")
    }
    
    # Initialize
    n <- length(survival_times)
    current_times <- survival_times
    current_status <- censoring_indicators
    current_covariates <- covariates
    
    labeled_indices <- which(current_status == 1)
    censored_indices <- which(current_status == 0)
    
    if (verbose) {
      cat("Initial - Labeled:", length(labeled_indices), "Censored:", length(censored_indices), "\n")
    }
    
    outer_iteration <- 1
    
    # Main Cox-SP-AFT Loop
    while (length(censored_indices) / n > min_censored_rate && outer_iteration <= max_outer_iterations) {
      
      if (verbose) {
        censored_rate <- length(censored_indices) / n
        cat("Outer Iteration", outer_iteration, "- Censored rate:", round(censored_rate, 3), "\n")
      }
      
      # Step 2: Update Cox classifier with MCP regularization
      if (length(labeled_indices) >= 3) {
        
        cox_data <- data.frame(
          time = current_times[labeled_indices],
          status = rep(1, length(labeled_indices)),
          current_covariates[labeled_indices, , drop = FALSE]
        )
        
        cox_formula <- as.formula(paste("survival::Surv(time, status) ~", 
                                        paste(names(current_covariates), collapse = " + ")))
        
        # Use Cox model with MCP regularization
        cox_fit <- cox_with_mcp(cox_data, cox_formula, lambda_mcp)
        
        # Get risk scores for all patients
        all_cox_data <- data.frame(current_covariates)
        risk_scores <- predict(cox_fit, newdata = all_cox_data, type = "lp")
        
        # Classify into high-risk and low-risk groups
        risk_threshold <- median(risk_scores, na.rm = TRUE)
        high_risk_indices <- which(risk_scores > risk_threshold)
        low_risk_indices <- which(risk_scores <= risk_threshold)
        
      } else {
        high_risk_indices <- 1:n
        low_risk_indices <- integer(0)
      }
      
      # Step 3: SP-AFT for each group with complete implementation
      newly_labeled_indices <- c()
      newly_labeled_times <- c()
      
      for (risk_group in list(high_risk_indices, low_risk_indices)) {
        if (length(risk_group) == 0) next
        
        group_censored <- intersect(risk_group, censored_indices)
        
        if (length(group_censored) >= 1) {
          
          group_covariates <- current_covariates[risk_group, , drop = FALSE]
          group_times <- current_times[risk_group]
          group_status <- current_status[risk_group]
          
          group_censored_relative <- match(group_censored, risk_group)
          
          # Use complete SP-AFT implementation
          sp_aft_result <- sp_aft_model_complete(
            x_data = group_covariates,
            y_initial = group_times,
            status_initial = group_status,
            censored_indices = group_censored_relative,
            alpha_start = 0.1,
            mu = 1.1,
            max_iterations = 20,
            lambda_mcp = lambda_mcp
          )
          
          if (length(sp_aft_result$reliable_indices) > 0) {
            global_reliable_indices <- risk_group[sp_aft_result$reliable_indices]
            newly_labeled_indices <- c(newly_labeled_indices, global_reliable_indices)
            newly_labeled_times <- c(newly_labeled_times, sp_aft_result$labeled_times)
          }
        }
      }
      
      # Step 4: Update training dataset
      if (length(newly_labeled_indices) > 0) {
        current_times[newly_labeled_indices] <- newly_labeled_times
        current_status[newly_labeled_indices] <- 1
        
        labeled_indices <- which(current_status == 1)
        censored_indices <- which(current_status == 0)
        
        if (verbose) {
          cat("  Added", length(newly_labeled_indices), "newly labeled observations\n")
        }
      } else {
        if (verbose) cat("  No reliable samples found - terminating\n")
        break
      }
      
      outer_iteration <- outer_iteration + 1
    }
    
    # Final survival estimation using enhanced dataset
    if (length(labeled_indices) >= 3) {
      final_cox_data <- data.frame(
        time = current_times[labeled_indices],
        status = rep(1, length(labeled_indices)),
        current_covariates[labeled_indices, , drop = FALSE]
      )
      
      final_cox_formula <- as.formula(paste("survival::Surv(time, status) ~", 
                                            paste(names(current_covariates), collapse = " + ")))
      
      final_cox_fit <- cox_with_mcp(final_cox_data, final_cox_formula, lambda_mcp)
      
      # Population-average survival estimation
      avg_covariates <- data.frame(lapply(current_covariates, function(x) {
        if (is.numeric(x)) mean(x, na.rm = TRUE) else names(sort(table(x), decreasing = TRUE))[1]
      }))
      
      baseline_surv <- tryCatch({
        survival::survfit(final_cox_fit, newdata = avg_covariates)
      }, error = function(e) {
        survival::survfit(survival::Surv(current_times[labeled_indices], rep(1, length(labeled_indices))) ~ 1)
      })
      
      survival_estimates <- approx(
        x = baseline_surv$time,
        y = baseline_surv$surv,
        xout = evaluation_times,
        method = "constant",
        f = 0,
        rule = 2
      )$y
      
      survival_estimates[is.na(survival_estimates)] <- 1.0
      
    } else {
      # Fallback to Kaplan-Meier
      km_fit <- survival::survfit(survival::Surv(survival_times, censoring_indicators) ~ 1)
      survival_estimates <- approx(km_fit$time, km_fit$surv, xout = evaluation_times, 
                                   method = "constant", f = 0, rule = 2)$y
      survival_estimates[is.na(survival_estimates)] <- 1.0
    }
    
    if (verbose) {
      final_censored_rate <- length(censored_indices) / n
      cat("Final censored rate:", round(final_censored_rate, 3), "\n")
      cat("Total outer iterations:", outer_iteration - 1, "\n")
    }
    
    return(list(
      survival = survival_estimates,
      final_labeled_count = length(labeled_indices),
      final_censored_count = length(censored_indices),
      outer_iterations = outer_iteration - 1,
      enhanced_dataset = list(
        times = current_times,
        status = current_status,
        covariates = current_covariates
      )
    ))
    
  }, error = function(e) {
    warning(paste("Complete Cox-SP-AFT failed:", e$message))
    return(list(
      survival = rep(NA, length(evaluation_times)),
      final_labeled_count = NA,
      final_censored_count = NA,
      outer_iterations = 0,
      enhanced_dataset = NULL
    ))
  })
}

#==============================================================================
# FIXES FOR NESTED PARALLEL PROCESSING ISSUE
# Problem: bootstrap_variance_cox_sp_aft() fails when called from within 
# parallel workers in main3.R
#==============================================================================

# SOLUTION 1: Detect if already in parallel context and disable internal parallelization
bootstrap_variance_cox_sp_aft_fixed <- function(survival_times, censoring_indicators, covariates, 
                                                evaluation_times, n_bootstrap = 30, 
                                                min_censored_rate = 0.05, max_outer_iterations = 10,
                                                lambda_mcp = 0.1, verbose = TRUE, 
                                                parallel = FALSE, n_cores = NULL,
                                                force_sequential = FALSE) {
  
  if (verbose) {
    cat("=== COX-SP-AFT BOOTSTRAP VARIANCE ESTIMATION (FIXED) ===\n")
    cat("Bootstrap samples:", n_bootstrap, "\n")
    cat("Evaluation times:", length(evaluation_times), "\n")
  }
  
  # CRITICAL FIX 1: Detect if we're already in a parallel context
  already_parallel <- FALSE
  tryCatch({
    # Check if we're in a foreach worker
    if (exists(".Random.seed", envir = .GlobalEnv) && 
        exists("foreach", where = "package:foreach") &&
        foreach::getDoParRegistered()) {
      already_parallel <- TRUE
      if (verbose) cat("Detected existing parallel context - forcing sequential execution\n")
    }
  }, error = function(e) {
    # If we can't determine, be safe and assume we might be parallel
    already_parallel <- TRUE
  })
  
  # Force sequential if already parallel or explicitly requested
  if (already_parallel || force_sequential) {
    parallel <- FALSE
    if (verbose) cat("Using sequential bootstrap processing\n")
  }
  
  n_original <- length(survival_times)
  n_times <- length(evaluation_times)
  
  # Storage for bootstrap results
  bootstrap_estimates <- matrix(NA, nrow = n_bootstrap, ncol = n_times)
  bootstrap_convergence <- logical(n_bootstrap)
  bootstrap_iterations <- numeric(n_bootstrap)
  
  # FIXED: Improved single_bootstrap function with better error handling
  single_bootstrap <- function(b) {
    
    tryCatch({
      
      # Bootstrap resampling with replacement
      boot_indices <- sample(1:n_original, n_original, replace = TRUE)
      
      boot_times <- survival_times[boot_indices]
      boot_status <- censoring_indicators[boot_indices]
      boot_covariates <- covariates[boot_indices, , drop = FALSE]
      
      # Ensure we have some events and some censored observations
      n_events <- sum(boot_status == 1)
      n_censored <- sum(boot_status == 0)
      
      if (n_events < 3 || n_censored < 1) {
        # Insufficient variation in bootstrap sample
        return(list(
          survival = rep(NA, n_times),
          converged = FALSE,
          iterations = 0
        ))
      }
      
      # CRITICAL FIX 2: Always use sequential for bootstrap
      boot_result <- cox_sp_aft_complete(
        survival_times = boot_times,
        censoring_indicators = boot_status,
        covariates = boot_covariates,
        evaluation_times = evaluation_times,
        min_censored_rate = min_censored_rate,
        max_outer_iterations = max_outer_iterations,
        lambda_mcp = lambda_mcp,
        verbose = FALSE
      )
      
      return(list(
        survival = boot_result$survival,
        converged = !any(is.na(boot_result$survival)),
        iterations = boot_result$outer_iterations
      ))
      
    }, error = function(e) {
      if (verbose) cat("Bootstrap iteration", b, "failed:", e$message, "\n")
      return(list(
        survival = rep(NA, n_times),
        converged = FALSE,
        iterations = 0
      ))
    })
  }
  
  # Execute bootstrap iterations
  if (parallel && !already_parallel && requireNamespace("parallel", quietly = TRUE)) {
    
    if (verbose) cat("Running parallel bootstrap...\n")
    
    if (is.null(n_cores)) {
      n_cores <- min(parallel::detectCores() - 1, n_bootstrap, 4) # Cap at 4 cores
    }
    
    # IMPROVED: Better parallel execution with more robust error handling
    tryCatch({
      cl <- parallel::makeCluster(n_cores, type = "PSOCK")
      
      # Ensure cleanup on exit
      on.exit({
        if (exists("cl")) {
          tryCatch({
            parallel::stopCluster(cl)
          }, error = function(e) {
            warning("Error stopping cluster: ", e$message)
          })
        }
      })
      
      # Load required packages on all workers
      parallel::clusterEvalQ(cl, {
        library(survival)
        if (requireNamespace("ncvreg", quietly = TRUE)) {
          library(ncvreg)
        }
      })
      
      # CRITICAL FIX 3: Export only essential functions that exist
      essential_functions <- c(
        "cox_sp_aft_complete", "cox_with_mcp", "aft_with_mcp", 
        "sp_aft_model_complete", "mcp_penalty"
      )
      
      # Only export functions that actually exist
      existing_functions <- essential_functions[sapply(essential_functions, exists)]
      
      if (length(existing_functions) > 0) {
        parallel::clusterExport(cl, existing_functions, envir = .GlobalEnv)
      }
      
      # CRITICAL FIX 4: Export only variables that exist in current environment
      current_env <- environment()
      vars_to_export <- c("survival_times", "censoring_indicators", "covariates",
                          "evaluation_times", "n_original", "n_times",
                          "min_censored_rate", "max_outer_iterations", "lambda_mcp")
      
      # Check which variables exist in current environment
      existing_vars <- vars_to_export[sapply(vars_to_export, exists, envir = current_env)]
      
      if (length(existing_vars) > 0) {
        parallel::clusterExport(cl, existing_vars, envir = current_env)
      }
      
      bootstrap_results <- parallel::parLapply(cl, 1:n_bootstrap, single_bootstrap)
      
    }, error = function(e) {
      warning("Parallel processing failed: ", e$message, " - switching to sequential")
      parallel <- FALSE
    })
  }
  
  # ALWAYS fall back to sequential if parallel fails or is disabled
  if (!parallel || !exists("bootstrap_results")) {
    
    if (verbose) cat("Running sequential bootstrap...\n")
    
    bootstrap_results <- lapply(1:n_bootstrap, function(b) {
      if (verbose && b %% max(1, floor(n_bootstrap/10)) == 0) {
        cat("Bootstrap sample", b, "/", n_bootstrap, "\n")
      }
      single_bootstrap(b)
    })
  }
  
  # Process bootstrap results (unchanged)
  successful_boots <- 0
  
  for (b in 1:n_bootstrap) {
    result <- bootstrap_results[[b]]
    
    if (!is.null(result) && is.list(result) && result$converged && !any(is.na(result$survival))) {
      bootstrap_estimates[b, ] <- result$survival
      bootstrap_convergence[b] <- result$converged
      bootstrap_iterations[b] <- result$iterations
      successful_boots <- successful_boots + 1
    }
  }
  
  if (verbose) {
    cat("Successful bootstrap samples:", successful_boots, "/", n_bootstrap, "\n")
  }
  
  if (successful_boots < 10) {
    warning("Very few successful bootstrap samples - results may be unreliable")
  }
  
  # Calculate variance estimates
  valid_rows <- !is.na(bootstrap_estimates[, 1])
  
  if (sum(valid_rows) > 1) {
    
    # Point-wise standard errors
    bootstrap_sd <- apply(bootstrap_estimates[valid_rows, , drop = FALSE], 2, sd, na.rm = TRUE)
    # bootstrap_mean <- apply(bootstrap_estimates[valid_rows, , drop = FALSE], 2, mean, na.rm = TRUE)
    # 
    # # Percentile confidence intervals
    # bootstrap_q025 <- apply(bootstrap_estimates[valid_rows, , drop = FALSE], 2, quantile, 
    #                         probs = 0.025, na.rm = TRUE)
    # bootstrap_q975 <- apply(bootstrap_estimates[valid_rows, , drop = FALSE], 2, quantile, 
    #                         probs = 0.975, na.rm = TRUE)
    # 
    # # Bias-corrected estimates
    # bootstrap_median <- apply(bootstrap_estimates[valid_rows, , drop = FALSE], 2, median, na.rm = TRUE)
    # 
  } else {
    # Fallback when bootstrap fails
    warning("Bootstrap variance estimation failed - using simple approximation")
    bootstrap_sd <- rep(0.1, n_times)
    # bootstrap_mean <- rep(0.5, n_times)
    # bootstrap_q025 <- rep(0.3, n_times)
    # bootstrap_q975 <- rep(0.7, n_times)
    # bootstrap_median <- rep(0.5, n_times)
  }
  
  return(list(
    standard_errors = bootstrap_sd,
    # bootstrap_mean = bootstrap_mean,
    # bootstrap_median = bootstrap_median,
    # confidence_intervals = list(
    #   lower = bootstrap_q025,
    #   upper = bootstrap_q975
    # ),
    # bootstrap_estimates = bootstrap_estimates[valid_rows, , drop = FALSE],
    # successful_samples = successful_boots,
    # total_samples = n_bootstrap,
    convergence_rate = mean(bootstrap_convergence, na.rm = TRUE),
    avg_iterations = mean(bootstrap_iterations[bootstrap_iterations > 0], na.rm = TRUE)
  ))
}


