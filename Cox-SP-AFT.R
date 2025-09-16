
#==============================================================================

# MCP penalty function
mcp_penalty <- function(beta, lambda_mcp, gamma = 3) {
  abs_beta <- abs(beta)
  penalty <- ifelse(abs_beta <= gamma * lambda_mcp,
                    lambda_mcp * abs_beta - beta^2 / (2 * gamma),
                    gamma * lambda_mcp^2 / 2)
  return(sum(penalty))
}

# Cox model with MCP regularization
cox_with_mcp <- function(cox_data, cox_formula, lambda_mcp = 0.1) {
  tryCatch({
    if (requireNamespace("ncvreg", quietly = TRUE)) {
      X <- model.matrix(cox_formula, data = cox_data)[, -1, drop = FALSE]
      y <- with(cox_data, survival::Surv(time, status))
      
      # Fit Cox model with MCP penalty
      mcp_fit <- ncvreg::ncvsurv(X = X, y = y, penalty = "MCP")
      
      # Select lambda via cross-validation
      cv_fit <- ncvreg::cv.ncvsurv(X = X, y = y, penalty = "MCP")
      best_coefs <- coef(mcp_fit, lambda = cv_fit$lambda.min)
      
      # Create coxph-like object
      regular_fit <- survival::coxph(cox_formula, data = cox_data)
      regular_fit$coefficients <- best_coefs
      
      return(regular_fit)
    } else {
      warning("ncvreg package not available - using standard Cox model")
      return(survival::coxph(cox_formula, data = cox_data))
    }
  }, error = function(e) {
    return(survival::coxph(cox_formula, data = cox_data))
  })
}

# AFT model with MCP regularization
aft_with_mcp <- function(aft_data, aft_formula, weight, lambda_mcp = 0.1) {
  tryCatch({
    if (requireNamespace("ncvreg", quietly = TRUE)) {
      X <- model.matrix(aft_formula, data = aft_data)[, -1, drop = FALSE]
      y <- aft_data$log_time
      
      if (length(weight) == nrow(X)) {
        # Fit AFT model with MCP penalty and weights
        mcp_fit <- ncvreg::ncvreg(X = X, y = y, penalty = "MCP", weights = weight)
        
        # Select lambda via cross-validation
        cv_fit <- ncvreg::cv.ncvreg(X = X, y = y, penalty = "MCP", weights = weight)
        best_lambda <- cv_fit$lambda.min
        
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

# Kaplan-Meier initialization for Cox model scenarios
enhanced_km_initialization <- function(y_initial, status_initial, censored_indices, 
                                       x_data, model_type = "cox") {
  m <- length(censored_indices)
  if (m == 0) return(numeric(0))
  
  Y <- numeric(m)
  
  tryCatch({
    # Fit Kaplan-Meier
    km_fit <- survival::survfit(survival::Surv(y_initial, status_initial) ~ 1)
    
    # For Cox model scenarios, use more conservative estimation
    if (model_type == "cox") {
      # Use more aggressive extrapolation for Cox models
      multiplier_base <- 2.0
      min_extension <- 1.5
    } else {
      # Original approach for logistic models
      multiplier_base <- 1.5
      min_extension <- 1.2
    }
    
    for (i in 1:m) {
      idx <- censored_indices[i]
      censor_time <- y_initial[idx]
      
      # Get survival probability at censoring time
      surv_prob <- approx(km_fit$time, km_fit$surv, xout = censor_time, 
                          method = "constant", f = 0, rule = 2)$y
      
      if (is.na(surv_prob)) surv_prob <- 0.5
      
      # More sophisticated time estimation based on survival curve
      if (surv_prob > 0.2) {
        # Find time where survival drops to a lower level
        target_surv <- max(0.1, surv_prob * 0.7)
        est_time <- approx(km_fit$surv, km_fit$time, xout = target_surv, 
                           method = "linear", rule = 2)$y
        
        if (!is.na(est_time) && est_time > censor_time) {
          Y[i] <- est_time
        } else {
          # Use survival probability to determine multiplier
          multiplier <- multiplier_base + (1 - surv_prob) * 2
          Y[i] <- censor_time * multiplier
        }
      } else {
        # Low survival probability - use larger multiplier
        multiplier <- multiplier_base + 3
        Y[i] <- censor_time * multiplier
      }
      
      # Ensure minimum extension
      Y[i] <- max(Y[i], censor_time * min_extension)
    }
    
  }, error = function(e) {
    # Fallback initialization with model-specific multipliers
    if (model_type == "cox") {
      Y <- y_initial[censored_indices] * 3  # More aggressive for Cox
    } else {
      Y <- y_initial[censored_indices] * 2  # Original for logistic
    }
  })
  
  return(Y)
}

# SP-AFT model with enhanced self-paced learning (FIXED VERSION)
sp_aft_model_complete_fixed <- function(x_data, y_initial, status_initial, censored_indices, 
                                        alpha_start = 0.05, mu = 1.3, max_iterations = 30, 
                                        tolerance = 1e-4, lambda_mcp = 0.1, model_type = "cox") {
  
  m <- length(censored_indices)
  if (m == 0) return(list(beta = NULL, labeled_times = NULL, reliable_indices = NULL))
  
  # Step 1: Initialize weights and age parameter
  V <- rep(1, m)  # Initial weights (all reliable)
  lambda_param <- alpha_start
  
  # Step 1: Enhanced initialization with model-specific approach
  Y <- enhanced_km_initialization(y_initial, status_initial, censored_indices, 
                                  x_data, model_type)
  
  # Track convergence history
  convergence_history <- list()
  
  ################################
  # Main SP-AFT iteration loop
  for (iteration in 1:max_iterations) {
    # Step 3: Update β using AFT model with MCP regularization
    complete_indices <- setdiff(1:length(y_initial), censored_indices)
    
    # Combine complete data with current estimates
    combined_x <- rbind(x_data[complete_indices, , drop = FALSE], 
                        x_data[censored_indices, , drop = FALSE])
    combined_y <- c(log(pmax(y_initial[complete_indices], 0.001)), 
                    log(pmax(Y, 0.001)))
    combined_weights <- c(rep(1, length(complete_indices)), V)
    
    # Fit AFT model
    aft_data <- data.frame(log_time = combined_y, combined_x)
    aft_formula <- as.formula(paste("log_time ~", paste(names(x_data), collapse = " + ")))
    beta_fit <- aft_with_mcp(aft_data = aft_data, aft_formula = aft_formula, 
                             weight = combined_weights, lambda_mcp = lambda_mcp)
    
    # Step 4: Update V according to Eq.11
    prev_V <- V
    reliable_count <- 0
    
    for (j in 1:m) {
      idx <- censored_indices[j]
      censor_time <- y_initial[idx]
      
      # Predict survival time using current β
      newx <- data.frame(x_data[idx, , drop = FALSE])
      colnames(newx) <- names(x_data)
      pred_log <- if ("ncvreg" %in% class(beta_fit$fit)) {
        predict(beta_fit$fit, X = as.matrix(newx), lambda = beta_fit$lambda, type = "response")[1]
      } else {
        predict(beta_fit$fit, newdata = newx, type = "response")[1]
      }
      pred_time <- exp(pred_log)
      
      # Calculate loss according to paper Eq.8 (FIXED)
      if (pred_time < censor_time) {
        # Constraint violation - assign large finite loss
        loss_j <- lambda_param * 10  # Large but finite penalty
      } else {
        # Original time scale loss (not log scale)
        loss_j <- (Y[j] - pred_time)^2
      }
      
      # Update weight V_j according to Eq.11 (CORRECTED)
      if (loss_j <= lambda_param) {
        V[j] <- 1  # Reliable sample
        reliable_count <- reliable_count + 1
      } else {
        V[j] <- 0  # Unreliable sample
      }
    }
    
    # Step 5: Update Y according to Eq.12 (ENHANCED)
    for (j in 1:m) {
      if (V[j] == 1) {  # Only update reliable samples
        idx <- censored_indices[j]
        newx <- data.frame(x_data[idx, , drop = FALSE])
        colnames(newx) <- names(x_data)
        
        pred_log <- if ("ncvreg" %in% class(beta_fit$fit)) {
          predict(beta_fit$fit, X = as.matrix(newx), lambda = beta_fit$lambda, type = "response")[1]
        } else {
          predict(beta_fit$fit, newdata = newx, type = "response")[1]
        }
        pred_time <- exp(pred_log)
        
        # Ensure constraint satisfaction
        if (pred_time >= y_initial[idx]) {
          Y[j] <- pred_time
        } else {
          # Keep current estimate if constraint would be violated
          V[j] <- 0  # Mark as unreliable
        }
      }
      # Unreliable samples keep their current Y[j]
    }
    
    # Step 6: Update age parameter λ (CRITICAL FIX)
    lambda_param <- mu * lambda_param
    
    # Enhanced convergence check
    weight_change <- mean(abs(V - prev_V))
    
    # Store convergence info
    convergence_history[[iteration]] <- list(
      weight_change = weight_change,
      reliable_count = reliable_count,
      lambda_param = lambda_param
    )
    
    # Multiple stopping criteria
    if (weight_change < tolerance) {
      break
    }
    
    # Stop if no reliable samples
    if (reliable_count == 0) {
      break
    }
    
    # Stop if all samples are reliable and stable
    if (reliable_count == m && weight_change < tolerance * 10) {
      break
    }
  }
  
  # Return reliable labeled samples
  reliable_indices <- censored_indices[V > 0]
  reliable_times <- Y[V > 0]
  
  return(list(
    beta_fit = beta_fit,
    labeled_times = reliable_times,
    reliable_indices = reliable_indices,
    final_weights = V,
    convergence_history = convergence_history,
    final_iteration = iteration
  ))
}

# Main Cox-SP-AFT algorithm (FIXED VERSION)
cox_sp_aft_complete <- function(survival_times, censoring_indicators, covariates, 
                                      evaluation_times, min_censored_rate = 0.05, 
                                      max_outer_iterations = 10, lambda_mcp = 0.1, 
                                      model_type = "cox", verbose = FALSE) {
  tryCatch({
    if (verbose) cat("=== FIXED COX-SP-AFT MODEL ===\n")
    
    # Check for required packages
    if (!requireNamespace("ncvreg", quietly = TRUE) && verbose) {
      cat("Warning: ncvreg package not available - using approximations\n")
    }
    
    # Initialize
    n <- length(survival_times)
    current_times <- survival_times
    current_status <- censoring_indicators
    current_covariates <- covariates
    
    labeled_indices <- which(current_status == 1)
    censored_indices <- which(current_status == 0)
    
    # Determine model-specific parameters
    if (model_type == "cox") {
      # Parameters optimized for Cox model scenarios
      sp_aft_alpha <- 0.05
      #sp_aft_mu <- 1.3
      sp_aft_mu <- 1.1 #bias smaller
      # sp_aft_iterations <- 30
      sp_aft_iterations <- 20
      
    } else {
      # Parameters for logistic model scenarios
      #sp_aft_alpha <- 0.1
      sp_aft_alpha <- 0.05 #bias smaller
      sp_aft_mu <- 1.1
      sp_aft_iterations <- 20
    }
    
    if (verbose) {
      cat("Initial - Labeled:", length(labeled_indices), "Censored:", length(censored_indices), "\n")
      cat("Using model-specific parameters for:", model_type, "\n")
    }
    
    outer_iteration <- 1
    
    # Main Cox-SP-AFT Loop
    while (length(censored_indices) / n > min_censored_rate && 
           outer_iteration <= max_outer_iterations) {
      
      censored_rate <- length(censored_indices) / n
      if (verbose) {
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
        
        # Fit Cox model with MCP
        cox_fit <- cox_with_mcp(cox_data, cox_formula, lambda_mcp)
        
        # Get risk scores for all patients
        all_cox_data <- data.frame(current_covariates)
        risk_scores <- predict(cox_fit, newdata = all_cox_data, type = "lp")
        
        # Use more robust risk group classification
        risk_threshold <- median(risk_scores, na.rm = TRUE)
        high_risk_indices <- which(risk_scores > risk_threshold)
        low_risk_indices <- which(risk_scores <= risk_threshold)
        
        if (verbose) {
          cat("Risk groups - High:", length(high_risk_indices), "Low:", length(low_risk_indices), "\n")
        }
      } else {
        high_risk_indices <- 1:n
        low_risk_indices <- integer(0)
      }
      
      # Step 3: SP-AFT for each risk group (FIXED PARAMETERS)
      newly_labeled_indices <- integer(0)
      newly_labeled_times <- numeric(0)
      
      # Process high-risk group
      if (length(high_risk_indices) > 0) {
        group_censored <- intersect(high_risk_indices, censored_indices)
        
        if (length(group_censored) >= 1) {
          group_covariates <- current_covariates[high_risk_indices, , drop = FALSE]
          group_times <- current_times[high_risk_indices]
          group_status <- current_status[high_risk_indices]
          group_censored_relative <- match(group_censored, high_risk_indices)
          
          # Run SP-AFT with model-specific parameters
          sp_aft_result <- sp_aft_model_complete_fixed(
            x_data = group_covariates,
            y_initial = group_times,
            status_initial = group_status,
            censored_indices = group_censored_relative,
            alpha_start = sp_aft_alpha,
            mu = sp_aft_mu,
            max_iterations = sp_aft_iterations,
            lambda_mcp = lambda_mcp,
            model_type = model_type
          )
          
          if (length(sp_aft_result$reliable_indices) > 0) {
            global_reliable_indices <- high_risk_indices[sp_aft_result$reliable_indices]
            newly_labeled_indices <- c(newly_labeled_indices, global_reliable_indices)
            newly_labeled_times <- c(newly_labeled_times, sp_aft_result$labeled_times)
            
            if (verbose) {
              cat("    High-risk group: added", length(sp_aft_result$reliable_indices), "samples\n")
            }
          }
        }
      }
      
      # Process low-risk group
      if (length(low_risk_indices) > 0) {
        group_censored <- intersect(low_risk_indices, censored_indices)
        
        if (length(group_censored) >= 1) {
          group_covariates <- current_covariates[low_risk_indices, , drop = FALSE]
          group_times <- current_times[low_risk_indices]
          group_status <- current_status[low_risk_indices]
          group_censored_relative <- match(group_censored, low_risk_indices)
          
          # Run SP-AFT with model-specific parameters
          sp_aft_result <- sp_aft_model_complete_fixed(
            x_data = group_covariates,
            y_initial = group_times,
            status_initial = group_status,
            censored_indices = group_censored_relative,
            alpha_start = sp_aft_alpha,
            mu = sp_aft_mu,
            max_iterations = sp_aft_iterations,
            lambda_mcp = lambda_mcp,
            model_type = model_type
          )
          
          if (length(sp_aft_result$reliable_indices) > 0) {
            global_reliable_indices <- low_risk_indices[sp_aft_result$reliable_indices]
            newly_labeled_indices <- c(newly_labeled_indices, global_reliable_indices)
            newly_labeled_times <- c(newly_labeled_times, sp_aft_result$labeled_times)
            
            if (verbose) {
              cat("    Low-risk group: added", length(sp_aft_result$reliable_indices), "samples\n")
            }
          }
        }
      }
      
      # Step 4: Update training dataset with newly labeled samples
      if (length(newly_labeled_indices) > 0) {
        current_times[newly_labeled_indices] <- newly_labeled_times
        current_status[newly_labeled_indices] <- 1
        
        labeled_indices <- which(current_status == 1)
        censored_indices <- which(current_status == 0)
        
        if (verbose) {
          cat("  Total newly labeled:", length(newly_labeled_indices), "\n")
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
      
      # Interpolate survival estimates at evaluation times
      survival0_estimates <- approx(
        x = baseline_surv$time,
        y = baseline_surv$surv,
        xout = evaluation_times,
        method = "constant",
        f = 0,
        rule = 2
      )$y
      
      survival0_estimates[is.na(survival0_estimates)] <- 1.0
      
      # Get linear predictors (risk scores) for all patients
      risk_scores <- predict(final_cox_fit, type = "lp")
      
      # Apply individual risk scores (vectorized operation)
      n_patients <- length(risk_scores)
      n_times <- length(evaluation_times)
      
      # Create matrix: rows = patients, cols = time points
      survival_matrix <- matrix(0, nrow = n_patients, ncol = n_times)
      
      for (i in 1:n_patients) {
        # S(t|x_i) = S0(t)^exp(risk_score_i)
        survival_matrix[i, ] <- survival0_estimates^exp(risk_scores[i])
      }
      survival_estimates <- colMeans(survival_matrix)
      
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
    warning(paste("Cox-SP-AFT failed:", e$message))
    return(list(
      survival = rep(NA, length(evaluation_times)),
      final_labeled_count = NA,
      final_censored_count = NA,
      outer_iterations = 0,
      enhanced_dataset = NULL
    ))
  })
}

