# ==============================================================================
# Classical Supervised Learning for Doubly-Censored Survival Data
# ==============================================================================
# 
# This file implements classical supervised learning methods for survival analysis
# with doubly-censored data. The approach uses only labeled
# data and creates separate estimators for different types of observations:
# exact observed, left-censored, and right-censored.
#
# The method works by:
# 1. Creating separate survival estimators for each observation type
# 2. Using kernel smoothing for censored observations
# 3. Combining estimators optimally using covariance information
# 4. Providing variance estimation and confidence intervals
#
# Functions included:
# - Supervised_est(): Main classical supervised learning function
#
# ==============================================================================

#' Classical Supervised Learning for Doubly-Censored Survival Data
#' 
#' This function implements classical supervised learning estimation for survival
#' analysis with doubly-censored data. It creates separate estimators for exact,
#' left-censored, and right-censored observations, then combines them optimally
#' using their estimated covariance structure.
#' 
#' @param time Numeric vector, evaluation time points for survival estimation
#' @param obse Numeric vector, observed exact survival times
#' @param ldelt Numeric vector, left censoring indicators (1 = left-censored, 0 = not)
#' @param rdelt Numeric vector, right censoring indicators (1 = right-censored, 0 = not)
#' @param lcen Numeric vector, left censoring times
#' @param rcen Numeric vector, right censoring times
#' @param lh Numeric scalar, bandwidth for left censoring kernel
#' @param rh Numeric scalar, bandwidth for right censoring kernel
#' 
#' @return Matrix with 8 rows and length(time) columns containing:
#'   \item{Row 1}{Direct estimator survival probabilities}
#'   \item{Row 2}{Direct estimator standard deviations}
#'   \item{Row 3}{Left estimator survival probabilities}
#'   \item{Row 4}{Left estimator standard deviations}
#'   \item{Row 5}{Right estimator survival probabilities}
#'   \item{Row 6}{Right estimator standard deviations}
#'   \item{Row 7}{Combined estimator survival probabilities}
#'   \item{Row 8}{Combined estimator standard deviations}
#'   
#' @details
#' The classical supervised learning approach creates three types of estimators:
#' 
#' **1. Direct Estimator (for exact observations)**
#' Uses only exact observations within the observation window:
#' \deqn{S_D(t) = \frac{\sum_{i} I(L_i < t \leq X_i)}{\sum_{i} I(L_i < t \leq R_i)}}
#' where the numerator counts exact deaths by time t and the denominator 
#' counts individuals at risk at time t.
#' 
#' **2. Left Estimator (for left-censored observations)**
#' Uses kernel smoothing around left censoring times:
#' \deqn{S_L(t) = \frac{\sum_{i} w_L(L_i - t) (1 - \delta_{L,i})}{\sum_{i} w_L(L_i - t)}}
#' where \eqn{w_L(x)} is a kernel function and \eqn{\delta_{L,i}} indicates left censoring.
#' 
#' **3. Right Estimator (for right-censored observations)**
#' Uses kernel smoothing around right censoring times:
#' \deqn{S_R(t) = \frac{\sum_{i} w_R(R_i - t) (1 - \delta_{R,i})}{\sum_{i} w_R(R_i - t)}}
#' 
#' **4. Combined Estimator**
#' Optimally combines the three estimators using inverse covariance weighting:
#' \deqn{S_C(t) = w_D S_D(t) + w_L S_L(t) + w_R S_R(t)}
#' where weights satisfy \eqn{w_D + w_L + w_R = 1} and minimize the variance
#' \eqn{w^T \Sigma w} subject to the constraint.
#' 
#' **Variance Estimation**
#' Uses the influence function approach and kernel variance formulas:
#' - Direct: \eqn{Var(S_D) = Var(I(exact)) / n}
#' - Left/Right: \eqn{Var(S_{L/R}) = \int K^2(u) du \cdot Var(residuals) / (nh)}
#' 
#' **Covariance Estimation** 
#' Estimates cross-covariances between different estimator types using
#' the sample covariance of their influence functions.
#' 
#' @references
#' Wang, Y., Zhou, Q., Cai, T. & Wang, X. (2023). Semi-supervised Estimation 
#' of Event Rate with Doubly-censored Survival Data.
#' 
#' Turnbull, B. W. (1976). The empirical distribution function with arbitrarily
#' grouped, censored and truncated data. Journal of the Royal Statistical Society, 38, 290-295.
#' 
#' @seealso \code{\link{IntrSSL_est}} for semi-supervised approach
Supervised_est <- function(time, obse, ldelt, rdelt, lcen, rcen, lh, rh) {
  
  # ==============================================================================
  # INPUT VALIDATION
  # ==============================================================================
  
  if (length(time) == 0) stop("time vector cannot be empty")
  
  # Check dimensions match
  n_obs <- length(obse)
  if (length(ldelt) != n_obs || length(rdelt) != n_obs || 
      length(lcen) != n_obs || length(rcen) != n_obs) {
    stop("All input vectors must have the same length")
  }
  
  # Check censoring indicators are binary
  if (!all(ldelt %in% c(0, 1)) || !all(rdelt %in% c(0, 1))) {
    stop("Censoring indicators must be binary (0/1)")
  }
  
  # Check bandwidth parameters
  if (lh <= 0 || rh <= 0) {
    stop("Bandwidth parameters must be positive")
  }
  
  # Check time ordering
  if (any(lcen > rcen)) {
    stop("Left censoring times must be <= right censoring times")
  }
  
  if (any(obse < lcen | obse > rcen)) {
    stop("Observed times must be within [lcen, rcen] intervals")
  }
  
  n_times <- length(time)
  cat(sprintf("Classical Supervised Learning: %d observations, %d time points\n", 
              n_obs, n_times))
  cat(sprintf("Left-censored: %d, Right-censored: %d, Exact: %d\n",
              sum(ldelt), sum(rdelt), sum(1 - ldelt - rdelt + ldelt*rdelt)))
  
  # ==============================================================================
  # MAIN ESTIMATION LOOP OVER TIME POINTS
  # ==============================================================================
  
  result_matrix <- sapply(time, function(t) {
    
    # cat(sprintf("  Processing time point %.3f\n", t))
    
    # ==========================================================================
    # STEP 1: DIRECT ESTIMATOR (Exact Observations)
    # ==========================================================================
    
    # Indicator for being at risk at time t: I(L < t <= R)
    Inter_obse <- as.numeric(t <= obse & t > lcen)
    
    # Indicator for observation window at time t: I(L < t <= R) 
    Inter_cen <- as.numeric(t <= rcen & t > lcen)
    
    # Direct estimator and variance
    if (sum(Inter_cen) == 0) {
      # No one at risk at this time point
      DSt_sl <- 0
      DSt_sl_var <- 0
    } else {
      # Survival probability = (# exact observations) / (# at risk)
      DSt_sl <- sum(Inter_obse) / sum(Inter_cen)
      
      # Variance using Kaplan-Meier type formula
      # Var(S) = Var(# events) / (# at risk)^2
      DSt_sl_var <- mean(Inter_cen^2 * (Inter_obse - DSt_sl)^2) / 
        (mean(Inter_cen)^2 * n_obs)
    }
    
    # ==========================================================================
    # STEP 2: LEFT ESTIMATOR (Left-Censored Information)
    # ==========================================================================
    
    # Compute kernel weights centered at time t
    lw <- ker(lcen - t, lh)
    
    # Normalize weights to have mean 1 (if not all zero)
    if (mean(lw) == 0) {
      lw <- rep(0, length(lw))
    } else {
      lw <- lw / mean(lw)
    }
    
    # Left estimator and variance
    if (sum(ldelt) == 0) {
      # No left-censored observations
      LSt_sl <- 1  # All survive past left censoring
      LSt_sl_var <- 0
    } else {
      # Survival probability = weighted average of survival indicators
      # (1 - ldelt) indicates survival past left censoring time
      LSt_sl <- sum(lw * (1 - ldelt)) / sum(lw)
      
      # Kernel-based variance estimate
      # Formula accounts for kernel smoothing: Var = integral(K^2) * Var(residuals) / (nh)
      LSt_sl_var <- lh * mean(lw^2 * ((1 - ldelt) - LSt_sl)^2) / 
        (mean(lw)^2 * lh * n_obs)
    }
    
    # ==========================================================================
    # STEP 3: RIGHT ESTIMATOR (Right-Censored Information)
    # ==========================================================================
    
    # Compute kernel weights centered at time t
    rw <- ker(rcen - t, rh)
    
    # Normalize weights
    if (mean(rw) == 0) {
      rw <- rep(0, length(rw))
    } else {
      rw <- rw / mean(rw)
    }
    
    # Right estimator and variance
    if (sum(rdelt) == 0) {
      # No right-censored observations
      RSt_sl <- 1  # All survive past right censoring
      RSt_sl_var <- 0
    } else {
      # Survival probability using right-censored information
      RSt_sl <- sum(rw * (1 - rdelt)) / sum(rw)
      
      # Kernel-based variance estimate
      RSt_sl_var <- rh * mean(rw^2 * ((1 - rdelt) - RSt_sl)^2) / 
        (mean(rw)^2 * rh * n_obs)
    }
    
    # ==========================================================================
    # STEP 4: COVARIANCE ESTIMATION
    # ==========================================================================
    
    # Covariance between direct and left estimators
    Cov_dl <- mean((Inter_cen * (Inter_obse - DSt_sl)) * 
                     (lw * ((1 - ldelt) - LSt_sl))) / 
      (mean(Inter_cen) * mean(lw) * n_obs)
    
    # Covariance between direct and right estimators
    Cov_dr <- mean((Inter_cen * (Inter_obse - DSt_sl)) * 
                     (rw * ((1 - rdelt) - RSt_sl))) / 
      (mean(Inter_cen) * mean(rw) * n_obs)
    
    # Covariance between left and right estimators
    Cov_lr <- mean((lw * ((1 - ldelt) - LSt_sl)) * 
                     (rw * ((1 - rdelt) - RSt_sl))) / 
      (mean(lw) * mean(rw) * n_obs)
    
    # Compile variance-covariance matrix
    St_sl_vec <- c(DSt_sl, LSt_sl, RSt_sl)
    Cov_sl <- matrix(c(DSt_sl_var, Cov_dl, Cov_dr, 
                       Cov_dl, LSt_sl_var, Cov_lr, 
                       Cov_dr, Cov_lr, RSt_sl_var), nrow = 3)
    
    # ==========================================================================
    # STEP 5: HANDLE MISSING ESTIMATORS
    # ==========================================================================
    
    # Check which estimator types have observations
    response_counts <- c(sum(Inter_obse), sum(ldelt), sum(rdelt))
    missing_estimators <- which(response_counts == 0)
    
    if (length(missing_estimators) > 0) {
      # Remove estimators with no information
      St_sl_vec <- St_sl_vec[-missing_estimators]
      Cov_sl <- as.matrix(Cov_sl[-missing_estimators, -missing_estimators])
    }
    
    # ==========================================================================
    # STEP 6: OPTIMAL COMBINATION
    # ==========================================================================
    
    # Try multiple regularization levels to find stable solution
    CSt_wall <- list()
    Diff <- numeric()
    lambda_vals <- round(seq(0.2, 1, by = 0.15), 3)
    
    for (i in seq_along(lambda_vals)) {
      tryCatch({
        # Regularization parameter
        epsilon <- n_obs^(-lambda_vals[i])
        
        # Regularized inverse covariance matrix
        Cov_inv <- solve(Cov_sl + epsilon * diag(ncol(Cov_sl)))
        
        # Optimal weights under sum-to-one constraint
        # Minimize w^T Sigma w subject to sum(w) = 1
        ones <- rep(1, ncol(Cov_sl))
        denominator <- as.vector(t(ones) %*% Cov_inv %*% ones)
        
        if (abs(denominator) > 1e-10) {
          weights <- t(ones) %*% Cov_inv / denominator
          
          # Objective function for weight selection
          # Prefer weights that deviate more from uniform (better discrimination)
          if (length(weights) >= 2) {
            weight_diff <- (weights[1] - weights[2])^2
            if (length(weights) >= 3) {
              weight_diff <- weight_diff + (weights[3] - weights[2])^2
            }
          } else {
            weight_diff <- 0
          }
          
          Diff <- c(Diff, weight_diff)
          CSt_wall <- c(CSt_wall, list(as.vector(weights)))
        } else {
          # Fallback to uniform weights
          Diff <- c(Diff, 0)
          CSt_wall <- c(CSt_wall, list(rep(1/length(St_sl_vec), length(St_sl_vec))))
        }
        
      }, error = function(e) {
        # Fallback to uniform weights on error
        Diff <- c(Diff, 0)
        CSt_wall <- c(CSt_wall, list(rep(1/length(St_sl_vec), length(St_sl_vec))))
      })
    }
    
    # Select best regularization parameter
    if (length(Diff) > 0 && max(Diff, na.rm = TRUE) > 0) {
      best_idx <- which.max(Diff)
      CSt_w <- CSt_wall[[best_idx]]
    } else {
      # Default to uniform weights
      CSt_w <- rep(1/length(St_sl_vec), length(St_sl_vec))
    }
    
    # Compute combined estimator
    CSt_sl <- as.vector(CSt_w %*% St_sl_vec)
    
    # ==========================================================================
    # STEP 7: COMBINED VARIANCE ESTIMATION
    # ==========================================================================
    
    # Map weights back to full 3-element vector
    CSt_wei <- rep(0, 3)
    if (length(missing_estimators) == 0) {
      CSt_wei <- CSt_w
    } else {
      CSt_wei[-missing_estimators] <- CSt_w
    }
    
    # Combined variance using influence function approach
    combined_influence <- (CSt_wei[1] * mean(Inter_cen)^(-1) * Inter_cen * (Inter_obse - DSt_sl) + 
                             CSt_wei[2] * mean(lw)^(-1) * lw * ((1 - ldelt) - LSt_sl) +
                             CSt_wei[3] * mean(rw)^(-1) * rw * ((1 - rdelt) - RSt_sl))
    
    CSt_sl_var <- mean(combined_influence^2) / n_obs
    
    # ==========================================================================
    # STEP 8: RETURN RESULTS FOR THIS TIME POINT
    # ==========================================================================
    
    return(c(
      'DSt_sl' = DSt_sl,                    # Direct survival estimate
      'DSt_sl_sd' = sqrt(max(0, DSt_sl_var)), # Direct standard deviation
      'LSt_sl' = LSt_sl,                    # Left survival estimate  
      'LSt_sl_sd' = sqrt(max(0, LSt_sl_var)), # Left standard deviation
      'RSt_sl' = RSt_sl,                    # Right survival estimate
      'RSt_sl_sd' = sqrt(max(0, RSt_sl_var)), # Right standard deviation
      'CSt_sl' = CSt_sl,                    # Combined survival estimate
      'CSt_sl_sd' = sqrt(max(0, CSt_sl_var))  # Combined standard deviation
    ))
    
  })  # End of sapply over time points
  
  # ==============================================================================
  # FINAL PROCESSING AND VALIDATION
  # ==============================================================================
  
  # Ensure result is a matrix
  if (!is.matrix(result_matrix)) {
    result_matrix <- matrix(result_matrix, nrow = 8)
  }
  
  # Add descriptive row names
  rownames(result_matrix) <- c(
    "Direct_Survival", "Direct_SD",
    "Left_Survival", "Left_SD", 
    "Right_Survival", "Right_SD",
    "Combined_Survival", "Combined_SD"
  )
  
  # Add column names with time points
  colnames(result_matrix) <- paste0("t_", round(time, 3))
  
  # Validation checks
  n_missing <- sum(is.na(result_matrix))
  if (n_missing > 0) {
    warning(sprintf("Classical supervised learning completed with %d missing values", n_missing))
  }
  
  # Check survival probability bounds
  surv_rows <- c(1, 3, 5, 7)  # Survival estimate rows
  surv_values <- result_matrix[surv_rows, ]
  bound_violations <- sum(surv_values < 0 | surv_values > 1, na.rm = TRUE)
  if (bound_violations > 0) {
    warning(sprintf("Survival probability bounds violated %d times", bound_violations))
  }
  
  # Check for negative variances
  sd_rows <- c(2, 4, 6, 8)  # Standard deviation rows
  sd_values <- result_matrix[sd_rows, ]
  negative_vars <- sum(sd_values < 0, na.rm = TRUE)
  if (negative_vars > 0) {
    warning(sprintf("Negative standard deviations detected %d times", negative_vars))
  }
  
  cat("Classical supervised learning estimation completed successfully\n")
  
  return(result_matrix)
}

