# ==============================================================================
# Cross-Validation for Semi-Supervised Survival Analysis
# ==============================================================================
# 
# This file implements cross-validation procedures for semi-supervised survival
# analysis. The main function estimates three types of survival estimators
# (Exact observed, Left, Right) and combines them optimally using cross-validation
# to minimize prediction error.
#
# The CV procedure:
# 1. Estimates parameters for each estimator type using labeled data
# 2. Computes survival estimates using both labeled and unlabeled data
# 3. Estimates variance-covariance matrix via cross-validation
# 4. Finds optimal linear combination weights
# 5. Provides final combined estimator with uncertainty quantification
#
# Functions included:
# - CV_function(): Main cross-validation function for combining estimators
#
# ==============================================================================
#' Cross-Validation for Semi-Supervised Survival Estimator Combination
#' 
#' This function implements cross-validation to optimally combine three types
#' of survival estimators (Exact observed, Left, Right) in a semi-supervised learning
#' framework. It estimates parameters, computes individual estimators, estimates
#' their variance-covariance structure, and finds optimal combination weights.
#' 
#' @param ddat_unlabel Numeric matrix, design matrix for exact observed estimator (unlabeled data)
#' @param ddat_label Numeric matrix, design matrix for exact observed estimator (labeled data)
#' @param ldat_unlabel Numeric matrix, design matrix for left estimator (unlabeled data)
#' @param ldat_label Numeric matrix, design matrix for left estimator (labeled data)
#' @param rdat_unlabel Numeric matrix, design matrix for right estimator (unlabeled data)
#' @param rdat_label Numeric matrix, design matrix for right estimator (labeled data)
#' @param dw_label Numeric vector, weights for exact observed estimator (labeled data)
#' @param dw_unlabel Numeric vector, weights for exact observed estimator (unlabeled data)
#' @param dy Numeric vector, binary response for exact observed estimator
#' @param dh Numeric scalar, bandwidth for exact observed estimator
#' @param lw_label Numeric vector, weights for left estimator (labeled data)
#' @param lw_unlabel Numeric vector, weights for left estimator (unlabeled data)
#' @param ly Numeric vector, binary response for left estimator
#' @param lh Numeric scalar, bandwidth for left estimator
#' @param rw_label Numeric vector, weights for right estimator (labeled data)
#' @param rw_unlabel Numeric vector, weights for right estimator (unlabeled data)
#' @param ry Numeric vector, binary response for right estimator
#' @param rh Numeric scalar, bandwidth for right estimator
#' @param num_folds Integer, number of cross-validation folds
#' @param a Integer, estimation type (0 = semi-supervised, 1 = intrinsic)
#' 
#' @return List containing:
#'   \item{DSt}{Exact observed survival estimate}
#'   \item{LSt}{Left survival estimate}
#'   \item{RSt}{Right survival estimate}
#'   \item{CSt}{Combined survival estimate}
#'   \item{D_sd}{Exact observed estimator standard deviation}
#'   \item{L_sd}{Left estimator standard deviation}
#'   \item{R_sd}{Right estimator standard deviation}
#'   \item{C_sd}{Combined estimator standard deviation}
#'   \item{wd}{Weight estimator for exact observed estimator in combination}
#'   \item{wl}{Weight estimator for left estimator in combination}
#'   \item{wr}{Weight estimator for right estimator in combination}
#'   
#' @details
#' The function implements a comprehensive cross-validation procedure:
#' 
#' **Step 1: Individual Estimator Computation**
#' For each estimator type (Exact observed, Left, Right):
#' \deqn{S_k = \frac{\sum_{i \in U} w_i^{(k)} \sigma(X_i^T \hat{\beta}_k)}{\sum_{i \in U} w_i^{(k)}}}
#' where k ∈ {D,L,R}, U is unlabeled data, and \eqn{\hat{\beta}_k} is estimated from labeled data.
#' 
#' **Step 2: Variance-Covariance Estimation**
#' Uses K-fold cross-validation to estimate:
#' \deqn{Var(S_k) = \frac{1}{K} \sum_{j=1}^K Var_j(S_k)}
#' \deqn{Cov(S_k, S_l) = \frac{1}{K} \sum_{j=1}^K Cov_j(S_k, S_l)}
#' 
#' **Step 3: Optimal Weight Selection**
#' Solves the constrained optimization:
#' \deqn{\min_w w^T \Sigma w \text{ subject to } \mathbf{1}^T w = 1}
#' where \eqn{\Sigma} is the estimated covariance matrix.
#' 
#' **Step 4: Combined Estimator**
#' \deqn{S_{combined} = w_D S_D + w_L S_L + w_R S_R}
#' 
#' The function handles missing estimators (when response sums to zero) by
#' excluding them from the combination and adjusting weights accordingly.
#' 
#' 
#' @seealso \code{\link{Beta_estimate}}, \code{\link{IntrSSL_est}}
#' 
CV_function <- function(ddat_unlabel, ddat_label, ldat_unlabel, ldat_label, 
                        rdat_unlabel, rdat_label, dw_label, dw_unlabel, dy, dh, 
                        lw_label, lw_unlabel, ly, lh, rw_label, rw_unlabel,
                        ry, rh, num_folds, a) {
  
  # ==============================================================================
  # INPUT VALIDATION
  # ==============================================================================
  
  if (num_folds < 2) stop("num_folds must be at least 2")
  if (!a %in% c(0, 1)) stop("a must be 0 (semi-supervised) or 1 (intrinsic)")
  
  # Get dimensions
  p_basis <- ncol(ddat_label)        # Number of covariates
  n_label <- nrow(ddat_label)        # Number of labeled observations
  n_unlabel <- nrow(ddat_unlabel)    # Number of unlabeled observations
  
  if (n_label < num_folds) {
    warning(sprintf("Fewer labeled observations (%d) than folds (%d), reducing folds", 
                    n_label, num_folds))
    num_folds <- max(2, n_label %/% 2)
  }
  
  cat(sprintf("CV: %d labeled, %d unlabeled, %d folds, method=%s\n", 
              n_label, n_unlabel, num_folds, ifelse(a == 0, "SSL", "Intrinsic")))
  
  # ==============================================================================
  # STEP 1: COMPUTE INDIVIDUAL SURVIVAL ESTIMATORS
  # ==============================================================================
  
  #' Internal function to compute survival estimator
  #' @param dat_unlabel Design matrix for unlabeled data
  #' @param dat_label Design matrix for labeled data  
  #' @param w_label Weights for labeled data
  #' @param w_unlabel Weights for unlabeled data
  #' @param y Response vector
  #' @param h Bandwidth parameter
  #' @return Survival probability estimate
  St_est <- function(dat_unlabel, dat_label, w_label, w_unlabel, y, h) {
    
    tryCatch({
      # Estimate parameters using labeled data
      Beta_n <- Beta_estimate(dat_label, y, h, w_label, w_unlabel)
      
      # Extract appropriate coefficients based on estimation type
      coef_start <- 1 + a * (p_basis + 1)
      coef_end <- (p_basis + 1) + a * (p_basis + 1)
      beta_n <- Beta_n[coef_start:coef_end]
      
      # Compute survival estimate using unlabeled data
      # Add intercept column to design matrix
      X_unlabel_aug <- cbind(1, dat_unlabel)
      
      # Compute predicted probabilities
      linear_pred <- as.vector(X_unlabel_aug %*% beta_n)
      prob_pred <- Expit(linear_pred)
      
      # Weighted average survival estimate
      St <- sum(w_unlabel * prob_pred) / sum(w_unlabel)
      
      return(St)
      
    }, error = function(e) {
      warning(paste("Survival estimation failed:", e$message))
      return(NA)
    })
  }
  
  # Compute individual estimators
  cat("  Computing individual estimators...\n")
  
  DSt <- St_est(ddat_unlabel, ddat_label, dw_label, dw_unlabel, dy, dh)
  LSt <- St_est(ldat_unlabel, ldat_label, lw_label, lw_unlabel, ly, lh)
  RSt <- St_est(rdat_unlabel, rdat_label, rw_label, rw_unlabel, ry, rh)
  
  cat(sprintf("  Individual estimates: D=%.4f, L=%.4f, R=%.4f\n", DSt, LSt, RSt))
  
  # ==============================================================================
  # STEP 2: CROSS-VALIDATION FOR VARIANCE ESTIMATION
  # ==============================================================================
  
  cat("  Running cross-validation for variance estimation...\n")
  
  # Create cross-validation folds
  ind_cv <- cv_split(n_label, num_folds)
  
  # Cross-validation loop
  Var_cv <- sapply(1:num_folds, function(k) {
    
    tryCatch({
      # Split data into train and validation
      inds_fold <- as.vector(unlist(ind_cv[-k]))  # Training indices
      inds_foldk <- as.vector(unlist(ind_cv[k]))  # Validation indices
      
      # === Exact observed Estimator ===
      if (sum(dy) > 0) {  # Only if there are positive responses
        Dbeta_K <- Beta_estimate(ddat_label[inds_fold, , drop = FALSE], 
                                 dy[inds_fold], dh, dw_label[inds_fold], dw_unlabel)
        coef_start <- 1 + a * (p_basis + 1)
        coef_end <- (p_basis + 1) + a * (p_basis + 1)
        dbeta_K <- Dbeta_K[coef_start:coef_end]
        
        # Predict on validation fold
        X_val_d <- cbind(1, ddat_label[inds_foldk, , drop = FALSE])
        DImp_k <- as.vector(Expit(X_val_d %*% dbeta_K))
        DVar_k <- dh * mean(dw_label[inds_foldk]^2 * (dy[inds_foldk] - DImp_k)^2) / (mean(dw_unlabel))^2
      } else {
        DImp_k <- rep(0, length(inds_foldk))
        DVar_k <- 0
      }
      
      # === Left Estimator ===
      if (sum(ly) > 0) {
        Lbeta_K <- Beta_estimate(ldat_label[inds_fold, , drop = FALSE], 
                                 ly[inds_fold], lh, lw_label[inds_fold], lw_unlabel)
        lbeta_K <- Lbeta_K[coef_start:coef_end]
        
        X_val_l <- cbind(1, ldat_label[inds_foldk, , drop = FALSE])
        LImp_k <- as.vector(Expit(X_val_l %*% lbeta_K))
        LVar_k <- lh * mean(lw_label[inds_foldk]^2 * (ly[inds_foldk] - LImp_k)^2) / (mean(lw_unlabel))^2
      } else {
        LImp_k <- rep(0, length(inds_foldk))
        LVar_k <- 0
      }
      
      # === Right Estimator ===
      if (sum(ry) > 0) {
        Rbeta_K <- Beta_estimate(rdat_label[inds_fold, , drop = FALSE], 
                                 ry[inds_fold], rh, rw_label[inds_fold], rw_unlabel)
        rbeta_K <- Rbeta_K[coef_start:coef_end]
        
        X_val_r <- cbind(1, rdat_label[inds_foldk, , drop = FALSE])
        RImp_k <- as.vector(Expit(X_val_r %*% rbeta_K))
        RVar_k <- rh * mean(rw_label[inds_foldk]^2 * (ry[inds_foldk] - RImp_k)^2) / (mean(rw_unlabel))^2
      } else {
        RImp_k <- rep(0, length(inds_foldk))
        RVar_k <- 0
      }
      
      # === Covariance Terms ===
      DLCov_k <- mean(dw_label[inds_foldk] * (dy[inds_foldk] - DImp_k) * 
                        lw_label[inds_foldk] * (ly[inds_foldk] - LImp_k)) / 
        (mean(dw_unlabel) * mean(lw_unlabel))
      
      DRCov_k <- mean(dw_label[inds_foldk] * (dy[inds_foldk] - DImp_k) * 
                        rw_label[inds_foldk] * (ry[inds_foldk] - RImp_k)) / 
        (mean(dw_unlabel) * mean(rw_unlabel))
      
      LRCov_k <- mean(lw_label[inds_foldk] * (ly[inds_foldk] - LImp_k) * 
                        rw_label[inds_foldk] * (ry[inds_foldk] - RImp_k)) / 
        (mean(lw_unlabel) * mean(rw_unlabel))
      
      return(c(DVar_k, LVar_k, RVar_k, DLCov_k, DRCov_k, LRCov_k))
      
    }, error = function(e) {
      warning(paste("CV fold", k, "failed:", e$message))
      return(rep(NA, 6))
    })
  })
  
  # Average across folds
  Var_mean <- apply(Var_cv, 1, mean, na.rm = TRUE)
  
  # Extract variance components
  DVar <- Var_mean[1] / (dh * n_label)
  LVar <- Var_mean[2] / (lh * n_label)
  RVar <- Var_mean[3] / (rh * n_label)
  DLCov <- Var_mean[4] / n_label
  DRCov <- Var_mean[5] / n_label
  LRCov <- Var_mean[6] / n_label
  
  # ==============================================================================
  # STEP 3: OPTIMAL COMBINATION
  # ==============================================================================
  
  cat("  Computing optimal combination...\n")
  
  # Compile estimators and covariance matrix
  VS <- c(DSt, LSt, RSt)
  Mcov <- matrix(c(DVar, DLCov, DRCov, 
                   DLCov, LVar, LRCov, 
                   DRCov, LRCov, RVar), ncol = 3, byrow = TRUE)
  
  # Check for missing estimators (zero response sums)
  response_sums <- c(sum(dy), sum(ly), sum(ry))
  missing_estimators <- which(response_sums == 0)
  
  if (length(missing_estimators) > 0) {
    cat(sprintf("  Excluding %d estimators with zero responses\n", length(missing_estimators)))
    VS <- VS[-missing_estimators]
    Mcov <- Mcov[-missing_estimators, -missing_estimators, drop = FALSE]
  }
  
  # Regularized optimization for combination weights
  if (length(VS) > 1) {
    
    # Try multiple regularization parameters
    Comb_wall <- list()
    Diff <- numeric()
    lambda_vals <- round(seq(0.2, 1, by = 0.15), 3)
    
    for (i in seq_along(lambda_vals)) {
      lambda_reg <- lambda_vals[i]
      epsilon <- n_label^(-lambda_reg)
      
      tryCatch({
        # Regularized covariance matrix
        Mcov_reg <- Mcov + epsilon * diag(ncol(Mcov))
        Mcov_inv <- solve(Mcov_reg)
        
        # Optimal weights under sum-to-one constraint
        ones <- rep(1, ncol(Mcov))
        denominator <- as.vector(t(ones) %*% Mcov_inv %*% ones)
        
        if (abs(denominator) > 1e-10) {
          weights <- as.vector(t(ones) %*% Mcov_inv) / denominator
          
          # Objective function for weight selection
          # Prefer weights that are more different from uniform
          if (length(weights) >= 2) {
            diff_score <- (weights[1] - weights[2])^2
            if (length(weights) >= 3) {
              diff_score <- diff_score + (weights[3] - weights[2])^2
            }
          } else {
            diff_score <- 0
          }
          
          Diff <- c(Diff, diff_score)
          Comb_wall <- c(Comb_wall, list(weights))
        } else {
          # Fallback to uniform weights
          Diff <- c(Diff, 0)
          Comb_wall <- c(Comb_wall, list(rep(1/length(VS), length(VS))))
        }
        
      }, error = function(e) {
        # Fallback to uniform weights
        Diff <- c(Diff, 0)
        Comb_wall <- c(Comb_wall, list(rep(1/length(VS), length(VS))))
      })
    }
    
    # Select best regularization
    if (length(Diff) > 0 && max(Diff, na.rm = TRUE) > 0) {
      best_idx <- which.max(Diff)
      Comb_w <- Comb_wall[[best_idx]]
    } else {
      # Fallback to uniform weights
      Comb_w <- rep(1/length(VS), length(VS))
    }
    
  } else {
    # Only one valid estimator
    Comb_w <- 1
  }
  
  # Compute combined estimator
  CSt <- as.vector(Comb_w %*% VS)
  
  # Map weights back to full 3-element vector
  New_Comb_w <- rep(0, 3)
  if (length(missing_estimators) == 0) {
    New_Comb_w <- Comb_w
  } else {
    New_Comb_w[-missing_estimators] <- Comb_w
  }
  
  cat(sprintf("  Combination weights: D=%.3f, L=%.3f, R=%.3f\n", 
              New_Comb_w[1], New_Comb_w[2], New_Comb_w[3]))
  
  # ==============================================================================
  # STEP 4: COMBINED ESTIMATOR VARIANCE
  # ==============================================================================
  
  cat("  Estimating combined variance...\n")
  
  # Cross-validation for combined estimator variance
  CVar_cv <- sapply(1:num_folds, function(k) {
    
    tryCatch({
      inds_fold <- as.vector(unlist(ind_cv[-k]))
      inds_foldk <- as.vector(unlist(ind_cv[k]))
      
      # Re-estimate parameters on training fold
      if (sum(dy) > 0) {
        Dbeta_K <- Beta_estimate(ddat_label[inds_fold, , drop = FALSE], 
                                 dy[inds_fold], dh, dw_label[inds_fold], dw_unlabel)
        coef_start <- 1 + a * (p_basis + 1)
        coef_end <- (p_basis + 1) + a * (p_basis + 1)
        dbeta_K <- Dbeta_K[coef_start:coef_end]
        X_val_d <- cbind(1, ddat_label[inds_foldk, , drop = FALSE])
        DImp_k <- as.vector(Expit(X_val_d %*% dbeta_K))
      } else {
        DImp_k <- rep(0, length(inds_foldk))
      }
      
      if (sum(ly) > 0) {
        Lbeta_K <- Beta_estimate(ldat_label[inds_fold, , drop = FALSE], 
                                 ly[inds_fold], lh, lw_label[inds_fold], lw_unlabel)
        lbeta_K <- Lbeta_K[coef_start:coef_end]
        X_val_l <- cbind(1, ldat_label[inds_foldk, , drop = FALSE])
        LImp_k <- as.vector(Expit(X_val_l %*% lbeta_K))
      } else {
        LImp_k <- rep(0, length(inds_foldk))
      }
      
      if (sum(ry) > 0) {
        Rbeta_K <- Beta_estimate(rdat_label[inds_fold, , drop = FALSE], 
                                 ry[inds_fold], rh, rw_label[inds_fold], rw_unlabel)
        rbeta_K <- Rbeta_K[coef_start:coef_end]
        X_val_r <- cbind(1, rdat_label[inds_foldk, , drop = FALSE])
        RImp_k <- as.vector(Expit(X_val_r %*% rbeta_K))
      } else {
        RImp_k <- rep(0, length(inds_foldk))
      }
      
      # Combined variance using optimal weights
      combined_residual <- (New_Comb_w[1] * mean(dw_unlabel)^(-1) * dw_label[inds_foldk] * (dy[inds_foldk] - DImp_k) +
                              New_Comb_w[2] * mean(lw_unlabel)^(-1) * lw_label[inds_foldk] * (ly[inds_foldk] - LImp_k) +
                              New_Comb_w[3] * mean(rw_unlabel)^(-1) * rw_label[inds_foldk] * (ry[inds_foldk] - RImp_k))
      
      CVar_k <- mean(combined_residual^2)
      return(CVar_k)
      
    }, error = function(e) {
      warning(paste("Combined variance estimation failed for fold", k))
      return(NA)
    })
  })
  
  CVarn <- mean(CVar_cv, na.rm = TRUE)
  CVar <- CVarn / n_label
  
  # ==============================================================================
  # STEP 5: COMPILE RESULTS
  # ==============================================================================
  
  # Compute standard deviations
  D_sd <- sqrt(max(0, DVar))
  L_sd <- sqrt(max(0, LVar))
  R_sd <- sqrt(max(0, RVar))
  C_sd <- sqrt(max(0, CVar))
  
  # Final validation
  if (any(is.na(c(DSt, LSt, RSt, CSt)))) {
    warning("Some survival estimates are NA")
  }
  
  if (any(c(D_sd, L_sd, R_sd, C_sd) < 0)) {
    warning("Negative variance detected")
  }
  
  # Ensure weights sum to 1 (approximately)
  weight_sum <- sum(New_Comb_w)
  if (abs(weight_sum - 1) > 0.01 && weight_sum > 0) {
    New_Comb_w <- New_Comb_w / weight_sum
  }
  
  cat("  Cross-validation completed successfully\n")
  
  return(list(
    DSt = DSt,           # Exact observed survival estimate
    LSt = LSt,           # Left survival estimate  
    RSt = RSt,           # Right survival estimate
    CSt = CSt,           # Combined survival estimate
    D_sd = D_sd,         # Exact observed standard deviation
    L_sd = L_sd,         # Left standard deviation
    R_sd = R_sd,         # Right standard deviation
    C_sd = C_sd,         # Combined standard deviation
    wd = New_Comb_w[1],  # Exact observed weight
    wl = New_Comb_w[2],  # Left weight  
    wr = New_Comb_w[3]   # Right weight
  ))
}





