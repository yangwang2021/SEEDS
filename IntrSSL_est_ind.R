# ==============================================================================
# SEEDS: Semi-supervised Estimation of Event rate with Doubly-censored Survival data
# ==============================================================================
# 
# This file implements the SEEDS (Semi-supervised Estimation of Event rate with
# Doubly-censored Survival data) method. SEEDS creates three specific types of
# estimators from doubly-censored data and combines them optimally using
# cross-validation.
#
# The method works by:
# 1. Creating three types of survival estimators (Direct, Left, Right)
# 2. Computing both semi-supervised and intrinsic versions
# 3. Combining estimators optimally using cross-validation
# 4. Providing variance estimation and optimal weights
#
# Functions included:
# - IntrSSL_est(): Main SEEDS estimation function
# ==============================================================================

#' SEEDS: Semi-supervised Estimation of Event rate with Doubly-censored Survival data
#' 
#' This function implements the SEEDS method for estimating survival rates using
#' both labeled and unlabeled data. The method creates three specific estimators
#' based on different types of censoring information and combines them optimally.
#' 
#' @param time Numeric vector, evaluation time points for survival estimation
#' @param base_cov Numeric matrix, baseline time-independent covariates (optional)
#' @param cova_tim List of vectors, observation times for time-dependent covariates
#' @param cova_ct List of vectors, time-dependent covariate values
#' @param lcen_ct Numeric vector, left censoring covariate values (optional)
#' @param rcen_ct Numeric vector, right censoring covariate values (optional)
#' @param xstar_all Numeric vector, surrogate variable observations for all subjects
#' @param deltastar_indv Numeric matrix, surrogate censoring indicators (2 columns)
#' @param label_id Numeric vector, indices of labeled observations
#' @param lcen_all Numeric vector, left censoring times for all subjects
#' @param rcen_all Numeric vector, right censoring times for all subjects
#' @param obse Numeric vector, observed survival times for labeled subjects
#' @param ldelt Numeric vector, left censoring indicators for labeled subjects
#' @param rdelt Numeric vector, right censoring indicators for labeled subjects
#' @param lh_labeled Numeric scalar, left bandwidth for labeled data
#' @param lh_unlabeled Numeric scalar, left bandwidth for unlabeled data  
#' @param rh_labeled Numeric scalar, right bandwidth for labeled data
#' @param rh_unlabeled Numeric scalar, right bandwidth for unlabeled data
#' @param num_folds Integer, number of cross-validation folds (default: 10)
#' 
#' @return Matrix with the following rows (19 total):
#'   \item{Rows 1-2}{Direct estimator: semi-supervised and intrinsic}
#'   \item{Rows 3-4}{Left estimator: semi-supervised and intrinsic}
#'   \item{Rows 5-6}{Right estimator: semi-supervised and intrinsic}
#'   \item{Rows 7-8}{Combined estimator: semi-supervised and intrinsic}
#'   \item{Rows 9-16}{Standard deviations for above estimators}
#'   \item{Rows 17-19}{Optimal combination weights for Direct, Left, Right}
#'   
#' @details
#' The SEEDS method works in two main steps:
#' 
#' **Step 1: Create Three Types of Estimators**
#' 
#' 1. **Direct Estimator**: Uses exact observations
#'    \deqn{S_D(t) = E[I(T \leq t, \text{exact observation}) | X]}
#'    
#' 2. **Left Estimator**: Uses left-censored information  
#'    \deqn{S_L(t) = E[I(\text{not left-censored at } t) | X]}
#'    
#' 3. **Right Estimator**: Uses right-censored information
#'    \deqn{S_R(t) = E[I(\text{not right-censored at } t) | X]}
#' 
#' **Step 2: Optimal Combination**
#' Uses cross-validation to find optimal weights w_D, w_L, w_R:
#' \deqn{S_{SEEDS}(t) = w_D S_D(t) + w_L S_L(t) + w_R S_R(t)}
#' subject to \eqn{w_D + w_L + w_R = 1}
#' 
#' Each estimator comes in two versions:
#' - **Semi-supervised**: Uses labeled data with imputation model
#' - **Intrinsic**: Uses both labeled and unlabeled data with constraints
#' 
#' @examples
#' # This function is typically called from the main simulation pipeline
#' # See basic_example.R for complete usage
#' 
#' # Example structure (not runnable without full data setup):
#' # results <- IntrSSL_est(
#' #   time = seq(1, 5, length=20),
#' #   base_cov = matrix(rnorm(1000), ncol=2),
#' #   cova_tim = list(...),  # Time-dependent covariate times
#' #   cova_ct = list(...),   # Time-dependent covariate values
#' #   # ... other parameters
#' # )
#' 
#' 
#' @seealso \code{\link{CV_function}} for cross-validation implementation

IntrSSL_est <- function(time, base_cov, cova_tim, cova_ct, lcen_ct, rcen_ct, 
                        xstar_all, deltastar_indv, label_id, lcen_all, rcen_all, 
                        obse, ldelt, rdelt, lh_labeled, lh_unlabeled, rh_labeled, 
                        rh_unlabeled, num_folds) {
  
  # ==============================================================================
  # INPUT VALIDATION
  # ==============================================================================
  
  if (length(time) == 0) stop("time vector cannot be empty")
  if (length(label_id) == 0) stop("No labeled observations provided")
  if (num_folds < 2) stop("num_folds must be at least 2")
  
  n_times <- length(time)
  n_total <- length(xstar_all)
  n_labeled <- length(label_id)
  
  cat(sprintf("SEEDS estimation for %d time points, %d total obs (%d labeled)\n",
              n_times, n_total, n_labeled))
  
  # ==============================================================================
  # MAIN ESTIMATION LOOP OVER TIME POINTS
  # ==============================================================================
  
  # Process each time point independently
  result_matrix <- sapply(time, function(t) {
    
    # ============================================================================
    # STEP 1: CONSTRUCT COVARIATE MATRICES FOR CURRENT TIME POINT
    # ============================================================================
    
    # Create surrogate indicator: I(X* <= t)
    Ind_star <- as.numeric(xstar_all <= t)
    
    # --- Direct Estimator Covariates ---
    # Combines surrogate censoring indicators, time indicator, and covariates
    if (is.null(cova_ct)) {
      dcova_all <- NULL
    } else {
      # Extract time-dependent covariate values at time t
      dcova_all <- TimeCovar(
        time_point = t, 
        List_time = cova_tim, 
        List_count = cova_ct, 
        leng = n_total
      )
    }
    
    # Build design matrix for direct estimator
    Dbasis_all <- cbind(deltastar_indv, Ind_star, base_cov, dcova_all)
    Dbasis_labeled <- Dbasis_all[label_id, , drop = FALSE]
    Dbasis_unlabeled <- Dbasis_all[-label_id, , drop = FALSE]
    
    # --- Left Estimator Covariates ---
    # Uses left censoring covariate information
    if (is.null(lcen_ct)) {
      lcova_all <- NULL
    } else {
      lcova_all <- lcen_ct
    }
    
    Lbasis_all <- cbind(deltastar_indv, Ind_star, base_cov, lcova_all)
    Lbasis_labeled <- Lbasis_all[label_id, , drop = FALSE]
    Lbasis_unlabeled <- Lbasis_all[-label_id, , drop = FALSE]
    
    # --- Right Estimator Covariates ---
    # Uses right censoring covariate information
    if (is.null(rcen_ct)) {
      rcova_all <- NULL
    } else {
      rcova_all <- rcen_ct
    }
    
    Rbasis_all <- cbind(deltastar_indv, Ind_star, base_cov, rcova_all)
    Rbasis_labeled <- Rbasis_all[label_id, , drop = FALSE]
    Rbasis_unlabeled <- Rbasis_all[-label_id, , drop = FALSE]
    
    # ============================================================================
    # STEP 2: PREPARE RESPONSE VARIABLES AND WEIGHTS
    # ============================================================================
    
    # Extract censoring times for labeled observations
    lcen_labeled <- lcen_all[label_id]
    rcen_labeled <- rcen_all[label_id]
    
    # Extract censoring times for unlabeled observations
    lcen_unlabeled <- lcen_all[-label_id]
    rcen_unlabeled <- rcen_all[-label_id]
    
    # --- Direct Estimator Response ---
    # Indicator for being in observation window at time t
    Inter_labeled <- as.numeric(t <= rcen_labeled & t > lcen_labeled)
    Inter_unlabeled <- as.numeric(t <= rcen_unlabeled & t > lcen_unlabeled)
    
    # Response: indicator for exact observation at time t
    dy_obse <- as.numeric(t <= obse & t > lcen_labeled)
    
    # --- Left Estimator Response and Weights ---
    # Kernel weights centered at time t for left censoring
    lw_labeled <- Normalize_fun(
      data = lcen_labeled, 
      pt = t, 
      bd = lh_labeled
    )
    lw_unlabeled <- Normalize_fun(
      data = lcen_unlabeled, 
      pt = t, 
      bd = lh_unlabeled
    )
    
    # --- Right Estimator Response and Weights ---
    # Kernel weights centered at time t for right censoring
    rw_labeled <- Normalize_fun(
      data = rcen_labeled, 
      pt = t, 
      bd = rh_labeled
    )
    rw_unlabeled <- Normalize_fun(
      data = rcen_unlabeled, 
      pt = t, 
      bd = rh_unlabeled
    )
    
    # ============================================================================
    # STEP 3: SEMI-SUPERVISED ESTIMATION (a = 0)
    # ============================================================================
    
    cat(sprintf("  Time %.3f: Semi-supervised estimation\n", t))
    
    SSL_est <- tryCatch({
      CV_function(
        ddat_unlabel = Dbasis_unlabeled,     # Direct: unlabeled covariates
        ddat_label = Dbasis_labeled,         # Direct: labeled covariates
        ldat_unlabel = Lbasis_unlabeled,     # Left: unlabeled covariates
        ldat_label = Lbasis_labeled,         # Left: labeled covariates
        rdat_unlabel = Rbasis_unlabeled,     # Right: unlabeled covariates
        rdat_label = Rbasis_labeled,         # Right: labeled covariates
        dw_label = Inter_labeled,            # Direct: labeled weights
        dw_unlabel = Inter_unlabeled,        # Direct: unlabeled weights
        dy = dy_obse,                        # Direct: response
        dh = 1,                              # Direct: bandwidth (fixed)
        lw_label = lw_labeled,               # Left: labeled weights
        lw_unlabel = lw_unlabeled,          # Left: unlabeled weights
        ly = 1 - ldelt,                      # Left: response (survival indicator)
        lh = lh_labeled,                     # Left: bandwidth
        rw_label = rw_labeled,               # Right: labeled weights
        rw_unlabel = rw_unlabeled,          # Right: unlabeled weights
        ry = 1 - rdelt,                      # Right: response (survival indicator)
        rh = rh_labeled,                     # Right: bandwidth
        num_folds = num_folds,               # Cross-validation folds
        a = 0                                # Semi-supervised flag
      )
    }, error = function(e) {
      warning(paste("Semi-supervised estimation failed at time", t, ":", e$message))
      return(list(
        DSt = NA, LSt = NA, RSt = NA, CSt = NA,
        D_sd = NA, L_sd = NA, R_sd = NA, C_sd = NA,
        wd = NA, wl = NA, wr = NA
      ))
    })
    
    # Extract semi-supervised results
    DSt_ssl <- SSL_est$DSt      # Direct survival estimate
    LSt_ssl <- SSL_est$LSt      # Left survival estimate  
    RSt_ssl <- SSL_est$RSt      # Right survival estimate
    CSt_ssl <- SSL_est$CSt      # Combined survival estimate
    Dsd_ssl <- SSL_est$D_sd     # Direct standard deviation
    Lsd_ssl <- SSL_est$L_sd     # Left standard deviation
    Rsd_ssl <- SSL_est$R_sd     # Right standard deviation
    Csd_ssl <- SSL_est$C_sd     # Combined standard deviation
    
    # ============================================================================
    # STEP 4: INTRINSIC ESTIMATION (a = 1)
    # ============================================================================
    
    cat(sprintf("  Time %.3f: Intrinsic estimation\n", t))
    
    Intr_est <- tryCatch({
      CV_function(
        ddat_unlabel = Dbasis_unlabeled,     # Same covariate setup
        ddat_label = Dbasis_labeled,
        ldat_unlabel = Lbasis_unlabeled,
        ldat_label = Lbasis_labeled,
        rdat_unlabel = Rbasis_unlabeled,
        rdat_label = Rbasis_labeled,
        dw_label = Inter_labeled,            # Same weight setup
        dw_unlabel = Inter_unlabeled,
        dy = dy_obse,                        # Same response setup
        dh = 1,
        lw_label = lw_labeled,
        lw_unlabel = lw_unlabeled,
        ly = 1 - ldelt,
        lh = lh_labeled,
        rw_label = rw_labeled,
        rw_unlabel = rw_unlabeled,
        ry = 1 - rdelt,
        rh = rh_labeled,
        num_folds = num_folds,
        a = 1                                # Intrinsic flag
      )
    }, error = function(e) {
      warning(paste("Intrinsic estimation failed at time", t, ":", e$message))
      return(list(
        DSt = NA, LSt = NA, RSt = NA, CSt = NA,
        D_sd = NA, L_sd = NA, R_sd = NA, C_sd = NA,
        wd = NA, wl = NA, wr = NA
      ))
    })
    
    # Extract intrinsic results
    DSt_intr <- Intr_est$DSt    # Direct survival estimate
    LSt_intr <- Intr_est$LSt    # Left survival estimate
    RSt_intr <- Intr_est$RSt    # Right survival estimate
    CSt_intr <- Intr_est$CSt    # Combined survival estimate
    Dsd_intr <- Intr_est$D_sd   # Direct standard deviation
    Lsd_intr <- Intr_est$L_sd   # Left standard deviation
    Rsd_intr <- Intr_est$R_sd   # Right standard deviation
    Csd_intr <- Intr_est$C_sd   # Combined standard deviation
    
    # Optimal combination weights
    Weig_d <- Intr_est$wd       # Weight for direct estimator
    Weig_l <- Intr_est$wl       # Weight for left estimator
    Weig_r <- Intr_est$wr       # Weight for right estimator
    
    # ============================================================================
    # STEP 5: RETURN RESULTS FOR THIS TIME POINT
    # ============================================================================
    
    return(c(
      # Row 1-2: Direct estimator (semi-supervised, intrinsic)
      'DSt_ssl' = DSt_ssl, 'DSt_intr' = DSt_intr,
      
      # Row 3-4: Left estimator (semi-supervised, intrinsic)
      'LSt_ssl' = LSt_ssl, 'LSt_intr' = LSt_intr,
      
      # Row 5-6: Right estimator (semi-supervised, intrinsic)
      'RSt_ssl' = RSt_ssl, 'RSt_intr' = RSt_intr,
      
      # Row 7-8: Combined estimator (semi-supervised, intrinsic)
      'CSt_ssl' = CSt_ssl, 'CSt_intr' = CSt_intr,
      
      # Row 9-12: Standard deviations for semi-supervised estimators
      'Dsd_ssl' = Dsd_ssl, 'Dsd_intr' = Dsd_intr,
      
      # Row 13-16: Standard deviations for intrinsic estimators
      'Lsd_ssl' = Lsd_ssl, 'Lsd_intr' = Lsd_intr,
      
      # Row 17-19: Standard deviations continued
      'Rsd_ssl' = Rsd_ssl, 'Rsd_intr' = Rsd_intr,
      'Csd_ssl' = Csd_ssl, 'Csd_intr' = Csd_intr,
      
      # Row 17-19: Optimal combination weights
      'Weig_d' = Weig_d, 'Weig_l' = Weig_l, 'Weig_r' = Weig_r
    ))
    
  })  # End of sapply over time points
  
  # ==============================================================================
  # FINAL PROCESSING AND VALIDATION
  # ==============================================================================
  
  # Ensure result is a matrix
  if (!is.matrix(result_matrix)) {
    result_matrix <- matrix(result_matrix, nrow = 19)
  }
  
  # Add descriptive row names
  rownames(result_matrix) <- c(
    "DSt_ssl", "DSt_intr",           # Direct estimators
    "LSt_ssl", "LSt_intr",           # Left estimators  
    "RSt_ssl", "RSt_intr",           # Right estimators
    "CSt_ssl", "CSt_intr",           # Combined estimators
    "Dsd_ssl", "Dsd_intr",           # Direct standard deviations
    "Lsd_ssl", "Lsd_intr",           # Left standard deviations
    "Rsd_ssl", "Rsd_intr",           # Right standard deviations
    "Csd_ssl", "Csd_intr",           # Combined standard deviations
    "Weig_d", "Weig_l", "Weig_r"     # Combination weights
  )
  
  # Add column names with time points
  colnames(result_matrix) <- paste0("t_", round(time, 3))
  
  # Validation checks
  n_missing <- sum(is.na(result_matrix))
  if (n_missing > 0) {
    warning(sprintf("SEEDS estimation completed with %d missing values", n_missing))
  }
  
  # Check weight constraints (should sum to 1)
  weight_sums <- result_matrix["Weig_d", ] + result_matrix["Weig_l", ] + result_matrix["Weig_r", ]
  weight_violations <- sum(abs(weight_sums - 1) > 0.01, na.rm = TRUE)
  if (weight_violations > 0) {
    warning(sprintf("Weight constraint violated at %d time points", weight_violations))
  }
  
  # Check survival probability bounds (should be in [0,1])
  surv_rows <- c("DSt_ssl", "DSt_intr", "LSt_ssl", "LSt_intr", 
                 "RSt_ssl", "RSt_intr", "CSt_ssl", "CSt_intr")
  surv_values <- result_matrix[surv_rows, ]
  bound_violations <- sum(surv_values < 0 | surv_values > 1, na.rm = TRUE)
  if (bound_violations > 0) {
    warning(sprintf("Survival probability bounds violated %d times", bound_violations))
  }
  
  cat(sprintf("SEEDS estimation completed for %d time points\n", n_times))
  cat(sprintf("Average combination weights: D=%.3f, L=%.3f, R=%.3f\n",
              mean(result_matrix["Weig_d", ], na.rm = TRUE),
              mean(result_matrix["Weig_l", ], na.rm = TRUE),
              mean(result_matrix["Weig_r", ], na.rm = TRUE)))
  
  return(result_matrix)
}

