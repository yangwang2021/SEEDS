# ==============================================================================
# Parameter Estimation for Semi-Supervised Survival Analysis
# ==============================================================================
# 
# This file implements parameter estimation methods for semi-supervised learning
# in survival analysis. It provides both conventional maximum likelihood estimation
# and intrinsic semi-supervised estimation approaches.
#
# The estimation process includes:
# 1. Semi-supervised estimator using labeled data with logistic regression
# 2. Intrinsic estimator using Newton-Raphson optimization with constraints
# 3. Proper handling of both labeled and unlabeled data
#
# Functions included:
# - Beta_estimate(): Main parameter estimation function
# ==============================================================================

#' Parameter Estimation for Semi-Supervised Survival Models
#' 
#' This function estimates regression parameters using both semi-supervised
#' and intrinsic estimation approaches. The semi-supervised estimator uses
#' standard logistic regression on labeled data, while the intrinsic estimator
#' incorporates information from unlabeled data through constrained optimization.
#' 
#' @param dat_label Numeric matrix, covariate matrix for labeled observations
#' @param y Numeric vector, binary response variable for labeled observations (0/1)
#' @param h Numeric scalar, bandwidth parameter for kernel smoothing
#' @param w_label Numeric vector, weights for labeled observations
#' @param w_unlabel Numeric vector, weights for unlabeled observations
#' 
#' @return Numeric vector containing:
#'   - First (p+1) elements: Semi-supervised estimator coefficients
#'   - Next (p+1) elements: Intrinsic estimator coefficients
#'   where p is the number of covariates in dat_label
#'   
#' @details
#' The function implements two estimation approaches:
#' 
#' **Step 1: Semi-supervised Estimator**
#' Uses standard logistic regression with weights:
#' \deqn{\beta_{ssl} = \arg\min_\beta \sum_{i \in L} w_i \ell(y_i, X_i^T \beta)}
#' where L is the set of labeled observations and \ell is the logistic loss.
#' 
#' **Step 2: Intrinsic Estimator**  
#' Uses constrained optimization that incorporates unlabeled data:
#' \deqn{\beta_{intr} = \arg\min_\beta h \sum_{i \in L} w_i^2 (y_i - \sigma(X_i^T \beta))^2}
#' subject to normalization constraints involving unlabeled data weights.
#' 
#' The intrinsic estimator is computed using Newton-Raphson optimization
#' implemented in NewtonCLseh() function.
#' 
#' @seealso \code{\link{NewtonCLseh}} for the optimization algorithm
Beta_estimate <- function(dat_label, y, h, w_label, w_unlabel) {
  
  # ==============================================================================
  # INPUT VALIDATION
  # ==============================================================================
  
  # Check dimensions
  if (!is.matrix(dat_label) && !is.data.frame(dat_label)) {
    stop("dat_label must be a matrix or data frame")
  }
  
  if (nrow(dat_label) != length(y)) {
    stop("Number of rows in dat_label must equal length of y")
  }
  
  if (nrow(dat_label) != length(w_label)) {
    stop("Number of rows in dat_label must equal length of w_label")
  }
  
  if (!all(y %in% c(0, 1))) {
    stop("y must be binary (0/1)")
  }
  
  if (h <= 0) {
    stop("Bandwidth h must be positive")
  }
  
  if (any(w_label < 0) || any(w_unlabel < 0)) {
    stop("Weights must be non-negative")
  }
  
  # Convert to matrix if needed
  dat_label <- as.matrix(dat_label)
  
  # Get dimensions
  n_labeled <- nrow(dat_label)
  p_basis <- ncol(dat_label)  # Number of covariates
  
  cat(sprintf("Estimating parameters for %d labeled observations with %d covariates\n", 
              n_labeled, p_basis))
  
  # ==============================================================================
  # STEP 1: SEMI-SUPERVISED ESTIMATOR
  # ==============================================================================
  
  cat("Computing semi-supervised estimator...\n")
  
  # Fit weighted logistic regression using labeled data only
  # This serves as the baseline estimator and initialization for intrinsic method
  beta_ssl <- tryCatch({
    
    # Create data frame for glm
    glm_data <- data.frame(
      y = y,
      dat_label
    )
    
    # Fit weighted logistic regression
    glm_fit <- glm(y ~ ., 
                   data = glm_data, 
                   family = binomial(link = "logit"), 
                   weights = w_label)
    
    # Extract coefficients
    coef_ssl <- glm_fit$coefficients
    
    # Handle potential convergence issues
    if (any(is.na(coef_ssl))) {
      warning("Some coefficients are NA in semi-supervised estimation")
      coef_ssl[is.na(coef_ssl)] <- 0
    }
    
    # Remove names for consistency
    unname(coef_ssl)
    
  }, error = function(e) {
    warning(paste("Semi-supervised estimation failed:", e$message))
    warning("Using zero initialization")
    return(rep(0, p_basis + 1))  # Intercept + covariates
  })
  
  # Ensure we have the right number of coefficients
  if (length(beta_ssl) != (p_basis + 1)) {
    warning("Coefficient dimension mismatch in semi-supervised estimation")
    beta_ssl <- rep(0, p_basis + 1)
  }
  
  cat(sprintf("Semi-supervised estimation completed. Coefficients: %s\n", 
              paste(round(beta_ssl, 4), collapse = ", ")))
  
  # ==============================================================================
  # STEP 2: INTRINSIC ESTIMATOR
  # ==============================================================================
  
  cat("Computing intrinsic estimator...\n")
  
  # The intrinsic estimator uses Newton-Raphson optimization with constraints
  # that incorporate information from unlabeled data weights
  beta_intr <- tryCatch({
    
    # Call Newton-Raphson optimization
    # This implements the constrained optimization problem that leverages
    # both labeled and unlabeled data
    NewtonCLseh(
      X = dat_label,          # Covariate matrix (without intercept)
      y = y,                  # Binary response 
      h = h,                  # Bandwidth parameter
      w_label = w_label,      # Labeled data weights
      w_unlabel = w_unlabel,  # Unlabeled data weights
      max_iter = 100,         # Maximum iterations
      tol = 1e-4,            # Convergence tolerance
      uper = 1e+15,          # Upper bound for condition number
      initial = beta_ssl      # Use semi-supervised estimate as initialization
    )
    
  }, error = function(e) {
    warning(paste("Intrinsic estimation failed:", e$message))
    warning("Falling back to semi-supervised estimate")
    return(beta_ssl)
  })
  
  # Validate intrinsic estimator results
  if (length(beta_intr) != (p_basis + 1)) {
    warning("Coefficient dimension mismatch in intrinsic estimation")
    beta_intr <- beta_ssl  # Fallback to semi-supervised
  }
  
  if (any(is.na(beta_intr)) || any(is.infinite(beta_intr))) {
    warning("Invalid coefficients in intrinsic estimation")
    beta_intr <- beta_ssl  # Fallback to semi-supervised
  }
  
  cat(sprintf("Intrinsic estimation completed. Coefficients: %s\n", 
              paste(round(beta_intr, 4), collapse = ", ")))
  
  # ==============================================================================
  # STEP 3: COMBINE AND RETURN RESULTS
  # ==============================================================================
  
  # Combine both estimators into a single vector
  # First (p+1) elements: semi-supervised estimator
  # Next (p+1) elements: intrinsic estimator
  result <- c(
    beta_ssl = beta_ssl,    # Semi-supervised coefficients
    beta_intr = beta_intr   # Intrinsic coefficients
  )
  
  # Add names for clarity (optional)
  covariate_names <- if (is.null(colnames(dat_label))) {
    paste0("X", 1:p_basis)
  } else {
    colnames(dat_label)
  }
  
  ssl_names <- c("(Intercept)_ssl", paste0(covariate_names, "_ssl"))
  intr_names <- c("(Intercept)_intr", paste0(covariate_names, "_intr"))
  
  names(result) <- c(ssl_names, intr_names)
  
  
  return(result)
}

