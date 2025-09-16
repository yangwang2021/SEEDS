# ==============================================================================
# Helper Functions for Semi-Supervised Survival Analysis
# ==============================================================================
# 
# This file contains utility functions used throughout the survival analysis
# methods. These include mathematical functions, kernel functions, data
# processing utilities, and cross-validation helpers.
#
# Functions included:
# - Expit and derivative functions
# - Kernel functions for smoothing
# - Time-dependent covariate processing
# - Cross-validation utilities  
# - Spline basis functions
# - Data normalization functions
# ==============================================================================

# ==============================================================================
# MATHEMATICAL FUNCTIONS
# ==============================================================================

#' Expit (logistic) function
#' 
#' Computes the expit function: expit(x) = 1 / (1 + exp(-x))
#' This is the inverse of the logit function and gives probabilities.
#' 
#' @param x Numeric vector or scalar
#' @return Numeric vector with values between 0 and 1
#' @examples
#' Expit(0)     # Returns 0.5
#' Expit(c(-2, 0, 2))  # Returns vector of probabilities
Expit <- function(x) {
  return(1 / (1 + exp(-x)))
}

#' Derivative of the expit function
#' 
#' Computes the derivative of the expit function: d/dx expit(x) = exp(x)/(1+exp(x))^2
#' This is useful for gradient computations in optimization.
#' 
#' @param x Numeric vector or scalar
#' @param na_correction Logical, whether to replace NA/NaN values with 0
#' @return Numeric vector of derivatives
#' @examples
#' ExpitDerivative(0)   # Returns 0.25
#' ExpitDerivative(c(-10, 0, 10))  # Returns derivatives
ExpitDerivative <- function(x, na_correction = TRUE) {
  expit_deriv <- exp(x) / (1 + exp(x))^2
  
  if (na_correction) {
    expit_deriv[which(is.na(expit_deriv))] <- 0
  }
  
  return(expit_deriv)
}

# ==============================================================================
# KERNEL FUNCTIONS
# ==============================================================================

#' Gaussian (Normal) kernel function
#' 
#' Computes the Gaussian kernel K(x/h) = exp(-(x/h)^2/2) / (h * sqrt(2*pi))
#' Used for kernel smoothing in non-parametric estimation.
#' 
#' @param x Numeric vector, distances from the center
#' @param h Numeric scalar, bandwidth parameter
#' @return Numeric vector of kernel weights
#' @examples
#' ker(c(-1, 0, 1), h = 1)  # Gaussian weights at -1, 0, 1
ker <- function(x, h) {
  return(exp(-(x/h)^2/2) / (h * sqrt(2*pi)))
}

# Alternative: Epanechnikov kernel (commented out)
# #' Epanechnikov kernel function
# #' 
# #' Computes the Epanechnikov kernel: K(u) = 0.75*(1-u^2)*I(|u|<1)/h
# #' This kernel has compact support and is sometimes preferred.
# #' 
# #' @param x Numeric vector, distances from the center
# #' @param h Numeric scalar, bandwidth parameter  
# #' @return Numeric vector of kernel weights
# ker_epanechnikov <- function(x, h){
#   u <- x / h
#   return(0.75 * (1 - u^2) * (abs(u) < 1) / h)
# }

# ==============================================================================
# TIME-DEPENDENT COVARIATE FUNCTIONS
# ==============================================================================

#' Extract time-dependent covariate values at specific time point
#' 
#' For each subject, finds the covariate value that corresponds to a specific
#' time point by looking up the most recent observation before that time.
#' 
#' @param time_point Numeric scalar, the time point of interest
#' @param List_time List of vectors, each containing observation times for one subject
#' @param List_count List of vectors, each containing covariate values for one subject
#' @param leng Integer, number of subjects
#' @return Numeric vector of covariate values at time_point for each subject
#' @examples
#' # Subject 1 has observations at times c(0, 1, 2) with values c(0, 3, 5)
#' # TimeCovar(1.5, list(c(0,1,2)), list(c(0,3,5)), 1) returns 3
TimeCovar <- function(time_point, List_time, List_count, leng) {
  
  result <- sapply(1:leng, function(i) {
    # Find the last observation time <= time_point
    valid_times <- which(List_time[[i]] <= time_point)
    
    if (length(valid_times) == 0) {
      # No observations before time_point, return 0
      return(0)
    } else {
      # Return the count at the last valid time
      last_time_index <- max(valid_times)
      return(List_count[[i]][last_time_index])
    }
  })
  
  return(result)
}

# ==============================================================================
# CROSS-VALIDATION UTILITIES
# ==============================================================================

#' Create cross-validation folds
#' 
#' Randomly splits n observations into num_folds approximately equal groups
#' for cross-validation.
#' 
#' @param n Integer, total number of observations
#' @param num_folds Integer, number of cross-validation folds
#' @return List of vectors, each containing indices for one fold
#' @examples
#' cv_split(100, 5)  # Creates 5 folds of ~20 observations each
cv_split <- function(n, num_folds) {
  
  # Create fold assignments
  fold_assignments <- sample(rep(1:num_folds, length.out = n))
  
  # Split indices into folds
  folds <- split(1:n, fold_assignments)
  
  return(folds)
}

# ==============================================================================
# DATA NORMALIZATION FUNCTIONS
# ==============================================================================

#' Normalize kernel weights
#' 
#' Computes kernel weights at a point and normalizes them to have mean 1.
#' This ensures proper weighting in kernel regression.
#' 
#' @param data Numeric vector, data points
#' @param pt Numeric scalar, evaluation point
#' @param bd Numeric scalar, bandwidth
#' @return Numeric vector of normalized weights
#' @examples
#' Normalize_fun(c(1,2,3,4,5), pt=3, bd=1)  # Weights centered at 3
Normalize_fun <- function(data, pt, bd) {
  
  # Compute raw kernel weights
  weights <- ker(data - pt, bd)
  
  # Normalize weights
  if (mean(weights) == 0) {
    # If all weights are zero, return zero weights
    weights <- rep(0, length(weights))
  } else {
    # Normalize to have mean 1
    weights <- weights / mean(weights)
  }
  
  return(weights)
}

# ==============================================================================
# SPLINE BASIS FUNCTIONS
# ==============================================================================

#' Truncated cubic function
#' 
#' Computes truncated cubic basis function: max(0, x - knot)^3
#' Used in constructing natural spline bases.
#' 
#' @param x Numeric vector, evaluation points
#' @param knot_location Numeric scalar, knot location
#' @return Numeric vector of truncated cubic values
TruncatedCubic <- function(x, knot_location) {
  return(pmax(0, x - knot_location)^3)
}

#' Natural spline basis matrix
#' 
#' Constructs a natural spline basis matrix with specified number of knots.
#' Natural splines are linear beyond the boundary knots.
#' 
#' @param X Numeric matrix, predictor variables
#' @param num_knots Integer, number of knots to use
#' @return Numeric matrix of natural spline basis functions
#' @details 
#' For each column of X, creates a natural spline basis with num_knots knots
#' placed at equally spaced quantiles. The basis includes the original variable
#' plus truncated cubic terms.
NaturalSplineBasis <- function(X, num_knots) {
  
  X <- as.matrix(X)
  basis.X <- NULL
  
  for (i in 1:ncol(X)) {
    X_i <- X[, i]
    
    # Determine knot locations using quantiles
    knots <- quantile(X_i, seq(0, 1, length = num_knots))
    
    # Handle case where there aren't enough unique values
    j <- 0
    while (length(unique(knots)) != num_knots) {
      j <- j + 1
      knots <- unique(quantile(X_i, seq(0, 1, length = num_knots + j)))
      
      # Prevent infinite loop
      if (j > 10) {
        warning("Could not find enough unique knots for natural spline")
        break
      }
    }
    
    # Compute natural spline basis
    if (num_knots >= 2) {
      # Reference knot (last knot)
      d_k <- TruncatedCubic(X_i, knots[num_knots - 1]) / 
        (knots[num_knots] - knots[num_knots - 1])
      
      # Compute basis functions for interior knots
      if (num_knots > 2) {
        evals <- sapply(1:(num_knots - 2), function(ii) {
          d_i <- TruncatedCubic(X_i, knots[ii]) / (knots[num_knots] - knots[ii])
          return(d_i - d_k)
        })
        
        # Combine original variable with spline terms
        basis.X <- cbind(basis.X, X_i, evals)
      } else {
        # Only two knots - just use original variable
        basis.X <- cbind(basis.X, X_i)
      }
    } else {
      # Less than 2 knots - just use original variable
      basis.X <- cbind(basis.X, X_i)
    }
  }
  
  return(basis.X)
}

# ==============================================================================
# ADDITIONAL UTILITY FUNCTIONS
# ==============================================================================

#' Check if a value is effectively zero
#' 
#' @param x Numeric value to check
#' @param tol Tolerance for zero comparison
#' @return Logical indicating if |x| < tol
is_effectively_zero <- function(x, tol = 1e-10) {
  return(abs(x) < tol)
}

#' Safely compute log probabilities
#' 
#' @param p Numeric vector of probabilities
#' @param eps Small value to add to avoid log(0)
#' @return Numeric vector of log probabilities
safe_log <- function(p, eps = 1e-10) {
  return(log(pmax(p, eps)))
}

#' Print summary statistics for a vector
#' 
#' @param x Numeric vector
#' @param name Character string, name for the variable
print_summary <- function(x, name = "Variable") {
  cat(paste(name, "summary:\n"))
  cat(paste("  Mean:", round(mean(x, na.rm = TRUE), 4), "\n"))
  cat(paste("  SD:", round(sd(x, na.rm = TRUE), 4), "\n"))
  cat(paste("  Min:", round(min(x, na.rm = TRUE), 4), "\n"))
  cat(paste("  Max:", round(max(x, na.rm = TRUE), 4), "\n"))
  cat(paste("  NA count:", sum(is.na(x)), "\n"))
}

# ==============================================================================
# END OF HELPER FUNCTIONS
# ==============================================================================