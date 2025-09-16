# ==============================================================================
# Newton-Raphson Optimization with Constraints for Semi-Supervised Learning
# ==============================================================================
# 
# This file implements Newton-Raphson optimization algorithm with linear constraints
# for semi-supervised parameter estimation. The method solves a constrained
# optimization problem that incorporates both labeled and unlabeled data.
#
# The optimization problem is:
# minimize: h * E[w_label^2 * (y - expit(X*beta))^2] / E[w_unlabel]^2
# subject to: E[w_label * expit(X*beta)] = E[w_label * y]
#
# Functions included:
# - NewtonCLseh(): Main Newton-Raphson optimization with constraints
#
# ==============================================================================

#' Newton-Raphson Optimization with Linear Constraints
#' 
#' This function implements Newton-Raphson algorithm for solving constrained
#' optimization problems in semi-supervised learning. The algorithm minimizes
#' a weighted squared loss while satisfying linear constraints that ensure
#' consistency between labeled and unlabeled data.
#' 
#' @param X Numeric matrix, design matrix of covariates (n x p)
#' @param y Numeric vector, binary response variable (0/1)
#' @param h Numeric scalar, bandwidth parameter for kernel smoothing
#' @param w_label Numeric vector, weights for labeled observations
#' @param w_unlabel Numeric vector, weights for unlabeled observations  
#' @param max_iter Integer, maximum number of iterations (default: 100)
#' @param tol Numeric, convergence tolerance (default: 1e-4)
#' @param uper Numeric, upper bound for condition number (default: 1e+15)
#' @param initial Numeric vector, initial parameter values (default: zeros)
#' 
#' @return Numeric vector of optimized parameters (length p+1 including intercept)
#' 
#' @details
#' The algorithm solves the constrained optimization problem:
#' 
#' **Objective Function:**
#' \deqn{L(\beta) = h \cdot \frac{1}{n} \sum_{i=1}^n w_i^2 (y_i - \sigma(x_i^T \beta))^2 / \bar{w}_{unlabel}^2}
#' 
#' **Constraint:**
#' \deqn{\frac{1}{n} \sum_{i=1}^n w_i \sigma(x_i^T \beta) = \frac{1}{n} \sum_{i=1}^n w_i y_i}
#' 
#' where \eqn{\sigma(z) = 1/(1+e^{-z})} is the expit (logistic) function.
#' 
#' **Algorithm Steps:**
#' 1. Compute gradient and Hessian of objective function
#' 2. Compute constraint gradient  
#' 3. Solve the constrained Newton system:
#'    \deqn{\begin{pmatrix} H & C^T \\ C & 0 \end{pmatrix} \begin{pmatrix} \Delta\beta \\ \lambda \end{pmatrix} = \begin{pmatrix} -g \\ -c \end{pmatrix}}
#' 4. Update parameters: \eqn{\beta^{(t+1)} = \beta^{(t)} + \Delta\beta}
#' 5. Check convergence based on parameter change
#' 
#' @examples
#' # Generate example data
#' n <- 100
#' p <- 3
#' X <- matrix(rnorm(n * p), nrow = n)
#' y <- rbinom(n, 1, 0.5)
#' w_label <- rep(1, n)
#' w_unlabel <- rep(1, 200)
#' 
#' # Run optimization
#' beta_opt <- NewtonCLseh(X, y, h = 0.1, w_label, w_unlabel)
#' 
#' # Check convergence
#' if (any(is.na(beta_opt))) {
#'   cat("Optimization failed to converge\n")
#' } else {
#'   cat("Optimized coefficients:", round(beta_opt, 4), "\n")
#' }
#' 
#' @references
#' Nocedal, J. and Wright, S. J. (2006). Numerical Optimization, 2nd edition.
#' Springer-Verlag, New York.
#' 
#' @seealso \code{\link{Beta_estimate}} for the main parameter estimation function
NewtonCLseh <- function(X, y, h, w_label, w_unlabel, max_iter = 100, tol = 1e-4, 
                        uper = 1e+15, initial = rep(0, 1 + ncol(X))) {
  
  # ==============================================================================
  # INPUT VALIDATION AND SETUP
  # ==============================================================================
  
  # Validate inputs
  if (!is.matrix(X)) X <- as.matrix(X)
  if (nrow(X) != length(y)) stop("X and y dimensions don't match")
  if (length(y) != length(w_label)) stop("y and w_label dimensions don't match")
  if (!all(y %in% c(0, 1))) stop("y must be binary (0/1)")
  if (h <= 0) stop("Bandwidth h must be positive")
  if (any(w_label < 0) || any(w_unlabel < 0)) stop("Weights must be non-negative")
  
  # Setup variables
  n_dat <- nrow(X)
  p_coef <- ncol(X) + 1  # Number of coefficients (including intercept)
  
  # Add intercept column to design matrix
  X_design <- cbind(1, X)
  
  # Initialize parameters
  beta <- initial
  if (length(beta) != p_coef) {
    warning("Initial parameter vector length mismatch, using zeros")
    beta <- rep(0, p_coef)
  }
  
  # Initialize convergence tracking
  error <- Inf
  iter <- 0
  conds <- 0  # Condition number
  
  # Constraint setup: E[w_label * expit(X*beta)] = E[w_label * y]
  Xc <- as.vector(rep(1, n_dat))  # Constraint vector (constant constraint)
  
  # Compute weight normalization
  weight <- mean(w_unlabel)
  if (weight <= 0) {
    warning("Mean unlabeled weight is non-positive")
    weight <- 1
  }
  
  # Initial objective function value
  z <- as.vector(X_design %*% beta)
  prob_pred <- Expit(z)
  sqloss <- h * mean(w_label^2 * (y - prob_pred)^2) / weight^2
  
  cat(sprintf("Newton optimization starting: %d parameters, initial loss = %.6f\n", 
              p_coef, sqloss))
  
  # ==============================================================================
  # MAIN NEWTON-RAPHSON ITERATION LOOP  
  # ==============================================================================
  
  while (iter < max_iter && error > tol) {
    
    iter <- iter + 1
    beta_old <- beta
    sqloss_old <- sqloss
    
    # ==========================================================================
    # STEP 1: COMPUTE CURRENT PREDICTIONS AND DERIVATIVES
    # ==========================================================================
    
    # Linear combination
    z <- as.vector(X_design %*% beta)
    
    # Predicted probabilities and derivatives
    prob_pred <- Expit(z)
    prob_deriv <- ExpitDerivative(z)
    
    # ==========================================================================
    # STEP 2: CONSTRUCT LINEARIZED PROBLEM 
    # ==========================================================================
    
    # For Newton's method, we solve a linearized version:
    # y_tilde = y - expit(z) + expit'(z) * z
    # X_tilde = expit'(z) * X
    y_tilde <- y - prob_pred + prob_deriv * z
    X_tilde <- prob_deriv * X_design  # Element-wise multiplication
    
    # ==========================================================================
    # STEP 3: COMPUTE HESSIAN AND GRADIENT
    # ==========================================================================
    
    # Weighted Hessian: X_tilde^T * diag(h * w_label^2 / weight^2) * X_tilde
    weight_factor <- h * w_label^2 / weight^2
    XTX <- crossprod(X_tilde, weight_factor * X_tilde) / n_dat
    
    # Gradient: X_tilde^T * (y_tilde * h * w_label^2 / weight^2)
    XTy <- t(X_tilde) %*% (y_tilde * weight_factor) / n_dat
    
    # ==========================================================================
    # STEP 4: CONSTRAINT HANDLING
    # ==========================================================================
    
    # Constraint matrix: C = Xc^T * (w_label * X_tilde)
    C <- crossprod(Xc, w_label * X_tilde) / n_dat
    
    # Constraint violation: b = Xc^T * (y_tilde * w_label) 
    b <- as.vector(t(Xc) %*% (y_tilde * w_label)) / n_dat
    
    # ==========================================================================
    # STEP 5: SOLVE CONSTRAINED NEWTON SYSTEM
    # ==========================================================================
    
    # Construct the augmented system matrix:
    # [H    C^T] [delta_beta]   [-g]
    # [C    0  ] [lambda   ] = [-b]
    
    n_constraints <- nrow(C)
    augmented_matrix <- rbind(
      cbind(XTX, t(C)),
      cbind(C, matrix(0, nrow = n_constraints, ncol = n_constraints))
    )
    
    augmented_rhs <- c(XTy, b)
    
    # Check condition number
    conds <- kappa(augmented_matrix)
    if (conds > uper) {
      warning(sprintf("Condition number %.2e exceeds threshold %.2e at iteration %d", 
                      conds, uper, iter))
      break
    }
    
    # Solve the linear system
    tryCatch({
      solution <- solve(augmented_matrix, augmented_rhs)
      delta_beta <- solution[1:p_coef]
      lambda <- solution[(p_coef + 1):length(solution)]
    }, error = function(e) {
      warning(paste("Linear system solution failed at iteration", iter, ":", e$message))
      delta_beta <<- rep(0, p_coef)
    })
    
    # ==========================================================================
    # STEP 6: UPDATE PARAMETERS
    # ==========================================================================
    
    beta <- beta_old + delta_beta
    
    # ==========================================================================
    # STEP 7: EVALUATE NEW OBJECTIVE FUNCTION
    # ==========================================================================
    
    z_new <- as.vector(X_design %*% beta)
    prob_pred_new <- Expit(z_new)
    sqloss <- h * mean(w_label^2 * (y - prob_pred_new)^2) / weight^2
    
    # Line search: if objective increased, reduce step size
    step_size <- 1.0
    max_backtrack <- 5
    backtrack_iter <- 0
    
    while (sqloss > sqloss_old && backtrack_iter < max_backtrack) {
      step_size <- step_size * 0.5
      beta <- beta_old + step_size * delta_beta
      
      z_new <- as.vector(X_design %*% beta)
      prob_pred_new <- Expit(z_new)
      sqloss <- h * mean(w_label^2 * (y - prob_pred_new)^2) / weight^2
      
      backtrack_iter <- backtrack_iter + 1
    }
    
    if (backtrack_iter > 0) {
      cat(sprintf("  Iteration %d: backtracking %d steps, step_size = %.3f\n", 
                  iter, backtrack_iter, step_size))
    }
    
    # ==========================================================================
    # STEP 8: CHECK CONVERGENCE
    # ==========================================================================
    
    error <- sqrt(mean((beta - beta_old)^2))
    
    if (iter %% 10 == 0 || iter <= 5) {
      cat(sprintf("  Iteration %d: error = %.6f, loss = %.6f, cond = %.2e\n", 
                  iter, error, sqloss, conds))
    }
    
    # Additional convergence checks
    if (step_size < 1e-8) {
      cat("Convergence: step size too small\n")
      break
    }
    
    if (abs(sqloss - sqloss_old) < tol * abs(sqloss_old)) {
      cat("Convergence: objective change too small\n")
      break
    }
  }
  
  # ==============================================================================
  # FINAL VALIDATION AND CLEANUP
  # ==============================================================================
  
  # Check convergence status
  if (iter >= max_iter) {
    warning(sprintf("Maximum iterations (%d) reached without convergence", max_iter))
  } else {
    cat(sprintf("Converged in %d iterations with error %.6f\n", iter, error))
  }
  
  # Validate final parameters
  if (any(is.na(beta)) || any(is.infinite(beta))) {
    warning("Final parameters contain NA or infinite values")
    return(rep(NA, p_coef))
  }
  
  # Check parameter magnitudes
  max_coef <- max(abs(beta))
  if (max_coef > 50) {
    warning(sprintf("Large coefficient detected: %.2f", max_coef))
  }
  
  # Final constraint violation check
  final_z <- as.vector(X_design %*% beta)
  final_prob <- Expit(final_z)
  constraint_violation <- abs(mean(w_label * final_prob) - mean(w_label * y))
  
  if (constraint_violation > 0.01) {
    warning(sprintf("Large constraint violation: %.6f", constraint_violation))
  }
  
  cat(sprintf("Final: loss = %.6f, constraint violation = %.6f\n", 
              sqloss, constraint_violation))
  
  return(beta)
}

