# ==============================================================================
# Cox Model Data Generation with Time-Dependent Covariates
# ==============================================================================
# 
# This file generates simulation data for Cox proportional hazards model
# with doubly-censored observations and time-dependent covariates.
#
# The data generation process follows these steps:
# 1. Generate surrogate variable T* ~ Uniform(0, 1/2)
# 2. Generate survival time T from Cox model conditional on T*
# 3. Generate left and right censoring times
# 4. Generate time-dependent covariates using Poisson process
# 5. Create appropriate censoring indicators
# 6. Randomly assign labeled/unlabeled status
#
# Functions included:
# - GenerateCoxtcor(): Main data generation function
# ==============================================================================

#' Generate Cox Model Data with Time-dependent Covariates and Correlation
#' 
#' This function generates simulation data for survival analysis following
#' a Cox proportional hazards model with time-dependent covariates. The data
#' includes doubly-censored observations and a correlation structure between
#' the true survival time and a surrogate variable.
#' 
#' @param N Integer, number of unlabeled observations
#' @param n Integer, number of labeled observations  
#' @param phi01 Numeric, intercept parameter in Cox model (default: -12)
#' @param M Numeric, upper bound for right censoring interval (default: 1.6)
#' @param a Numeric, power parameter in survival time generation (default: 1)
#' @param rate Numeric, rate parameter for Poisson process generating covariates (default: 0.1)
#' 
#' @return List containing the following elements:
#'   \item{X}{Vector of observed survival times (exact, left, or right censored)}
#'   \item{Delta}{Vector of censoring indicators (1=exact, 2=right, 3=left)}
#'   \item{Ltime}{Vector of left censoring times}
#'   \item{Rtime}{Vector of right censoring times}
#'   \item{Xstar}{Vector of observed surrogate times}
#'   \item{Deltastar}{Vector of surrogate censoring indicators}
#'   \item{Label_id}{Vector indicating labeled observations (1) vs unlabeled (NA)}
#'   \item{Lcount}{Vector of time-dependent covariate values at left censoring}
#'   \item{Rcount}{Vector of time-dependent covariate values at right censoring}
#'   \item{Z_time}{List of observation time sequences for each subject}
#'   \item{Z_count}{List of cumulative count sequences for each subject}
#'   
#' @details
#' The data generation follows this process:
#' 
#' 1. **Surrogate Generation**: T* ~ Uniform(0, 1/2)
#' 
#' 2. **Survival Time Generation**: 
#'    T = (-a * log(U) * exp(-phi01 * T*))^(1/4)
#'    where U ~ Uniform(0,1)
#'    
#' 3. **Censoring Times**:
#'    - Left: L ~ Weibull(shape=1.5, scale=1.35)  
#'    - Right: R ~ L + Uniform(0, M)
#'    
#' 4. **Observed Times**:
#'    - X = max(min(T, R), L) (observed survival time)
#'    - X* = max(min(T*, R), L) (observed surrogate time)
#'    
#' 5. **Time-dependent Covariates**:
#'    Generated via Poisson process with rate λ(T) = T
#'    
#' 6. **Censoring Indicators**:
#'    - δ = 1 if L < T < R (exact observation)
#'    - δ = 2 if T ≥ R (right-censored)  
#'    - δ = 3 if T ≤ L (left-censored)
#'    
#' @examples
#' # Generate data with default parameters
#' data <- GenerateCoxtcor(N = 5000, n = 250)
#' 
#' # Check censoring rates
#' table(data$Delta) / length(data$Delta)
#' 
#' # Generate data with custom parameters
#' data <- GenerateCoxtcor(N = 1000, n = 100, phi01 = -10, M = 2.0, rate = 0.05)
#' 
GenerateCoxtcor <- function(N, n, phi01 = -12, M = 1.6, a = 1, rate = 0.1) {
  
  # Input validation
  if (N <= 0 || n <= 0) {
    stop("N and n must be positive integers")
  }
  if (n > N + n) {
    stop("Number of labeled observations (n) cannot exceed total sample size")
  }
  if (rate <= 0) {
    stop("Rate parameter must be positive")
  }
  
  total_n <- N + n
  
  # ==============================================================================
  # STEP 1: Generate Surrogate Variable
  # ==============================================================================
  
  # Surrogate variable T* ~ Uniform(0, 1/2)
  Tstar <- runif(total_n, 0, 1/2)
  
  # ==============================================================================
  # STEP 2: Generate True Survival Times  
  # ==============================================================================
  
  # Random variable for survival time generation
  u <- runif(total_n, 0, 1)
  
  # Generate survival time T from Cox model
  # Formula: T = (-a * log(U) * exp(-phi01 * T*))^(1/4)
  # This creates correlation between T and T*
  Tim <- (-a * log(u) * exp(-(phi01 * Tstar)))^{1/4}
  
  # ==============================================================================
  # STEP 3: Generate Censoring Times
  # ==============================================================================
  
  # Left censoring times from Weibull distribution
  CL <- rweibull(total_n, shape = 1.5, scale = 1.35)
  
  # Right censoring times: L + Uniform(0, M)  
  CR <- runif(total_n, 0, M) + CL
  
  # ==============================================================================
  # STEP 4: Create Observed Times
  # ==============================================================================
  
  # Observed survival time: X = max(min(T, CR), CL)
  X <- pmax(pmin(Tim, CR), CL)
  
  # Observed surrogate time: X* = max(min(T*, CR), CL)  
  Xstar <- pmax(pmin(Tstar, CR), CL)
  
  # ==============================================================================
  # STEP 5: Generate Time-Dependent Covariates
  # ==============================================================================
  
  # Poisson process rate depends on true survival time
  lambda_t <- 1 * Tim
  
  # Generate time-dependent covariate trajectories for each subject
  cova_tim <- sapply(1:total_n, function(k) {
    
    # Number of observations in the interval [CL[k], CR[k]]
    interval_length <- CR[k] - CL[k]
    obs_n <- ceiling(interval_length / rate)
    
    # Generate observation times within the censoring interval
    # Include boundary times and random times within interval
    z_tim <- unique(sort(c(
      CL[k],                                          # Left boundary
      X[k],                                          # Observed time  
      CR[k],                                         # Right boundary
      runif(obs_n, min = CL[k], max = CR[k])       # Random times
    )))
    
    # Generate cumulative counts via Poisson process
    z_ct <- numeric(length(z_tim))
    z_ct[1] <- rpois(1, lambda_t[k] * z_tim[1])
    
    # Incremental Poisson process
    for (j in 2:length(z_tim)) {
      increment <- rpois(1, lambda_t[k] * (z_tim[j] - z_tim[j-1]))
      z_ct[j] <- z_ct[j-1] + increment
    }
    
    return(list(z_tim, z_ct))
  })
  
  # Extract time and count sequences
  Z_tim <- cova_tim[1, ]  # Time sequences
  Z_ct <- cova_tim[2, ]   # Count sequences
  
  # Extract covariate values at censoring boundaries  
  CL_ct <- sapply(Z_ct, function(L) L[[1]])           # Count at left boundary
  CR_ct <- sapply(Z_ct, function(L) L[[length(L)]])   # Count at right boundary
  
  # ==============================================================================
  # STEP 6: Create Censoring Indicators
  # ==============================================================================
  
  # Initialize all as exact observations (δ = 1)
  Delta <- rep(1, total_n)
  
  # Left-censored: T ≤ L, so δ = 3
  Delta[which(CL > Tim)] <- 3
  
  # Right-censored: T ≥ R, so δ = 2  
  Delta[which(CR < Tim)] <- 2
  
  # Same for surrogate variable
  Deltastar <- rep(1, total_n)
  Deltastar[which(CL > Tstar)] <- 3
  Deltastar[which(CR < Tstar)] <- 2
  
  # ==============================================================================
  # STEP 7: Assign Labeled/Unlabeled Status
  # ==============================================================================
  
  # Initialize all as unlabeled
  Label_id <- rep(NA, total_n)
  
  # Randomly select n observations to be labeled
  labeled_indices <- sample.int(total_n, n)
  Label_id[labeled_indices] <- 1
  
  # ==============================================================================
  # STEP 8: Compile Dataset
  # ==============================================================================
  
  dataset <- list(
    # Core survival data
    'X' = X,                    # Observed survival times
    'Delta' = Delta,            # Survival censoring indicators
    'Ltime' = CL,              # Left censoring times
    'Rtime' = CR,              # Right censoring times
    
    # Surrogate data
    'Xstar' = Xstar,           # Observed surrogate times  
    'Deltastar' = Deltastar,   # Surrogate censoring indicators
    
    # Labeling information
    'Label_id' = Label_id,     # Labeled vs unlabeled indicator
    
    # Time-dependent covariates
    'Lcount' = CL_ct,          # Covariate values at left boundary
    'Rcount' = CR_ct,          # Covariate values at right boundary  
    'Z_time' = Z_tim,          # Time sequence lists
    'Z_count' = Z_ct           # Count sequence lists
  )
  
  # ==============================================================================
  # STEP 9: Data Quality Checks
  # ==============================================================================
  
  # Verify censoring consistency
  stopifnot(all(CL <= CR))                    # Left ≤ Right censoring
  stopifnot(all(X >= CL & X <= CR))          # Observed within interval
  stopifnot(all(Xstar >= CL & Xstar <= CR))  # Surrogate within interval
  
  # Check censoring indicator consistency
  exact_obs <- which(Delta == 1)
  left_cens <- which(Delta == 3)  
  right_cens <- which(Delta == 2)
  
  if (length(exact_obs) > 0) {
    stopifnot(all(CL[exact_obs] < X[exact_obs] & X[exact_obs] < CR[exact_obs]))
  }
  if (length(left_cens) > 0) {
    stopifnot(all(X[left_cens] == CL[left_cens]))
  }
  if (length(right_cens) > 0) {
    stopifnot(all(X[right_cens] == CR[right_cens]))
  }
  
  return(dataset)
}

# ==============================================================================
# UTILITY FUNCTIONS FOR COX DATA GENERATION
# ==============================================================================

#' Print summary of generated Cox data
#' 
#' @param data List returned by GenerateCoxtcor()
print_cox_data_summary <- function(data) {
  
  cat("=== Cox Model Data Generation Summary ===\n")
  cat("Total observations:", length(data$X), "\n")
  cat("Labeled observations:", sum(data$Label_id == 1, na.rm = TRUE), "\n")
  cat("Unlabeled observations:", sum(is.na(data$Label_id)), "\n\n")
  
  cat("Censoring distribution:\n")
  cens_table <- table(data$Delta)
  cens_props <- prop.table(cens_table)
  for (i in 1:length(cens_table)) {
    status <- switch(names(cens_table)[i],
                     "1" = "Exact",
                     "2" = "Right-censored", 
                     "3" = "Left-censored")
    cat(sprintf("  %s: %d (%.1f%%)\n", status, cens_table[i], 100*cens_props[i]))
  }
  
  cat("\nSurvival time summary:\n")
  cat("  Min:", round(min(data$X), 3), "\n")
  cat("  Max:", round(max(data$X), 3), "\n") 
  cat("  Mean:", round(mean(data$X), 3), "\n")
  cat("  Median:", round(median(data$X), 3), "\n")
  
  cat("\nSurrogate correlation:", round(cor(data$X, data$Xstar), 3), "\n")
}

# ==============================================================================
# END OF COX DATA GENERATION
# ==============================================================================