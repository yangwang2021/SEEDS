# ==============================================================================
# Logistic Model Data Generation with Time-Dependent Covariates
# ==============================================================================
# 
# This file generates simulation data for logistic regression model with 
# doubly-censored observations and time-dependent covariates. The logistic
# model is used as an alternative to the Cox model for survival analysis.
#
# The data generation process follows these steps:
# 1. Generate surrogate variable T* ~ Uniform(-1, 1)
# 2. Generate baseline covariate Z ~ Normal(5, 1)
# 3. Generate survival time T from logistic model conditional on T* and Z
# 4. Generate left and right censoring times
# 5. Generate time-dependent covariates using Poisson process
# 6. Create appropriate censoring indicators
# 7. Randomly assign labeled/unlabeled status
#
# Functions included:
# - Generatelogbtcor(): Main data generation function
# ==============================================================================

#' Generate Logistic Model Data with Time-dependent Covariates and Correlation
#' 
#' This function generates simulation data for survival analysis following
#' a logistic regression model with time-dependent covariates. The data includes
#' doubly-censored observations and correlation between the true survival time
#' and a surrogate variable.
#' 
#' @param N Integer, number of unlabeled observations
#' @param n Integer, number of labeled observations
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
#'   \item{Z}{Vector of baseline covariates}
#'   \item{Lcount}{Vector of time-dependent covariate values at left censoring}
#'   \item{Rcount}{Vector of time-dependent covariate values at right censoring}
#'   \item{Z_time}{List of observation time sequences for each subject}
#'   \item{Z_count}{List of cumulative count sequences for each subject}
#'   
#' @details
#' The data generation follows this process:
#' 
#' 1. **Surrogate Generation**: T* ~ Uniform(-1, 1)
#' 
#' 2. **Baseline Covariate**: Z ~ Normal(μ=5, σ=1)
#' 
#' 3. **Survival Time Generation**: 
#'    T ~ Logistic(location = 2 + 1.28*T* + 0.1*Z, scale = 0.3)
#'    This creates correlation ≈ 0.8 between T and T*
#'    
#' 4. **Censoring Times**:
#'    - Left: L ~ Weibull(shape=1.9, scale=1.9)
#'    - Right: R ~ L + Uniform(0, 2.3)
#'    
#' 5. **Observed Times**:
#'    - X = max(min(T, R), L) (observed survival time)
#'    - X* = max(min(T*, R), L) (observed surrogate time)
#'    
#' 6. **Time-dependent Covariates**:
#'    Generated via Poisson process with rate λ(T) = |T|
#'    
#' 7. **Censoring Indicators**:
#'    - δ = 1 if L < T < R (exact observation)
#'    - δ = 2 if T ≥ R (right-censored)
#'    - δ = 3 if T ≤ L (left-censored)
#'    
#' @examples
#' # Generate data with default parameters
#' data <- Generatelogbtcor(N = 5000, n = 250)
#' 
#' # Check censoring rates
#' table(data$Delta) / length(data$Delta)
#' 
#' # Check correlation between T and T*
#' cor(data$X, data$Xstar)
#' 
Generatelogbtcor <- function(N, n, rate = 0.1) {
  
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
  
  # Surrogate variable T* ~ Uniform(-1, 1)
  # Wider range than Cox model to create different correlation structure
  Tstar <- runif(total_n, -1, 1)
  
  # ==============================================================================
  # STEP 2: Generate Baseline Covariate
  # ==============================================================================
  
  # Baseline covariate Z ~ Normal(5, 1)
  # This represents a continuous baseline characteristic (e.g., age, biomarker level)
  Z <- rnorm(total_n, mean = 5, sd = 1)
  
  # ==============================================================================
  # STEP 3: Generate True Survival Times
  # ==============================================================================
  
  # Generate survival time T from logistic distribution
  # Location parameter depends on surrogate T* and baseline covariate Z
  # Formula: T ~ Logistic(2 + 1.28*T* + 0.1*Z, scale = 0.3)
  # Coefficients chosen to achieve approximately 0.8 correlation with T*
  location_param <- 2 + 1.28 * Tstar + 0.1 * Z
  scale_param <- 0.3
  
  Tim <- rlogis(total_n, location = location_param, scale = scale_param)
  
  # ==============================================================================
  # STEP 4: Generate Censoring Times
  # ==============================================================================
  
  # Left censoring times from Weibull distribution
  # Different parameters than Cox model to create appropriate censoring rates
  CL <- rweibull(total_n, shape = 1.9, scale = 1.9)
  
  # Right censoring times: L + Uniform(0, 2.3)
  # Interval width chosen to balance exact vs censored observations
  CR <- runif(total_n, 0, 2.3) + CL
  
  # ==============================================================================
  # STEP 5: Create Observed Times
  # ==============================================================================
  
  # Observed survival time: X = max(min(T, CR), CL)
  # This ensures X is always within the censoring interval [CL, CR]
  X <- pmax(pmin(Tim, CR), CL)
  
  # Observed surrogate time: X* = max(min(T*, CR), CL)
  Xstar <- pmax(pmin(Tstar, CR), CL)
  
  # ==============================================================================
  # STEP 6: Generate Time-Dependent Covariates
  # ==============================================================================
  
  # Poisson process rate depends on absolute value of true survival time
  # This creates a reasonable non-negative rate
  lambda_t <- abs(Tim)
  
  # Generate time-dependent covariate trajectories for each subject
  cova_tim <- sapply(1:total_n, function(k) {
    
    # Number of observations in the interval [CL[k], CR[k]]
    # More observations for longer intervals
    interval_length <- CR[k] - CL[k]
    obs_n <- ceiling(interval_length / rate)
    
    # Generate observation times within the censoring interval
    # Include key boundary times plus random intermediate times
    z_tim <- unique(sort(c(
      CL[k],                                          # Left boundary
      X[k],                                          # Observed survival time
      CR[k],                                         # Right boundary
      runif(obs_n, min = CL[k], max = CR[k])       # Random observation times
    )))
    
    # Generate cumulative counts via Poisson process
    # Start with initial count at first time point
    z_ct <- numeric(length(z_tim))
    z_ct[1] <- rpois(1, lambda_t[k] * z_tim[1])
    
    # Generate incremental counts for subsequent time points
    for (j in 2:length(z_tim)) {
      # Increment based on time difference and rate
      time_increment <- z_tim[j] - z_tim[j-1]
      count_increment <- rpois(1, lambda_t[k] * time_increment)
      z_ct[j] <- z_ct[j-1] + count_increment
    }
    
    return(list(z_tim, z_ct))
  })
  
  # Extract time and count sequences
  Z_tim <- cova_tim[1, ]  # List of time observation vectors
  Z_ct <- cova_tim[2, ]   # List of cumulative count vectors
  
  # Extract covariate values at censoring boundaries
  CL_ct <- sapply(Z_ct, function(L) L[[1]])           # Count at left boundary
  CR_ct <- sapply(Z_ct, function(L) L[[length(L)]])   # Count at right boundary
  
  # ==============================================================================
  # STEP 7: Create Censoring Indicators
  # ==============================================================================
  
  # Initialize all observations as exact (δ = 1)
  Delta <- rep(1, total_n)
  
  # Left-censored: True time ≤ Left boundary, so δ = 3
  Delta[which(CL > Tim)] <- 3
  
  # Right-censored: True time ≥ Right boundary, so δ = 2
  Delta[which(CR < Tim)] <- 2
  
  # Same process for surrogate variable
  Deltastar <- rep(1, total_n)
  Deltastar[which(CL > Tstar)] <- 3
  Deltastar[which(CR < Tstar)] <- 2
  
  # ==============================================================================
  # STEP 8: Assign Labeled/Unlabeled Status
  # ==============================================================================
  
  # Initialize all observations as unlabeled
  Label_id <- rep(NA, total_n)
  
  # Randomly select n observations to be labeled
  labeled_indices <- sample.int(total_n, n)
  Label_id[labeled_indices] <- 1
  
  # ==============================================================================
  # STEP 9: Compile Dataset
  # ==============================================================================
  
  dataset <- list(
    # Core survival data
    'X' = X,                    # Observed survival times
    'Delta' = Delta,            # Survival censoring indicators (1/2/3)
    'Ltime' = CL,              # Left censoring times
    'Rtime' = CR,              # Right censoring times
    
    # Surrogate data
    'Xstar' = Xstar,           # Observed surrogate times
    'Deltastar' = Deltastar,   # Surrogate censoring indicators
    
    # Labeling information
    'Label_id' = Label_id,     # Labeled (1) vs unlabeled (NA)
    
    # Baseline covariates
    'Z' = Z,                   # Baseline continuous covariate
    
    # Time-dependent covariates
    'Lcount' = CL_ct,          # Covariate values at left boundary
    'Rcount' = CR_ct,          # Covariate values at right boundary
    'Z_time' = Z_tim,          # Lists of observation times
    'Z_count' = Z_ct           # Lists of cumulative counts
  )
  
  # ==============================================================================
  # STEP 10: Data Quality Checks
  # ==============================================================================
  
  # Verify basic consistency
  stopifnot(all(CL <= CR))                    # Left ≤ Right censoring
  stopifnot(all(X >= CL & X <= CR))          # Observed within interval
  stopifnot(all(Xstar >= CL & Xstar <= CR))  # Surrogate within interval
  
  # Check censoring indicator consistency
  exact_obs <- which(Delta == 1)
  left_cens <- which(Delta == 3)
  right_cens <- which(Delta == 2)
  
  # Exact observations should be strictly within censoring interval
  if (length(exact_obs) > 0) {
    stopifnot(all(CL[exact_obs] < X[exact_obs] & X[exact_obs] < CR[exact_obs]))
  }
  
  # Left-censored observations should be at left boundary
  if (length(left_cens) > 0) {
    stopifnot(all(abs(X[left_cens] - CL[left_cens]) < 1e-10))
  }
  
  # Right-censored observations should be at right boundary
  if (length(right_cens) > 0) {
    stopifnot(all(abs(X[right_cens] - CR[right_cens]) < 1e-10))
  }
  
  # Verify reasonable correlation between T and T*
  correlation <- cor(X, Xstar)
  if (correlation < 0.5) {
    warning(sprintf("Low correlation between survival and surrogate times: %.3f", correlation))
  }
  
  return(dataset)
}

# ==============================================================================
# UTILITY FUNCTIONS FOR LOGISTIC DATA GENERATION
# ==============================================================================

#' Print summary of generated logistic data
#' 
#' @param data List returned by Generatelogbtcor()
print_logistic_data_summary <- function(data) {
  
  cat("=== Logistic Model Data Generation Summary ===\n")
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
  
  cat("\nBaseline covariate (Z) summary:\n")
  cat("  Min:", round(min(data$Z), 3), "\n")
  cat("  Max:", round(max(data$Z), 3), "\n")
  cat("  Mean:", round(mean(data$Z), 3), "\n")
  cat("  SD:", round(sd(data$Z), 3), "\n")
  
  cat("\nCorrelations:\n")
  cat("  Survival-Surrogate:", round(cor(data$X, data$Xstar), 3), "\n")
  cat("  Survival-Baseline:", round(cor(data$X, data$Z), 3), "\n")
  cat("  Surrogate-Baseline:", round(cor(data$Xstar, data$Z), 3), "\n")
  
  # Time-dependent covariate summary
  if (!is.null(data$Z_count)) {
    all_counts <- unlist(data$Z_count)
    cat("\nTime-dependent covariate summary:\n")
    cat("  Min count:", min(all_counts), "\n")
    cat("  Max count:", max(all_counts), "\n")
    cat("  Mean count:", round(mean(all_counts), 2), "\n")
    
    # Average number of observations per subject
    obs_per_subject <- sapply(data$Z_time, length)
    cat("  Avg observations per subject:", round(mean(obs_per_subject), 1), "\n")
  }
}

#' Compare logistic and Cox data characteristics
#' 
#' @param logistic_data List from Generatelogbtcor()
#' @param cox_data List from GenerateCoxtcor() 
compare_data_characteristics <- function(logistic_data, cox_data) {
  
  cat("=== Data Generation Comparison ===\n")
  
  # Censoring rates
  log_cens <- table(logistic_data$Delta)
  cox_cens <- table(cox_data$Delta)
  
  cat("\nCensoring rates:\n")
  cat("Model      Exact   Right   Left\n")
  cat("---------- ------- ------- -------\n")
  cat(sprintf("Logistic   %5.1f%% %5.1f%% %5.1f%%\n", 
              100*log_cens[1]/sum(log_cens),
              100*log_cens[2]/sum(log_cens), 
              100*log_cens[3]/sum(log_cens)))
  cat(sprintf("Cox        %5.1f%% %5.1f%% %5.1f%%\n",
              100*cox_cens[1]/sum(cox_cens),
              100*cox_cens[2]/sum(cox_cens),
              100*cox_cens[3]/sum(cox_cens)))
  
  # Survival time characteristics
  cat("\nSurvival time characteristics:\n")
  cat("Model      Min     Max     Mean    SD\n")
  cat("---------- ------- ------- ------- -------\n")
  cat(sprintf("Logistic   %7.3f %7.3f %7.3f %7.3f\n",
              min(logistic_data$X), max(logistic_data$X),
              mean(logistic_data$X), sd(logistic_data$X)))
  cat(sprintf("Cox        %7.3f %7.3f %7.3f %7.3f\n",
              min(cox_data$X), max(cox_data$X),
              mean(cox_data$X), sd(cox_data$X)))
  
  # Correlations
  cat("\nSurvival-Surrogate correlations:\n")
  cat(sprintf("Logistic: %.3f\n", cor(logistic_data$X, logistic_data$Xstar)))
  cat(sprintf("Cox:      %.3f\n", cor(cox_data$X, cox_data$Xstar)))
}

# ==============================================================================
# END OF LOGISTIC DATA GENERATION
# ==============================================================================