#!/usr/bin/env Rscript
# ==============================================================================
# Semi-Supervised Survival Analysis: Main Simulation Script
# ==============================================================================
# 
# This script implements and compares multiple semi-supervised survival analysis
# methods for doubly-censored data:
# 
# Methods implemented:
# 1. SEEDS: Semi-supervised Estimation of Event rate with Doubly-censored Survival data  
# 2. CSL: Combined Supervised Learning
# 3. Self-consistency: Non-parametric method for doubly-censored data
#
# ==============================================================================

# Clean environment and set options
options(warn = -1)
rm(list = ls())

# Set working directory (modify as needed)
# setwd("/path/to/your/project")

# ==============================================================================
# PACKAGE MANAGEMENT
# ==============================================================================

# Define required packages
required_packages <- c(
  "survival", "dplyr", "ggplot2", "tidyr", "parallel", 
  "foreach", "doParallel", "knitr", "kableExtra", "gridExtra", 
  "viridis", "mgcv", "splines", "dblcens", "pracma", "MASS",
  "ncvreg"
)

# Function to install and load packages
install_and_load_packages <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cat("Installing package:", pkg, "\n")
      install.packages(pkg, dependencies = TRUE)
    }
    library(pkg, character.only = TRUE)
  }
}

# Install and load all required packages
install_and_load_packages(required_packages)

# ==============================================================================
# SOURCE HELPER FUNCTIONS
# ==============================================================================

# Core utility functions
source('src/utils/HelperFunctions.R')

# Parameter estimation functions  
source('src/estimation/Beta_estimate.R')
source('src/estimation/NewtonCLseh.R')

# Estimation methods
source('src/methods/Supervised_est.R')
source('src/utils/CV_function.R')
source('src/methods/IntrSSL_est_ind.R')
source('src/methods/d011_new.R')

# Data generation functions
source('src/simulation/Generatelogbtcor.R')
source('src/simulation/GenerateCoxtcor.R')

# ==============================================================================
# DATA PROCESSING FUNCTIONS
# ==============================================================================

#' Process simulation data for survival analysis
#' 
#' @param data List containing simulation data with required fields
#' @param adl Lower bandwidth adjustment parameter  
#' @param adh Upper bandwidth adjustment parameter
#' @param verbose Logical, whether to print processing information
#' @return List of processed data components
process_simulation_data <- function(data, adl, adh, verbose = FALSE) {
  
  tryCatch({
    
    # Validate required fields
    required_fields <- c("Label_id", "Delta", "X", "Rtime", "Ltime", "Xstar", "Deltastar")
    missing_fields <- setdiff(required_fields, names(data))
    if (length(missing_fields) > 0) {
      stop(paste("Missing required fields:", paste(missing_fields, collapse = ", ")))
    }
    
    # Extract labeled observations
    label_id <- which(data$Label_id == 1)
    if (length(label_id) == 0) {
      stop("No labeled observations found")
    }
    
    # Extract survival information from labeled data
    delt <- data$Delta[label_id]
    obse <- data$X[label_id]
    
    # Extract censoring times for labeled data
    rcen_labeled <- data$Rtime[label_id]
    lcen_labeled <- data$Ltime[label_id]
    n_l <- length(rcen_labeled)
    
    # Calculate bandwidths for kernel smoothing (labeled data)
    lcen_sd <- sd(lcen_labeled)
    rcen_sd <- sd(rcen_labeled)
    
    # Prevent zero standard deviations
    if (lcen_sd == 0) lcen_sd <- 0.1
    if (rcen_sd == 0) rcen_sd <- 0.1
    
    # Calculate rule-of-thumb bandwidths with adjustment
    lh_labeled <- 1.06 * lcen_sd * n_l^{-(0.2 + adl)}
    rh_labeled <- 1.06 * rcen_sd * n_l^{-(0.2 + adh)}
    
    # Handle unlabeled data
    rcen_unlabeled <- data$Rtime[-label_id]
    lcen_unlabeled <- data$Ltime[-label_id]
    N_unl <- length(rcen_unlabeled)
    
    if (N_unl > 0) {
      lcen_unl_sd <- sd(lcen_unlabeled)
      rcen_unl_sd <- sd(rcen_unlabeled)
      
      if (lcen_unl_sd == 0) lcen_unl_sd <- 0.1
      if (rcen_unl_sd == 0) rcen_unl_sd <- 0.1
      
      lh_unlabeled <- 1.06 * lcen_unl_sd * N_unl^{-(0.2 + adl)}
      rh_unlabeled <- 1.06 * rcen_unl_sd * N_unl^{-(0.2 + adh)}
    } else {
      lh_unlabeled <- lh_labeled
      rh_unlabeled <- rh_labeled
    }
    
    # Process surrogate information
    xstar_all <- data$Xstar
    deltastar_all <- data$Deltastar
    
    # Create indicator variables for censoring types
    deltastar_ind1 <- deltastar_ind2 <- rep(0, length(deltastar_all))
    deltastar_ind1[deltastar_all == 2] <- 1  # Right-censored
    deltastar_ind2[deltastar_all == 3] <- 1  # Left-censored
    deltastar_indv <- cbind(deltastar_ind1, deltastar_ind2)
    
    # Extract additional covariates if available
    base_cov <- if("Z" %in% names(data)) data$Z else NULL
    lcen_ct <- if ("Lcount" %in% names(data)) data$Lcount else NULL
    rcen_ct <- if ("Rcount" %in% names(data)) data$Rcount else NULL
    cova_tim <- if ("Z_time" %in% names(data)) data$Z_time else NULL
    cova_ct <- if ("Z_count" %in% names(data)) data$Z_count else NULL
    
    # Create binary indicators for different censoring types
    ldelt <- rep(0, n_l)  # Left-censored indicator
    ldelt[delt == 3] <- 1
    
    rdelt <- rep(0, n_l)  # Right-censored
    rdelt[delt %in% c(1, 3)] <- 1
    
    # Compile processed data
    processed_data <- list(
      label_id = label_id, 
      delt = delt, 
      obse = obse,
      rcen_labeled = rcen_labeled, 
      lcen_labeled = lcen_labeled,
      lh_labeled = lh_labeled, 
      rh_labeled = rh_labeled,
      lh_unlabeled = lh_unlabeled, 
      rh_unlabeled = rh_unlabeled,
      xstar_all = xstar_all, 
      deltastar_indv = deltastar_indv,
      base_cov = base_cov, 
      lcen_ct = lcen_ct, 
      rcen_ct = rcen_ct,
      cova_tim = cova_tim, 
      cova_ct = cova_ct,
      ldelt = ldelt, 
      rdelt = rdelt,
      lcen_all = data$Ltime, 
      rcen_all = data$Rtime,
      data = data
    )
    
    if (verbose) {
      cat("Data processing completed successfully:\n")
      cat("- Labeled observations:", length(label_id), "\n")
      cat("- Unlabeled observations:", length(xstar_all) - length(label_id), "\n")
      cat("- Total observations:", length(xstar_all), "\n")
      cat("- Left bandwidth (labeled):", round(lh_labeled, 4), "\n")
      cat("- Right bandwidth (labeled):", round(rh_labeled, 4), "\n")
    }
    
    return(processed_data)
    
  }, error = function(e) {
    warning(paste("Data processing failed:", e$message))
    return(NULL)
  })
}

# ==============================================================================
# ESTIMATION METHODS
# ==============================================================================

#' Bootstrap variance estimation for self-consistency method
#' 
#' @param delt Censoring indicators
#' @param obse Observed times
#' @param tim Evaluation time points
#' @param n_bootstrap Number of bootstrap samples
#' @return Vector of bootstrap standard deviations
bootstrap_variance_sc <- function(delt, obse, tim, n_bootstrap = 100) {
  
  bootstrap_estimates <- replicate(n_bootstrap, {
    tryCatch({
      n <- length(delt)
      boot_indices <- sample(1:n, n, replace = TRUE)
      boot_delt <- delt[boot_indices]
      boot_obse <- obse[boot_indices]
      
      sc_boot <- estimate_self_consistency(boot_delt, boot_obse, tim)
      return(sc_boot$survival)
      
    }, error = function(e) {
      return(rep(NA, length(tim)))
    })
  })
  
  bootstrap_sd <- apply(bootstrap_estimates, 1, sd, na.rm = TRUE)
  
  # Handle missing standard deviations
  missing_sd <- is.na(bootstrap_sd)
  if (any(missing_sd)) {
    avg_sd <- mean(bootstrap_sd, na.rm = TRUE)
    if (is.na(avg_sd)) avg_sd <- 0.1
    bootstrap_sd[missing_sd] <- avg_sd
  }
  
  return(bootstrap_sd)
}

#' Self-consistency estimation for doubly-censored data
#' 
#' @param delt Censoring indicators (1=exact, 2=right, 3=left)
#' @param obse Observed times
#' @param tim Evaluation time points
#' @return List containing survival estimates
estimate_self_consistency <- function(delt, obse, tim) {
  
  # Convert censoring indicators to d011 format
  d011_d <- rep(1, length(delt))
  d011_d[delt == 2] <- 0  # Right-censored
  d011_d[delt == 3] <- 2  # Left-censored
  
  tryCatch({
    # Use d011 algorithm for doubly-censored data
    D011 <- d011_new(z = obse, d = d011_d, influence.fun = TRUE)
    stimd011 <- D011$exttime
    surtd011 <- D011$extsurv.Sx
    
    # Interpolate survival estimates at evaluation times
    std <- sapply(tim, function(x) {
      if (x <= min(stimd011)) {
        return(max(surtd011))
      } else if (x >= max(stimd011)) {
        return(min(surtd011))
      } else {
        # Weighted interpolation
        index <- order(abs(stimd011 - x))[1:2]
        weights <- 1 / (abs(stimd011[index] - x) + 1e-8)
        return(sum(surtd011[index] * weights) / sum(weights))
      }
    })
    
    return(list(survival = std))
    
  }, error = function(e) {
    warning(paste("Self-consistency estimation failed:", e$message))
    return(list(survival = rep(NA, length(tim))))
  })
}

#' Wrapper function for supervised estimation
#' 
#' @param processed_data Processed simulation data
#' @param tim Evaluation time points
#' @return List containing survival estimates and standard deviations
estimate_supervised <- function(processed_data, tim) {
  
  tryCatch({
    
    SupSt <- Supervised_est(
      time = tim, 
      processed_data$obse, 
      processed_data$ldelt, 
      processed_data$rdelt,
      lcen = processed_data$lcen_labeled, 
      rcen = processed_data$rcen_labeled,
      lh = processed_data$lh_labeled, 
      rh = processed_data$rh_labeled
    )
    
    # Validate results
    if (is.null(SupSt) || !is.matrix(SupSt) || nrow(SupSt) < 8) {
      stop("Invalid supervised estimation results")
    }
    
    return(list(
      d_surv = SupSt[1, ], d_sd = SupSt[2, ],
      l_surv = SupSt[3, ], l_sd = SupSt[4, ],
      r_surv = SupSt[5, ], r_sd = SupSt[6, ],
      c_surv = SupSt[7, ], c_sd = SupSt[8, ]
    ))
    
  }, error = function(e) {
    warning(paste("Supervised estimation failed:", e$message))
    n_times <- length(tim)
    na_vec <- rep(NA, n_times)
    return(list(
      d_surv = na_vec, d_sd = na_vec,
      l_surv = na_vec, l_sd = na_vec,
      r_surv = na_vec, r_sd = na_vec,
      c_surv = na_vec, c_sd = na_vec
    ))
  })
}

#' Wrapper function for SEEDS estimation
#' 
#' @param processed_data Processed simulation data
#' @param tim Evaluation time points
#' @return List containing survival estimates, standard deviations, and weights
estimate_seeds <- function(processed_data, tim) {
  
  tryCatch({
    seeds_results <- IntrSSL_est(
      time = tim,
      processed_data$base_cov, 
      processed_data$cova_tim, 
      processed_data$cova_ct,
      processed_data$lcen_ct, 
      processed_data$rcen_ct, 
      processed_data$xstar_all,
      processed_data$deltastar_indv, 
      processed_data$label_id, 
      processed_data$lcen_all, 
      processed_data$rcen_all,
      processed_data$obse, 
      processed_data$ldelt, 
      processed_data$rdelt,
      processed_data$lh_labeled, 
      processed_data$lh_unlabeled, 
      processed_data$rh_labeled, 
      processed_data$rh_unlabeled, 
      num_folds = 10
    )
    
    # Validate results
    if (is.null(seeds_results) || !is.matrix(seeds_results) || nrow(seeds_results) < 19) {
      stop("Invalid SEEDS estimation results")
    }
    
    return(list(
      d_surv = seeds_results[1, ], d_surv_intr = seeds_results[2, ],
      l_surv = seeds_results[3, ], l_surv_intr = seeds_results[4, ],
      r_surv = seeds_results[5, ], r_surv_intr = seeds_results[6, ],
      c_surv = seeds_results[7, ], c_surv_intr = seeds_results[8, ],
      d_sd = seeds_results[9, ], d_sd_intr = seeds_results[10, ],
      l_sd = seeds_results[11, ], l_sd_intr = seeds_results[12, ],
      r_sd = seeds_results[13, ], r_sd_intr = seeds_results[14, ],
      c_sd = seeds_results[15, ], c_sd_intr = seeds_results[16, ],
      w_d = seeds_results[17, ], w_l = seeds_results[18, ], w_r = seeds_results[19, ]
    ))
    
  }, error = function(e) {
    warning(paste("SEEDS estimation failed:", e$message))
    n_times <- length(tim)
    na_vec <- rep(NA, n_times)
    return(list(
      d_surv = na_vec, d_surv_intr = na_vec,
      l_surv = na_vec, l_surv_intr = na_vec,
      r_surv = na_vec, r_surv_intr = na_vec,
      c_surv = na_vec, c_surv_intr = na_vec,
      d_sd = na_vec, d_sd_intr = na_vec,
      l_sd = na_vec, l_sd_intr = na_vec,
      r_sd = na_vec, r_sd_intr = na_vec,
      c_sd = na_vec, c_sd_intr = na_vec,
      w_d = na_vec, w_l = na_vec, w_r = na_vec 
    ))
  })
}

# ==============================================================================
# MAIN SIMULATION LOOP
# ==============================================================================

# Initialize simulation parameters
setting <- 1  # Change to 2 for logistic model
seed <- 1     # Random seed for reproducibility

# Initialize result storage
# Self-consistency results
sc_st_all <- sc_sd_all <- NULL  

# Supervised learning results
D_sl_all <- D_sl_sd_all <- NULL 
L_sl_all <- L_sl_sd_all <- NULL
R_sl_all <- R_sl_sd_all <- NULL
C_sl_all <- C_sl_sd_all <- NULL 

# SEEDS results
D_intr_all <- D_intr_sd_all <- NULL
L_intr_all <- L_intr_sd_all <- NULL
R_intr_all <- R_intr_sd_all <- NULL
C_intr_all <- C_intr_sd_all <- NULL
W_D_all <- W_L_all <- W_R_all <- NULL

# Main simulation loop
cat("Starting simulation with setting", setting, "and seed", seed, "\n")

for(i in (seed*100 + 1):(seed*100 + 500)){
  set.seed(i)
  cat(paste('Iteration', i, '\n'))
  
  # Generate data based on setting
  if(setting == 1){
    cat("Using Cox model (Setting 1)\n")
    data <- GenerateCoxtcor(N = 5000, n = 250, phi01 = -12, M = 1.6, a = 1, rate = 0.1)
    tim <- seq(0.82, 2.69, length=50)
    adl <- 0.1; adh <- 0.15
  }
  
  if(setting == 2){
    cat("Using Logistic model (Setting 2)\n")
    data <- Generatelogbtcor(N = 5000, n = 100, rate = 0.1)
    tim <- seq(1.3, 3.356, length=50)
    adl <- 0.1; adh <- 0.15
  }
  
  # ==============================================================================
  # PROCESS DATA AND RUN ESTIMATIONS
  # ==============================================================================
  
  # Process simulation data
  processed_data <- process_simulation_data(data, adl, adh, verbose = FALSE)
  if (is.null(processed_data)) {
    stop("Data processing failed at iteration ", i)
  }
  
  # 1. Self-consistency method
  sc_results <- list(
    survival = estimate_self_consistency(
      processed_data$delt, processed_data$obse, tim
    )$survival,
    sd = bootstrap_variance_sc(
      processed_data$delt, processed_data$obse, tim, n_bootstrap = 100 
    )
  )
  
  sc_st <- sc_results$survival
  sc_sd <- sc_results$sd
  sc_st_all <- rbind(sc_st_all, sc_st)
  sc_sd_all <- rbind(sc_sd_all, sc_sd) 

  
  # 2. CSL (Classical Supervised Learning) method
  csl_results <- estimate_supervised(processed_data, tim)
  
  DSt_sl <- csl_results$d_surv
  DSt_sl_sd <- csl_results$d_sd
  LSt_sl <- csl_results$l_surv
  LSt_sl_sd <- csl_results$l_sd
  RSt_sl <- csl_results$r_surv
  RSt_sl_sd <- csl_results$r_sd
  CSt_sl <- csl_results$c_surv
  CSt_sl_sd <- csl_results$c_sd
  
  # Store CSL results
  D_sl_all <- rbind(D_sl_all, DSt_sl)
  D_sl_sd_all <- rbind(D_sl_sd_all, DSt_sl_sd)
  L_sl_all <- rbind(L_sl_all, LSt_sl)
  L_sl_sd_all <- rbind(L_sl_sd_all, LSt_sl_sd)
  R_sl_all <- rbind(R_sl_all, RSt_sl)
  R_sl_sd_all <- rbind(R_sl_sd_all, RSt_sl_sd)
  C_sl_all <- rbind(C_sl_all, CSt_sl)
  C_sl_sd_all <- rbind(C_sl_sd_all, CSt_sl_sd)
  
  # 5. SEEDS method
  seeds_results <- estimate_seeds(processed_data, tim)
  
  D_intr <- seeds_results$d_surv_intr
  L_intr <- seeds_results$l_surv_intr
  R_intr <- seeds_results$r_surv_intr
  C_intr <- seeds_results$c_surv_intr
  D_intr_sd <- seeds_results$d_sd_intr
  L_intr_sd <- seeds_results$l_sd_intr
  R_intr_sd <- seeds_results$r_sd_intr
  C_intr_sd <- seeds_results$c_sd_intr
  W_D <- seeds_results$w_d
  W_L <- seeds_results$w_l
  W_R <- seeds_results$w_r
  
  # Store SEEDS results
  D_intr_all <- rbind(D_intr_all, D_intr)
  D_intr_sd_all <- rbind(D_intr_sd_all, D_intr_sd)
  L_intr_all <- rbind(L_intr_all, L_intr)
  L_intr_sd_all <- rbind(L_intr_sd_all, L_intr_sd)
  R_intr_all <- rbind(R_intr_all, R_intr)
  R_intr_sd_all <- rbind(R_intr_sd_all, R_intr_sd)
  C_intr_all <- rbind(C_intr_all, C_intr)
  C_intr_sd_all <- rbind(C_intr_sd_all, C_intr_sd)
  W_D_all <- rbind(W_D_all, W_D)
  W_L_all <- rbind(W_L_all, W_L)
  W_R_all <- rbind(W_R_all, W_R)
  
  # ==============================================================================
  # SAVE RESULTS
  # ==============================================================================
  
  # Create output directory if it doesn't exist
  output_dir <- 'data/simulation_results'
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Save all results
  save(D_sl_all, D_sl_sd_all, D_intr_all, D_intr_sd_all, 
       L_sl_all, L_sl_sd_all, L_intr_all, L_intr_sd_all,
       R_sl_all, R_sl_sd_all, R_intr_all, R_intr_sd_all,
       C_sl_all, C_sl_sd_all, C_intr_all, C_intr_sd_all,
       sc_st_all, sc_sd_all, 
       W_D_all, W_L_all, W_R_all,
       file = file.path(output_dir, paste0('setting', setting, 'seed', seed, '.RData')))
}

