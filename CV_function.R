
CV_function <- function(ddat_unlabel, ddat_label, ldat_unlabel, ldat_label, rdat_unlabel, rdat_label, 
                        dw_label, dw_unlabel, dy, dh, lw_label, lw_unlabel, ly, lh, rw_label, rw_unlabel,
                        ry, rh, num_folds, a){
  
  p_basis <- ncol(ddat_label)        # the dimension of the data set
  n_label <- nrow(ddat_label)        # the number of the labeled data
  n_unlabel <- nrow(ddat_unlabel)    # the number of the unlabeled data 
  
  # obtain the survival estimators for the three different types of data
St_est <- function(dat_unlabel, dat_label, w_label, w_unlabel, y, h){
    Beta_n <- Beta_estimate(dat_label, y, h, w_label, w_unlabel)  # use the labeled data to estimate beta for three different types of data
    beta_n <- Beta_n[(1 + a*(p_basis+1)) : (p_basis + 1 + a*(p_basis+1))]
    St <- mean(w_unlabel * as.vector(Expit(cbind(1, dat_unlabel) %*% beta_n))) / mean(w_unlabel) # use the unlabeled data to make marginalization for the survival estimator 
    return(St)
  }
  DSt <- St_est(ddat_unlabel, ddat_label, dw_label, dw_unlabel, dy, dh)  
  LSt <- St_est(ldat_unlabel, ldat_label, lw_label, lw_unlabel, ly, lh)
  RSt <- St_est(rdat_unlabel, rdat_label, rw_label, rw_unlabel, ry, rh)
  
  # randomly split the n_label data into K folds
  ind_cv <- cv_split(n_label, num_folds)
  
  Var_cv <- sapply(1:num_folds, function(k){
    # use the K-1 folds data to estimate beta of the three different type of data
    inds_fold <- as.vector(unlist(ind_cv[-k]))
    # obtain the bate estimator for the exact observed information
    Dbeta_K <- Beta_estimate(ddat_label[inds_fold, ], dy[inds_fold], dh, dw_label[inds_fold], dw_unlabel)
    dbeta_K <- Dbeta_K[(1 + a*(p_basis+1)) : (p_basis + 1 + a*(p_basis+1))]
    # obtain the bate estimator for the left censored data
    Lbeta_K <- Beta_estimate(ldat_label[inds_fold, ], ly[inds_fold], lh, lw_label[inds_fold], lw_unlabel)
    lbeta_K <- Lbeta_K[(1 + a*(p_basis+1)) : (p_basis + 1 + a*(p_basis+1))]
    # obtain the bate estimator for the right censored data
    Rbeta_K <- Beta_estimate(rdat_label[inds_fold, ], ry[inds_fold], rh, rw_label[inds_fold], rw_unlabel)
    rbeta_K <- Rbeta_K[(1 + a*(p_basis+1)) : (p_basis + 1 + a*(p_basis+1))]
    
    # use the k-th fold to estimate the survival estimators and variance and covariance estimator for the three different type of data
    inds_foldk <- as.vector(unlist(ind_cv[k]))
    DImp_k <- as.vector(Expit(cbind(1, ddat_label[inds_foldk, ]) %*% dbeta_K))
    LImp_k <- as.vector(Expit(cbind(1, ldat_label[inds_foldk, ]) %*% lbeta_K))
    RImp_k <- as.vector(Expit(cbind(1, rdat_label[inds_foldk, ]) %*% rbeta_K))
    # obtain the variance estimator for the exact observed information
    DVar_k <- dh * mean(dw_label[inds_foldk]^2 * (dy[inds_foldk] - DImp_k)^2) / (mean(dw_unlabel))^2
    # obtain the variance estimator for the left censored data
    LVar_k <- lh * mean(lw_label[inds_foldk]^2 * (ly[inds_foldk] - LImp_k)^2) / (mean(lw_unlabel))^2
    # obtain the variance estimator for the right censored data
    RVar_k <- rh * mean(rw_label[inds_foldk]^2 * (ry[inds_foldk] - RImp_k)^2) / (mean(rw_unlabel))^2
    # obtain the covariance estimator of the exact observed inf. and the left censored data
    DLCov_k <- mean(dw_label[inds_foldk] * (dy[inds_foldk] - DImp_k) * 
                      lw_label[inds_foldk] * (ly[inds_foldk] - LImp_k)) / (mean(dw_unlabel) * mean(lw_unlabel))
    # obtain the covariance estimator of the exact observed inf. and the right censored data
    DRCov_k <- mean(dw_label[inds_foldk] * (dy[inds_foldk] - DImp_k) * 
                      rw_label[inds_foldk] * (ry[inds_foldk] - RImp_k)) / (mean(dw_unlabel) * mean(rw_unlabel))
    # obtain the covariance estimator of the left and the right censored data
    LRCov_k <- mean(lw_label[inds_foldk] * (ly[inds_foldk] - LImp_k) * 
                      rw_label[inds_foldk] * (ry[inds_foldk] - RImp_k)) / (mean(lw_unlabel) * mean(rw_unlabel))
    
    return(c(DVar_k, LVar_k, RVar_k, DLCov_k, DRCov_k, LRCov_k))
  })
  
  Var_mean <- apply(Var_cv, 1, mean)
  DVar <- Var_mean[1] / (dh * n_label)
  LVar <- Var_mean[2] / (lh * n_label)
  RVar <- Var_mean[3] / (rh * n_label)
  DLCov <- Var_mean[4] / n_label
  DRCov <- Var_mean[5] / n_label
  LRCov <- Var_mean[6] / n_label
  
  # remove the survival estimator and standard deviation estimator of the exact observed information if there no observed time
  VS <- c(DSt, LSt, RSt)
  Mcov <- matrix(c(DVar, DLCov, DRCov, DLCov, LVar, LRCov, DRCov, LRCov, RVar), ncol=3, byrow=T)
  
  codi <- c(sum(dy), sum(ly), sum(ry))
  mis_id <- NULL
  if(0 %in% codi){
    mis_id <- which(codi == 0)
    VS <- VS[-mis_id]
    Mcov <- as.matrix(Mcov[-mis_id,-mis_id])
  }

  epsilon <- n_label^{-1/2}
  Mcov_inv <- solve(Mcov + diag(ncol(Mcov)) * epsilon)
  # obtain the combination weight
  ones <- rep(1, ncol(Mcov))
  Denom <- as.vector(t(ones) %*% Mcov_inv %*% ones)
  Comb_w <- (t(ones) %*% Mcov_inv) / Denom
  # obtain the combination estimator
  CSt <- as.vector(Comb_w %*% VS)
  
  # remove the weights of the exact observed information if there no observed time
  New_Comb_w <- rep(0, 3)
  if(length(mis_id) == 0){
    New_Comb_w <- Comb_w}else{New_Comb_w[-mis_id] <- Comb_w}
  
  # obtain the combination variance estimator by using CV method
  CVar_cv <- sapply(1:num_folds, function(k){
    # use the K-1 folds data to estimate the beta of the three different types of data
    inds_fold <- as.vector(unlist(ind_cv[-k]))
    # obtain the bate estimator for the exact observed information
    Dbeta_K <- Beta_estimate(ddat_label[inds_fold, ], dy[inds_fold], dh, dw_label[inds_fold], dw_unlabel)
    dbeta_K <- Dbeta_K[(1 + a*(p_basis+1)) : (p_basis + 1 + a*(p_basis+1))]
    # obtain the bate estimator for the left censored data
    Lbeta_K <- Beta_estimate(ldat_label[inds_fold, ], ly[inds_fold], lh, lw_label[inds_fold], lw_unlabel)
    lbeta_K <- Lbeta_K[(1 + a*(p_basis+1)) : (p_basis + 1 + a*(p_basis+1))]
    # obtain the bate estimator for the right censored data
    Rbeta_K <- Beta_estimate(rdat_label[inds_fold, ], ry[inds_fold], rh, rw_label[inds_fold], rw_unlabel)
    rbeta_K <- Rbeta_K[(1 + a*(p_basis+1)) : (p_basis + 1 + a*(p_basis+1))]
    
    # use the k-th fold to estimate the combination variance estimator
    inds_foldk <- as.vector(unlist(ind_cv[k]))
    DImp_k <- as.vector(Expit(cbind(1, ddat_label[inds_foldk, ]) %*% dbeta_K))
    LImp_k <- as.vector(Expit(cbind(1, ldat_label[inds_foldk, ]) %*% lbeta_K))
    RImp_k <- as.vector(Expit(cbind(1, rdat_label[inds_foldk, ]) %*% rbeta_K))
    CVar_k <- mean((New_Comb_w[1] * mean(dw_unlabel)^{-1} * dw_label[inds_foldk] * (dy[inds_foldk] - DImp_k) +
                    New_Comb_w[2] * mean(lw_unlabel)^{-1} * lw_label[inds_foldk] * (ly[inds_foldk] - LImp_k) +
                    New_Comb_w[3] * mean(rw_unlabel)^{-1} * rw_label[inds_foldk] * (ry[inds_foldk] - RImp_k))^2) 
    
    return(CVar_k)
  })
  CVarn <- mean(CVar_cv)
  CVar <- CVarn / n_label
  
  return(list(DSt = DSt, LSt = LSt, RSt = RSt, CSt = CSt, D_sd = sqrt(DVar), L_sd = sqrt(LVar),
              R_sd = sqrt(RVar), C_sd = sqrt(CVar)))
}



