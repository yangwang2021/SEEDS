Supervised_est <- function(time, obse, ldelt, rdelt, lcen, rcen, lh, rh){
  
  re <- sapply(time, function(x){
    # obtain I(L < t <= X) at fixed time point
    Inter_obse <- as.numeric(I(x <= obse & x > lcen))
    # obtain I(L < t <= U) at fixed time point
    Inter_cen <- as.numeric(I(x <= rcen & x > lcen))
    
    # obtain the survival estimator and standard deviation estimator for the exact observed information 
    if(sum(Inter_obse) == 0){DSt_sl <- 0; DSt_sl_var <- 0;}else{   #### the survival estimator and standard deviation equal to zero if there no observed time
      DSt_sl <- sum(Inter_obse) / sum(Inter_cen) 
      DSt_sl_var <- mean(Inter_cen^2 * (Inter_obse - DSt_sl)^2) / (mean(Inter_cen)^2 * length(obse))
    }
    DSt_sl; sqrt(DSt_sl_var)
    
    # normalized the kernel weights for the left censored data 
    lw <- ker(lcen - x, lh)
    if(mean(lw) == 0){ lw <- rep(0, length(lw)) }else{ lw <- lw / mean(lw)}
    # obtain the survival estimator and standard deviation estimator for the left censored data
    if(sum(ldelt) ==0 ){LSt_sl <- 0; LSt_sl_var <- 0;}else{
      LSt_sl <- sum(lw * (1 - ldelt)) / sum(lw)
      LSt_sl_var <- lh * mean(lw^2 * ((1 - ldelt) - LSt_sl)^2) / (mean(lw)^2 * lh * length(obse))
    }
    LSt_sl; sqrt(LSt_sl_var)
    
    # normalized the kernel weights for the right censored data 
    rw <- ker(rcen - x, rh)
    if(mean(rw) == 0){ rw <- rep(0, length(rw)) }else{ rw <- rw / mean(rw)}
    # obtain the survival estimator and standard deviation estimator for the right censored data
    if(sum(rdelt) == 0){RSt_sl <- 0; RSt_sl_var <- 0;}else{
      RSt_sl <- sum(rw * (1 - rdelt)) / sum(rw)
      RSt_sl_var <- rh * mean(rw^2 * ((1 - rdelt) - RSt_sl)^2) / (mean(rw)^2 * rh * length(obse))
    }
    RSt_sl; sqrt(RSt_sl_var)
    
    # obtain the covariance of estimators of the three different types of data
    St_sl_vec <- c(DSt_sl, LSt_sl, RSt_sl)
    Cov_dl <- mean((Inter_cen * (Inter_obse - DSt_sl)) * (lw * ((1 - ldelt) - LSt_sl))) / (mean(Inter_cen) *  mean(lw) * length(obse))
    Cov_dr <- mean((Inter_cen * (Inter_obse - DSt_sl)) * (rw * ((1 - rdelt) - RSt_sl))) / (mean(Inter_cen) * mean(rw)  * length(obse))
    Cov_lr <- mean((lw * ((1 - ldelt) - LSt_sl)) * (rw * ((1 - rdelt) - RSt_sl))) / (mean(lw) * mean(rw)  * length(obse))
    Cov_sl <- matrix(c(DSt_sl_var, Cov_dl, Cov_dr, Cov_dl, LSt_sl_var, Cov_lr, Cov_dr, Cov_lr, RSt_sl_var), nrow = 3)
    
   # remove the survival estimator and standard deviation estimator of the exact observed information if there no observed time
    codi <- c(sum(Inter_obse), sum(ldelt), sum(rdelt))
    mis_id <- NULL
    if(0 %in% codi){
      mis_id <- which(codi == 0)
      St_sl_vec <- St_sl_vec[-mis_id]
      Cov_sl <- as.matrix(Cov_sl[-mis_id, -mis_id])
    }
    
    # add a penalty to the covariance matrix
     epsilon <- length(obse)^{-1} 
     Cov_inv <- solve(Cov_sl + epsilon * diag(ncol(Cov_sl)))
     
    # obtain the combination weight
     ones <- rep(1, ncol(Cov_sl))
     Denom <- as.vector(t(ones) %*% Cov_inv %*% ones)
     CSt_w <- t(ones) %*% Cov_inv / Denom
    
    # obtain the combination estimator
    CSt_sl <- as.vector(CSt_w %*% St_sl_vec) 
    
    # obtain the combination variance estimator
    CSt_wei <- rep(0, 3)
    if(length(mis_id) == 0){
      CSt_wei <- CSt_w}else{CSt_wei[-mis_id] <- CSt_w}
    
    
    CSt_sl_var <- mean((CSt_wei[1] * mean(Inter_cen)^{-1} * Inter_cen * (Inter_obse - DSt_sl) + 
                        CSt_wei[2] * mean(lw)^{-1} * lw * ((1 - ldelt) - LSt_sl) +
                        CSt_wei[3] * mean(rw)^{-1} * rw * ((1 - rdelt) - RSt_sl))^2) / length(obse)
  
    CSt_w; St_sl_vec; CSt_sl; sqrt(CSt_sl_var); sqrt(DSt_sl_var); sqrt(LSt_sl_var); sqrt(RSt_sl_var)
    
    return(c('DSt_sl' = DSt_sl, 'DSt_sl_sd' = sqrt(DSt_sl_var), 'LSt_sl' = LSt_sl, 
            'LSt_sl_sd' = sqrt(LSt_sl_var), 'RSt_sl' = RSt_sl, 'RSt_sl_sd' = sqrt(RSt_sl_var), 
            'CSt_sl' = CSt_sl, 'CSt_sl_sd' = sqrt(CSt_sl_var)))
  })
  return(re)
}
