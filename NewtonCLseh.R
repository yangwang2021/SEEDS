NewtonCLseh <- function(X, y, h, w_label, w_unlabel, max_iter = 100, tol = 1e-4, 
                        uper = 1e+15, initial = rep(0, 1 + ncol(X))){
  error <- Inf
  iter <- 0
  conds <- 0
  
  beta <- initial
  X <- cbind(1, X)
  n_dat <- nrow(X)
  ###to get the constrain condition
  Xc <- as.vector(rep(1, n_dat))
  
  weight <- mean(w_unlabel)
  # the original minimizaiton
  sqloss <- h * mean(w_label^2 * (y - as.vector(Expit(X %*% beta)))^2) / weight^2
  
  while(iter < max_iter & error > tol){

    iter <- iter + 1
    beta_old <- beta
    sqloss_old <- sqloss
    # Update the minimization by using newton algorithm
    z <- as.vector(X %*% beta)
    y_ <- y - Expit(z) + ExpitDerivative(z) * z
    x_ <- as.vector(ExpitDerivative(z)) * X
    xTx <- crossprod(x_, h * w_label^2 * weight^{-2} * x_) / n_dat
    xTy <- t(x_) %*% (as.vector(y_) * h * as.vector(w_label)^2 * weight^{-2}) / n_dat
    C <- crossprod(Xc, as.vector(w_label) * x_) / n_dat
    b <- as.vector(t(Xc) %*% (as.vector(y_) * as.vector(w_label))) / n_dat
    mat_bind <- cbind(rbind(xTx, C), rbind(t(C), matrix(0, nrow(C), nrow(C))))
    vec_bind <- c(xTy, b)
    conds <- cond(mat_bind)
    if(conds > uper){break}  # if cond(mat_bind) large than 1e+14, the matrix mat_bind is singularity
    solution_bind <- solve(mat_bind) %*% vec_bind
    beta <- solution_bind[1:ncol(X)]
    
    sqloss <- h * mean(w_label^2 * (y - as.vector(Expit(X %*% beta)))^2) / weight^2
  
    if (sqloss_old < sqloss){
      beta <- beta_old
    }
    error <- sqrt(mean((beta - beta_old)^2))
  }
  return(beta)
}

