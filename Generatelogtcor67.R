Generatelogtcor67 <- function(N, n, rate){

  Tstar <- runif(N+n, -1, 1)
  Tim <- rlogis(N+n, 3 + 1 * Tstar, 0.35)
  CL <- rweibull(N+n, 2.5, 2.45)
  CR <- runif(N+n, 0, 2.25) + CL
  X <- pmax(pmin(Tim, CR), CL)
  Xstar <- pmax(pmin(Tstar, CR), CL)
  
  lambda_t = abs(Tim) / 1
  ###generate time-dependent covariate by using poisson process
  cova_tim <- sapply(1 : (n+N), function(k){
    obs_n <- ceiling(CR[k] - CL[k]) / rate
    z_tim <- unique(sort(c(CL[k], X[k], CR[k], runif(obs_n, min = CL[k], max = CR[k]))))
    z_ct <- rep()
    z_ct[1] <- rpois(1, lambda_t[k] * z_tim[1])
    for (j in 2 : length(z_tim)) {
      z_ct[j] <- z_ct[j-1] + rpois(1, lambda_t[k] * (z_tim[j] - z_tim[j-1]))
    }
    return(list(z_tim, z_ct))
  })
  
  Z_tim <- cova_tim[1, ]
  Z_ct <- cova_tim[2, ]
  
  CL_ct <- sapply(Z_ct, function(L)L[[1]])
  CR_ct <- sapply(Z_ct, function(L)L[[length(L)]])

  Delta <- rep(1, N+n) 
  Delta[which(CL > Tim)] <- 3
  Delta[which(CR < Tim)] <- 2
  
  Deltastar <- rep(1, N+n) 
  Deltastar[which(CL>Tstar)] <- 3
  Deltastar[which(CR<Tstar)] <- 2
  
  Label_id <- rep(NA, N+n)
  ii <- sample.int(N+n, n)
  Label_id[ii] <- 1
  
dataset <- list('X' = X, 'Delta' = Delta, 'Ltime' = CL, 'Rtime' = CR, 'Xstar' = Xstar, 
                'Deltastar' = Deltastar, 'Label_id' = Label_id, 'Lcount' = CL_ct,
                'Rcount' = CR_ct, 'Z_time' = Z_tim, 'Z_count' = Z_ct)
  return(dataset)
}

