Generatelogweiminic <- function(N, n){

  Tstar <- rweibull(N+n, 3, 4)
  Tim <- rlogis(N+n, 1 + 1 * Tstar, 0.2)
  CL <- rweibull(N+n, 2.3, 1.92)
  CR <- runif(N+n, 0, 2.1) + CL
  X <- pmax(pmin(Tim, CR), CL)
  Xstar <- pmax(pmin(Tstar, CR), CL)
  
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
                'Deltastar' = Deltastar, 'Label_id' = Label_id)
  return(dataset)
}

