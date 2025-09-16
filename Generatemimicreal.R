Generatelognormminic <- function(N, n){
 
  # N=113623; n=1613
  # set.seed(10000)
  
  Label_id <- rep(NA, N+n)
  ii <- sample.int(N+n, n)
  Label_id[ii] <- 1
  
  #Tstar <- rnorm(N+n, 72.5, 17.2)
  Tstar <- rnorm(N+n, 75, 15)
  #Tstar <- rnorm(N+n, 75, 1)
  #mean(Tstar); sd(Tstar); range(Tstar); hist(Tstar)
  
  Tim <- rlogis(N+n,  1.25 * Tstar + 80 - mean(1.25 * Tstar), 5)
  #mean(Tim); sd(Tim); range(Tim); cor(Tim, Tstar)
  
  CL <- rnorm(N+n, 33, 14.5)
  CL <- pmax(0, pmin(99, CL)) 
  # real data: 42.8; 16.7; c(0, 99); 
  #mean(CL); sd(CL);range(CL); hist(CL, freq = F); 
  
  #CR <- runif(N+n, 0, 30) + CL
  CR <- rnorm(N+n, 0, 16) + 15 + CL
  #CR <- rnorm(N+n, 15, 16) + CL
  CR <- pmax(0, pmin(105, CR))
  # real data: 55.4; 17.7; c(0, 105); 
  # mean(CR); sd(CR); range(CR); cor(CL, CR); hist(CR, freq = F)
  
  X <- pmax(pmin(Tim, CR), CL)
  Delta <- rep(1, N+n) 
  Delta[which(CL > Tim)] <- 3
  Delta[which(CR < Tim)] <- 2
  # real data: 10.4%; 86.86%; 2.73%
  #mean(Delta[which(Label_id==1)] ==1 ); mean(Delta[which(Label_id==1)] == 2); mean(Delta[which(Label_id==1)] == 3)
  # real data: 53.74; 17.19; (6, 96.75)
  # mean(X[which(Label_id==1)]); sd(X[which(Label_id==1)]); range(X[which(Label_id==1)]); hist(X[which(Label_id==1)],freq = F)
  

  Xstar <- pmax(pmin(Tstar, CR), CL)
  Deltastar <- rep(1, N+n) 
  Deltastar[which(CL>Tstar)] <- 3
  Deltastar[which(CR<Tstar)] <- 2
  # real data: 14.7%; 83%; 2.3%
  #mean(Deltastar == 1);mean(Deltastar == 2);mean(Deltastar == 3);
  # 53.89; 17.23; c(0, 105)
  #mean(Xstar); sd(Xstar); range(Xstar); hist(Xstar, freq = F)
  
  
  # par(mfrow=c(2,2))
  # hist(CL, freq = F); hist(CR, freq = F)
  # hist(X[which(Label_id==1)],freq = F)
  # hist(Xstar, freq = F)
  # 
  # par(mfrow=c(2,2))
  # hist(CL[which(Label_id==1)], freq = F);
  # hist(CR[which(Label_id==1)], freq = F)
  # hist(X[which(Label_id==1)],freq = F)
  # hist(Xstar[which(Label_id==1)], freq = F)
  
  
dataset <- list('X' = X, 'Delta' = Delta, 'Ltime' = CL, 'Rtime' = CR, 'Xstar' = Xstar, 
                'Deltastar' = Deltastar, 'Label_id' = Label_id)
  return(dataset)
}


# set.seed(100)
# data = Generatelognormminic(N=113623, n=1613)
# quantile(data$X[which(data$Label_id==1)], probs = c(0.1, 0.9))

# True_fun <- function(time, N){
#   set.seed(1000)
#   Tstar <- rnorm(N, 75, 15)
#   #Tim <- rlogis(N, 7 * Tstar + 80 - mean(7 * Tstar), 5)
#   Tim <- rlogis(N, 1.25 * Tstar + 80 - mean(1.25 * Tstar), 5)
#   #mean(Tim); sd(Tim); cor(Tim, Tstar)
#   # hist(Tstar); 
#   # hist(Tim, freq = F)
#   # lines(density(Tim))
#   St <- sapply(time, function(t) mean(Tim > t))
#   return(St)
# }
# 
# tim = seq(30, 70, length=100) 
# X_quant = quantile(X[which(Label_id==1)], prob = c(0.1, 0.9))
# tim = seq(X_quant[1], X_quant[2], length=100)
# St = True_fun(time = tim, N=1000000)
# plot(tim, St, type = "l", ylim = c(0.5, 1))


# par(mfrow=c(2,2))
# hist(CL, freq = F)
# lines(density(CL))
# 
# hist(CR, freq = F)
# lines(density(CR))
# 
# hist(X, freq = F, na.rm = T)
# lines(density(X, na.rm=T))
# 
# hist(Xstar, freq = F)
# lines(density(Xstar))
# cor(X, Xstar)
# # 0.9423
# 
# 
# par(mfrow=c(2,2))
# hist(CL[Label_id==1], freq = F)
# lines(density(CL[which(Label_id==1)]))
# 
# hist(CR[Label_id==1], freq = F)
# lines(density(CR[which(Label_id==1)]))
# 
# hist(X[which(Label_id==1)], freq = F)
# lines(density(X[which(Label_id==1)]))
# 
# hist(Xstar[Label_id==1], freq = F)
# lines(density(Xstar[which(Label_id==1)]))
# cor(X[which(Label_id==1)], Xstar[which(Label_id==1)])
# # 0.9411

