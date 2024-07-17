
library(matrixStats); library(stepPlr); library(evd); library(methods); library(pracma);
library(MASS); library(survival); library(plyr); library(splines); library(glmnet); 
library(matrixcalc);

source('HelperFunctions.R');
source('Beta_estimate.R');
source('NewtonCLseh.R');
source('Supervised_est.R');
source('CV_function.R')
source('IntrSSL_est_ind.R')
source('GenerateCoxbtcor67.R')
source('GenerateCoxbtcor8.R')
source('Generatelogbtcor8.R')
source('d011.R')


setting = 1; seed = 1; 

std0_all = sdd0_all = NULL 
D_sl_all = D_sl_sd_all = NULL 
D_intr_all = D_intr_sd_all = NULL
L_sl_all = L_sl_sd_all = NULL
L_intr_all = L_intr_sd_all = NULL
R_sl_all = R_sl_sd_all = NULL
R_intr_all = R_intr_sd_all = NULL
C_sl_all = C_sl_sd_all = NULL 
C_intr_all = C_intr_sd_all = NULL

for(i in (seed*1 + 1):(seed*1 + 500)){
  set.seed(i);
  print(paste('iter',i));

  ####cox model
  if(setting == 1){ 
    print(setting)
    data <- GenerateCoxbtcor67(N = 5000, n = 250, phi1 = -7.6, phi2 = -0.15, M = 3.3, a = 0.6, rate = 0.1);
    tim <- seq(0.9, 3.02, length=50); 
    adl = 0.13; adh = 0.2;
  }
  
  if(setting == 2){ 
    print(setting)
    data <- GenerateCoxbtcor8(N = 5000, n = 250, phi01 = -12.5, phi02 = -0.2, M = 3, a = 0.3, rate = 0.1);
    tim <- seq(1, 3.08, length=50); 
    adl = 0.1; adh = 0.15;
  }
  
  ####log model
  if(setting == 3){ 
    print(setting)
    data <- Generatelogbtcor8(N = 5000, n = 250, rate = 0.1);
    tim <- seq(1.3, 3.356, length=50); 
    adl = 0.1; adh = 0.15;
  }
  
  #########################################################################
  ###for labeled data id
  label_id <- which(data$Label_id == 1); 
  #### for double censor ######
  delt <- data$Delta[label_id];
  obse <- data$X[label_id];
  
  rcen_labeled <- data$Rtime[label_id];
  lcen_labeled <- data$Ltime[label_id];
  n_l <- length(rcen_labeled)
  lh_labeled <- 1.06 * sd(lcen_labeled) * n_l^{-(0.2 + adl)};
  rh_labeled <- 1.06 * sd(rcen_labeled) * n_l^{-(0.2 + adh)};
  
  rcen_unlabeled <- data$Rtime[-label_id];
  lcen_unlabeled <- data$Ltime[-label_id];
  N_unl <- length(rcen_unlabeled)
  lh_unlabeled <- 1.06 * sd(lcen_unlabeled) * N_unl^{-(0.2 + adl)};
  rh_unlabeled <- 1.06 * sd(rcen_unlabeled) * N_unl^{-(0.2 + adh)};
  
  rcen_all <- data$Rtime;
  lcen_all <- data$Ltime;
  
  ##########################################################################
  #### for surrgate and covariates ###### 
  xstar_all <- data$Xstar;
  deltastar_all <- data$Deltastar;
  #### transfer the catergorival to indicator
  deltastar_all <- data$Deltastar;
  deltastar_ind1 <- rep(); deltastar_ind2 <- rep();
  deltastar_ind1[which(deltastar_all==1)] <- 0;
  deltastar_ind2[which(deltastar_all==1)] <- 0;
  deltastar_ind1[which(deltastar_all==2)] <- 1;
  deltastar_ind2[which(deltastar_all==2)] <- 0;
  deltastar_ind1[which(deltastar_all==3)] <- 0;
  deltastar_ind2[which(deltastar_all==3)] <- 1;
  deltastar_indv <- cbind(deltastar_ind1, deltastar_ind2) 
  #### time-indepnedent covariate
  base_cov <- data$Z
  #### time-depnedent covariate
  lcen_ct <- data$Lcount
  rcen_ct <- data$Rcount
  cova_tim <- data$Z_time
  cova_ct <- data$Z_count
  
  ##########################################################################
  ###### for right censor and left censor #########
  ldelt <- rep(0, n_l);
  ldelt[which(delt == 3)] = 1;
  
  rdelt <- rep(0, n_l);
  rdelt[which(delt==1|delt==3)] = 1;
  
  ############################################################################################
  ################# For survival function estimate ##########
  ###### for self-consistent method #####
  d011_d <- rep(1, length(delt))
  d011_d[which(delt==2)] <- 0
  d011_d[which(delt==3)] <- 2
  
  D011<-d011_new(z=obse, d=d011_d, influence.fun = T)
  stimd011 <- D011$exttime
  surtd011 <- D011$extsurv.Sx
  vtimd011 <- D011$Nodes
  sdd011 <- sqrt(D011$VarFt)
  
  std0 <- sapply(tim,function(x){
    index<-order(abs(stimd011-x))[c(1,2,3,4)]
    mst<-mean(surtd011[index])
    return(mst)
  })
  
  sdd0 <- sapply(tim,function(x){
    index<-order(abs(vtimd011-x))[c(1,2)]
    msd<-mean(sdd011[index])
    return(msd)
  })  
 
  std0_all = rbind(std0_all, std0);
  sdd0_all = rbind(sdd0_all, sdd0);
  
######################################################################################################
  ############################################################################################
  ################# For survival function estimate ##########
  ################## for supervised estimate #########    
  SupSt <- Supervised_est(time = tim, obse, ldelt, rdelt, lcen = lcen_labeled, rcen = rcen_labeled, 
                          lh = lh_labeled, rh = rh_labeled);
  DSt_sl <- SupSt[1, ];
  DSt_sl_sd <- SupSt[2, ];
  LSt_sl <- SupSt[3, ];
  LSt_sl_sd <- SupSt[4, ];
  RSt_sl <- SupSt[5, ];
  RSt_sl_sd <- SupSt[6, ];
  CSt_sl <- SupSt[7, ];
  CSt_sl_sd <- SupSt[8, ];
  D_sl_all = rbind(D_sl_all, DSt_sl);
  D_sl_sd_all = rbind(D_sl_sd_all, DSt_sl_sd);
  L_sl_all = rbind(L_sl_all, LSt_sl);
  L_sl_sd_all = rbind(L_sl_sd_all, LSt_sl_sd);
  R_sl_all = rbind(R_sl_all, RSt_sl);
  R_sl_sd_all = rbind(R_sl_sd_all, RSt_sl_sd);
  C_sl_all = rbind(C_sl_all, CSt_sl);
  C_sl_sd_all = rbind(C_sl_sd_all, CSt_sl_sd);
  ####### for Intrinsic and semi-supervised estimate #############
  Intr_SSL <- IntrSSL_est(time = tim, base_cov, cova_tim, cova_ct, lcen_ct, rcen_ct, xstar_all, 
                          deltastar_indv, label_id, lcen_all, rcen_all, obse, ldelt, rdelt, 
                          lh_labeled, lh_unlabeled, rh_labeled, rh_unlabeled, num_folds = 10)
  D_intr <- Intr_SSL[1, ];
  L_intr <- Intr_SSL[2, ];
  R_intr <- Intr_SSL[3, ];
  C_intr <- Intr_SSL[4, ];
  D_intr_sd <- Intr_SSL[5, ];
  L_intr_sd <- Intr_SSL[6, ];
  R_intr_sd <- Intr_SSL[7, ];
  C_intr_sd <- Intr_SSL[8, ];
  D_intr_all <- rbind(D_intr_all, D_intr);
  D_intr_sd_all <- rbind(D_intr_sd_all, D_intr_sd);
  L_intr_all <- rbind(L_intr_all, L_intr);
  L_intr_sd_all <- rbind(L_intr_sd_all, L_intr_sd);
  R_intr_all <- rbind(R_intr_all, R_intr);
  R_intr_sd_all <- rbind(R_intr_sd_all, R_intr_sd);
  C_intr_all <- rbind(C_intr_all, C_intr);
  C_intr_sd_all <- rbind(C_intr_sd_all, C_intr_sd);
  ######################################################################################################
  
  save(std0_all, sdd0_all, D_sl_all, D_sl_sd_all, D_intr_all, D_intr_sd_all, 
       L_sl_all, L_sl_sd_all, L_intr_all, L_intr_sd_all,
       R_sl_all, R_sl_sd_all, R_intr_all, R_intr_sd_all,
       C_sl_all, C_sl_sd_all, C_intr_all, C_intr_sd_all,
       file=paste0('setting',setting,'seed',seed,'.RData'))
}
