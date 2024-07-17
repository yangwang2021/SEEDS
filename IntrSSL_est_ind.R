
IntrSSL_est <- function(time, base_cov, cova_tim, cova_ct, lcen_ct, rcen_ct, 
                        xstar_all, deltastar_indv, label_id, lcen_all, rcen_all, 
                        obse, ldelt, rdelt, lh_labeled, lh_unlabeled, rh_labeled, 
                        rh_unlabeled, num_folds){

  re <- sapply(time, function(t){
    # obtain the indicate of surrogate
    Ind_star <- as.numeric(I(xstar_all <= t))
    
    # construct the covariate matrix for the exact observed information
    if(is.null(cova_ct)){ dcova_all <- NULL }else{
      dcova_all <- TimeCovar(time_point = t, List_time = cova_tim, List_count = cova_ct, 
                             leng = length(xstar_all))}
    Dbasis_all <- cbind(deltastar_indv, Ind_star, base_cov, dcova_all);
    Dbasis_labeled <- Dbasis_all[label_id, ];
    Dbasis_unlabeled <- Dbasis_all[-label_id, ];
    
    # construct the covariate matrix for the left censored data
    if(is.null(lcen_ct)){ lcova_all <- NULL }else{ lcova_all <- lcen_ct }
    Lbasis_all <- cbind(deltastar_indv, Ind_star, base_cov, lcova_all);
    Lbasis_labeled <- Lbasis_all[label_id, ];
    Lbasis_unlabeled <- Lbasis_all[-label_id, ];
    
    # construct the covariate matrix for the right censored data
    if(is.null(rcen_ct)){ rcova_all <- NULL }else{ rcova_all <- rcen_ct }
    Rbasis_all <- cbind(deltastar_indv, Ind_star, base_cov, rcova_all);
    Rbasis_labeled <- Rbasis_all[label_id, ];
    Rbasis_unlabeled <- Rbasis_all[-label_id, ];
    
    # obtain the response variable and weight of the three different types of data
    lcen_labeled <- lcen_all[label_id];
    rcen_labeled <- rcen_all[label_id];
    
    lcen_unlabeled <- lcen_all[-label_id];
    rcen_unlabeled <- rcen_all[-label_id];
    
    Inter_labeled <- as.numeric(I(t <= rcen_labeled & t > lcen_labeled))
    Inter_unlabeled <- as.numeric(I(t <= rcen_unlabeled & t > lcen_unlabeled))
    dy_obse <- as.numeric(I(t <= obse & t > lcen_labeled))
    
    lw_labeled <- Normalize_fun(data = lcen_labeled, pt = t, bd = lh_labeled)
    lw_unlabeled <- Normalize_fun(data = lcen_unlabeled, pt = t, bd = lh_unlabeled)

    rw_labeled <- Normalize_fun(data = rcen_labeled, pt = t, bd = rh_labeled)
    rw_unlabeled <- Normalize_fun(data = rcen_unlabeled, pt = t, bd = rh_unlabeled)
    # the intrinsic estimation
    Intr_est <- CV_function(ddat_unlabel = Dbasis_unlabeled, ddat_label = Dbasis_labeled, 
                            ldat_unlabel = Lbasis_unlabeled, ldat_label = Lbasis_labeled, 
                            rdat_unlabel = Rbasis_unlabeled, rdat_label = Rbasis_labeled,
                            dw_label = Inter_labeled, dw_unlabel = Inter_unlabeled, 
                            dy = dy_obse, dh = 1, lw_label = lw_labeled, lw_unlabel = lw_unlabeled, 
                            ly = 1 - ldelt, lh = lh_labeled, rw_label = rw_labeled, rw_unlabel = rw_unlabeled, 
                            ry = 1 - rdelt, rh = rh_labeled, num_folds, a = 1)
    DSt_intr <- Intr_est$DSt
    LSt_intr <- Intr_est$LSt
    RSt_intr <- Intr_est$RSt
    CSt_intr <- Intr_est$CSt
    Dsd_intr <- Intr_est$D_sd
    Lsd_intr <- Intr_est$L_sd
    Rsd_intr <- Intr_est$R_sd
    Csd_intr <- Intr_est$C_sd
    
    return(c( 'DSt_intr' = DSt_intr,  'LSt_intr' = LSt_intr, 
              'RSt_intr' = RSt_intr,  'CSt_intr' = CSt_intr,  
              'Dsd_intr' = Dsd_intr,  'Lsd_intr' = Lsd_intr, 
              'Rsd_intr' = Rsd_intr,  'Csd_intr' = Csd_intr))
    })
  return(re)
}
