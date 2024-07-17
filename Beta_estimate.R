Beta_estimate <- function(dat_label, y, h, w_label, w_unlabel){
  
  p_basis <- ncol(dat_label) # the dimension the data set
  
  # Step 1: calculate the semi-supervised estimator
  # Basis function regression for imputation
  beta_imp <- tryCatch(glm(y ~ dat_label, family = 'binomial', weights = w_label)$coeff,
                  error = function(e) rep(NA, p_basis + 1))
  beta_imp[which(is.na(beta_imp))] <- 0
  beta_imp <- unname(beta_imp)
  
  ###Step 2: calculate the intrinsic estimator
  beta_intr <- NewtonCLseh(X = dat_label, y, h, w_label, w_unlabel, max_iter = 100, 
                         tol = 1e-4, uper = 1e+15, initial = beta_imp)

  return(c(beta_ssl = beta_imp, beta_intr = beta_intr))
}
