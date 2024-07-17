# the expit function
Expit <- function(x) {
  return(1 / (1 + exp(-x)))
}

# the deviation of the expit function
ExpitDerivative <- function(x, na_correction = T) {
  expit_deriv <- exp(x)/(1+exp(x))^2
  if(na_correction){
    expit_deriv[which(is.na(expit_deriv))] = 0
  }
  return(expit_deriv)
}

# the normal kernel function
ker <- function(x,h) {
  re <- exp(-(x/h)^2/2)/(h*sqrt(2*pi))
}


# the time-dependent count number
TimeCovar <- function(time_point, List_time, List_count, leng){
  re<-sapply(1:leng, function(i){
    lt <- length(which(List_time[[i]] <= time_point))
    if(lt == 0){lc <- 0}else{lc <- List_count[[i]][lt]}
    return(lc)
  })
  return(re)
}

# the cross-validation split function
cv_split <- function(n, num_folds){
  re <- split(1:n, sample(rep(1:num_folds, floor(n / num_folds))))
  return(re)
} 

# the normalization function for weight
Normalize_fun <- function(data, pt, bd){
  weigt <- ker(data - pt, bd)
  if(mean(weigt) == 0){weigt <- rep(0, length(weigt))}else{
    weigt <- weigt / mean(weigt)}
  return(weigt)
}





