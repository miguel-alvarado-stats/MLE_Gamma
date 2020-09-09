# MAXIMUM LIKELIHOOD ESTIMATION: FISHER-SCORING
# GAMMA DISTRIBUTION

MLE_FS_Gamma <- function(y, shape = 2, maxiter = 100, epsilon = 0.000000001, stop_criteria = 10000){

  n <- length(y)
  ysum <- sum(y)

  U <- function(n, ysum, scale){
    -((n * shape)/scale) + (1/scale^2) * ysum
  }

  J <- function(n, scale){
    (n*shape)/(scale^2)
  }

  # first guess:
  scale0 <- min(y)
  Estimator <- matrix(NA, ncol = 1)
  Estimator[1, 1] <- scale0

  iter <- 1
  while( (stop_criteria > epsilon) & (iter <= maxiter) ){

    num <- U(n, ysum, Estimator[iter, 1])
    den <- J(n, Estimator[iter, 1])

    UPD <- (num/den)
    Estimator_iter <- as.matrix(Estimator[iter, 1]) + UPD
    Estimator <- rbind(Estimator, Estimator_iter)
    stop_criteria <- UPD^2
    iter <- iter + 1
  }

  Likelihood <- matrix(NA, ncol = 1)
  LogLikelihood <- matrix(NA, ncol = 1)

  for(i in 1:length(Estimator)){
    Likelihood[i] <- prod( (1/(gamma(shape)*(Estimator[i,1]^shape))) * (y^(shape - 1)) * exp(-(y/Estimator[i,1])) )
  }

  for(i in 1:length(Estimator)){
    LogLikelihood[i] <- sum( -log(gamma(shape)) - (shape*log(Estimator[i,1])) + ((shape - 1)*log(y)) - (y/Estimator[i,1]) )
  }

  results <- cbind(format(Estimator, nsmall = 6), Likelihood, LogLikelihood)
  colnames(results) <- c("ML Estimator", "Likelihood", "Log-Likelihood")
  return(results)

}

