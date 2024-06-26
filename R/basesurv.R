basesurv <- function(Time, Status, X, beta, Lambda, w) {
   kappaa <- uniroot(function(t) { sum(w * log(Time) * (Status - Time^t* exp(beta %*% t(X))) + Status/t  ) }, c(0.1, 3))$root
  
   # kappaa <- -sum(Status)/sum(w * log(Time) * (Status - Lambda * exp(beta %*% t(X))) )
  bcumhaz <- Time^(kappaa)
  bsurv <- exp(-bcumhaz)
  uncuresurv <- exp(-bcumhaz * exp(beta %*% t(X)))
  list(uncuresurv = uncuresurv, bcumhaz = bcumhaz, bsurv = bsurv, kappa = kappaa)
}