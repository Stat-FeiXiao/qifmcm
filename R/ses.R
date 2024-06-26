ses <- function(Time,Status,id,Z,X,corstr,invfun,ResultPara,arrayRM,eps){
  ir <- 0       
  d <- length(ResultPara)
  M <- M1 <- matrix(0,d,d)
  eps1 <- 0.01
  repeat{
    ir <- ir + 1
    RM <- arrayRM[,,ir]
    for(i in 1:d){
      for(j in 1:d){
        thetai <- RM[i,]
        Map <- ESMap(Time,Status,id,Z,X,thetai,corstr,invfun,eps1)$Map
        M[i,j] <-  (Map[j] - ResultPara[j])/ (thetai[i] - ResultPara[i])             
        
      }
    }

    if ( max(abs(M - M1)) < eps1 ) break
    M1 <- M
    
  }
  
  list(ir=ir,M=M) #M是ES收敛率矩阵J的估计
}

ESMap <- function(Time,Status,id,Z,X,thetai,corstr,invfun,eps){
  Lgamma <- dim(Z)[2]
  Lbeta <- dim(X)[2]
  gamma1 <- thetai[1:Lgamma]
  beta1 <- thetai[(Lgamma+1):(Lgamma+Lbeta)]
  kappaa <- thetai[(Lgamma+Lbeta+1) : length(thetai)]
  #
  Lambda <- Time^(kappaa)
  survival <- exp(-Lambda * exp(beta1 %*% t(X)))
 
  #
  ####ES -- one iteration
  ##E
  w <- Status + ((1 - Status) * exp(gamma1 %*% t(Z)) * survival)/(1 + exp(gamma1 %*% t(Z)) * survival)
  w <- as.vector(w)
  ##S
#  Gammadata <- data.frame(w=w,zzz=Z[,-1])
#  Gamma1 <- eval(parse(text = paste("qif", "(", "w~zzz",",data=Gammadata", ",id=", "id", ",family = binomial", ",corstr='", 
#                                    corstr, "'", ",invfun=","invfun",",tol=","eps", ")", sep = "") ))  
#  gamma1 <- as.double( (Gamma1$parameter)[,1] )
  gamma2 <-   gamma1        #        rep(0, dim(Z)[2])
  Gamma1 <- geega(w,Z,gamma=gamma2,id,corstr)
  gamma1 <- Gamma1$gamma
  beta2 <-    beta1        #       rep(0, dim(X)[2])
  alpha2 <- kappaa
  SK2 <- 0
 repeat{ 
  Beta1 <- geebt(Status, Lambda, X, beta = beta2, w, id, corstr,invfun,eps)
  beta1 <- as.double(Beta1$beta)
  Basesurv <- basesurv(Time, Status, X, beta = beta1, Lambda, w)
  gSS3 <- Basesurv$bcumhaz
  alpha1 <- Basesurv$kappa

  if ((any(abs(c(beta1, alpha1) - c(beta2,alpha2)) > eps)) & (SK2 <= 200)) {
    beta2 <- beta1
    alpha2 <- alpha1
    Lambda <- gSS3
    SK2 <- SK2 + 1
  } else break
}
  ####
  Map <- c(gamma1,beta1,alpha1)
  list(Map=Map)
}