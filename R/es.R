es <- function(Time, Status, X, Z, id,  corstr, stdz, esmax, eps,invfun) {
    Kn <- length(id)
    K <- length(unique(id))
    n <- as.vector(table(id))
    Z1 <- Z
    X1 <- X
    Map1 <- list()
    if (stdz) {
        for (i in 2:ncol(Z1)) {
            Z[, i] <- (Z[, i] - mean(Z[, i]))/sd(Z[, i])
        }
     
      for (i in 2:ncol(X1)) {
        X[, i] <- (X[, i] - mean(X[, i]))/sd(X[, i])
      }
        
    }
    
    w <- Status
    Lambda <- initial_Lambda(Time, Status, X, Z, id, model, corstr = "independence",invfun,eps)$Lambda
    KK1 <- -1
    gamma2<-rep(0,dim(Z)[2])
    beta2 <- rep(0, dim(X)[2])
    alpha2 <- 0
    ppmt2 <- c(gamma2, beta2,alpha2)
    repeat {
    #  Gammadata <- data.frame(w=w,zzz=Z[,-1])
    #    Gamma1 <- eval(parse(text = paste("qif", "(", "w~zzz",",data=Gammadata", ",id=", "id", ",family = binomial", ",corstr='", corstr, "'", ",invfun=","invfun",",tol=","eps", ")", sep = "") ))  
    #    gamma1 <- as.double( (Gamma1$parameter)[,1] ) 
       
      #independence exchangeable AR-1
       Gamma1 <- geega(w,Z,gamma=gamma2,id,corstr)
        gamma1 <- as.double(Gamma1$gamma)

        SK2 <- 1
        Lambda <- Lambda
 
        repeat {         
       
            Beta1 <- geebt(Status, Lambda, X, beta = beta2, w, id, corstr,invfun,eps)
      
            beta1 <- as.double(Beta1$beta) 
            Basesurv <- basesurv(Time, Status, X, beta = beta1, Lambda, w)
            gSS3 <- Basesurv$bcumhaz
            alpha1 <- Basesurv$kappa
            if ((any(abs(c(beta1, alpha1) - c(beta2,alpha2)) > eps)) & (SK2 <= esmax)) {
                beta2 <- beta1
                alpha2 <- alpha1
                Lambda <- gSS3
                SK2 <- SK2 + 1
            } else {
                ppmt1 <- c(gamma1, beta1,alpha1)
                survival <- Basesurv$uncuresurv
                alpha1 <- Basesurv$kappa
                w <- Status + ((1 - Status) * exp(gamma1 %*% t(Z)) * survival)/(1 + exp(gamma1 %*% t(Z)) * survival)
                w <- as.vector(w)
                break
            }
        }
        if (any(abs(ppmt1 - ppmt2) > eps) && (KK1 < esmax)) {
            ppmt2 <- ppmt1
            gamma2 <- gamma1
            w1 <- w
            KK1 <- KK1 + 1
            if(KK1 > 0){
            Map1 [[KK1]] <- c(gamma1, beta1, alpha1)
            }
        } else {
          KK1 <- KK1 + 1
          if(KK1 > 0){
            Map1 [[KK1]] <- c(gamma1, beta1, alpha1)
          }
          break
        }       
    }
    gSS3 <- Basesurv$bcumhaz
    basesurv <- Basesurv$bsurv
    if (stdz == "TRUE") {
        gamma1 <- c(gamma1[1] - sum((gamma1[-1] * apply(Z1[, -1, drop = FALSE], 2, mean)/apply(Z1[, -1, drop = FALSE], 2, sd))), gamma1[-1]/apply(Z1[, -1, drop = FALSE], 2, sd))
        beta1 <- c(beta1[1] - sum((beta1[-1] * apply(X1[, -1, drop = FALSE], 2, mean)/apply(X1[, -1, drop = FALSE], 2, sd))),beta1[-1]/apply(X1[, -1, drop = FALSE], 2, sd))
    }
  
    ResultPara <- c(gamma1, beta1, alpha1)
    lRP <- length(ResultPara)

    arrayRM <- array(0,c(lRP,lRP,KK1))
    resP11 <- rep(0,lRP)
    for(ii in 1:KK1){
      resP <- Map1 [[ii]]
      for(ii1 in 1:lRP){
        
         if(ii1 == 1 ){
           resP11[ii1] <- resP[ii1]
           resP11[(ii1+1):lRP] <- ResultPara[(ii1+1):lRP]
           arrayRM[ii1,,ii] <- resP11
         }else if(ii1 == lRP){
           resP11[1:(ii1-1)] <- ResultPara[1:(ii1-1)]
           resP11[ii1] <- resP[ii1]
           arrayRM[ii1,,ii] <- resP11
         }else{
           resP11[1:(ii1-1)] <- ResultPara[1:(ii1-1)]
           resP11[ii1] <- resP[ii1]
           resP11[(ii1+1):lRP] <- ResultPara[(ii1+1):lRP]
           arrayRM[ii1,,ii] <- resP11
         }
      }
    }
  
    uncureprob <- as.vector(exp(Z%*%gamma1)/(1 + exp(Z%*%gamma1)))
    convergence <- sum((ppmt1 - ppmt2)^2)
    gamma_Q2 <- Gamma1$Q2
    beta_Q2 <- Beta1$Q2
    gamma_seQ2 <- Gamma1$seQ2
    beta_seQ2 <- Beta1$seQ2
   list(gamma = gamma1, beta = beta1, kappaa = alpha1, Lambda = gSS3, 
      w = w1, survival = survival,basesurv=basesurv, Uncureprob = uncureprob, tau = convergence,gS=alpha1,uncuresurv=survival, gamma_Q2 =  gamma_Q2,beta_Q2 = beta_Q2,gamma_seQ2 = gamma_seQ2 ,beta_seQ2 =beta_seQ2 ,arrayRM=arrayRM)
}