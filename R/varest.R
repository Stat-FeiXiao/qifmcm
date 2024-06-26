varest <- function(Time, Status, X, Z, id, gamma, beta, Lambda, w,kappaa,gSfordev,uncuresurv,baseuncuresurv,gamma_Q2 ,beta_Q2 ,gamma_seQ2,beta_seQ2,invfun,corstr,arrayRM,eps) { 
  K <- length(unique(id))
  n <- as.vector(table(id))
  sdm <- matrix(0, dim(Z)[2] + dim(X)[2] + 1, dim(Z)[2] + dim(X)[2] + 1)
  Time <- as.double(Time)

    t2 <- as.double(Time)
    c1 <- Status
    t11 <- sort(Time)
    c11 <- Status[order(Time)]
    tt1 <- unique(t11[c11 == 1])
    kk <- length(table(t11[c11 == 1]))
    newppmt2c <- gamma
    newppmt2s <- beta
    ResultPara <- c(gamma,beta,kappaa)
    gg1 <- w
    z1 <- Z
    p1 <- exp(z1 %*% newppmt2c)/(1 + exp(z1 %*% newppmt2c))
    g1 <- gg1
    sdm_gamma <- gamma_Q2 
    be <- as.vector(newppmt2s)
    z2 <- X
    c2 <- Status
    mu2 <- exp(z2 %*% be)
    Lambda0 <- Lambda
    VA1 <- beta_Q2 
    F1 <- (t2^kappaa) * mu2 * log(t2)
    CCC <- -sum((g1 * ((log(t2))^2) * (t2^kappaa) * mu2 + c2/(kappaa^2)))  
    BBC <- t(g1 * F1) %*% z2 
    
    BC1=matrix(0, dim(z2)[2], 1)
    BC2=matrix(0, dim(z2)[2], 1)
    
    
    OABC1 <- rep(0, K)
    for (v in 1:dim(z2)[2]) {
      for (i in 1:K) {
        M1<-diag(n[i])
        z22 <- matrix(z2[id == i, ], nrow = n[i], )
        A2 <- t(z22[, v])
        g22 <- g1[id == i]
        mu22 <- mu2[id == i]
        D2 <- (mu22^(1/2)) %*% ((t(mu22))^(-1/2)) * M1
        G2 <- diag(g22, n[i], n[i])
        F2 <- F1[id == i]
        OABC1[i] <- A2 %*% D2 %*% G2 %*% F2
      }
      BC1[v, 1] <- -sum(OABC1) 
    }
    
    OABC2 <- rep(0, K)
    for (v in 1:dim(z2)[2]) {
      for (i in 1:K) {
        M2<-matrix(1,n[i],n[i])
        diag(M2)<-0
        z22 <- matrix(z2[id == i, ], nrow = n[i], )
        A2 <- t(z22[, v])
        g22 <- g1[id == i]
        mu22 <- mu2[id == i]
        D2 <- (mu22^(1/2)) %*% ((t(mu22))^(-1/2)) * M2
        G2 <- diag(g22, n[i], n[i])
        F2 <- F1[id == i]
        OABC1[i] <- A2 %*% D2 %*% G2 %*% F2
      }
      BC2[v, 1] <- -sum(OABC1) 
    }
    
    BC=beta_seQ2%*%rbind((BC1),(BC2))
   
    sdm_betalpha <- rbind(cbind(VA1, BC), cbind(-BBC, CCC ))
    sdm <- bdiag(sdm_gamma, sdm_betalpha)
    sdm <- as.matrix(sdm)
    
  
    ses <- ses(Time,Status,id,Z,X,corstr,invfun,ResultPara,arrayRM,eps)

    estJ <- ses$Map  
    lsdm <- -sdm %*% ( diag(dim(z1)[2] + dim(z2)[2] + 1) - estJ )
    fdv <- rep(0, dim(z1)[2] + dim(z2)[2] + 1)
    fdm <- matrix(0, dim(z1)[2] + dim(z2)[2] + 1, dim(z1)[2] + dim(z2)[2] + 1)
    for (i in 1:K) {
      z11 <- z1[id == i, ]
      p11 <- p1[id == i]
      g11 <- g1[id == i]
      t22 <- t2[id == i]
      pp11 <- p11 * (1 - p11)
      pp11m <- diag(pp11)
      C1 <- g11 - p11
      R0 <- diag(n[i])
      
      
      fdv10<- t(pp11m %*% z11) %*% ginv(sqrt(pp11m) )%*% R0 %*% ginv(sqrt(pp11m) ) %*% C1
      R1 <- matrix(1, n[i], n[i])
      diag(R1) <- 0
      fdv11<- t(pp11m %*% z11) %*% ginv(sqrt(pp11m) )%*% R1 %*% ginv(sqrt(pp11m) ) %*% C1
      
      
      
      fdv[1:dim(z1)[2]]<- gamma_seQ2%*%rbind(fdv10,fdv11)
      
      ##############
      
      
      z22 <- z2[id == i, ]
      mu22 <- mu2[id == i]
      mu22m <- diag(mu22, n[i], n[i])
      Lambda01 <- Lambda0[id == i]
      G2 <- diag(g11 * (t22^kappaa), n[i], n[i])
      c22 <- c2[id == i]
      C2 <- (c22/(t22^kappaa)) - mu22
  
      R0 <- diag(n[i])
      
      fdv20<- t(mu22m %*% z22) %*% ginv(sqrt(mu22m)) %*% R0 %*% ginv(sqrt(mu22m) ) %*% G2 %*% C2
      R1 <- matrix(1, n[i], n[i])
      diag(R1) <- 0
      fdv21<- t(mu22m %*% z22) %*% ginv(sqrt(mu22m)) %*% R1 %*% ginv(sqrt(mu22m) ) %*% G2 %*% C2
      
      
      
      fdv[(dim(z1)[2]+1):(dim(z1)[2]+dim(z2)[2])]<- beta_seQ2%*%rbind(fdv20,fdv21)
      
      ######################
      
     
      fdv[(dim(z1)[2] + dim(z2)[2] + 1)]<- sum(g11 * (log(t22)) * (c22 - (t22^kappaa) * mu22) + c22/kappaa)
      fdm1 <- fdv %*% t(fdv)
      fdm <- fdm + fdm1
    }
    if (invfun=="finv") {
      a11 <- lsdm
      storage.mode(a11)<-"double"
      lda<-as.integer(nrow(a11))
      n11<-as.integer(ncol(a11))
      z <- .Fortran("finv", a=a11, lda, n11)
      invlsdm <- z$a
    }
    else invlsdm <- ginv(lsdm)
    vcm <- invlsdm %*% fdm %*% t(invlsdm)
    var_gamma <- diag(vcm)[1:dim(z1)[2]]
    var_beta <- diag(vcm)[(dim(z1)[2] + 1):(dim(z1)[2] + dim(z2)[2])]
    var_kappa <- diag(vcm)[dim(z1)[2] + dim(z2)[2] + 1]
    sd_gamma <- sqrt(var_gamma)
    sd_beta <- sqrt(var_beta)
    sd_kappa <- sqrt(var_kappa)
  
  list(varga = var_gamma, varbe = var_beta, varka = var_kappa, sdga = sd_gamma, sdbe = sd_beta, sdka = sd_kappa)
}