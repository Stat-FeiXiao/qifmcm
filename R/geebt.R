geebt<-function(Status, Lambda, X, beta, w, id, corstr,invfun,eps) {
  K <- length(unique(id))
  n <- as.vector(table(id))
  gbeta <- rep(0, length(beta))
  newY1 <- Status/Lambda
  mu <- exp(X%*%beta)
  SK1 <- 1
  corstr<- corstr
  SK <- 1
  repeat {
    S1=newY1-mu
    M1 <- diag(n[1])
    M2 <-	matrix(1,n[1],n[1])
    diag(M2) <- 0
    
    W1=diag(w[id==1]*Lambda[id==1])
    D1=diag(mu[id==1])%*%diag(rep(1,n[1]))%*%(X[id==1,])
    
    V1=solve( sqrt(diag(mu[id==1])) )%*%M1%*%solve( sqrt(diag(mu[id==1])) ) 
    V2=solve( sqrt(diag(mu[id==1])) )%*%M2%*%solve( sqrt(diag(mu[id==1])) )      
    G=t(D1)%*%V1%*%W1%*%S1[id==1]
    G=rbind(G,t(D1)%*%V2%*%W1%*%S1[id==1])
    C=G%*%t(G)
    
    G1=-t(D1)%*%V1%*%W1%*%D1
    G1=rbind(G1,-t(D1)%*%V2%*%W1%*%D1)
    
    for(i in 2:K)
    {
      M1 <- diag(n[i])
      M2 <-	matrix(1,n[i],n[i])
      diag(M2) <- 0
      W1=diag(w[id==i]*Lambda[id==i])
      D1=diag(mu[id==i])%*%diag(rep(1,n[i]))%*%(X[id==i,])
      V1=solve( sqrt(diag(mu[id==i])) )%*%M1%*%solve( sqrt(diag(mu[id==i])) ) 
      V2=solve( sqrt(diag(mu[id==i])) )%*%M2%*%solve( sqrt(diag(mu[id==i])) )
      
      G_1=t(D1)%*%V1%*%W1%*%S1[id==i]
      G_1=rbind(G_1,t(D1)%*%V2%*%W1%*%S1[id==i])
      G=G+G_1
      C=C+G_1%*%t(G_1)
     
      G11=-t(D1)%*%V1%*%W1%*%D1
      G11=rbind(G11,-t(D1)%*%V2%*%W1%*%D1)
      G1=G1+G11
    }

    if (invfun=="finv") {
      a11 <- C
      storage.mode(a11)<-"double"
      lda<-as.integer(nrow(a11))
      n11<-as.integer(ncol(a11))
      z <- .Fortran("finv", a=a11, lda, n11)
      invC <- z$a
    }
    else invC <- ginv(C)
    
    Q =t(G)%*%invC%*%G
    Q1=t(G1)%*%invC%*%G
    Q2=t(G1)%*%invC%*%G1

    if (invfun=="finv") {
      a11 <- Q2
      storage.mode(a11)<-"double"
      lda<-as.integer(nrow(a11))
      n11<-as.integer(ncol(a11))
      z <- .Fortran("finv", a=a11, lda, n11)
      invQ2 <- z$a
    }
    else invQ2 <- ginv(Q2)
    
    newbeta=beta-invQ2%*%Q1
    
    if (any(abs(newbeta - beta) > eps) && (SK <= 500)) {
      beta <- newbeta
      mu <- exp(X%*%beta)
      SK <- SK + 1
    } else break
  }
  seQ2 = t(G1)%*%invC
  list(beta = newbeta,Q2=Q2,seQ2 = seQ2, C=C )
}
