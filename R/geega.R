geega<-function(w, Z, gamma, id, corstr) {
  gamma1 <- rep(0, length(gamma))
  id <- id
  K <- length(unique(id))
  n <- as.vector(table(id))
  alpha <- 0
  phi <- 1
  gamma<-as.matrix(gamma)
  P <- exp(Z%*%gamma )/(1 + exp(Z%*%gamma) )
  PP <- P * (1 - P)
  z1<-Z
  p1<-P
  g1<-w
  repeat {
    S1 <-  w - P
    M1 <- diag(n[1])
    M2 <-	matrix(1,n[1],n[1])
    diag(M2) <- 0
    D1= diag(c(PP[id == 1]))%*% diag(rep(1, n[1])) %*% Z[id == 1, ]
    
    V1=ginv(sqrt(diag(c(PP[id == 1])))) %*% M1 %*% ginv(sqrt(diag(c(PP[id == 1]))))
    V2=ginv(sqrt(diag(c(PP[id == 1])))) %*% M2 %*% ginv(sqrt(diag(c(PP[id == 1]))))   
    G=t(D1)%*%V1%*%S1[id==1]  
    G=rbind(G,t(D1)%*%V2%*%S1[id==1]) 
    C=G%*%t(G)
    
    G1=-t(D1)%*%V1%*%D1
    G1=rbind(G1,-t(D1)%*%V2%*%D1)
    for(i in 2:K)
    {
      M1 <- diag(n[i])
      M2 <-	matrix(1,n[i],n[i])
      diag(M2) <- 0
      D1= diag(c(PP[id == i]))%*% diag(rep(1, n[i])) %*% Z[id == i, ]
      V1=ginv(sqrt(diag(c(PP[id == i])))) %*% M1 %*% ginv(sqrt(diag(c(PP[id == i]))))
      V2=ginv(sqrt(diag(c(PP[id == i])))) %*% M2 %*% ginv(sqrt(diag(c(PP[id == i]))))    
      G_1=t(D1)%*%V1%*%S1[id==i]  
      G_1=rbind(G_1,t(D1)%*%V2%*%S1[id==i]) 
      G=G+G_1
      C=C+G_1%*%t(G_1)
      
      G11=-t(D1)%*%V1%*%D1
      G11=rbind(G11,-t(D1)%*%V2%*%D1)
      G1=G1+G11
      
      
    }
    ################
    ABC1 <- rep(0, K)       
    VA1 <- matrix(0, dim(z1)[2], dim(z1)[2])
    for (v in 1:(dim(z1)[2])) {
      for (w1 in 1:(dim(z1)[2])) {
        for (i in 1:K) {
          IR1 <- diag(n[i])
          
          B1 <- matrix(0, n[i], n[i])
          z11 <- matrix(z1[id == i, ], nrow = n[i], )
          A1 <- t(z11[, v])
          p11 <- p1[id == i]
          g11 <- g1[id == i]
          pp11 <- p11 * (1 - p11)
          BB <- (pp11^(1/2)) %*% ((t(pp11))^(-1/2)) * IR1
          for (s in 1:n[i]) {
            for (l in 1:n[i]) {
              B1[s, l] <- (1/2) * (z11[s, w1] * (1 - 2 * p11[s]) - z11[l, w1] * (1 - 2 * p11[l])) * BB[s, l]
            }
          }
          C1 <- g11 - p11
          D1 <- BB
          E1 <- z11[, w1] * pp11
          ABC1[i] <- A1 %*% (B1 %*% C1 - D1 %*% E1)
        }
        VA1[v, w1] <- sum(ABC1) 
        ABC1 <- rep(0, K)
      }
    }
    ###
    ABC2 <- rep(0, K)
    VA2 <- matrix(0, dim(z1)[2], dim(z1)[2])
    for (v in 1:(dim(z1)[2])) {
      for (w1 in 1:(dim(z1)[2])) {
        for (i in 1:K) {
          IR1 <- matrix(1, n[i], n[i])
          diag(IR1) <- 0
          
          B1 <- matrix(0, n[i], n[i])
          z11 <- matrix(z1[id == i, ], nrow = n[i], )
          A1 <- t(z11[, v])
          p11 <- p1[id == i]
          g11 <- g1[id == i]
          pp11 <- p11 * (1 - p11)
          BB <- (pp11^(1/2)) %*% ((t(pp11))^(-1/2)) * IR1
          for (s in 1:n[i]) {
            for (l in 1:n[i]) {
              B1[s, l] <- (1/2) * (z11[s, w1] * (1 - 2 * p11[s]) - z11[l, w1] * (1 - 2 * p11[l])) * BB[s, l]
            }
          }
          C1 <- g11 - p11
          D1 <- BB
          E1 <- z11[, w1] * pp11
          ABC2[i] <- A1 %*% (B1 %*% C1 - D1 %*% E1)
        }
        VA2[v, w1] <- sum(ABC2) 
        ABC2 <- rep(0, K)
      }
    }
    #############
 #   G1<-rbind(VA1,VA2)     #G1_ga=\sum_i g关于beta的导
    
    Q =t(G)%*%MASS::ginv(C)%*%G
    Q1=t(G1)%*%MASS::ginv(C)%*%G
    Q2=t(G1)%*%MASS::ginv(C)%*%G1
    newgamma=gamma-MASS::ginv(Q2)%*%Q1
    
    
    if (any(abs(newgamma - gamma) > 1e-06)) {
      gamma <- newgamma
      P <- exp(Z%*%gamma)/(1 + exp(Z%*%gamma) )
      PP <- P * (1 - P)
    } else break
  }

  seQ2 <- t(G1)%*%MASS::ginv(C)
  
  list(gamma = newgamma,Q2=Q2,seQ2=seQ2)
}
