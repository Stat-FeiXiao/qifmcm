qifcure <- function(formula, cureform, data, id, corstr, Var = TRUE, stdz = FALSE, boots = FALSE, invfun="finv", nboot = 100, esmax = 100, eps = 1e-06) {
    call <- match.call()
    data <- data
    id <- id
    uid <- sort(unique(id))
    newid <- rep(0, length(id))
    for (i in 1:length(id)) {
        j <- 1
        repeat {
            if (id[i] != uid[j]) 
                j <- j + 1 else {
                newid[i] <- j
                break
            }
        }
    }
    data$id <- newid    
    data1 <- data[data$id == 1, ]
    for (i in 2:length(uid)) {
        data1 <- rbind(data1, data[data$id == i, ])
    }    
    data <- data1
    id <- data$id    
    Kn <- length(id)
    K <- length(unique(id))
    n <- as.vector(table(id))
    mf <- model.frame(formula, data)
    mf1 <- model.frame(cureform, data)
    Z <- model.matrix(attr(mf1, "terms"), mf1)
    
    X <- model.matrix(attr(mf, "terms"), mf)
        
    gamma_name <- colnames(Z)  
    gamma_length <- ncol(Z)  
    beta_name <- colnames(X) 
    beta_length <- ncol(X)
    Y <- model.extract(mf, "response")
    if (!inherits(Y, "Surv")) 
        stop("Response must be a survival object")
    Time <- Y[, 1]      
    Status <- Y[, 2]    
    stdz <- stdz
    esmax <- esmax
    eps <- eps
    esfit <- es(Time, Status, X, Z, id,  corstr, stdz, esmax, eps,invfun)
    gamma <- esfit$gamma
    names(gamma) <- gamma_name
    beta <- esfit$beta
    w <- esfit$w
    kappaa <- esfit$kappaa
    alpha2 <- exp(beta[1])
    Lambda <- esfit$Lambda
    uncuresurv <- esfit$uncuresurv
    baseuncuresurv <- esfit$basesurv  
    gamma_Q2 <- esfit$gamma_Q2 
    beta_Q2  <- esfit$beta_Q2 
    gamma_seQ2 <- esfit$gamma_seQ2
    beta_seQ2 <- esfit$beta_seQ2
    arrayRM <- esfit$arrayRM
    gS<-esfit$gS
    if (Var) {
        varfit <- varest(Time, Status, X, Z, id, gamma, beta, Lambda, w,kappaa,gSfordev,uncuresurv,baseuncuresurv,gamma_Q2 ,beta_Q2 ,gamma_seQ2,beta_seQ2,invfun,corstr,arrayRM,eps)
        var_gamma <- varfit$varga
        var_beta <- varfit$varbe
        var_kappa <- varfit$varka
        var_alpha2 <- (exp(beta[1]))^2 * var_beta[1]
        sd_gamma <- varfit$sdga
        sd_beta <- varfit$sdbe
        sd_kappa <- varfit$sdka
        sd_alpha2 <- sqrt(var_alpha2)
        fdm <- varfit$fdm
        lsdm <- varfit$lsdm
        LV <- varfit$LV
    }    
    if (boots) {
        Bootsample <- nboot
        stdz <- stdz
        corstr <- corstr
       
        BMc <- matrix(0, Bootsample, ncol(Z))
        BMs <- matrix(0, Bootsample, ncol(X))
        BMnus <- matrix(0, Bootsample, 4)
            
        for (rrn in 1:Bootsample) {
        repeat{  bootid<- sample((1:K), replace = TRUE)
                 bootdata <- data[id == bootid[1], ]
                 bootdata$id <- rep(1, sum(id == bootid[1]))
                 for (ll in 2:K) {
                 bootdata1 <- data[id == bootid[ll], ]
                 bootdata1$id <- rep(ll, sum(id == bootid[ll]))
                 bootdata <- rbind(bootdata, bootdata1)
                 }            
                 id_boot <- bootdata$id
                 Kn_boot <- length(id_boot)
                 K_boot <- length(unique(id_boot))
                 n_boot <- as.vector(table(id_boot))
                 mf <- model.frame(formula, bootdata)
                 mf1 <- model.frame(cureform, bootdata)
                 Z <- model.matrix(attr(mf1, "terms"), mf1)
                 
              
                 X <- model.matrix(attr(mf, "terms"), mf)[, -1]
                 X <- as.matrix(X)
                 colnames(X) <- colnames(model.matrix(attr(mf, "terms"), mf))[-1]
                             
                 Y <- model.extract(mf, "response")
                 if (!inherits(Y, "Surv")) 
                 stop("Response must be a survival object")
                 Time <- Y[, 1]
                 Status <- Y[, 2]  
                 tryboot <- try(es(Time, Status, X, Z, id = id_boot, corstr, stdz, esmax = 100, eps,invfun), silent = TRUE)  
                 if(is(tryboot,"try-error") == FALSE)
                 break
              }
            esfitboot <- tryboot
            BMc[rrn, ] <- esfitboot$gamma
            BMs[rrn, ] <- esfitboot$beta
           
            
        }
        var_gamma_boots <- apply(BMc, 2, var)
        sd_gamma_boots <- sqrt(var_gamma_boots)
        var_beta_boots <- apply(BMs, 2, var)
        sd_beta_boots <- sqrt(var_beta_boots)
        var_gcor_boots <- var(BMnus[, 1])
        sd_gcor_boots <- sqrt(var_gcor_boots)
        var_bcor_boots <- var(BMnus[, 2])
        sd_bcor_boots <- sqrt(var_bcor_boots)
        
    }    
    
    fit <- list()
    class(fit) <- c("geecure")    
    fit$gamma <- gamma
    if (Var) {
        fit$gamma_var <- var_gamma
        fit$gamma_sd <- sd_gamma
        fit$gamma_zvalue <- gamma/sd_gamma
        fit$gamma_pvalue <- (1 - pnorm(abs(fit$gamma_zvalue))) * 2
    }
    fit$beta <- beta
    if (Var) {
        fit$beta_var <- var_beta
        fit$beta_sd <- sd_beta
        fit$beta_zvalue <- beta/sd_beta
        fit$beta_pvalue <- (1 - pnorm(abs(fit$beta_zvalue))) * 2
        fit$LV<-LV
    }    
   
    fit$kappa <- kappaa
    fit$alpha2 <- alpha2
    if (Var) {
      fit$kappa_var <- var_kappa
      fit$kappa_sd <- sd_kappa
      fit$kappa_zvalue <- kappaa/sd_kappa
      fit$kappa_pvalue <- (1 - pnorm(abs(fit$kappa_zvalue))) * 2
      fit$alpha2_var <- var_alpha2
      fit$alpha2_sd <- sd_alpha2
      fit$alpha2_zvalue <- alpha2/sd_alpha2
      fit$alpha2_pvalue <- (1 - pnorm(abs(fit$alpha2_zvalue))) * 2
    }
    
    fit$num_of_clusters <- K
    fit$max_cluster_size <- max(n)   

     
    if (boots) {
        fit$boots_gamma_sd <- sd_gamma_boots
        fit$boots_beta_sd <- sd_beta_boots
        fit$boots_gcor_sd <- sd_gcor_boots
        fit$boots_bcor_sd <- sd_bcor_boots
        fit$boots_gamma_zvalue <- gamma/sd_gamma_boots
        fit$boots_beta_zvalue <- beta/sd_beta_boots
        fit$boots_gcor_zvalue <- gcor/sd_gcor_boots
        fit$boots_bcor_zvalue <- bcor/sd_bcor_boots
        fit$boots_gamma_pvalue <- (1 - pnorm(abs(fit$boots_gamma_zvalue))) * 2
        fit$boots_beta_pvalue <- (1 - pnorm(abs(fit$boots_beta_zvalue))) * 2
        fit$boots_gcor_pvalue <- (1 - pnorm(abs(fit$boots_gcor_zvalue))) * 2
        fit$boots_bcor_pvalue <- (1 - pnorm(abs(fit$boots_bcor_zvalue))) * 2
        fit$boots_kappa_sd <- sd_kappa_boots
        fit$boots_kappa_zvalue <- kappa/sd_kappa_boots
        fit$boots_kappa_pvalue <- (1 - pnorm(abs(fit$boots_kappa_zvalue))) * 2
        fit$boots_alpha2_sd <- sd_alpha2_boots
        fit$boots_alpha2_zvalue <- exp(beta[1])/sd_alpha2_boots
        fit$boots_alpha2_pvalue <- (1 - pnorm(abs(fit$boots_alpha2_zvalue))) * 2
    }    
    fit$call <- call   
    fit$gamma_name <- gamma_name
    fit$beta_name <- beta_name    
    fit$Time <- Time
    fit$Var <- Var
    fit$boots <- boots
    class(fit) = "geecure"
    fit
}

print.geecure <- function(x, ...) {
    if (!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl)
    }
    cat("\nCure Probability Model:\n")
    if (x$Var) {
        if (x$boots) {
            gm <- array(x$gamma, c(length(x$gamma), 4))
            rownames(gm) <- x$gamma_name
            colnames(gm) <- c("Estimate", "Std.Error", "Z value", "Pr(>|Z|)")
            gm[, 2] <- x$boots_gamma_sd
            gm[, 3] <- x$boots_gamma_zvalue
            gm[, 4] <- x$boots_gamma_pvalue
        }
		else {
            gm <- array(x$gamma, c(length(x$gamma), 4))
            rownames(gm) <- x$gamma_name
            colnames(gm) <- c("Estimate", "Std.Error", "Z value", "Pr(>|Z|)")
            gm[, 2] <- x$gamma_sd
            gm[, 3] <- x$gamma_zvalue
            gm[, 4] <- x$gamma_pvalue
        }
    }
	else {
        if (x$boots) {
            gm <- array(x$gamma, c(length(x$gamma), 4))
            rownames(gm) <- x$gamma_name
            colnames(gm) <- c("Estimate", "Std.Error", "Z value", "Pr(>|Z|)")
            gm[, 2] <- x$boots_gamma_sd
            gm[, 3] <- x$boots_gamma_zvalue
            gm[, 4] <- x$boots_gamma_pvalue
        }
		else {
            gm <- array(x$gamma, c(length(x$gamma), 1))
            rownames(gm) <- x$gamma_name
            colnames(gm) <- "Estimate"
        }
    }
    print(gm)
    cat("\n")
    cat("\nFailure Time Distribution Model:\n")
    if (x$Var) {

            if (x$boots) {
                bt <- array(x$beta, c(length(x$beta), 4))
                rownames(bt) <- x$beta_name
                colnames(bt) <- c("Estimate", "Std.Error", "Z value", "Pr(>|Z|)")
                bt[, 2] <- x$boots_beta_sd
                bt[, 3] <- x$boots_beta_zvalue
                bt[, 4] <- x$boots_beta_pvalue
            }
			else {
                bt <- array(x$beta, c(length(x$beta), 4))
                rownames(bt) <- x$beta_name
                colnames(bt) <- c("Estimate", "Std.Error", "Z value", "Pr(>|Z|)")
                bt[, 2] <- x$beta_sd
                bt[, 3] <- x$beta_zvalue
                bt[, 4] <- x$beta_pvalue
            }
        
	
    }
	else {

            if (x$boots) {
                bt <- array(x$beta, c(length(x$beta), 4))
                rownames(bt) <- x$beta_name
                colnames(bt) <- c("Estimate", "Std.Error", "Z value", "Pr(>|Z|)")
                bt[, 2] <- x$boots_beta_sd
                bt[, 3] <- x$boots_beta_zvalue
                bt[, 4] <- x$boots_beta_pvalue
            }
			else {
                bt <- array(x$beta, c(length(x$beta), 1))
                rownames(bt) <- x$beta_name
                colnames(bt) <- "Estimate"
            }
        
	
    }
    print(bt)

    
    cat("Number of clusters:", x$num_of_clusters)
    cat("       Maximum cluster size:", x$max_cluster_size)
    cat("\n")
    invisible(x)
}
