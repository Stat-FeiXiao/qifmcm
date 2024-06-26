#' Print Function for qif Object
#'
#' Print a \code{qif} model object.
#'
#' @param x the \code{qif} model object.
#' @param digits number of digits to print.
#' @param quote logical, indicating whether or not strings should be printed with
#' surrounding quotes.
#' @param prefix string, only \code{""} is implemented.
#' @param ... further arguments passed to or from other methods.
#'
#' @return The invisible object from the arguments.
#'
#' @author Zhichang Jiang, Alberta Health Services, and Peter X.K. Song, University
#' of Michigan.
#'
#'
#' @seealso \code{\link[base]{print}} \code{\link[qif]{qif}}
#'
#' @export
#'
print.qif <- function(x, digits = NULL, quote = FALSE, prefix = "", ...)
{
    if(is.null(digits)) digits <- options()$digits else options(digits =
                                                                digits)
    cat("\n", x$title)
    cat("\n", x$version, "\n")
    cat("\nModel:\n")
    cat(" Link:                     ", x$model$link, "\n")
    cat(" Variance to Mean Relation:", x$model$varfun, "\n")
    if(!is.null(x$model$M))
        cat(" Correlation Structure:    ", x$model$corstr, ", M =", x$
            model$M, "\n")
    else cat(" Correlation Structure:    ", x$model$corstr, "\n")
    cat("\nCall:\n")
    dput(x$call)                        #       cat("\nTerms:\n")
    cat("\nNumber of observations : ", length(x$y), "\n")
    cat("\nNumber of clusters : ", length(x$nobs), "\n")
    cat("\nMaximum cluster size   : ", x$max.id, "\n")

    cat("\n\nStatistics:\n")
    print(x$statistics, digits = digits)

    cat("\n\nCoefficients:\n")
    print(x$coefficients, digits = digits)
    cat("\nEstimated Scale Parameter: ", x$scale, "\n")
    cat("\nNumber of Iterations: ", x$iteration,"\n")
    invisible(x)
}

#' Print the Summary of qif Object
#'
#' Prints the summary of a \code{qif} object.
#'
#' @param x the \code{qif} model object.
#' @param digits number of digits to print.
#' @param quote logical, indicating whether or not strings should be printed with
#' surrounding quotes.
#' @param prefix string, only \code{""} is implemented.
#' @param ... further arguments passed to or from other methods.
#'
#' @return The invisible object from the arguments.
#'
#' @author Zhichang Jiang, Alberta Health Services, and Peter X.K. Song, University
#' of Michigan.
#'
#' @return The invisible object from the arguments.
#'
#'
#' @seealso \code{\link[base]{summary}} \code{\link[qif]{qif}}
#'
#' @export
#'
print.summary.qif <- function(x, digits = NULL, quote = FALSE, prefix = "", ... )
{
    if(is.null(digits))
        digits <- options()$digits
    else options(digits = digits)
    cat("\n",x$title)
    cat("\n",x$version,"\n")
    cat("\nModel:\n")
    cat(" Link:                     ",x$model$link,"\n")
    cat(" Variance to Mean Relation:",x$model$varfun,"\n")
    if(!is.null(x$model$M))
        cat(" Correlation Structure:    ",x$model$corstr,", M =",x$model$M,"\n")
    else 	cat(" Correlation Structure:    ",x$model$corstr,"\n")
    cat("\nCall:\n")
    dput(x$call)
    cat("\nSummary of Residuals:\n")
    print(x$residual.summary, digits = digits)
    cat("\n\nCoefficients:\n")
    print(x$coefficients, digits = digits)
    cat("\nEstimated Scale Parameter: ", x$scale,"\n")
    cat("\nNumber of Iterations: ", x$iteration,"\n")
    invisible(x)
}

#' Summary of a qif Object
#'
#' Procuce a summary of a \code{qif} object.
#'
#' @param object an object for which a summary is desired.
#' @param correlation binary, include correlation.
#' @param ... additional arguements to be passed to \code{summary}.
#'
#' @return The \code{summary.qif} object.
#'
#' @author Zhichang Jiang, Alberta Health Services, and Peter X.K. Song, University
#' of Michigan.
#'
#' @importFrom stats quantile
#'
#' @seealso \code{\link[qif]{qif}}
#'
#'
#' @export
#'
summary.qif <- function(object, correlation = TRUE, ...)
{
    coef <- object$coefficients
    resid <- object$residuals
    n <- length(resid)
    summary <- list()
    summary$call <- object$call
    summary$version <- object$version
    summary$nobs <- length(object$nobs)
    summary$residual.summary <- quantile(as.vector(object$residuals))
    names(summary$residual.summary) <- c("Min", "1Q", "Median", "3Q", "Max")
    summary$model<- object$model
    summary$title <- object$title
    summary$coefficients <- object$parameter
    summary$scale <- object$scale
    summary$iteration <- object$iteration
    attr(summary,"class") <- "summary.qif"
    summary
}


#' Function to Solve a Quadratic Inference Function Model
#'
#' Produces an object of class "\code{qif}" which is a Quadratic Inference Function fit
#' of the balanced longitudinal data.
#'
#' \code{qif} provides two options of computing matrix inverses. The default
#' is from Fortran math library, and the other one is generalized inverse "\code{ginv}"
#' given in R package \code{MASS}. You can call option "\code{ginv}" through argument "\code{invfun}"
#' in "\code{qif()}".
#'
#' @param formula a formula expression as for other regression models, of the form
#' \code{response ~ predictors}. See the documentation of \code{\link[stats]{lm}}
#' and \code{\link[stats]{formula}} for details.
#' @param id a vector which identifies the clusters. The length of \code{id} should be
#' the same as the number of observations. Data are assumed to be sorted so that
#' observations on a cluster are contiguous rows for all entities in the formula.
#' @param data an optional data frame in which to interpret the variables occurring
#' in the \code{formula}, along with the \code{id} variables.
#' @param b an initial estimate for the parameters.
#' @param tol the tolerance used in the fitting algorithm.
#' @param maxiter the maximum number of iterations.
#' @param family a \code{family} object: a list of functions and expressions for defining
#' canonical link and variance functions. Families supported in \code{qif} are \code{gaussian},
#' \code{binomial}, \code{poisson}, and \code{gamma}; see the \code{\link[stats]{glm}} and \code{\link[stats]{formula}} documentation. Some links
#' are not currently available: probit link for binomial family and log link for
#' gamma family.
#' @param corstr a character string specifying the correlation structure. The
#' following are permitted: \code{"independence"}, \code{"exchangeable"}, \code{"AR-1"} and \code{"unstructured"}.
#' @param invfun a character string specifying the matrix inverse function. The
#' following are permitted: \code{"finv"} and \code{"ginv"}.
#'
#' @return A list containing:
#'
#' \itemize{
#'   \item{\code{title}: }{name of qif}
#'   \item{\code{version}: }{the current version of qif}
#'   \item{\code{model}: }{analysis model for link function, variance function and correlation struture}
#'   \item{\code{terms}: }{analysis model for link function, variance function and correlation struture}
#'   \item{\code{iteration}: }{the number of iterations}
#'   \item{\code{coefficients}: }{beta esitmates value}
#'   \item{\code{linear.perdictors}: }{linear predictor value}
#'   \item{\code{fitted.value}: }{fitted value of y}
#'   \item{\code{x}: }{the perdicted matrix}
#'   \item{\code{y}: }{the response}
#'   \item{\code{residuals}: }{y-mu}
#'   \item{\code{pearson.resi}: }{pearson residuals}
#'   \item{\code{scale}: }{the scale of fitted model}
#'   \item{\code{family}: }{the type of distribution}
#'   \item{\code{id}: }{model fitted value}
#'   \item{\code{max.id}: }{max number of each steps}
#'   \item{\code{xnames}: }{the values are X name of qif}
#'   \item{\code{statistics}: }{The qif statistics}
#'   \item{\code{Xnames}: }{the name X matrix in qif}
#'   \item{\code{parameter}: }{parameter estimates}
#'   \item{\code{covariance}: }{Covariance of coefficients}
#' }
#'
#' @note This R package is created by transplanting a SAS macro QIF developed
#' originally by Peter Song and Zhichang Jiang (2006). This is version 1.5 of
#' this user documentation file, revised 2019-07-02.
#'
#' @author Zhichang Jiang, Alberta Health Services, and Peter X.K. Song, University
#' of Michigan.
#'
#' @references
#' Qu A, Lindsay BG, Li B. Improving generalized estimating equations using quadratic
#' inference functions. Biometrika 2000, 87 823-836.
#'
#' Qu A, Song P X-K. Assessing robustness of generalised estimating equations and
#' quadratic inference functions. Biometrika 2004, 91 447-459.
#'
#' Qu A, Lindsay BG. Building adaptive estimating equations when inverse of covariance
#' estimation is difficult. J. Roy. Statist. Soc. B 2003, 65, 127-142.
#'
#' @examples
#' ## Marginal log-linear model for the epileptic seizures count data
#' ## (Diggle et al., 2002, Analysis of Longitudinal Data, 2nd Ed., Oxford Press).
#'
#' # Read in the epilepsy data set:
#' data(epil)
#'
#' # Fit the QIF model:
#' fit <- qif(y ~ base + trt + lage + V4, id=subject, data=epil,
#'                                        family=poisson, corstr="AR-1")
#' \donttest{
#' # Alternately, use ginv() from package MASS
#' fit <- qif(y ~ base + trt + lage + V4, id=subject, data=epil,
#'                       family=poisson, corstr="AR-1", invfun = "ginv")
#' }
#' # Print summary of QIF fit:
#' summary(fit)
#' \donttest{
#' ## Second example: MS study
#' data(exacerb)
#'
#' qif_BIN_IND<-qif(exacerbation ~ treatment + time + duration + time2, id=id,
#'                         data=exacerb, family=binomial, corstr="independence")
#' qif_BIN_AR1<-qif(exacerbation ~ treatment + time + duration + time2, id=id,
#'                         data=exacerb, family=binomial, corstr="AR-1")
#' qif_BIN_CS<-qif(exacerbation ~ treatment + time + duration + time2, id=id,
#'                         data=exacerb, family=binomial, corstr="exchangeable")
#' qif_BIN_UN<-qif(exacerbation ~ treatment + time + duration + time2, id=id,
#'                         data=exacerb, family=binomial, corstr="unstructured")
#'
#' summary(qif_BIN_CS)
#'
#' qif_BIN_CS$statistics
#'
#' qif_BIN_CS$covariance
#' }
#'
#' @seealso glm, lm, formula.
#' @importFrom MASS ginv
#' @importFrom stats model.matrix
#' @importFrom stats gaussian
#' @importFrom stats model.extract
#' @importFrom stats pchisq
#' @importFrom stats pnorm
#' @useDynLib qif
#' @export
#'
qif <- function (formula = formula(data), id = id, data = parent.frame(), b = NULL,
	tol = 1e-8, maxiter = 1000, family = gaussian, corstr = "independence", invfun="finv")
{

  if ((invfun!="finv")&&(invfun!="ginv")){
    stop("Unknown inverse function. Only finv or ginv.")
  }

  #if (invfun=="ginv") && (length(findFunction("ginv", generic = TRUE))==0) {
  #	stop("Please install library MASS first to apply the ginv function.", call. = FALSE )
  #}

  #if (invfun=="ginv") # library(MASS) # moved imported function to @importFrom

  # message("Beginning QIF function")

  call <- match.call()
  m <- match.call(expand.dots = FALSE)
  m$b <- m$tol <- m$maxiter <- m$link <- m$varfun <- m$corstr <- m$family <- m$invfun <- NULL

  if (is.null(m$id))
    m$id <- as.name("id")

  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  Terms <- attr(m, "terms")
  y <- as.matrix(model.extract(m, "response"))
  x <- model.matrix(Terms, m)

  QR <- qr(x)
  if (QR$rank < ncol(x))
    stop("rank-deficient model matrix")

  #    N <- rep(1, length(y))

  #    if (dim(y)[2] == 2) {
  #        N <- as.vector(y %*% c(1, 1))
  #        y <- y[, 1]
  #    }
  #    else {
  #        if (dim(y)[2] > 1)
  #            stop("Only one response (1 column)")
  #    }

  # only one response
  #message("\n","dim(y)[2]:")
  #print(dim(y)[2])

  if (dim(y)[2] > 1)
    stop("Only one response (1 column)")

  # ? get offset for?????
  #    offset <- model.extract(m, offset)
  #message("\n","offset:")
  #print(offset)

  # get class id
  id <- model.extract(m, id)
  if (is.null(id)) {
    stop("Id variable not found")
  }

  # get number of observe persons
  nobs <- nrow(x)

  # get number of parameters, include (intercept)
  np <- ncol(x)

  # get the name of X matrix colomn (variables)
  xnames <- dimnames(x)[[2]]

  # sign colomn name if not name for colomn
  if (is.null(xnames)) {
    xnames <- paste("x", 1:np, sep = "")
    dimnames(x) <- list(NULL, xnames)
  }

  if (is.character(family))
    family <- get(family)
  if (is.function(family))
    family <- family()

  # if b is not NULL then sign b to beta; otherwise, use gml to get initial beta
  if (!is.null(b)) {
    beta <- matrix(as.double(b), ncol = 1)
    if (nrow(beta) != np) {
      stop("Dim beta != ncol(x)")
    }
    # message("user's initial regression estimate")
    # print(beta) # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  }
  else {
    # message("\n","running glm to get initial regression estimate")
    mm <- match.call(expand.dots = FALSE)
    mm$b <- mm$tol <- mm$maxiter <- mm$link <- mm$varfun <- mm$corstr <- mm$id <- mm$invfun <- NULL
    mm[[1]] <- as.name("glm")
    beta <- eval(mm, parent.frame())$coef
    # print(beta) #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    beta <- as.numeric(beta)
  }
  if (length(id) != length(y))
    stop("Id and y not same length")

  #get the maximum number of iteration
  maxiter <- as.integer(maxiter)

  # Linkage function
  links <- c("identity", "log", "logit", "inverse")

  # Distribution family must be one of GAUSSIAN, POISSON, GAMMA or BIN
  fams <- c("gaussian", "poisson", "binomial", "Gamma")

  # Variance function
  varfuns <- c("constant", "mu", "mu(1-mu)", "mu^2")

  # Correlation structure variable CORR must be one of AR1, EXCH, IND or UNST
  corstrs <- c("independence", "exchangeable", "AR-1", "unstructured")

  famv <- match(family$family, fams, -1)

  dist <- family$family
  #cat("dist=",dist,"\n")
  #cat("corstr=",corstr,"\n")

  linkv <- as.integer(match(c(family$link), links, -1))

  # check the input parameter of family, corstr correct or not
  if (famv < 1)
    stop("unknown family")
  if (famv <= 4)
    varfunv <- famv
  else varfunv <- match(family$varfun, varfuns, -1)
  varfunv <- as.integer(varfunv)

  corstrv <- as.integer(match(corstr, corstrs, -1))

  if (linkv < 1)
    stop("unknown link.")
  if (varfunv < 1)
    stop("unknown varfun.")
  if (corstrv < 1)
    stop("unknown corstr.")

  ### start sign value to calculation ###
  #### define:
  # x- X matrix,
  # y- response matrix
  # beta - coefficients of X varibles
  # id- id list
  # uid- unique id
  # nobs- number of rep time
  # nsub- number of unique id
  # np - number of parameter
  # m0 - m0 matrix, nobs*nobs matrix, depend on correlation structure
  # m1 - m1 matrix, nobs*nobs matrix, depend on correlation structure

  # sumg - np*1 matrix for g
  # sumc - np*np matrix for c
  # arsumg - np*1 matrix
  # arsumc - np*np matrix
  # gi - np*1 matrix for i-th id
  # arsumgfirstdev - np*np matrix for first derivative of g
  # firstdev - np*np matrix of first derivative

  # for each id:
  # xi- i-th x
  # yi - i-th y
  # ni- number of row of xi

  # ui- i-th mu
  # fui - link function
  # fui_dev - devirate of link function
  # vui - variance of ui

  ### Calculation

  y <- as.matrix(y)
  x <- as.matrix(x)

  obs <- lapply(split(id, id), "length")
  nobs <- as.numeric(obs)
  nsub <- length(nobs)
  np <- dim(x)[[2]]

  ################ 	QIF iteration #######

  time1 <- date()

  betadiff <- 1
  iteration <- 0
  betanew <- beta

  # qif iteration
  while(betadiff > tol && iteration < maxiter)
  {

    # initial value
    beta <- betanew
    if (corstr == "independence") {
      sumg <- matrix(rep(0,np),nrow=np)
      sumc <- matrix(rep(0,np*np),nrow=np)
      arsumg <- matrix(rep(0,np),nrow=np)
      arsumc <- matrix(rep(0,np*np),nrow=np)
      gi <- matrix(rep(0,np),nrow=np)
      arsumgfirstdev <- matrix(rep(0,np*np),nrow=np)
      firstdev <- matrix(rep(0,np*np),nrow=np)
    }
    else {
      sumg <- matrix(rep(0,2*np),nrow=2*np)
      sumc <- matrix(rep(0,2*np*2*np),nrow=2*np)
      arsumg <- matrix(rep(0,2*np),nrow=2*np)
      arsumc <- matrix(rep(0,2*np*2*np),nrow=2*np)
      gi <- matrix(rep(0,2*np),nrow=2*np)
      arsumgfirstdev <- matrix(rep(0,2*np*np),nrow=2*np)
      firstdev <- matrix(rep(0,2*np*np),nrow=2*np)
    }
    # one iteration

    # set unstructured correlation, m0 and m1, only good for balance cluster
    if (corstr == "unstructured") {
      n1 <- nobs[1]
      m0 <- diag(n1)
      m1 <- matrix(rep(0,n1*n1),n1)

      loc1 <- 0
      loc2 <- 0
      for (i in 1:nsub) {
        # set start location for next xi
        loc1 <- loc2+1
        loc2 <- loc1+nobs[i]-1
        yi <- as.matrix(y[loc1:loc2,])
        xi <- x[loc1:loc2,]
        ni <- nrow(yi)

        if (dist == "gaussian") {
          ui <- xi %*% beta
        }
        else if (dist == "poisson") {
          ui <- exp(xi %*% beta)
        }
        else if (dist == "Gamma") {
          ui <- 1 / (xi %*% beta)
        }
        else if (dist == "binomial") {
          ui <- 1 /(1 + exp(- xi %*% beta))
        }
        m1 <- m1 + (yi-ui) %*% t(yi-ui)
      }
      m1 <- 1/nsub * m1
    }


    # loc1 - start location for row in xi
    # loc2 - end location for row in xi
    loc1 <- 0
    loc2 <- 0
    for (i in 1:nsub)
    {
      # set start location for next xi
      loc1 <- loc2+1
      loc2 <- loc1+nobs[i]-1
      #cat("loc1=",loc1, "  loc2=",loc2,"\n")


      yi <- as.matrix(y[loc1:loc2,])
      xi <- x[loc1:loc2,]
      ni <- nrow(yi)

      #print(xi)
      #cat("yi=",yi, "  xi=",xi,"   ni=", ni, "\n")

      # set m0, m1
      m0 <- diag(ni)
      # set m1 by corr structure
      if (corstr == "independence") {
        m1 <- matrix(rep(0,ni*ni),ni)
      }
      else if (corstr == "exchangeable") {
        m1 <- matrix(rep(1,ni*ni),ni) - m0
      }
      else if (corstr == "AR-1") {
        m1 <- matrix(rep(0,ni*ni),ni)
        for (k in 1:ni) {
          for (l in 1:ni) {
            if (abs(k-l)==1) m1[k,l] <-1
          }
        }
      }

      # change ui, fui, fui_dev, vui depending on distribution
      if (dist == "gaussian") {
        ui <- xi %*% beta
        fui <- ui
        fui_dev <- diag(ni)
        vui <- diag(ni)
        #cat("ui=",ui, " @fui=",fui, " @fui_dev=",fui_dev, " @vui=",vui," @i=",i,"\n")

      }
      else if (dist == "poisson") {
        ui <- exp(xi %*% beta)
        fui <- log(ui)
        fui_dev <- diag(as.vector(ui))
        vui <- diag(as.vector(sqrt(1/ui)))
        #cat("ui=",ui, " @fui=",fui, " @fui_dev=",fui_dev, " @vui=",vui," @i=",i,"\n")

      }
      else if (dist == "Gamma") {
        ui <- 1 / (xi %*% beta)
        fui <- 1 / ui
        fui_dev <- -diag(as.vector(ui)) %*% diag(as.vector(ui))
        vui <- diag(as.vector(1/ui))
        #cat("ui=",ui, " @fui=",fui, " @fui_dev=",fui_dev, " @vui=",vui," @i=",i,"\n")

      }
      else if (dist == "binomial") {
        ui <- 1 /(1 + exp(- xi %*% beta))
        fui <- log(ui) - log(1-ui)
        fui_dev <- diag(as.vector(ui)) %*% diag(as.vector(1-ui))
        vui <- diag(as.vector(sqrt(1/ui))) %*% diag(as.vector(sqrt(1/(1-ui))))
      }

      # calculate gi, wi, zi, c, arsumc, arsumg, di, firstdev, arsumgfirstdev
      # depending on corr structure
      if (corstr == "independence") {
        wi <- t(xi) %*% fui_dev %*% vui %*% m0 %*% vui
        gi0 <- (1/nsub)*wi %*% (yi-ui)
        gi[1:np,] <- gi0
        arsumc <- arsumc + gi %*% t(gi)
        arsumg <- arsumg + gi

        di0 <- -(1/nsub) * wi %*% fui_dev %*% xi
        firstdev[1:np,] <- di0
        arsumgfirstdev <- arsumgfirstdev + firstdev
      }

      else if (corstr == "unstructured") {
        wi <- t(xi) %*% fui_dev %*% m0
        zi <- t(xi) %*% fui_dev %*% m1

        gi0 <- (1/nsub)*wi %*% (yi-ui)
        gi1 <- (1/nsub)*zi %*% (yi-ui)

        gi[1:np,] <- gi0
        gi[(np+1):(2*np),] <- gi1

        arsumc <- arsumc + gi %*% t(gi)
        arsumg <- arsumg + gi

        if (is.na(arsumc[1,1])) {
          # print(iteration) # Commented out printing outside of print() methods
          # print(gi)
          # print(arsumc)
        }

        di0 <- -(1/nsub) * wi %*% fui_dev %*% xi
        di1 <- -(1/nsub) * zi %*% fui_dev %*% xi

        firstdev[1:np,] <- di0
        firstdev[(np+1):(2*np),] <- di1
        arsumgfirstdev <- arsumgfirstdev + firstdev
      }

      else {
        wi <- t(xi) %*% fui_dev %*% vui %*% m0 %*% vui
        zi <- t(xi) %*% fui_dev %*% vui %*% m1 %*% vui

        gi0 <- (1/nsub)*wi %*% (yi-ui)
        gi1 <- (1/nsub)*zi %*% (yi-ui)

        gi[1:np,] <- gi0
        gi[(np+1):(2*np),] <- gi1

        arsumc <- arsumc + gi %*% t(gi)
        arsumg <- arsumg + gi

        if (is.na(arsumc[1,1])) {
          # print(iteration) # Commented out printing outside of print() methods
          # print(gi)
          # print(arsumc)
        }

        di0 <- -(1/nsub) * wi %*% fui_dev %*% xi
        di1 <- -(1/nsub) * zi %*% fui_dev %*% xi

        firstdev[1:np,] <- di0
        firstdev[(np+1):(2*np),] <- di1
        arsumgfirstdev <- arsumgfirstdev + firstdev
      }
    }

    # after calculating all persons and sum them,
    # calculate Q, betanew,

    if (invfun=="finv") {
      a <- arsumc
      storage.mode(a)<-"double"
      lda<-as.integer(nrow(a))
      n<-as.integer(ncol(a))
      #print("a:")
      #print(a)
      z <- .Fortran("finv", a=a, lda, n)
      arcinv <- z$a
      #print("a inv:")
      #print(arcinv)
    }
    else arcinv=ginv(arsumc)

    Q <- t(arsumg) %*% arcinv %*% arsumg

    arqif1dev <- t(arsumgfirstdev) %*% arcinv %*% arsumg
    arqif2dev <- t(arsumgfirstdev) %*% arcinv %*% arsumgfirstdev


    if (invfun=="finv") {
      a <- arqif2dev
      storage.mode(a)<-"double"
      lda<-as.integer(nrow(a))
      n<-as.integer(ncol(a))
      z <- .Fortran("finv", a=a, lda, n)
      invarqif2dev <- z$a
    }
    else invarqif2dev <- ginv(arqif2dev)

    betanew <- beta - invarqif2dev %*% arqif1dev
    betadiff <- abs(sum(betanew - beta))
    iteration <- iteration +1
    #cat("Q=",Q, " @betanew=",betanew, " @betadiff=",betadiff, " @iteration=",iteration,"\n")

  }


  time2 <- date()
  #runtime <- time2-time1
  #cat("Runtime=",runtime,"\n")
  #cat("time1=",time1, " @time2=",time2,"\n")

  ################  QIF end   ###########

  ################ Output QIF parameter estimate and feacture  #########
  fit <- list()
  attr(fit, "class") <- c("qif", "glm")
  fit$title <- "QIF: Quadratic Inference Function"
  fit$version <- "QIF function, version 1.5 modified 2019/07/2"
  links <- c("Identity", "Logarithm", "Logit", "Reciprocal")
  varfuns <- c("Gaussian", "Poisson", "Binomial", "Gamma")
  corstrs <- c("Independent", "Exchangeable", "AR-1","Unstructured")
  fit$model <- list()
  fit$model$link <- links[linkv]
  fit$model$varfun <- varfuns[varfunv]
  fit$model$corstr <- corstrs[corstrv]
  fit$terms <- Terms
  fit$formula <- as.vector(attr(Terms, "formula"))
  fit$call <- call
  fit$nobs <- nobs
  fit$iteration <- iteration
  fit$coefficients <- as.vector(beta)
  names(fit$coefficients) <- xnames
  fit$linear.predictors <- as.matrix(x %*% beta)
  # for mu
  if (dist == "gaussian") {
    mu <- x %*% beta
    pearson <- y - mu
  }
  else if (dist == "poisson") {
    mu <- exp(x %*% beta)
    pearson <- (y - mu) / sqrt(mu)
  }
  else if (dist == "Gamma") {
    mu <- 1 / (x %*% beta)
    pearson <- (y - mu) / mu
  }
  else if (dist == "binomial") {
    mu <- 1 /(1 + exp(- x %*% beta))
    pearson <- (y - mu) / sqrt(mu * (1-mu))
  }
  fit$fitted.value <- as.matrix(mu)
  fit$residuals <- y - mu
  fit$pearson.resi <- pearson

  # scale value estimate
  fit$scale <- sum(pearson**2)/(length(y)-length(beta))

  fit$family <- family
  fit$y <- y
  fit$x <- x
  fit$id <- unique(id)
  fit$max.id <- max(nobs)
  fit$xnames <- xnames

  # for QIF goodness of fit test

  # TGI: trace of inverse matrix of variance matrix
  #tgi <- sum(diag(arqif2dev))
  #print(tgi)

  pvalue <- 1 - pchisq(Q,np)
  if (corstr == "indenpendence") {
    AIC <- Q
    BIC <- Q
  }
  else {
    AIC <- Q + 2*np
    BIC <- Q + np*log(nsub)
  }
  statistics <- c(Q, np, pvalue, AIC, BIC)
  names(statistics) <- c("Q", "D.F.", "pvalue","AIC","BIC")
  fit$statistics <- statistics

  # for beta, standarderror, Z, pvalue
  betase <- sqrt(diag(invarqif2dev))
  Z <- as.vector(beta)/betase
  betapvalue <- 2*(1-pnorm(abs(Z)))
  parameter <- cbind(beta, betase, Z, betapvalue)
  dimnames(parameter)[[2]] <- c("estimate", "stderr", "Z", "pvalue")
  dimnames(parameter)[[1]] <- xnames
  fit$parameter <- parameter

  # covariance matrix
  dimnames(invarqif2dev)[[1]] <- xnames
  dimnames(invarqif2dev)[[2]] <- xnames
  fit$Q2 <- arqif2dev
  seQ2 <- t(arsumgfirstdev) %*% arcinv
  fit$seQ2 <- seQ2
  fit
}
