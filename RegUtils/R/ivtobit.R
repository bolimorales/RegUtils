###   $Id: ivtobit.R 1126 2015-09-29  $
###
###            tobit with endogenous variables for R
###
###
### This file is part of the regUtils library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 2, or at your option, any later version,
### incorporated herein by reference.
###
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
###
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA

ivtobit <- function(formula, instruments, data, subset, na.action = NULL, weights, offset,
                    contrasts = NULL, model = TRUE, y = TRUE, x = FALSE, left=0,
                    right = Inf, method = "BHHH", ...)
{
  ## set up model.frame() call
  cl <- match.call()
  if(missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", "weights", "offset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE

  ## handle instruments for backward compatibility
  if(!missing(instruments)) {
    formula <- as.Formula(formula, instruments)
    cl$instruments <- NULL
    cl$formula <- formula(formula)
  } else {
    formula <- as.Formula(formula)
  }
  stopifnot(length(formula)[1] == 1L, length(formula)[2] %in% 1:2)

  ## try to handle dots in formula
  has_dot <- function(formula) inherits(try(terms(formula), silent = TRUE), "try-error")
  if(has_dot(formula)) {
    f1 <- formula(formula, rhs = 1)
    f2 <- formula(formula, lhs = 0, rhs = 2)
    if(!has_dot(f1) & has_dot(f2)) formula <- as.Formula(f1,
                                                         update(formula(formula, lhs = 0, rhs = 1), f2))
  }

  ## call model.frame()
  mf$formula <- formula
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  mtw <- delete.response(terms(formula, data = data, rhs = 1))
  W <- model.matrix(mtw, mf, contrasts)
  mt <- terms(formula, data = data)
  Y <- model.response(mf, "numeric")

  if (method == "twostep") {
    ## extract response, terms, model matrices
    y1 <- Y
    mtZ <- delete.response(terms(formula, data = data, rhs = 2))
    Z <- model.matrix(mtZ, mf, contrasts)
    rval = ivtobit_twostep.fit(y1,Z,W,left,right,...)
    rval$fitted.values = y1 - rval$residuals
  } else {   #Use maximun likelihood as default
    ## extract response, terms, model matrices
    mtX <- terms(formula, data = data, rhs = 1)
    X <- model.matrix(mtX, mf, contrasts)
    #Check if there are instruments
    if(length(formula)[2] < 2L) {
      mtZ <- NULL
      Z <- NULL
    } else {
      mtZ <- delete.response(terms(formula, data = data, rhs = 2))
      Z <- model.matrix(mtZ, mf, contrasts)
    }

    ## weights and offset
    weights <- model.weights(mf)
    offset <- model.offset(mf)
    if(is.null(offset)) offset <- 0
    if(length(offset) == 1) offset <- rep(offset, NROW(Y))
    offset <- as.vector(offset)

    ## call default interface
    rval <- ivtobit.fit(X, Y, Z, weights, offset, left, right, method, ...)

    # return starting values
    rval$start <- start
    rval$residuals = Y - W%*%rval$estimate[colnames(W)]
    rval$fitted.values = Y - rval$residuals
  }
  # save and return the call
  rval$call <- match.call()
  rval$method <- method
  # return the model terms
  rval$terms <- mt

  # save and return the number of oservations (in each category)
  obsBelow <- Y <= left
  obsAbove <- Y >= right
  obsBetween <- !obsBelow & !obsAbove
  rval$nObs <- c( sum( obsBelow ), sum( obsBetween ), sum( obsAbove ) )
  rval$nObs <- c( sum( rval$nObs ), rval$nObs )
  names(rval$nObs) <- c( "Total", "Left-censored", "Uncensored",
                         "Right-censored" )

  # return the degrees of freedom of the residuals
  rval$df.residual <- unname( rval$nObs[ 1 ] - length(coef(rval)))

  if (model)
    rval$model <- mf
  rval$real.v = length(mf)-1

  # censoring points
  rval$left <- left
  rval$right <- right
  rval$na.action <- na.action
  class(rval) <- c( "ivtobit", class(rval) )

  return (rval)
}

##Two Step estimation
ivtobit_twostep.fit <- function(y1,Z,W,left,right,...) {
  y2 <- W[,!colnames(W)%in%colnames(Z)]
  x1 = W[,colnames(W)%in%colnames(Z)]
  idx = !colnames(Z)%in%colnames(W)
  idx[1] = TRUE
  x2 = Z[,idx]
  X <- cbind(x1, x2)
  X <- X[, unique(colnames(X))]
  J1 <- matrix(0, nrow = ncol(X), ncol = ncol(x1))
  J2 <- matrix(0, nrow = ncol(X), ncol = ncol(x2))
  for (i in 1:ncol(x1)) J1[match(colnames(x1)[i], colnames(X)), i] <- 1
  for (i in 1:ncol(x2)) J2[match(colnames(x2)[i], colnames(X)), i] <- 1

  m1 <- lm(y2 ~ X-1)
  zi = cbind(x1,y2)
  m2 <- censReg(y1 ~ zi+m1$residuals-1, left=left, right = right)   #Aqui estan los betas
  PI1 <- m1$coefficients
  PI2 <- coef(m2)[-length(coef(m2))]
  m3 <- censReg(y1 ~ X+m1$residuals-1, left=left, right = right)
  PI3 <- coef(m3)[-length(coef(m3))]

  D = cbind(PI1, J1)
  n = ncol(X)
  alpha = PI3[1:n]
  lambda = PI3[n+1]
  Jaa = vcov(m3)[1:n,1:n]

  y2_hat = (lambda-PI2[n])*y2
  m1 <- lm(y2_hat ~ X-1)
  omega_hat = Jaa + vcov(m1)
  Var_delta = solve(t(D)%*%solve(omega_hat)%*%D)
  rval = list()
  rval$coefficients = PI2[1:n]
  names(rval$coefficients) = c(colnames(x1), colnames(W)[!colnames(W) %in% colnames(x1)])
  nn = c(colnames(W)[!colnames(W) %in% colnames(x1)], colnames(x1))
  rval$.vcov = Var_delta
  colnames(rval$.vcov) = row.names(rval$.vcov) = nn
  rval$.vcov = rval$.vcov[names(rval$coefficients),names(rval$coefficients)]
  rval$nObs = nrow(X)
  rval$residuals = y1 - W%*%rval$coefficients[colnames(W)]
  rval$df = length(y1) - ncol(X)
  return(rval)
}

#Prepare data for Maximun likelihood estimation
ivtobit.fit <- function(x, y, z, weights, offset, left, right, method, ...)
{
  ## model dimensions
  n <- NROW(y)
  p <- ncol(x)

  ## defaults
  if(missing(z)) z <- NULL
  if(missing(weights)) weights <- NULL
  if(missing(offset)) offset <- rep(0, n)

  ## sanity checks
  stopifnot(n == nrow(x))
  if(!is.null(z)) stopifnot(n == nrow(z))
  if(!is.null(weights)) stopifnot(n == NROW(weights))
  stopifnot(n == NROW(offset))

  ## project regressors x on image of instruments z
  if(!is.null(z)) {
    if(ncol(z) < ncol(x)) warning("more regressors than instruments")
    instruments = which(colnames(x)!=colnames(z))
    auxreg <- if(is.null(weights)) lm.fit(z, x, ...) else lm.wfit(z, x, weights, ...)
    xz <- as.matrix(auxreg$fitted.values)
    # pz <- z %*% chol2inv(auxreg$qr$qr) %*% t(z)
    colnames(xz) <- colnames(x)
    rdf <- auxreg$df.residual
    r <- auxreg$residuals

    w <- auxreg$weights
    if (is.null(w)) {
      rss <- sum(r^2)
    }
    else {
      rss <- sum(w * r^2)
    }
    resvar <- rss/rdf
    v_s = sqrt(resvar)

  } else {
    xz <- x
    # pz <- diag(NROW(x))
    # colnames(pz) <- rownames(pz) <- rownames(x)
  }

  ## main regression
  fit <- if(is.null(weights)) .fit.tobit(x, xz, y, v_s, instruments = instruments,
                                         offset = offset, left=left, right= right,
                                         method = method, ...)
  else
    stop("Error, weights with Maximun Likelihood are not implemented")
  rval <- fit
  rval$ninst = length(instruments)

  return(rval)
}

#Maximun likelihood estimation
.fit.tobit <- function(x, xz, y, v_s, instruments, offset, left, right,
                       method, ...) {
  vi = x[,instruments] - xz[,instruments]
  ni = length(instruments)
  b_i = (ni+1)*(ni+2)/2-2

  #Get start values
  start1 = censReg(y~xz-1, left = left, right = right)
  start=c(coef(start1)[-length(coef(start1))], exp(coef(start1)[length(coef(start1))]),
          rep(0,b_i),v_s)
  names(start) = c(colnames(x), "sigma_u", paste0("s", 1:(b_i+1)))

  #Optimize using maxLik package
  res = maxLik(.loglik.ivtobit, start=start, method=method, y = y, z = x,
               vi = vi, left = left, right = right, ni = ni)
  #Compute alpha
  k = ncol(x)
  last = length(start)
  s = coef(res)[(k+1):last]
  S = diag(ni+1)
  S[!lower.tri(S)] = s
  Sigma = crossprod(S)
  Sigma21 = Sigma[2:(ni+1),1]
  res$alpha = solve(Sigma[2:(ni+1),2:(ni+1)])%*%Sigma21
  res$sigma = Sigma
  return(res)
}

#Loglikehood function
.loglik.ivtobit = function(beta, y ,z, vi, left, right, ni) {
  k = ncol(z)
  last = length(beta)
  delta = beta[1:k]
  s = beta[(k+1):last]
  S = diag(ni+1)
  S[!lower.tri(S)] = s
  Sigma = crossprod(S)

  Sigma22 = Sigma[2:(ni+1),2:(ni+1)]
  Sigma21 = Sigma[2:(ni+1),1]
  sigmau.sq = diag(Sigma)[1]

  mi <- z%*%delta + vi%*%solve(Sigma22)%*%t(Sigma21)

  lnfy2i = -0.5*(log(2*pi)+log(det(as.matrix(Sigma22)))+diag(vi%*%(solve(Sigma22))%*%t(vi)))

  ## classify observations
  obsBelow <- y <= left
  obsAbove <- y >= right
  obsBetween <- !obsBelow & !obsAbove
  sigma_uv.sq = (sigmau.sq - t(Sigma21)%*%solve(Sigma22)%*%Sigma21)[1,1]
  sigma_uv = sqrt(sigma_uv.sq)

  ll <- rep(NA, length(y))
  ll[ obsBelow ] <-
    pnorm((left - mi[ obsBelow ]) / sigma_uv, log.p = TRUE)
  ll[obsBetween] <- -0.5*(log(2*pi)+log(sigma_uv.sq)+(y[obsBetween] - mi[obsBetween])^2/sigma_uv.sq)
  ll[obsAbove] <-
    pnorm((mi[obsAbove] - right) / sigma_uv, log.p = TRUE)

  return (ll + lnfy2i)
}

activePar.default <- function(x, ...) {
  if( !is.null( x$fixed ) ) {
    result <- !x$fixed
  } else {
    result <- x$activePar
  }
  if( is.null( result ) ) {
    result <- rep( TRUE, length( coef( x ) ) )
  }
  return( result )
}

summary.ivtobit <- function(object, robust=FALSE, ...) {
  if(!inherits(object, "ivtobit"))
    stop("'summary.ivtobit' called on a non-'ivtobit' object")
  if (object$method == "twostep") {
    if (robust) {
      warning("robust argument ignored on twostep method")
    }
    result = list()
    s <- sqrt(diag(object$.vcov))
    names(s) <- names(object$coefficients)
    t <- object$coefficients/s
    p <- 2*pnorm( -abs( t))
    results <- cbind("Estimate"=object$coefficients,
                     "Std. error"=s,
                     "t value"=t, "Pr(> t)"=p)
    result$estimate = results
    result$call <- object$call
    result$nObs <- object$nObs
    result$method <- object$method
    result$s <- sqrt(sum(object$residuals^2)/object$df)
    result$df <- object$df
    n = object$real.v
    wt = Wald.Test(b = coef(object)[1:n], Sigma = vcov(object)[1:n,1:n], Terms = 2:n)
    result$wald.test = wt
    class(result) <- c( "summary.ivtobit", class( result ) )
    return(result)
  }
  result <- NextMethod(summary, object)
  if (robust) {
    nParam <- length(object$estimate)
    activePar <- activePar(object)
    if((object$code < 100) & !is.null(object$estimate)) {
      # in case of infinity at initial values, the coefs are not provided
      if(!is.null(vc <- sandwich(object,...))) {
        s <- sqrt(diag(vc))
        names(s) <- names(object$estimate)
      } else {s=NULL}
      t <- object$estimate/s
      p <- 2*pnorm( -abs( t))
      t[!activePar(object)] <- NA
      p[!activePar(object)] <- NA
      results <- cbind("Estimate"=object$estimate,
                       "Std. error"=s,
                       "t value"=t, "Pr(> t)"=p)
      result$estimate = results
    }
  }
  lastRow = which(row.names(result$estimate) == "sigma_u")
  result$estimate = result$estimate[1:lastRow,]
  result$alpha = object$alpha
  result$sigma = object$sigma
  result$call <- object$call
  result$nObs <- object$nObs
  result$method <- object$method
  n = object$real.v
  wt = Wald.Test(b = coef(object)[1:n], Sigma = vcov(object)[1:n,1:n], Terms = 2:n)
  result$wald.test = wt
  class(result) <- c( "summary.ivtobit", class( result ) )
  return(result)
}

coef.ivtobit <- function( object, ... ) {
  if (object$method == "twostep") {
    return ( object$coefficients)
  } else {
    return( object$estimate )
  }
}

coef.summary.ivtobit <- function( object, ... ) {
  result <- object$estimate
  return( result )
}


print.summary.ivtobit <- function( x, logSigma = TRUE, digits = 4, ... ) {
  cat( "\n" )
  cat( "Call:\n" )
  cat( paste( deparse( x$call ), sep = "\n", collapse = "\n" ) )
  cat( "\n\n" )
  cat( "Observations:\n" )
  print( x$nObs )
  cat( "\n" )
  cat( "Coefficients:\n" )
  printCoefmat( coef( x, logSigma = logSigma ), digits = digits )
  if (x$method == "twostep") {
    cat( "\n" )
    cat(paste("\nResidual standard error:", round(x$s, digits),
              "on", x$df, "degrees of freedom\n"))
    v = x$wald.test[["result"]][["chi2"]]
    cat("Wald test:", formatC(v["chi2"], digits = digits),
        "on", v["df"],
        "DF,  p-value:", format.pval(v["P"], digits = digits), "\n\n")
  } else {
    cat( "\n" )
    cat( maximType( x ), ", ", nIter( x ), " iterations\n", sep = "" )
    cat( "Return code ", returnCode( x ), ": ", returnMessage( x ),
         "\n", sep = "" )
    cat( "Log-likelihood:", x$loglik, "on", sum( activePar( x ) ), "Df\n" )
    v = x$wald.test[["result"]][["chi2"]]

    cat("\nWald test:", formatC(v["chi2"], digits = digits),
        "on", v["df"],
        "DF,  p-value:", format.pval(v["P"], digits = digits), "\n\n")
    cat( "alpha: ",x$alpha)
    cat( "\n" )
    cat( "Sigma:\n" )
    print(x$sigma)
  }
  invisible( x )
}

print.ivtobit <- function( x, logSigma = TRUE, digits = 4, ... ) {

  cat( "\n" )
  cat( "Call:\n" )
  cat( paste( deparse( x$call ), sep = "\n", collapse = "\n" ) )
  cat( "\n\n" )
  cat( "Coefficients:\n" )
  print( coef( x, logSigma = logSigma ), digits = digits )
  cat( "\n" )
  invisible( x )
}

vcov.ivtobit <- function(object, eigentol=1e-12, ...) {
  if (object$method == "twostep") {
    return (object$.vcov)
  } else {
    return (NextMethod(vcov, object))
  }
}

fitted.ivtobit <- function(object, ...){
  if (is.null(object$na.action))
    object$fitted.values
  else napredict(object$na.action, object$fitted.values)
}
