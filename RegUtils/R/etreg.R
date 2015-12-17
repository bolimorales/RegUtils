###   $Id: etreg.R 1126 2015-09-29  $
###
### Linear regression with endogenous treatment effects for R
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

#Main function
etreg <- function(formula, treat, data, subset, na.action,
                  contrasts = NULL, model = TRUE, y = TRUE, x = FALSE, method = "NR", ...)
{
  ret.x <- x
  ret.y <- y
  ## set up model.frame() call
  cl <- match.call()
  if(missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0)
  mf <- mf[c(1, m)]

  mf$drop.unused.levels <- TRUE

  ## call model.frame()
  c.formula = as.Formula(formula, treat)
  #Check correct formulas
  stopifnot(length(c.formula) == c(2L,2L))

  mf$formula <- c.formula
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  ## extract response, terms, model matrices
  y <- as.matrix(model.part(c.formula, data = mf, lhs=1))
  t <- as.matrix(model.part(c.formula, data = mf, lhs=2))
  mt <- terms(c.formula, data = data)

  mtX <- terms(c.formula, data = data, lhs=0,rhs = 1)
  x <- model.matrix(mtX, mf, contrasts)

  mtW <- terms(c.formula, data = data, lhs=0,rhs = 2)
  w <- model.matrix(mtW, mf, contrasts)

  ##Compute start values:
  t = as.integer(t)
  fit1 = lm(t~w-1)
  z = cbind(x,fit1$fitted.values)
  fit2 = lm(y~z-1)
  fit2s = summary(fit2)
  rhos = cov(fit2$residuals, fit1$residuals)/fit2s$sigma
  lnsigma = if (fit2s$sigma > 0) log(fit2s$sigma) else 1
  start = c(fit2$coefficients, fit1$coefficients, atanh(rhos),lnsigma)
  names(start) = c(colnames(x), attr(terms(c.formula, data=data), "term.labels")[2],
                   colnames(t), colnames(w), "athRho", "lnsigma")

  ## call default interface
  rval <- etreg.fit(y, t, x, w, start, method, ...)

  # save and return the call
  rval$call <- cl

  # return the model terms
  rval$terms <- mt

  # return the degrees of freedom of the residuals
  rval$df.residual <- unname(rval$nObs[1] - length(coef(rval)))

  # return starting values
  rval$start <- start
  rval$xlevels <- .getXlevels(mt, mf)
  rval$contrasts <- attr(x, "contrasts")

  if (model)
    rval$model <- mf
  if (ret.x)
    rval$x <- x
  if (ret.y)
    rval$y <- y
  rval$real.v = length(mf)-1
  class(rval) <- c("etreg", class(rval))
  return (rval)
}

#Fit function: Prepara data for maximun likelihood estimation
etreg.fit <- function(y, t, x, w, start, method, ...)
{
  ## model dimensions
  n <- NROW(y)
  p <- ncol(x)

  ## sanity checks
  stopifnot(n == nrow(x))
  if(!is.null(w)) stopifnot(n == nrow(w))

  ## main regression
  rval = maxLik(.loglik.etreg, start=start, method=method, y = y, tt=t, x= x,
                w = w)
  ## Save NObs
  rval$nObs = n

  return(rval)
}

#loglikelihood function
.loglik.etreg = function(params, y, tt, x, w) {
  p = ncol(x)
  k = ncol(w)
  betas = params[1:p]
  gamas = params[(p+2):(p+k+1)]
  delta = params[p+1]
  rho = tanh(params[p+k+2])
  sigma.sq = exp(2*params[p+k+3])

  #Classify observations
  obst1 = tt==1
  obst0 = !obst1

  ll <- rep(NA, length(y))
  A = (y[obst1] - x[obst1,] %*% betas - delta)/sqrt(sigma.sq)
  ll[obst1] <-
    pnorm((w[obst1,] %*% gamas + A * rho)/sqrt(1 - rho^2), log.p = TRUE) -
    0.5*A^2 - log(sqrt(2*pi*sigma.sq))
  B = (y[obst0] - x[obst0,] %*% betas)/sqrt(sigma.sq)
  ll[obst0] <-
    pnorm((-w[obst0,] %*% gamas - B * rho)/sqrt(1 - rho^2), log.p = TRUE) -
    0.5*B^2 - log(sqrt(2*pi*sigma.sq))
  return (ll)
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

summary.etreg <- function(object, robust=FALSE, eigentol=1e-12,... ) {
  if(!inherits(object, "etreg"))
    stop("'summary.etreg' called on a non-'etreg' object")
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
  } else {
    vc = vcov(object)
  }

  ### Delta method para Rho, sigma and lambda
  object$rho = deltaMethod(object, "tanh(athRho)", vcov. = vc)
  object$sigma = deltaMethod(object, "exp(lnsigma)", vcov. = vc)
  object$lambda = deltaMethod(object, "tanh(athRho)*exp(lnsigma)", vcov.=vc)

  v = rbind(object$rho,object$sigma,object$lambda)
  t = unname(unlist(c(object$rho[1]/object$rho[2],
                      object$sigma[1]/object$sigma[2],
                      object$lambda[1]/object$lambda[2])))
  pv = 2*pnorm( -abs(t))
  estim = cbind(v,t,pv)
  row.names(estim) = c("rho", "sigma","lambda")
  colnames(estim) = c("Estimate", "Std. error", "t value", "Pr(> t)")
  result$estimate = rbind(result$estimate,
                          estim)
  class(result) <- c( "summary.etreg", class( result ) )
  result$call <- object$call
  result$nObs <- object$nObs
  n = object$real.v
  wt = Wald.Test(b = coef(object)[1:n], Sigma = vcov(object)[1:n,1:n], Terms = 2:n)
  result$wald.test = wt
  return(result)
}

print.summary.etreg <- function( x, logSigma = TRUE, digits = 4, ... ) {

  cat( "\n" )
  cat( "Call:\n" )
  cat( paste( deparse( x$call ), sep = "\n", collapse = "\n" ) )
  cat( "\n\n" )
  cat( "Observations:\n" )
  print( x$nObs )
  cat( "\n" )
  cat( "Coefficients:\n" )
  printCoefmat( coef( x, logSigma = logSigma ), digits = digits )
  cat( "\n" )
  cat( maximType( x ), ", ", nIter( x ), " iterations\n", sep = "" )
  cat( "Return code ", returnCode( x ), ": ", returnMessage( x ),
       "\n", sep = "" )
  cat( "Log-likelihood:", x$loglik, "on", sum( activePar( x ) ), "Df" )
  v = x$wald.test[["result"]][["chi2"]]

  cat("\nWald test:", formatC(v["chi2"], digits = digits),
      "on", v["df"],
      "DF,  p-value:", format.pval(v["P"], digits = digits), "\n\n")
  cat( "\n" )
  invisible( x )
}

print.etreg <- function( x, logSigma = TRUE, digits = 4, ... ) {
  rho = unname(tanh(coef(x)["athRho"]))
  sigma = unname(exp(coef(x)["lnsigma"]))
  lambda = rho*sigma
  cat( "\n" )
  cat( "Call:\n" )
  cat( paste( deparse( x$call ), sep = "\n", collapse = "\n" ) )
  cat( "\n\n" )
  cat( "Coefficients:\n" )
  print( coef( x, logSigma = logSigma ), digits = digits )
  print(list(rho = rho, sigma = sigma,
             lambda = lambda))
  cat( "\n" )
  invisible( x )
}
