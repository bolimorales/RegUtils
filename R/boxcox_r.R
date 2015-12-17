###   $Id: boxcox_r.R 1126 2015-09-29  $
###
###     Box-Cox regression for R
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

#Box Cox regression main function
boxcox.r <- function(formula, data, subset, na.action,
                     x = FALSE, y = FALSE, noTrans = NULL, optimize.bounds=c(-2,2),
                     model = "theta",test.params = NULL, ...)
{
  ret.x <- x
  ret.y <- y
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  formula_a = NULL
  model.o = model
  z = NULL

  #Check noTransformation argument
  if (!is.null(noTrans) && model=="lhs")
    warning(gettextf("noTrans = '%s' is ignored since method is 'lhs'",
                     noTrans), domain = NA)

  ######### Divide formula in two sections splitting transformed variables
  if (!is.null(noTrans) && model!="lhs") {
    if (length(noTrans)>1) {
      formula_a <- eval(substitute(update(b,~+a), list(b=formula,a=as.name(noTrans[1]))))
      for (abs_v in noTrans) {
        formula <- eval(substitute(update(b,~.-a), list(b=formula,a=as.name(abs_v))))
        formula_a <- eval(substitute(update(b,~.+a), list(b=formula_a,a=as.name(abs_v))))
      }
    }
    else {
      formula <- eval(substitute(update(b,~.-a), list(b=formula,a=as.name(noTrans))))
      formula_a <- eval(substitute(update(b,~+a), list(b=formula,a=as.name(noTrans))))
    }
    formula_a <- eval(substitute(update(b,~.-1), list(b=formula_a)))
  }
  #########
  m <- match(c("formula","data", "subset", "na.action"), names(mf), 0L)
  mf <- mf[c(1L,m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  if (!is.null(formula_a)) {
    m_y <- eval(mf, parent.frame())
    mf$formula <- formula_a
    mfa <- eval(mf, parent.frame())
    mf$formula <- formula
  }

  mf <- eval(mf, parent.frame())

  mt <- attr(mf, "terms")

  if (!is.null(formula_a)) {
    mta <- attr(mfa, "terms")
    y <- model.response(m_y, "numeric")
    index_set = names(y)
  } else {
    y <- model.response(mf, "numeric")
  }
  x <- model.matrix(mt, mf, contrasts)

  if (!is.null(formula_a)) {
    z <- model.matrix(mta, mfa, contrasts)
  }
  if (model == "lhs") {
    z <- x
    x = NULL
    model = "lambda"
  }

  #Check model and initilalize start values
  has.intercept = (attr(terms(formula), "intercept") == 1L)
  start_value = NULL
  if (model == "theta") {
    fcn.to.min <- .min.theta
    start_value = c(1,1)
  } else if (model == "lambda") {
    fcn.to.min <- .min.lambda
    start_value = 1
  } else if (model == "rhs") {
    fcn.to.min <- .min.rhs
    start_value = 1
  } else {
    stop("Model must be one of the following: theta, lambda, lhs or rhs")
  }

  #Check test parameters argument (no optimization just likelihood estimation)
  min.list = NULL
  if (!is.null(test.params)) {
    if (length(test.params) != 2){
      stop("Incorrect number of test parameters,
           input must be a vector: c(theta, lambda)")
    }
    else
      min.list$par = test.params
    }
  else
    min.list = nlminb(start = start_value, objective = fcn.to.min,
                      lower = optimize.bounds[1], upper = optimize.bounds[2],
                      data=list(x = x, y = y, z=z), has.intercept = has.intercept)

  res = list()
  class(res) <- c(if (is.matrix(y)) "mlm", "lm")
  params = NULL
  model.df = NULL
  model.r.df = NULL
  if (model == "theta") {
    res$lambda = min.list$par[2]
    res$theta = min.list$par[1]
    model.df = 1
  } else if (model == "lambda") {
    res$lambda = res$theta = min.list$par
    model.df = 0
  } else if (model == "rhs") {
    res$theta = 1
    res$lambda = min.list$par
    model.df = 1
  }

  outp = .model.parameters(c(res$theta, res$lambda), data=list(x = x, y = y, z=z), has.intercept=has.intercept)
  if (ret.x)
    res$x = outp$x
  if (ret.y)
    res$y = out$y_t
  res$model = model.o
  res$na.action <- attr(mf, "na.action")
  res$contrasts <- attr(z, "contrasts")
  res$xlevels <- .getXlevels(mt, mf)
  res$call <- cl
  res$terms <- mt
  res$residuals <- outp$residuals
  res$fitted.values <- outp$fitted
  res$coefficients = outp$coefs
  res$log_likelihood = outp$log_lik
  res$vcov = outp$vcov
  row.names(res$vcov) = row.names(res$coefficients)
  colnames(res$vcov) = row.names(res$coefficients)
  res$df.residual = length(y)
  res$observations = length(y)
  res$rank = outp$qr$rank
  res$qr = outp$qr
  res$df.model = model.df + res$rank - 1


  #Likelihood test for regression parameters
  min.list.g = nlminb(start = start_value, objective = fcn.to.min,
                      lower = optimize.bounds[1], upper = optimize.bounds[2],
                      data=list(x = x, y = y, z=z),
                      has.intercept = has.intercept, omit.x=2:res$rank)
  res$chi.sq = 2*(res$log_likelihood + min.list.g$objective)
  res$chi.sq.x = NULL
  for (i in 2:res$rank){
    min.list.x = nlminb(start = start_value, objective = fcn.to.min,
                        lower = optimize.bounds[1], upper = optimize.bounds[2],
                        data=list(x = x, y = y, z=z), has.intercept = has.intercept, omit.x=i)
    res$chi.sq.x = c(res$chi.sq.x,2*(res$log_likelihood + min.list.x$objective))
  }

  if (ret.x)
    res$x <- x
  if (ret.y)
    res$y <- y
  class(res) <- c("boxcox_r", class(res))
  return (res)
}

#Box-Cox transformation
.box_cox <- function(x, lambda) {
  if (!is.finite(lambda) || length(lambda) != 1)
    stop("'lambda' must be a non-missing, finite numeric scalar")
  if (any(x[!is.na(x)] <= 0))
    stop("All non-missing values of 'x' must be positive")
  return (as.matrix(if (abs(lambda)>.Machine$double.eps) (x^lambda-1)/lambda else log(x)))
}

#Loglikelihood function for theta model
.min.theta <- function(params, data, has.intercept, omit.x=NULL) {
  theta = params[1]
  lambda = params[2]
  y = data$y; x = data$x; z = data$z
  y_t = .box_cox(y, theta)
  N = length(y_t)
  if (!is.null(x)) {
    W_lambda = if (has.intercept) cbind(1,.box_cox(x[,2:ncol(x)], lambda), z) else cbind(.box_cox(x, lambda), z)
  } else {
    W_lambda = z
  }
  if (!is.null(omit.x)) W_lambda = W_lambda[,-omit.x]
  qw <- qr(W_lambda)
  d_hat <- solve.qr(qw, y_t)
  sigma_sq_hat = 1/N*crossprod((y_t-W_lambda%*%d_hat))
  log_lik = (-N/2)*(log(2*pi)+1+log(sigma_sq_hat))+(theta-1)*sum(log(y))
  -1*log_lik
}

#Loglikelihood function for lambda model
.min.lambda <- function(lambda, data, has.intercept, omit.x=NULL) {
  .min.theta(c(lambda, lambda), data=data, has.intercept = has.intercept,
            omit.x=omit.x)
}

#loglikehood functions for restricted models used on Wald Test
.min.rhs <- function(lambda, data, has.intercept, omit.x=NULL) {
  .min.theta(c(1, lambda), data=data, has.intercept = has.intercept,
            omit.x=omit.x)
}

#Estimation of residuals, fitted values, vcov
.model.parameters <- function(params, data, has.intercept) {
  theta = params[1];  lambda = params[2]
  y = data$y; x = data$x; z = data$z
  if (!is.null(x)) {
    W_lambda = if (has.intercept) cbind(1,.box_cox(x[,2:ncol(x)], lambda), z) else cbind(.box_cox(x, lambda), z)
  } else {
    W_lambda = z
  }
  y_t = .box_cox(y, theta)
  N = length(y_t)
  qw <- qr(W_lambda)
  d_hat <- solve.qr(qw, y_t)
  sigma_sq_hat = 1/N*crossprod((y_t-W_lambda%*%d_hat))
  row.names(d_hat)[1] = "(intercept)"
  residuals = .box_cox(y, theta) - W_lambda %*% d_hat
  fitted = .box_cox(y, theta) - residuals
  vcov = sigma_sq_hat[1,1] * chol2inv(qw$qr)
  log_lik = (-N/2)*(log(2*pi)+1+log(sigma_sq_hat))+(theta-1)*sum(log(y))
  return (list(y_t = y_t, coefs = d_hat, residuals = residuals, fitted = fitted,
               vcov=vcov, qr=qw, log_lik = log_lik, x = W_lambda))
}

summary.boxcox_r <- function (object, correlation = FALSE, symbolic.cor = FALSE,...)
{
  # construct the table of coefficients
  if (!is.null(object$vcov)){
    std.err <- sqrt(diag(object$vcov))
  }
  else{
    std.err <- sqrt(diag(vcov(object)))
  }
  b <- coefficients(object)
  chi.sq = c(NA,object$chi.sq.x)
  z = c(NA,rep(1,length(b)-1))
  p <- 2 * pchisq(abs(chi.sq), df = z, lower.tail = FALSE)
  object$coefficients <- cbind("Estimate"   = b,
                               "Chi.sq" = chi.sq,
                               "df Chi.sq"    = z,
                               "Pr(>|t|)"   = p)
  object$p.val = 2 * pchisq(abs(object$chi.sq), df = object$rank, lower.tail = FALSE)
  class(object) <- c("summary.boxcox_r", class(object))
  object
}

print.summary.boxcox_r <- function(x, digits = 4,...) {
  cat("Call:\n")
  print(x$call)
  cat("\nBox-Cox Estimates :\n")
  if (x$model == "theta") {
    cat("Theta = ", format(x$theta, digits = digits, nsmall = 4),
        "\nLambda = ", format(x$lambda, digits = digits, nsmall = 4), "\n", sep = "")
  } else if (x$model == "lambda" | x$model == "lhs" | x$model == "rhs") {
    cat("Lambda = ", format(x$lambda, digits = digits, nsmall = 4),"\n", sep = "")
  }
  cat("\nLR Test:\n")
  cat("X2 = ", format(x$chi.sq, digits = digits, nsmall = 1),
      ", df = ", x$df.model, ", P(> X2) = ", format(x$p.val, digits = digits,
                                                    nsmall = 1), "\n", sep = "")
  cat("\nCoefficients :\n")
  printCoefmat(x$coefficients)
  cat("---\nlog likehood = ", format(x$log_likelihood, digits = digits, nsmall = 1), "\n", sep = "")
}
