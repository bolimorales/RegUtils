###   $Id: eivlm.R 1126 2015-09-29  $
###
###     Errors in variables regression for R
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
eivlm <- function(formula, data, subset, weights, na.action,
                  model = TRUE, x = FALSE, method = NULL, y = FALSE,
                  singular.ok = TRUE, contrasts = NULL, rel = NULL, offset, ...) {
  ret.x <- x
  ret.y <- y
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action",
               "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  if (!is.null(method) && method == "model.frame")
    return(mf)

  mt <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  w <- as.vector(model.weights(mf))
  if(!is.null(rel))
    rel = c(1,rel)
  if (!is.null(w) && !is.numeric(w))
    stop("'weights' must be a numeric vector")
  offset <- as.vector(model.offset(mf))
  if (!is.null(offset)) {
    if (length(offset) != NROW(y))
      stop(gettextf("number of offsets is %d, should equal %d (number of observations)",
                    length(offset), NROW(y)), domain = NA)
  }
  if (is.empty.model(mt)) {
    x <- NULL
    z <- list(coefficients = if (is.matrix(y)) matrix(, 0,
                                                      3) else numeric(), residuals = y, fitted.values = 0 *
                y, weights = w, rank = 0L, df.residual = if (!is.null(w)) sum(w !=
                                                                                0) else if (is.matrix(y)) nrow(y) else length(y))
    if (!is.null(offset)) {
      z$fitted.values <- offset
      z$residuals <- y - offset
    }
  }
  else {
    x <- model.matrix(mt, mf, contrasts)
    z <- if (is.null(w))
      eivlm.fit(x, y, offset = offset, singular.ok = singular.ok,
                rel=rel,...)
    else eivlm.wfit(x, y, w, offset = offset, singular.ok = singular.ok,
                    rel=rel,...)
  }
  class(z) <- c(if (is.matrix(y)) "mlm", "lm")
  z$na.action <- attr(mf, "na.action")
  z$offset <- offset
  z$contrasts <- attr(x, "contrasts")
  z$xlevels <- .getXlevels(mt, mf)
  z$call <- cl
  z$terms <- mt
  if (model)
    z$model <- mf
  if (ret.x)
    z$x <- x
  if (ret.y)
    z$y <- y
  z$qr <- NULL
  class(z) <- c("eivlm", class(z))
  z
}

#Estimation function
eivlm.wfit <- function (x, y, w, offset = NULL, tol = 1e-07,
                        singular.ok = TRUE, rel=NULL, ...)
{
  z = list()
  if (is.null(n <- nrow(x)))
    stop("'x' must be a matrix")
  if (n == 0)
    stop("0 (non-NA) cases")
  ny <- NCOL(y)
  if (is.matrix(y) && ny == 1L)
    y <- drop(y)
  if (!is.null(offset))
    y <- y - offset
  if (NROW(y) != n | length(w) != n)
    stop("incompatible dimensions")
  if (any(w < 0 | is.na(w)))
    stop("missing or negative weights not allowed")

  dots <- list(...)
  if (length(dots) > 1L)
    warning("extra arguments ", paste(sQuote(names(dots)),
                                      sep = ", "), " are disregarded.", domain = NA)
  else if (length(dots) == 1L)
    warning("extra argument ", sQuote(names(dots)), " is disregarded.",
            domain = NA)
  x.asgn <- attr(x, "assign")
  zero.weights <- any(w == 0)
  if (zero.weights) {
    save.r <- y
    save.f <- y
    save.w <- w
    ok <- w != 0
    nok <- !ok
    w <- w[ok]
    x0 <- x[!ok, , drop = FALSE]
    x <- x[ok, , drop = FALSE]
    n <- nrow(x)
    y0 <- if (ny > 1L)
      y[!ok, , drop = FALSE]
    else y[!ok]
    y <- if (ny > 1L)
      y[ok, , drop = FALSE]
    else y[ok]
  }
  p <- ncol(x)
  if (p == 0) {
    return(list(coefficients = numeric(), residuals = y,
                fitted.values = 0 * y, weights = w, rank = 0L, df.residual = length(y)))
  }
  if (n == 0) {
    return(list(coefficients = rep(NA_real_, p), residuals = y,
                fitted.values = 0 * y, weights = w, rank = 0L, df.residual = 0L))
  }

  if (is.null(rel) | all(rel==1)) {
    return (lm.wfit(x, y, w, offset = offset, singular.ok = singular.ok,
                    ...))
  }
  if (length(rel) != p)
    stop("Number of reliabilities is different of number of variables")

  #########################
  var_i = apply(x,2,var, na.rm=TRUE)*(n-1)/n
  S = (n*(1-rel)*var_i) * diag(p)
  W = diag(w)
  A = t(x)%*%W%*%x-S
  coef = solve(A)%*%t(x)%*%W%*%y
  z$s.sq = (t(y)%*%W%*%y - t(coef)%*%A%*%coef)
  z$vcov = (z$s.sq[1,1]/(n-p))*solve(A)%*%t(x)%*%W%*%x%*%solve(A)
  z$rank = p
  z$n = n
  z$residuals = y - x%*%coef
  I1 = rep(1,n)
  L = diag(n) - I1%*%solve(crossprod(I1))%*%t(I1)
  z$TSS = t(y)%*%L%*%y


  #########################
  if (!singular.ok && z$rank < p)
    stop("singular fit encountered")
  z$coefficients <- coef
  z$fitted.values <- y - z$residuals
  z$weights <- w
  if (zero.weights) {
    coef[is.na(coef)] <- 0
    f0 <- x0 %*% coef
    if (ny > 1) {
      save.r[ok, ] <- z$residuals
      save.r[nok, ] <- y0 - f0
      save.f[ok, ] <- z$fitted.values
      save.f[nok, ] <- f0
    }
    else {
      save.r[ok] <- z$residuals
      save.r[nok] <- y0 - f0
      save.f[ok] <- z$fitted.values
      save.f[nok] <- f0
    }
    z$residuals <- save.r
    z$fitted.values <- save.f
    z$weights <- save.w
  }
  if (!is.null(offset))
    z$fitted.values <- z$fitted.values + offset
  c(z[c("coefficients", "residuals", "fitted.values", "effects",
        "weights", "rank", "vcov","n","s.sq","TSS")], df.residual = n - z$rank)
}

eivlm.fit <- function (x, y, offset = NULL, tol = 1e-07,
                       singular.ok = TRUE, rel=NULL, ...)
{
  n <- nrow(x)
  w = rep(1,n)
  return (eivlm.wfit(x, y, w, offset = offset, singular.ok = singular.ok,
                     rel=rel,...))

}

summary.eivlm <- function (object, correlation = FALSE, symbolic.cor = FALSE,
                           ...)
{
  z <- object
  p <- z$rank
  rdf <- z$df.residual
  if (p == 0) {
    r <- z$residuals
    n <- length(r)
    w <- z$weights
    rss <- z$s.sq
    resvar <- rss/rdf
    ans <- z[c("call", "terms", if (!is.null(z$weights)) "weights")]
    class(ans) <- "summary.lm"
    ans$aliased <- is.na(coef(object))
    ans$residuals <- r
    ans$df <- c(0L, n, length(ans$aliased))
    ans$coefficients <- matrix(NA, 0L, 4L)
    dimnames(ans$coefficients) <- list(NULL, c("Estimate",
                                               "Std. Error", "t value", "Pr(>|t|)"))
    ans$sigma <- sqrt(resvar)
    ans$r.squared <- ans$adj.r.squared <- 0
    return(ans)
  }
  if (is.null(z$terms))
    stop("invalid 'lm' object:  no 'terms' component")
  if (!inherits(object, "lm"))
    warning("calling summary.lm(<fake-lm-object>) ...")
  r <- z$residuals
  f <- z$fitted.values
  w <- z$weights
  mss <- z$TSS - z$s.sq
  rss <- z$s.sq
  if (!is.null(w)) {
    r <- sqrt(w) * r
  }
  resvar <- rss/rdf
  if (is.finite(resvar) && resvar < (mean(f)^2 + var(f)) *
      1e-30)
    warning("essentially perfect fit: summary may be unreliable")
  R <- z$vcov
  se <- sqrt(diag(z$vcov))
  est <- z$coefficients
  tval <- est/se
  ans <- z[c("call", "terms", if (!is.null(z$weights)) "weights")]
  ans$residuals <- r
  ans$coefficients <- cbind(est, se, tval, 2 * pt(abs(tval),
                                                  rdf, lower.tail = FALSE))
  dimnames(ans$coefficients) <- list(row.names(z$coefficients),
                                     c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
  ans$aliased <- is.na(coef(object))
  ans$sigma <- sqrt(resvar)[1,1]
  ans$df <- c(p, rdf, p)
  if (p != attr(z$terms, "intercept")) {
    df.int <- if (attr(z$terms, "intercept"))
      1L
    else 0L
    ans$r.squared <- mss/(mss + rss)
    ans$adj.r.squared <- 1 - (1 - ans$r.squared) * ((z$n -
                                                       df.int)/rdf)
    ans$fstatistic <- c(value = (mss/(p - df.int))/resvar,
                        numdf = p - df.int, dendf = rdf)
  }
  else ans$r.squared <- ans$adj.r.squared <- 0
  ans$cov.unscaled <- R
  dimnames(ans$cov.unscaled) <- dimnames(ans$coefficients)[c(1,
                                                             1)]
  if (correlation) {
    ans$correlation <- (R * resvar)/outer(se, se)
    dimnames(ans$correlation) <- dimnames(ans$cov.unscaled)
    ans$symbolic.cor <- symbolic.cor
  }
  if (!is.null(z$na.action))
    ans$na.action <- z$na.action
  class(ans) <- "summary.lm"
  class(ans) <- c("summary.eivlm", class(ans))
  ans
}

vcov.eivlm <- function (object, ...)
{
  so <- summary.eivlm(object, corr = FALSE)
  so$sigma^2 * so$cov.unscaled
}
