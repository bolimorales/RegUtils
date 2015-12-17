###   $Id: alm.R 1126 2015-09-29  $
###
###     Linear regression with a large dummy-variable set for R
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

#alm main function
alm <- function (formula, data, subset, weights, na.action, method = "qr",
                 model = TRUE, x = FALSE, y = FALSE, qr = TRUE, singular.ok = TRUE,
                 contrasts = NULL, absorb = NULL, offset, ...)
{
  ret.x <- x
  ret.y <- y
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  formula_a = NULL

  ######### Divide formula in two sections splitting absorbing variables
  if (!is.null(absorb)) {
    if (length(absorb)>1) {
      formula_a <- eval(substitute(update(b,~+a), list(b=formula,a=as.name(absorb[1]))))
      for (abs_v in absorb) {
        formula <- eval(substitute(update(b,~.-a), list(b=formula,a=as.name(abs_v))))
        formula_a <- eval(substitute(update(b,~.+a), list(b=formula_a,a=as.name(abs_v))))
      }
    }
    else {
      formula <- eval(substitute(update(b,~.-a), list(b=formula,a=as.name(absorb))))
      formula_a <- eval(substitute(update(b,~+a), list(b=formula,a=as.name(absorb))))
    }
    formula_a <- eval(substitute(update(b,~.-1), list(b=formula_a)))
  }
  #########

  m <- match(c("formula","data", "subset", "weights", "na.action",
               "offset"), names(mf), 0L)
  mf <- mf[c(1L,m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  if (!is.null(absorb)) {
    datos_y <- eval(mf, parent.frame())
    mf[[2]] <- formula_a
    datos <- eval(mf, parent.frame())
    mf[[2]] <- formula
  }

  mf <- eval(mf, parent.frame())


  if (method == "model.frame")
    return(mf)
  else if (method != "qr")
    warning(gettextf("method = '%s' is not supported. Using 'qr'",
                     method), domain = NA)
  mt <- attr(mf, "terms")
  mta <- attr(datos, "terms")
  if (!is.null(absorb)) {
    # almacenar en index_set las filas sin NA's en los dos conjuntos
    # (variables que se van a absorber y las que no)
    y <- model.response(datos_y, "numeric")
    index_set = names(y)
  } else {
    y <- model.response(mf, "numeric")
  }
  w <- as.vector(model.weights(mf))
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
    XP <- x
    YP <- y
    #########################
    if (!is.null(absorb)) {
      # read categorical absorbed variables matrix
      xt <- model.matrix(mta, datos, contrasts)
      # remove last variable (set is a partition)
      xt <- xt[index_set,1:(dim(xt)[2]-1)]
      x <- x[index_set,]

      #save number of observations
      dN = dim(xt)[1]
      if (is.null(dN)) {
        stop("Absorbed variable must be categorical, with more than one category")
      }

      #find orthogonal projection of regressors
      mPX = diag(dN) - xt %*% solve(crossprod(xt)) %*% t(xt)

      #Project response and regressor matrices
      YP = mPX %*% y
      XP = mPX %*% x
      #########################
    }
    z <- if (is.null(w))
      lm.fit(XP, YP, offset = offset, singular.ok = singular.ok,
             ...)
    else lm.wfit(XP, YP, w, offset = offset, singular.ok = singular.ok,
                 ...)
  }
  class(z) <- c(if (is.matrix(y)) "mlm", "lm")
  z$na.action <- attr(mf, "na.action")
  z$offset <- offset
  z$contrasts <- attr(x, "contrasts")
  z$xlevels <- .getXlevels(mt, mf)
  z$call <- cl
  z$terms <- mt
  ##############
  ## Correct degrees of fredom by absorbed variables
  if (!is.null(absorb)) {
    z$df.residual <- z$df.residual - dim(xt)[2]
    z$residuals <- YP - as.vector(XP %*% z$coefficients)
    z$fitted.values <- y - z$residuals
    ii = rep(1,dim(x)[1])
    ii = matrix(ii/sum(ii),nrow = 1)
    z$coefficients[1] = mean(y) - (ii%*%x[,-1])%*%z$coefficients[-1]
  }
  ##################

  if (model)
    z$model <- mf
  if (ret.x)
    z$x <- x
  if (ret.y)
    z$y <- y
  if (!qr)
    z$qr <- NULL
  class(z) <- c("alm", class(z))
  return (z)
}

summary.alm <- function (object, correlation = FALSE, symbolic.cor = FALSE,
                         ...)
{
  ## call lm method
  ans <- NextMethod()
  wt = waldTest(b = coef(object), Sigma = vcov(object), Terms = 2:length(coef(object)))
  ans$wald.test = wt
  class(ans) <- c("summary.alm", class(ans))
  ans
}

## based on AER ivreg
print.summary.alm <- function(x, digits = max(3, getOption("digits") - 3),
                                 signif.stars = getOption("show.signif.stars"), ...)
{
  cat("\nCall:\n")
  cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  cat(if(!is.null(x$weights) && diff(range(x$weights))) "Weighted ", "Residuals:\n", sep = "")
  if(NROW(x$residuals) > 5L) {
    nam <- c("Min", "1Q", "Median", "3Q", "Max")
    rq <- if(length(dim(x$residuals)) == 2)
      structure(apply(t(x$residuals), 1, quantile), dimnames = list(nam, dimnames(x$residuals)[[2]]))
    else structure(quantile(x$residuals), names = nam)
    print(rq, digits = digits, ...)
  } else {
    print(x$residuals, digits = digits, ...)
  }

  cat("\nCoefficients:\n")
  printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars,
               signif.legend = signif.stars & is.null(x$diagnostics), na.print = "NA", ...)

  if(!is.null(x$diagnostics)) {
    cat("\nDiagnostic tests:\n")
    printCoefmat(x$diagnostics, cs.ind = 1L:2L, tst.ind = 3L,
                 has.Pvalue = TRUE, P.values = TRUE, digits = digits,
                 signif.stars = signif.stars, na.print = "NA", ...)
  }

  cat("\nResidual standard error:", format(signif(x$sigma, digits)),
      "on", x$df[2L], "degrees of freedom\n")

  cat("Multiple R-Squared:", formatC(x$r.squared, digits = digits))
  v = x$wald.test[["result"]][["chi2"]]

  cat(",\tAdjusted R-squared:", formatC(x$adj.r.squared, digits = digits),
      "\nWald test:", formatC(v["chi2"], digits = digits),
      "on", v["df"],
      "DF,  p-value:", format.pval(v["P"], digits = digits), "\n\n")

  invisible(x)
}

#taken from aod package, performs Wald-Test
waldTest = function (Sigma, b, Terms = NULL, L = NULL, H0 = NULL, df = NULL,
                     verbose = FALSE)
{
  if (is.null(Terms) & is.null(L))
    stop("One of the arguments Terms or L must be used.")
  if (!is.null(Terms) & !is.null(L))
    stop("Only one of the arguments Terms or L must be used.")
  if (is.null(Terms)) {
    w <- nrow(L)
    Terms <- seq(length(b))[colSums(L) > 0]
  }
  else w <- length(Terms)
  if (is.null(H0))
    H0 <- rep(0, w)
  if (w != length(H0))
    stop("Vectors of tested coefficients and of null hypothesis have different lengths\n")
  if (is.null(L)) {
    L <- matrix(rep(0, length(b) * w), ncol = length(b))
    for (i in 1:w) L[i, Terms[i]] <- 1
  }
  dimnames(L) <- list(paste("L", as.character(seq(NROW(L))),
                            sep = ""), names(b))
  f <- L %*% b
  V <- Sigma
  mat <- qr.solve(L %*% V %*% t(L))
  stat <- t(f - H0) %*% mat %*% (f - H0)
  p <- 1 - pchisq(stat, df = w)
  if (is.null(df))
    res <- list(chi2 = c(chi2 = stat, df = w, P = p))
  else {
    fstat <- stat/nrow(L)
    df1 <- nrow(L)
    df2 <- df
    res <- list(chi2 = c(chi2 = stat, df = w, P = p), Ftest = c(Fstat = fstat,
                                                                df1 = df1, df2 = df2, P = 1 - pf(fstat, df1, df2)))
  }
  return (list(Sigma = Sigma, b = b, Terms = Terms, H0 = H0,
               L = L, result = res, verbose = verbose, df = df))
}
