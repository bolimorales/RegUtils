\name{etreg}
\alias{etreg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
etreg (Endogenous treatment effect model)
}
\description{
\code{etreg} fits an average treatment effect (ATE) linear regression model augmented with an endogenous treatment variable.
}
\usage{
etreg(formula, treat, data, subset, na.action,
contrasts = NULL, model = TRUE, y = TRUE, x = FALSE,
method = "NR", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{a symbolic description for the model to be estimated,}
  \item{treat}{a symbolic description for the treatment model to be estimated,}
  \item{data}{a \code{data.frame},}
  \item{subset}{see \code{\link{lm}} for \code{"alm"}, a character or
    numeric vector indicaing asubset of the table of coefficient to be
    printed for \code{"print.summary.alm"},}
  \item{na.action}{
a function which indicates what should happen when the data contain NAs. The default is set by the na.action setting of options, and is na.fail if that is unset. The ‘factory-fresh’ default is na.omit. Another possible value is NULL, no action. Value na.exclude can be useful.
}
  \item{method}{maximisation method (from maxLik package), currently either "NR" (for Newton-Raphson), "BFGS" (for Broyden-Fletcher-Goldfarb-Shanno), "BFGSR" (for the BFGS algorithm implemented in R), "BHHH" (for Berndt-Hall-Hall-Hausman), "SANN" (for Simulated ANNealing), "CG" (for Conjugate Gradients), or "NM" (for Nelder-Mead). Lower-case letters (such as "nr" for Newton-Raphson) are allowed. If missing, a suitable method is selected automatically. See \code{\link{maxLik}}}
  \item{\dots}{further arguments.}
}
\details{
The endogenous treatment-regression model is composed of an equation for the
outcome:
yj = xj.Beta + delta.tj + uj
and an equation for the endogenous treatmen:
tj = 1 (if wj.gamma + ej > 0) | 0 (otherwise)
where xj are the covariates used to model the outcome, wj are the covariates used to model treatment assignment, and the error terms uj and ej are bivariate normal with mean zero and covariance matrix:
|sigma^2    rho.sigma |
|rho.sigma  1         |
}
\value{
  object of class 'etrg' which inherits from class 'maxLik'.
  Components are identical to those of class 'maxLik',
  see \code{\link{maxLik}}.
}
\references{
Wooldridge, J. M. (2010). Econometric Analysis of Cross Section and Panel Data. 2nd ed. Cambridge. MIT Press.
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See Also as \code{\link{maxLik}}
}
\examples{
data(Mroz, package="car")
Mroz$wcb = as.integer(Mroz$wc)
fit1 = etreg(inc~age+wcb+k618+k5+lfp,wcb~hc+age,data=Mroz)
summary(fit1, robust=TRUE)
}
\keyword{ etreg }
\keyword{ endogenous treatment }% __ONLY ONE__ keyword per line
