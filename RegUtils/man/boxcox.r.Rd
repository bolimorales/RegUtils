\name{boxcox.r}
\alias{boxcox.r}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
boxcox.r Fit Box-Cox regressions
}
\description{
boxcox.r finds the maximum likelihood estimates of the parameters of the Box–Cox transform and the coefficients on the independent variables.It can fit theta models, lambda models, and models where all but some group of independent variables are not transformed.
}
\usage{
boxcox.r(formula, data, subset, na.action, x = FALSE, y = FALSE,
noTrans = NULL, optimize.bounds = c(-2, 2), model = "theta", test.params = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under ‘Details’.
}
  \item{data}{
an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which boxcox.r is called.
}
  \item{subset}{
an optional vector specifying a subset of observations to be used in the fitting process.
}
  \item{na.action}{
a function which indicates what should happen when the data contain NAs. The default is set by the na.action setting of options, and is na.fail if that is unset. The ‘factory-fresh’ default is na.omit. Another possible value is NULL, no action. Value na.exclude can be useful.
}
  \item{noTrans}{
a set of independet variables which are not going to be transformed.
}
  \item{optimize.bounds}{
bounds for lambda and theta optimization
}
  \item{model}{
one of the following models: theta, lambda, lhs, rhs
}
  \item{test.params}{
a pair of c(theta, lambda) parameters. If test.params are passed, no optimization is performed.
}
  \item{\dots}{
additional arguments to be passed to the low level regression fitting functions (see below).
}
}
\details{
The following models could be fitted:

Theta: Different parameters for box-cox transformation are used for independent and dependent variables.

Lambda: the same parameter is used for box-cox transformation of independent and dependent variables.

lhs: only left side of 'formula' is transformed.

rhs: only right side of 'formula' is transformed.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{coefficients }{a named vector of coefficients}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Box, G. E. P., and D. R. Cox. 1964. An analysis of transformations. Journal of the Royal Statistical Society, Series B 26: 211 - 252.
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
fit1 = boxcox.r(Volume ~ Height + Girth, data = trees)
summary(fit1)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
}
\keyword{ boxcox }
\keyword{ likelihood }% __ONLY ONE__ keyword per line
