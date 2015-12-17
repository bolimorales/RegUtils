\name{eivlm}
\alias{eivlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
eivlm (Errors in variables regression model)
}
\description{
\code{etreg} fits a regression models when one or more of the independent variables are
measured with additive noise.
}
\usage{
eivlm <- function(formula, data, subset, weights, na.action,
                  model = TRUE, x = FALSE, method = NULL, y = FALSE,
                  singular.ok = TRUE, contrasts = NULL, rel = NULL, offset, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{a symbolic description for the model to be estimated,}
  \item{data}{a \code{data.frame},}
  \item{subset}{see \code{\link{lm}} for \code{"alm"}, a character or
    numeric vector indicaing asubset of the table of coefficient to be
    printed for \code{"print.summary.alm"},}
  \item{na.action}{
a function which indicates what should happen when the data contain NAs. The default is set by the na.action setting of options, and is na.fail if that is unset. The ‘factory-fresh’ default is na.omit. Another possible value is NULL, no action. Value na.exclude can be useful.
}
  \item{weights}{an optional vector of weights to be used in the fitting process.}
  \item{model, x, y}{logicals.  If \code{TRUE} the corresponding components of
    the fit (the model frame, the model matrices , the response) are returned.}
  \item{\dots}{further arguments.}
}
\details{
Standard regression would underestimate the effect of the variable if it is measure with aditive noise. The regression can be adjusted for the biases if the reliability is known:

reliability = 1 - (noise variance)/(total variance)
}
\value{
  object of class 'eivlm' which inherits from class 'lm'.
  Components are identical to those of class 'lm',
  see \code{\link{lm}}.
}
\references{
Draper, N., and H. Smith. (1998). Applied Regression Analysis. 3rd ed. New York. Wiley.
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See Also \code{\link{lm}}
}
\examples{
#Assuming speed information was measured with a reliability of 0.85
fit1 = eivlm(dist~speed, data=cars, rel = c(0.85))
summary(fit1)
}
\keyword{ eivlm }
\keyword{ errors in variables }% __ONLY ONE__ keyword per line
