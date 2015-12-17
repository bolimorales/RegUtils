\name{ivtobit}
\alias{ivtobit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ivtobit (tobit model with instrumental variables)
}
\description{
\code{ivtobit} fits a tobit model where one or more of the regressors is endogenous.
}
\usage{
ivtobit <- function(formula, instruments, data, subset, na.action = NULL,
weights, offset, contrasts = NULL, model = TRUE, y = TRUE, x = FALSE,
left=0, right = Inf, method = "BHHH", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{a symbolic description for the model to be estimated,}
  \item{instruments}{legacy option: a symbolic description for instrumental variables,}
  \item{data}{a \code{data.frame},}
  \item{subset}{see \code{\link{lm}} for \code{"alm"}, a character or
    numeric vector indicaing asubset of the table of coefficient to be
    printed for \code{"print.summary.alm"},}
  \item{na.action}{
a function which indicates what should happen when the data contain NAs. The default is set by the na.action setting of options, and is na.fail if that is unset. The ‘factory-fresh’ default is na.omit. Another possible value is NULL, no action. Value na.exclude can be useful.
}
 \item{weights}{an optional vector of weights to be used in the fitting process.}
  \item{method}{"twostep" for two step estimator or maximization method (from maxLik package), currently either "NR" (for Newton-Raphson), "BFGS" (for Broyden-Fletcher-Goldfarb-Shanno), "BFGSR" (for the BFGS algorithm implemented in R), "BHHH" (for Berndt-Hall-Hall-Hausman), "SANN" (for Simulated ANNealing), "CG" (for Conjugate Gradients), or "NM" (for Nelder-Mead). Lower-case letters (such as "nr" for Newton-Raphson) are allowed. If missing, a suitable method is selected automatically. See \code{\link{maxLik}}}
  \item{\dots}{further arguments.}
}
\details{

}
\value{
  object of class 'ivtobit' which inherits from class 'maxLik' if maximun likelihood     estimation is used.
  Components are identical to those of class 'maxLik' or 'lm' depending on method used,
  see \code{\link{maxLik}} or \code{\link{lm}}.
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

}
\keyword{ ivtobit }
\keyword{ instrumental variables }
\keyword{ tobit }
