\name{alm}
\alias{alm}
\title{ALM}
\usage{
alm(formula, data, subset, weights, absorb = NULL, na.action,
   method = "qr", model = TRUE, x = FALSE, y = FALSE, qr = TRUE,
   singular.ok = TRUE, contrasts = NULL, offset, ...)
}
\arguments{
  \item{formula}{a symbolic description for the model to be estimated,}
  \item{data}{a \code{data.frame},}
  \item{subset}{see \code{\link{lm}} for \code{"alm"}, a character or
    numeric vector indicaing asubset of the table of coefficient to be
    printed for \code{"print.summary.alm"},}
  \item{absorb}{a variable encode as factors.}
  \item{na.action}{see \code{\link{lm}},}
  \item{weights}{an optional vector of weights to be used in the fitting process.}
  \item{offset}{an optional offset that can be used to specify an a priori known
    component to be included during fitting.}
  \item{contrasts}{an optional list. See the \code{contrasts.arg} of
    \code{\link[stats:model.matrix]{model.matrix.default}}.}
  \item{model, x, y}{logicals.  If \code{TRUE} the corresponding components of
    the fit (the model frame, the model matrices , the response) are returned.}
  \item{\dots}{further arguments.}
}
\description{
\code{alm} performs linear regression by absorbing one categorical variable.
Model adjustment is reported with Wald Test.
}
\details{
\code{alm} fits a linear model, absorbing a set of k mutually exclusive and exhaustive binary variables, based on Frisch-Waugh-Lovell Theorem.
The intercept reported by \code{alm}, is calculated by choosing the intercept that
makes the prediction calculated at the means of the independent variables equal to the mean of the dependent variable.
}
\examples{
# Without absorb it performs like lm
## Annette Dobson (1990) "An Introduction to Generalized Linear Models".
## Page 9: Plant Weight Data.
ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
weight <- c(ctl, trt)
lm.D9 <- alm(weight ~ group)
lm.D90 <- lm(weight ~ group)
summary(lm.D9)
summary(lm.D90)

#Using Chile dataframe from car package, to absorb categorical region variable:
data(Chile, package="car")
fit1 = alm(formula = income ~ education + age + statusquo + region,
absorb="region", data = Chile)
summary(fit1)
}
