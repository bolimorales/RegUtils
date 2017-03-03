# RegUtils

This R package adds regression models common to Stata users:

- alm: Linear regression with a large dummy-variable set
- boxcox_r: finds the maximum likelihood estimates of the parameters of the Box–Cox transform
- eivlm: Errors-in-variables regression
- etreg: fits an average treatment effect (ATE) linear regression model augmented with an endogenous treatment variable.
- ivtobit: fits a tobit model where one or more of the regressors is endogenous.

## Examples:
### alm
```R
#Using Chile dataframe from car package, to absorb categorical region variable:
data(Chile, package="car")
fit1 = alm(formula = income ~ education + age + statusquo + region,
absorb="region", data = Chile)
summary(fit1)
```

### boxcox_r
```R
fit1 = boxcox.r(Volume ~ Height + Girth, data = trees)
summary(fit1)
```

### eivlm
```R
#Assuming speed information was measured with a reliability of 0.85
fit1 = eivlm(dist~speed, data=cars, rel = c(0.85))
summary(fit1)
```

### etreg
```R
data(Mroz, package="car")
Mroz$wcb = as.integer(Mroz$wc)
fit1 = etreg(inc~age+wcb+k618+k5+lfp,wcb~hc+age,data=Mroz)
summary(fit1, robust=TRUE)
```

### ivtobit
```R
```
