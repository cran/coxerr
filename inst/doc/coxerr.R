## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----install, eval=FALSE, message=FALSE, warning=FALSE------------------------
#  install.packages("coxerr")

## ----simulation, eval=TRUE, message=FALSE, warning=FALSE----------------------
size <- 300
bt0 <- 1
## true covariate
x <- rnorm(size)
## survival time, censoring time, follow-up time, censoring indicator
s <- rexp(size) * exp(-bt0 * x)
c <- runif(size) * ifelse(x <= 0, 4.3, 8.6)
t <- pmin(s, c)
dlt <- as.numeric(s <= c)
## mismeasured covariate with heterogeneous error, IV
w <- x + rnorm(size) * sqrt(pnorm(x) * 2) * 0.5 + 1
u <- x * 0.8 + rnorm(size) * 0.6
wuz <- cbind(w, u)

## ----coxerr, eval=TRUE, message=FALSE, warning=FALSE--------------------------
library(coxerr)
## estimation using PROP1
fit1 <- coxerr(t, dlt, wuz, 1)
fit1
## estimation using PROP2
fit2 <- coxerr(t, dlt, wuz, 2)
fit2

