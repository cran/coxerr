\name{coxerr}
\alias{coxerr}
\title{Cox regression with dependent error in covariates}
\description{Estimation methods of Huang and Wang (2018)}
\usage{coxerr(t,dlt,wuz,method,initbt=rep(0,dim(as.matrix(wuz))[2]-1),
       derr=1e-6)}
\arguments{
  \item{t}{follow-up time.}
  \item{dlt}{censoring indicator: 1 - event, 0 - censored.}
  \item{wuz}{covariate-related variables: wuz[,1] -  mismeasured,
             wuz[,2] - instrumental variable (IV),
             wuz[,-c(1,2)] - accurately measured.}
  \item{method}{estimation method: 1 - Prop1, 2 - Prop 2.}
  \item{initbt}{initial value for the estimate.}
  \item{derr}{error tolerance.}
}
\value{
  \item{bt}{point estimate.}
  \item{va}{estimated variance-covariance matrix.}
  \item{succ}{indicator for estimate-finding success.}
}
\examples{
## simulate a dataset following Scenario 1 of Table 1 in Huang and Wang (2018)
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

## estimation using PROP1
fit1 <- coxerr(t, dlt, wuz, 1)
## estimation using PROP2
fit2 <- coxerr(t, dlt, wuz, 2)
}

\references{Huang, Y. and Wang, C. Y. (2018) Cox Regression with dependent
error in covariates, \emph{Biometrics} 74, 118--126.}
\author{Yijian Huang}
