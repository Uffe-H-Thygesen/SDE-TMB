## Simulation and reestimation in a linear stochastic differential equation
## using TMB and the Laplace approximation
##
## Uffe Høgsbro Thygesen, uhth@dtu.dk, 2024

set.seed(123456)

## The model: dX = f(X)*dt + g(X)*dB
f <- function(x,lambda=lambda,mu=mu) lambda*(mu-x)
g <- function(x,sigma=sigma) sigma*rep(1,length(x))

## Simulation control
dt <- 0.1
t <- seq(0,1000,dt)
x0 <- 0

## System parameters
lambda <- 1
sigma <- 1
mu <- 2
s <- 0.5

## Simulate using Euler-Maruyama
X <- numeric(length(t))
X[1] <- x0
for(i in 2:length(t))
    X[i] <- rnorm(1,X[i-1] + f(X[i-1],lambda,mu)*dt,g(X[i-1],sigma)*sqrt(dt))

## Generate observations at regular sample points of time 
n <- length(X)
tsample <- 1
iobs <- seq(1,n,tsample/dt)
Y <- X[iobs] + rnorm(length(iobs),sd=s)

## Core function: The joint negative log density of all random variables
nloglik <- function(p)
{
    dX <- diff(p$X)
    Xi <- head(p$X,-1)
    
    ## Use the Euler approximation for the SDE
    nll <- - sum(dnorm(dX,
                       f(Xi,p$lambda,p$mu)*dt,
                       g(Xi,exp(p$logsigma))*sqrt(dt),
                       log=TRUE))

    ## Contribution from observations
    nll <- nll - sum(dnorm(Y,p$X[iobs],exp(p$logs),log=TRUE))

    return(nll)
}
    
## Initial guess on all parameters
p0 <- list(X=numeric(n),lambda=0,mu=0,logsigma=0,logs=log(s))

nloglik(p0)

## Construct likelihood function by integrating out X. Fix measurement noise logs
require(RTMB)
obj <- MakeADFun(nloglik,p0,random=c("X"),map=list(logs=factor(NA)))
    
## Fit parameters
comp.time <- system.time(fit <- nlminb(obj$par, obj$fn, obj$gr))

## Obtain random effects, variances, etc.
rep <- sdreport(obj)

## Extract states: Estimates and marginal posterior variances, confidence limits
Xhat <- rep$par.random
VX <- rep$diag.cov.random
Xu <- Xhat+sqrt(VX)
Xl <- Xhat-sqrt(VX)

rangeX <- range(c(X,Y,Xu,Xl))

## Plot confidence regions, "true" states, observations, and estimates
## Plot only the first part of the time series
tlim <- c(0,10)

pdf("../figures/sdelinear.pdf",width=8)
plot(tlim,rangeX,type="n",xlab="t",ylab="x")
polygon(c(t,rev(t)),c(Xu,rev(Xl)),col="lightgray",border=NA)
lines(t,X)
points(t[iobs],Y,pch=16)
lines(t,Xhat,lwd=2,lty="dashed")
legend("topleft",legend=c("True signal","Data","Estimated signal"),
       lty=c("solid",NA,"dashed"),lwd=c(1,NA,2),pch=c(NA,16,NA))
dev.off()

print(comp.time)

print(fit)
