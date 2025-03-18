## Simulation and reestimation in two stochastic differential equations
## using TMB and the Laplace approximation
##
## THe model is a stochastic version of the Rosenzweig-MacArthur model
##
## This version uses the "tiny" approach and both X and dB as latent variables
##
## uhth, 2025

require(SDEtools)
set.seed(123456)

## RmA/Bazykin-model in natural coordinates
fN <- function(n,p,par=p0) with(par,r*n*(1-n/K) - beta*n*p/(1+beta*n/Cmax))
fP <- function(n,p,par=p0) with(par,epsilon*beta*n*p/(1+beta*n/Cmax) - mu*p)

e2 <- 1e-2

gN <- function(n,par=p0) with(par,sN*n) # sqrt(sqrt(e2+n^2)))
gP <- function(p,par=p0) with(par,sP*p) # sqrt(sqrt(e2+p^2)))

## In transformed coordinates X=h(N), Y=h(P)
h <- log
hi <- exp
hd <- function(n) 1/n
hdd <- function(n) -1/n^2

fX <- function(x,y,par=p0)
    fN(hi(x),hi(y),par) * hd(hi(x)) + 0.5*gN(hi(x),par)^2*hdd(hi(x))

fY <- function(x,y,par=p0)
    fP(hi(x),hi(y),par) * hd(hi(y)) + 0.5*gP(hi(y),par)^2*hdd(hi(y))

gX <- function(x,par=p0) with(par,gN(hi(x),par)*hd(hi(x)))
gY <- function(y,par=p0) with(par,gP(hi(y),par)*hd(hi(y)))

## The model for simulation: dX = f(X)*dt + g(X)*dB, X=(log N,logP)
fsim <- function(x,par=p0) c(fX(x[1],x[2],par),fY(x[1],x[2],par))
gsim <- function(x,par=p0) diag(c(gX(x[1],par),gY(x[2],par)))

## Simulation control
dt <- 0.1
T <- 100
t <- seq(0,T,dt)

## System parameters
p0 <- list(
    r = 1,
    K = 1,
    epsilon = 3,
    beta = 3,
    Cmax = 1,
    mu = 1,
    sN = 0.2, #0.1
    sP = 0.1  #0.1
)

x0 <- log(c(0.1,0.1))

## Simulate using Euler-Maruyama
BN <- rBM(t)
BP <- rBM(t)
sim <- euler(function(x)fsim(x,p0),
             function(x)gsim(x,p0),
             t,x0,B=cbind(BN,BP))

Nsim <- hi(sim$X[,1])
Psim <- hi(sim$X[,2])

## Generate observations at these time indeces
n <- length(t)
tsample <- 1 # 5 
iobs <- round(seq(1,n,tsample/dt))
vP <- 1e-4
vN <- 10 # 8 # 10 
YN <- rpois(length(iobs),vN*Nsim[iobs])
YP <- rpois(length(iobs),vP*Psim[iobs])

## Model family for reestimation
fest <- function(x,lambda=lambda,mu=mu) lambda*(mu-x)
gest <- function(x,sigma=sigma) sigma*rep(1,length(x))

## Core function: The joint negative log density of all random variables
nloglik <- function(p)
{
    tiny <- 1e-3
    
    dX <- diff(p$X)
    dY <- diff(p$Y)

    Xi <- head(p$X,-1)
    Yi <- head(p$Y,-1)

    dXpred <- fX(Xi,Yi,p)*dt + gX(Xi,p)*p$dBN
    dYpred <- fY(Xi,Yi,p)*dt + gY(Yi,p)*p$dBP

    ## Use the Euler approximation for the SDE and the "tiny" approach
    nll <- - sum(dnorm(dX,
                       dXpred,
                       tiny,
                       log=TRUE))
    nll <- nll - sum(dnorm(dY,
                       dYpred,
                       tiny,
                       log=TRUE))
    print(nll)

    ## The contribution from the Brownian motion 
    nll <- nll - sum(dnorm(p$dBN,0,sqrt(dt),log=TRUE))
    print(nll)
    nll <- nll - sum(dnorm(p$dBP,0,sqrt(dt),log=TRUE))
    print(nll)
    
    ## Contribution from observations
    nll <- nll - sum(dpois(YN,hi(p$X[iobs])*vN,log=TRUE))
    print(nll)
    nll <- nll - sum(dpois(YP,hi(p$Y[iobs])*vP,log=TRUE))
    print(nll)

    return(nll)
}

## Initial guess on all parameters
pinit <- c(p0,
           list(X=0*sim$X[,1],
                dBN=0*diff(BN),
                Y=0*sim$X[,2],
                dBP=0*diff(BP)))

nloglik(pinit)

## Construct likelihood function by integrating out X
require(RTMB)

map <- list(epsilon=factor(NA),Cmax=factor(NA))

obj <- MakeADFun(nloglik,pinit,random=c("X","Y","dBN","dBP"),
                 map=map)

obj$fn()


## Fit parameters
comp.time <- system.time(fit <- nlminb(obj$par, obj$fn, obj$gr))

## Obtain random effects, variances, etc.
rep <- sdreport(obj)

## Extract states: Estimates and marginal posterior variances, confidence limits
match.name <- function(x,name) x[names(x)==name]
find.name <- function(x,name) names(x)==name

Xhat <- match.name(rep$par.random,"X")
Yhat <- match.name(rep$par.random,"Y")

Nhat <- hi(Xhat)
Phat <- hi(Yhat)

sd.X <- sqrt(abs(rep$diag.cov.random[find.name(rep$par.random,"X")]))
sd.Y <- sqrt(abs(rep$diag.cov.random[find.name(rep$par.random,"Y")]))

Nu <- hi(Xhat + sd.X)
Nl <- hi(Xhat - sd.X)

Pu <- hi(Yhat + sd.Y)
Pl <- hi(Yhat - sd.Y)

pdf(file="../figures/RmA-fit.pdf",width=4.5,height=4.5)
par(mfrow=c(2,1),mar=c(3,3,2,2))
plot(range(t),c(0,max(Nu)),type="n",main="Prey")
polygon(c(t,rev(t)),c(Nu,rev(Nl)),col="grey",border=NA)
lines(t,Nhat)
lines(t,Nsim,col="red")
points(t[iobs],YN/vN)

plot(range(t),c(0,quantile(Pu,0.99)),type="n",main="Predators")
polygon(c(t,rev(t)),c(Pu,rev(Pl)),col="grey",border=NA)
lines(t,Phat)
lines(t,Psim,col="red")
dev.off()

pdf(file="../figures/RmA-sim-states.pdf",width=5,height=5)
par(mfrow=c(2,1),mar=c(3,3,2,2))
plot(range(t),c(0,max(Nu)*vN),type="n",main="Prey")
lines(t,Nsim*vN,col="red")

plot(range(t),c(0,quantile(Pu,0.99)),type="n",main="Predators")
lines(t,Psim,col="red")
dev.off()

pdf(file="../figures/RmA-sim.pdf",width=5,height=5)
par(mfrow=c(2,1),mar=c(3,3,2,2))
plot(range(t),c(0,max(Nu)*vN),type="n",main="Prey")
lines(t,Nsim*vN,col="red")
points(t[iobs],YN)

plot(range(t),c(0,quantile(Pu,0.99)),type="n",main="Predators")
lines(t,Psim,col="red")
dev.off()

ptrue <- p0[!names(p0) %in% names(map)]
est.tab <- cbind(ptrue,rep$par.fixed,sqrt(diag(rep$cov.fixed)))
colnames(est.tab) <- c("True","Estimate","std.dev")
