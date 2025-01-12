## Transition probabilities in the CIR model given by the Stratonovich equation 
##
## dX = lambda(xi-X)*dt + gamma*sqrt(x) o dB
##
## Comparison of the true transition probability, the continuous Laplace
## approximation given by ODE's, and the time-discrete Laplace approximation
## for different time steps.

graphics.off()

require(latex2exp)
require(deSolve)
require(SDEtools)
require(RTMB)

###########################################################
## Model specifiation
###########################################################

## Parameters
lambda <- 1
xi <- 1
gamma <- 0.5

## Analytical transition probability
pA <- function(y,x,T) dCIR(y,x,lambda,xi,gamma,T,Stratonovich=TRUE)

## Noise term
g <- function(x) gamma*sqrt(x)
gp <- function(x) gamma/2/sqrt(x)
gpp <- function(x) -gamma/4*x^(-1.5)

## Ito and Stratonovich drift
fS <- function(x) lambda*(xi-x)
fSp <- function(x) 0*x-lambda
fSpp <- function(x) 0*x

fI <- function(x) fS(x) + 0.5*g(x)*gp(x)
fIp <- function(x) fSp(x) + 0.5*gp(x)^2 + 0.5*g(x)*gpp(x)

## Hamiltonian
H <- function(x,l) l*fS(x) - 0.5*g(x)^2*l^2

dHdx <- function(x,l) l*fSp(x) - g(x)*gp(x)*l^2
dHdl <- function(x,l) fS(x) - g(x)^2*l

###############################################################
## Basic Strang for CIR
###############################################################

dB.Strang <- function(x,y,h)
{
    ## Forward drift to from t=0 to h/2
    X1 <- xi + (x-xi)*exp(-lambda*h/2)
    
    ## Backward drift to from t=h to h/2
    X2 <- xi + (y-xi)*exp(lambda*h/2)

    ## Compute increment of Brownian motion 
    dB <- (sqrt(X2)-sqrt(X1))*2/gamma

    return(list(X1=X1,X2=X2,dB=dB))
}


################################################################
## The following is generic, i.e., should apply to any model f,g
################################################################

myapprox <- function(xi,yi,xo) approx(xi,yi,xo,yleft=head(yi,1),yright=tail(yi,1))$y


## Hamiltonian system for the most probable path and the integrated effort
## xli is: State x, co-state l, integral i of 0.5*|u|^2
dHam <- function(t,xli,p)
{
    x <- xli[1] 
    l <- xli[2]
    i <- xli[3]
    
    dx <- dHdl(x,l)
    dl <- -dHdx(x,l)
    di <- 0.5*g(x)^2*l^2
    
    return(list(c(dx,dl,di)))
}

## Find transition probabilities from this x0 using the following l0,
## which then determines the end point

x0 <- 0.75
l0 <- -2.106
T <- 1

## Simulation control
dt <- 1e-3
tv <- seq(0,T,dt)
nt <- length(tv)

## Find the end point corresponding to the IC by solving the Hamiltonian
solHam <- ode(c(x0,l0,0),tv,dHam,NULL)

tsol <- solHam[,1]
xsol <- solHam[,2]
lsol <- solHam[,3]
isol <- solHam[,4]
xT <- xsol[nt]
iT <- isol[nt]

## The Riccati equation governing the curvature of the value function 
dRiccati <- function(t,QI,p)
{
    Q <- QI[1]
    
    x <- myapprox(tsol,xsol,t)
    l <- myapprox(tsol,lsol,t)
    
    Hxx <- l*fSpp(x) - gp(x)^2*l^2 - g(x)*gpp(x)*l^2
    Hxl <- fSp(x) - 2*g(x)*gp(x)*l
    Hll <- -g(x)^2
    
    dQ <- -(Hxx+2*Hxl*Q+Q^2*Hll)
    dI <-  -0.5*(g(x)^2*Q + l*g(x)*gp(x))
    
    dQI <- c(dQ,dI)
    
    return(list(dQI))
}

solRic <- ode(c(0,0),rev(tsol),dRiccati,NULL)

Qsol <- rev(solRic[,2])

dLyap <- function(t,S,p)
{ 
    x <- myapprox(tsol,xsol,t)
    l <- myapprox(tsol,lsol,t)
    Q <- myapprox(tsol,Qsol,t)
    
    dS <- 2*(fSp(x)-g(x)^2*Q-2*l*g(x)*gp(x))*S + g(x)^2
    return(list(dS))
}

solLyap <- ode(0,tsol,dLyap,NULL)
ST <- solLyap[nt,2]

## Compute transtion probability from ODE's
pLT <- 1/sqrt(2*pi*ST)*exp(-iT-solRic[nt,3])

if(exists("pA")) pAT <- pA(xT,x0,T)

## Number of time steps for TMB and order analysis
Nts <- 2^(1:6)

pTMBs <- pSTRANG <- numeric(length(Nts))

## Centered Euler approximation
nloglik <- function(p)
{
    dX <- diff(c(p$x,p$X,p$y))
    xX <- c(p$x,p$X)
    Xy <- c(p$X,p$y)

    ## Stratonovich
    fM <- 0.5*(fS(xX) + fS(Xy))
    gM <- 0.5*(g(xX) + g(Xy))
    dB <- (dX - fM*p$dt)/gM
    
    nll <- - sum(dnorm(dB,
                       0,
                       sqrt(p$dt),
                       log=TRUE))

    return(nll)
}

for(i in 1:length(Nts))
{
    Nt <- Nts[i]
    tv <- seq(0,T,length=Nt+1)
    dt <- T / Nt
    
    ## Arbitrary initial guess on all parameters. First the state
    Xinit <- myapprox(c(0,T),c(x0,xT),tv)[2:Nt]

    p <- list(X=Xinit,x=x0,y=xT,dt=dt)
    p0 <- c(x=x0,y=xT,dt=dt)
    
    ## Construct likelihood function by integrating out X
    obj <- MakeADFun(nloglik,p,random=c("X"))

    logpdf <- -obj$fn(p0)
    rep <- sdreport(obj,par.fixed=p0)
    X <- rep$par.random
    xX <- c(x0,X)
    Xy <- c(X,xT)

    dX <- diff(c(x0,X,xT))
    dB <- 2*(dX - 0.5*(fS(xX) + fS(Xy))*dt)/(g(xX) + g(Xy))

    gbar <- 0.5*(g(xX)+g(Xy))
    Lambda <- -dB/dt/gbar 

    logpdf <- logpdf - sum(log(gbar))
    logpdf <- logpdf + sum(log(1 -  0.5*fSp(Xy)*dt - 0.5*gp(Xy)*dB))

    pTMBs[i] <- exp(logpdf)
}


########################################################################
## Strang splitting
##
## This uses code that is specific to the CIR model
########################################################################

## Core function for TMB: The joint negative log density of all random variables
nloglik <- function(p)
{
    dX <- diff(c(p$x,p$X,p$y))
    xX <- c(p$x,p$X)
    Xy <- c(p$X,p$y)

    strang <- dB.Strang(xX,Xy,p$dt)
    
    REPORT(strang)
    REPORT(p)
    
    nll <- -sum(dnorm(strang$dB,0,sqrt(p$dt),log=TRUE))
    return(nll)
}

for(i in 1:length(Nts))
{
    Nt <- Nts[i]
    tv <- seq(0,T,length=Nt+1)
    dt <- T / Nt
    
    ## Arbitrary initial guess on all parameters. First the state
    Xinit <- myapprox(c(0,T),c(x0,xT),tv)[2:Nt]

    p <- list(X=Xinit,x=x0,y=xT,dt=dt)
    p0 <- c(x=x0,y=xT,dt=dt)
    
    ## Construct likelihood function by integrating out X
    obj <- MakeADFun(nloglik,p,random=c("X"))

    logpdf <- -obj$fn(p0)
    rep <- obj$report()
    X <- rep$p$X
    xX <- c(x0,X)
    Xy <- c(X,xT)

    dX <- diff(c(x0,X,xT))
    dB <- rep$strang$dB

    logpdf <- logpdf  - sum(log(gamma*exp(-lambda*dt/2)*(sqrt(rep$strang$X1) + gamma/2*dB)))
    
    pSTRANG[i] <- exp(logpdf)
}


pdf(file="../figures/CIR-Order.pdf",width=6)

absErr <- abs(pTMBs - pLT)
absErrSTRANG <- abs(pSTRANG - pLT)

par(mar = c(5,5,4,2))

ErrRange <- range(c(absErr,absErrSTRANG))

plot(T/Nts,absErr,ylim=ErrRange,log="xy",xlab="Time step h",ylab = TeX(r'(Absolute error $|\hat{p}(\cdot,h) - \hat{p}(\cdot,0)|$)'),pch=16)
points(T/Nts,absErrSTRANG,pch=1)

lines(T/range(Nts),rep(abs(pLT-pAT),2),lty="dashed",lwd=2)

rat <- 5
lines(T/tail(Nts,1)*c(1,rat),tail(absErr,1)*c(1,rat)*1.5,lwd=3,lty=1)
lines(T/tail(Nts,1)*c(1,rat),tail(absErrSTRANG,1)*c(1,rat)^2*1.5,lty=2,lwd=3)

legend("bottomright",
       legend=c("Euler-Stratonovich","Strang",TeX(r'($|\hat{p}(\cdot,0) - p(\cdot)|$)'),'Order 1 scaling','Order 2 scaling'),
       lty=c(NA,NA,2,1,2),
       lwd=c(NA,NA,2,3,3),
       pch=c(16,1,NA,NA,NA))


dev.off()
