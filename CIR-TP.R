## Transition probabilities in the CIR model given by the Ito equation 
##
## dX = lambda(xi-X)*dt + gamma*sqrt(x) * dB
##
## or the equivalanet Stratonovich equation
##
## Using several variants  of the Laplace approximation in TMB

###########################################################
## Model specifiation
###########################################################

## Parameters
lambda <- 1
xi <- 1
gamma <- 0.5

## Analytical transition probability
pA <- function(y,x,T) dCIR(y,x,lambda,xi,gamma,T,Stratonovich=FALSE)

## Noise term
g <- function(x) gamma*sqrt(x)
gp <- function(x) gamma/2/sqrt(x)
gpp <- function(x) -gamma/4*x^(-1.5)
gppp <- function(x) gamma/4*1.5*x^(-2.5)

## Ito and Stratonovich drift
fI <- function(x) lambda*(xi-x)
fIp <- function(x) 0*x-lambda
fIpp <- function(x) 0*x

fS <- function(x) fI(x) - 0.5*g(x)*gp(x)
fSp <- function(x) fIp(x) - 0.5*gp(x)^2 - 0.5*g(x)*gpp(x)
fSpp <- function(x) fIpp(x) - 1.5*gp(x)*gpp(x) - 0.5*g(x)*gppp(x) 
    

#############################################################

graphics.off()

require(SDEtools)

myapprox <- function(xi,yi,xo) approx(xi,yi,xo,yleft=head(yi,1),yright=tail(yi,1))$y


x0 <- 0.5
T <- 1

dt <- 2^-10
tv <- seq(0,T,dt)
nt <- length(tv)

xTs <- seq(0.1,2.5,0.1)

#############################################################
## TMB: Stratonovich formulation

## Core function for TMB: The joint negative log density of all random variables
nloglik.S <- function(p)
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

## Arbitrary initial guess on all parameters. First the bridge
Xinit <- rep(1,length(tv)-2)

p0.S <- list(X=Xinit,x=1,y=1,dt=diff(tv)[1])

## Construct likelihood function by integrating out X
require(RTMB)
obj.S <- MakeADFun(nloglik.S,p0.S,random=c("X"))

dSDE.S <- function(y,x,T,return.list = FALSE,log=FALSE)
{
    ## Initial guess on all parameters
    t <- seq(0,T,length=length(tv))
    Xinit <- approx(c(0,T),c(x,y),t)$y[-c(1,length(t))]

    dt <- diff(t)[1]
    
    p0 <- list(X=Xinit,x=x,y=y,dt=dt)

    obj.S$env$last.par <- unlist(p0)

    p <- c(x=x,y=y,dt=dt)
    logpdf <- -obj.S$fn(p)
    
    rep <- sdreport(obj.S,par.fixed=p,skip.delta.method=TRUE)
    X <- rep$par.random

    xX <- c(x,X)
    Xy <- c(X,y)

    dX <- diff(c(x,X,y))
    dB <- 2*(dX - 0.5*(fS(xX) + fS(Xy))*dt)/(g(xX) + g(Xy))

    logpdf <- logpdf - sum(log(0.5*(g(xX) + g(Xy))))
    logpdf <- logpdf + sum(log(1 -  0.5*fSp(Xy)*dt - 0.5*gp(Xy)*dB))

    if(return.list)
    {
        l <- list(logpdf = logpdf,
                  pdf = exp(logpdf),
                  X = X
                  )
        return(l)
    }
        
    if(log) return(logpdf) else return(exp(logpdf))
}

#####################################################################################
## TMB: Ito formulation, dB as only root variables
nloglik.dB <- function(p)
{
    ## Generate X using Euler-Maruyama
    X <- numeric(length(tv))
    X[1] <- p$x

    for(i in 2:length(tv))
        X[i] <- X[i-1] + fI(X[i-1])*p$dt + g(X[i-1])*p$dB[i-1]

    nll <- - sum(dnorm(p$dB,
                       0,
                       sqrt(p$dt),
                       log=TRUE))

    nll <- nll - dnorm(X[length(tv)],p$y,p$epsilon,log=TRUE)
    return(nll)
}

## Arbitrary initial guess on all parameters. First the bridge
dB <- numeric(length(tv)-1)

p0.dB <- list(dB=dB,x=1,y=1,dt=diff(tv)[1],epsilon=1e-4)

## Construct likelihood function by integrating out dB
obj.dB <- MakeADFun(nloglik.dB,p0.dB,random=c("dB"))

dSDE.dB <- function(y,x,T,return.list = FALSE,log=FALSE)
{
    ## Initial guess on all parameters
    ## t <- seq(0,T,length(tv))
    dt <- T/(length(tv)-1)
    
    ## Xinit <- approx(c(0,T),c(x,y),t)$y
    ## xX <- head(Xinit,-1)
    ## dB <- (diff(Xinit) - fI(xX)*dt ) / g(xX)

    epsilon <- 1e-4

    ## p0 <- list(dB=dB,x=x,y=y,dt=dt,epsilon=epsilon)
    ## obj.dB$env$last.par <- unlist(p0)

    p <- c(x=x,y=y,dt,epsilon)
    logpdf <- -obj.dB$fn(p)
    
    if(return.list)
    {
        l <- list(logpdf = logpdf,
                  pdf = exp(logpdf)
                  )
        return(l)
    }
        
    if(log) return(logpdf) else return(exp(logpdf))
}

##############################################################
## TMB: Ito formulation, dB and X as root variables
nloglik.XdB <- function(p)
{
    dX <- diff(p$X)
    xX <- head(p$X,-1)
    XY <- p$X[c(1,length(p$X))]
    
    dXpred <- fI(xX)*p$dt + g(xX)*p$dB
    

    nll <- - sum(dnorm(dX,
                       dXpred,
                       p$epsilon,
                       log=TRUE))

    nll <- nll - sum(dnorm(p$dB,0,sqrt(p$dt),log=TRUE))
    nll <- nll - sum(dnorm(XY,c(p$x,p$y),p$epsilon,log=TRUE))

    return(nll)
}

## Arbitrary initial guess on all parameters. First the bridge
dB <- numeric(length(tv)-1)
X <- numeric(length(tv))+1

p0.XdB <- list(X=X,dB=dB,x=1,y=1,dt=diff(tv)[1],epsilon=1e-4)

## Construct likelihood function by integrating out dB
obj.XdB <- MakeADFun(nloglik.XdB,p0.XdB,random=c("X","dB"))

dSDE.XdB <- function(y,x,T,return.list = FALSE,log=FALSE)
{
    ## Initial guess on all parameters
    t <- seq(0,T,length(tv))
    dt <- T/(length(tv)-1)
    
    X <- approx(c(0,T),c(x,y),t)$y
    xX <- head(X,-1)
    dB <- (diff(X) - fI(xX)*dt ) / g(xX)

    epsilon <- 1e-4
    p0 <- list(X=X,dB=dB,x=x,y=y,dt=dt,epsilon=epsilon)
    obj.XdB$env$last.par <- unlist(p0)

    p <- c(x=x,y=y,dt,epsilon)
    logpdf <- -obj.XdB$fn(p)
    
    if(return.list)
    {
        l <- list(logpdf = logpdf,
                  pdf = exp(logpdf)
                  )
        return(l)
    }
        
    if(log) return(logpdf) else return(exp(logpdf))
}


##############################################################
## TMB: Ito formulation, X as root variables
nloglik.X <- function(p)
{
    xX <- c(p$x,p$X)
    xXy <- c(p$x,p$X,p$y)

    dX <- diff(xXy)
    dB <- (dX - fI(xX)*p$dt)/g(xX)
    
    nll <- - sum(dnorm(dB,0,sqrt(p$dt),log=TRUE))

    return(nll)
}

## Arbitrary initial guess on all parameters. 
X <- rep(1,length(tv)-2)

p0.X <- list(X=X,x=1,y=1,dt=diff(tv)[1])

## Construct likelihood function by integrating out dB
obj.X <- MakeADFun(nloglik.X,p0.X,random=c("X"))

dSDE.X <- function(y,x,T,return.list = FALSE,log=FALSE)
{
    ## Initial guess on all parameters
    t <- seq(0,T,length(tv))
    dt <- T/(length(tv)-1)
    
    xXy <- approx(c(0,T),c(x,y),t)$y
    xX <- head(X,-1)

    p0 <- list(X=tail(xX,-1),x=x,y=y,dt=dt)
    obj.X$env$last.par <- unlist(p0)

    p <- c(x=x,y=y,dt)
    logpdf <- -obj.X$fn(p)

    rep <- sdreport(obj.X,par.fixed=p,skip.delta.method=TRUE)
    X <- rep$par.random

    xX <- c(x,X)
    
    logpdf <- logpdf - sum(log(g(xX)))
    if(return.list)
    {
        l <- list(logpdf = logpdf,
                  pdf = exp(logpdf)
                  )
        return(l)
    }
        
    if(log) return(logpdf) else return(exp(logpdf))
}

####################################################

tAs <- system.time(pAs <- sapply(xTs,function(y)pA(y,x0,T)))
tBs <- system.time(pBs <- sapply(xTs,function(y)dSDE.dB(y,x0,T)))
tXdBs <- system.time(pXdBs <- sapply(xTs,function(y)dSDE.XdB(y,x0,T)))
tXs <- system.time(pXs <- sapply(xTs,function(y)dSDE.X(y,x0,T)))
tSs <- system.time(pSs <- sapply(xTs,function(y)dSDE.S(y,x0,T)))

pdf(file="CIR-TP.pdf")
par(mfrow=c(2,1))

plot(xTs,pAs,xlab="y",ylab="p(0,x,T,y)")
lines(xTs,pBs,lwd=2,col=1)
lines(xTs,pXdBs,lwd=2,col=2)
lines(xTs,pXs,lwd=2,col=3)
lines(xTs,pSs,lwd=2,col=4)

legend("topright",legend=c("True","dB","XdB","X","S"),pty=c(1,NA,NA,NA,NA),lty=c(NA,1,1,1,1),lwd=2,col=c(1,1:4))

plot(xTs,pAs,log="y",xlab="y",ylab="p(0,x,T,y)")
lines(xTs,pBs,lwd=2,col=1)
lines(xTs,pXdBs,lwd=2,col=2)
lines(xTs,pXs,lwd=2,col=3)
lines(xTs,pSs,lwd=2,col=4)

legend("topright",legend=c("True","dB","XdB","X","S"),pty=c(1,NA,NA,NA,NA),lty=c(NA,1,1,1,1),lty=1,lwd=2,col=c(1,1:4))

dev.off()

pdf(file="CIR-TP-ERR.pdf")

par(mfrow=c(2,1))

ps <- cbind(pBs,pXdBs,pXs,pSs)
errs <- ps - pAs
errr <- range(errs)

plot(xTs,pBs-pAs,lwd=2,type="l",ylim=errr,xlab="y",ylab=expression("Error p-p"[A]))

lines(xTs,pXdBs-pAs,lwd=2,col=2,lty="dashed")
lines(xTs,pXs-pAs,lwd=2,col=3)
lines(xTs,pSs-pAs,lwd=2,col=4)

legend("topright",legend=c("dB","XdB","X","S"),lty=c(2,2,1,1),lwd=2,col=1:4)

abline(a=0,b=0,lty="dashed")

rerr <- errs / pAs
rerrr <- range(rerr)

plot(xTs,pBs/pAs-1,lwd=2,type="l",ylim=rerrr,xlab="y",ylab=expression("Relative Error p/p"[A]*"-1"))
lines(xTs,pXdBs/pAs-1,lwd=2,col=2,lty="dashed")
lines(xTs,pXs/pAs-1,lwd=2,col=3)
lines(xTs,pSs/pAs-1,lwd=2,col=4)

abline(a=0,b=0,lty="dashed")

legend("topright",legend=c("dB","XdB","X","S"),lty=c(2,2,1,1),lwd=2,col=1:4)

print(cbind(tAs,tBs,tXdBs,tXs,tSs))
dev.off()
