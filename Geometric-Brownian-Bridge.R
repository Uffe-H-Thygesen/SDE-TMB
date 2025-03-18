## Examine the geometric Brownian bridge:
## Find the most probable path between two points
## (In the sense of the mode of the finite-dimensional distributions)
## The take-home message is that this mode is useless!
##
## uhth, 2025

require(RTMB)

## End point of the bridge

x0 <- 1
xT <- 1
T <- 1

## The governing equation:
## Stratonovich SDE dX = rX dt + sigma o dB with
r <- 1
sigma <- 1

## For some odd reason, RTMB requires a reimplementation of the densities in lnorm
mydlnorm <- function(x,logm,logs,log=FALSE)
{
    d <- dnorm(log(x),logm,logs,log=TRUE) - log(x)
    if(log) return(d) else return(exp(d))
}


## Joint p.d.f. of a number of intermediate points
p.d.f <- function(xv,do.X=TRUE,do.plot=FALSE)
{
    xe <- c(x0,xv,xT)
    xold <- c(x0,xv)
    xnew <- c(xv,xT)
    
    ## Report either the joint pdf of the states X (if TRUE), or the
    ## joint pdf of the dB's
    if(do.X){
        ans <- - sum(mydlnorm(xe[-1],
                              log(xe[-(n+1)])+r*diff(te),
                              sigma*sqrt(diff(te)),log=TRUE))
        } else {
            ans <- -sum(dnorm(log(xe[-(n+1)]),log(xe[-1])+r*diff(te),sd=sigma*sqrt(diff(te)),log=TRUE))
        }
    

    if(do.plot) points(te,xe)
    return(ans)
}

## Determine the bridge for different number of sample points
ns <- 2^(1:5)

fits <- list()
reps <- list()
ress <- list()

for(i in 1:length(ns))
{
    ## Define a skeleton
    n <- ns[i] ## Number of time steps

    te <- seq(0,T,length=n+1)

    ## Initial guess: linear interpolation in log domain
    xe <- x0^(te/T)*xT^(1-te/T)
    xv <- xe[-c(1,n+1)]

    ## Create a RTMB object so that we can optimize
    negloglik <- function(p) p.d.f(p$xv)
    p0 <- list(xv=xv)
    obj <- MakeADFun(negloglik,p0)
    
    ## Store the fits; most impotantly, the mode
    fits[[i]] <- c(nlminb(obj$par,obj$fn,obj$gr),list(te=te))

    reps[[i]] <- sdreport(obj)

    ress[[i]] <- c(exp(fits[[i]]$objective) / sqrt(det(obj$he()/2/pi)),
                   exp(fits[[i]]$objective) * sqrt(det(2*pi*reps[[i]]$cov.fixed)))
}

## Plot the bridges
pdf(file="../figures/Geometric-Brownian-Bridge.pdf",width=8)
plot(c(0,T),range(c(x0,xT,unlist(lapply(fits,function(f)f$par)))),
     type="n",xlab="t",ylab="x")
for(i in 1:length(ns))
    {
        lines(fits[[i]]$te,c(x0,fits[[i]]$par,xT),lty=i)
        points(fits[[i]]$te,c(x0,fits[[i]]$par,xT),lty=i,pch=i)
    }

legend("bottomright",legend=paste("n =",ns),pch=1:length(ns),lty=1:length(ns))
dev.off()
