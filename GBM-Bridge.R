## Plot of the conditional p.d.f. in the Geometric Brownian bridge,
## with two intermediate points

T <- 1
tv <- c(0,1/3,2/3,1) * T
dt <- diff(tv)
           
x0 <- 1
xT <- 1

r <- 1
sigma <- 1

xmin <- 0.1
xmax <- 2.1
nx <- 100

x1 <- seq(xmin,xmax,length=nx)
x2 <- seq(xmin,xmax,length=nx+1)

## Transition probabilities in the g.b.m.
p <- function(x0,xt,t) dlnorm(xt,log(x0)+r*t,sigma*sqrt(t))

f <- function(x1,x2) p(x0,x1,dt[1])*p(x1,x2,dt[2])*p(x2,xT,dt[3])/p(x0,xT,T)

P <- outer(x1,x2,f)

pdf(file="../figures/GBM-Bridge.pdf",width=7,height=7)
xlim <- c(0,xmax)
contour(x1,x2,P,xlab=expression(x[1/3]),ylab=expression(x[2/3]),asp=1,xlim=xlim,ylim=xlim)
dev.off()

require(numDeriv)

xinit <- approx(c(0,T),c(x0,xT),tv[2:(length(tv)-1)])$y
xhat <- nlminb(xinit,function(x)-f(x[1],x[2]))$par

H <- hessian(function(x)-log(f(x[1],x[2])),xhat)

f(xhat[1],xhat[2])/sqrt(det(H/2/pi))
