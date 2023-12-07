## -----------------------------------------------------------------------------
require(SDEtools)
require(Matrix)
set.seed(1)  ## Fix the seed so that I know what the plot looks like :) 

sigma <- 1
alpha <- 1
gamma <- 0.5
r <- 0.5

T <- 10

dt <- 0.01
tvec <- seq(0,T,dt)
x0 <- 0.1

## Value function
lambda <- (1-gamma)*(r + (alpha - r)^2/(2*gamma*sigma^2))
b <- function(t) exp(lambda * (t-T))
Vopt <- function(x,t) b(t)*x^(1-gamma)/(1-gamma)

p <- (alpha-r)/(gamma*sigma^2)
Vopt_stationary <- function(x) min(1,p)

fopt <- function(x,t) x*(r+Vopt(x,t)*(alpha-r))
fopt_stationary <- function(x) x*(r + Vopt_stationary(x) * (alpha-r))
fci1 <- function(x) x*(r + 0.1 * (alpha - r))
fci2 <- function(x) x*(r + 0.5 * (alpha - r))

g_stationnary <- function(x) sigma*Vopt_stationary(x)*x

## Number of realizations
M <- 100

## Generate noise for all realizations 
B <- rvBM(tvec,n=M)

sol.opt <- sapply(1:M,function(i)euler(fopt_stationary,g,tvec,x0,B=B[,i],p=abs)$X)
sol.ci1 <-  sapply(1:M, function(i)euler(fci1,g,tvec,x0,B=B[,i],p=abs)$X)
sol.ci2 <-  sapply(1:M, function(i)euler(fci2,g,tvec,x0,B=B[,i],p=abs)$X)

J.opt <- mean(sqrt(sol.opt^2))
J.ci1 <- mean(sqrt(0.5*sol.opt))
J.ci2 <- mean(sqrt(0.5*sol.opt))


## -----------------------------------------------------------------------------
x11()
plot(tvec, sol.opt[,2], type="l", xlab="t", ylab="X", col="blue", lty=1, 
     ylim=range(c(sol.opt[,2], sol.ch[,2])))
lines(tvec, sol.ci1[,2], col="red", lty=1)
lines(tvec, sol.ci2[,2], col="green", lty=1)
legend("topright", 
       legend=c("Optimal control", 
                "Constant investment 0.1", 
                "Constant investment 0.5"), 
       col=c("blue", "red",'green'), lty=c(1,1,1))

## -----------------------------------------------------------------------------
print(c(J.opt,J.ch))


## -----------------------------------------------------------------------------
require(SDEtools)

### Discretization of state space
Xmax <- 4
dx <- 0.01
xi <- seq(0,Xmax,dx)
xc <- xi[-1] - 0.5*diff(xi)

sigma <- 1

### Functions entering into the model

### Uncontrolled system:
D <- function(x) 1/2*sigma^2*x^2
dD <- function(x) sigma^2*x

f <- function(x) x*(1-x)
advection <- function(x) f(x) - dD(x)

G0 <- fvade(advection,D,xi,'r')

### Effect of the fishing: The "generator" d/dx
G1 <- fvade(function(x)-1,function(x)0,xi,'r')

ubound <- c(0,rep(100,nrow(G1)-1))

k <- function(u) sqrt(u)

vbar <- c(Inf,rep(0.01,nrow(G0)-1))
hack <- function(dV) pmax(-dV,vbar)
uopt <- function(dV) 1/4/hack(dV)^2

sol <- PolicyIterationSingular(G0,G1,k,uopt,do.minimize = FALSE)


## -----------------------------------------------------------------------------
x11()
par(mfrow=c(2,1))
plot(xc,sol$V - approx(xc,sol$V,1)$y,xlab="x",ylab="V")
lines(xc,0.5*log(xc))

plot(xc,sol$u,ylim=c(0,max(xc)^2),xlab="x",ylab="u")
lines(xc,xc^2)


## -----------------------------------------------------------------------------
pv <- c(0.5,1,2)

sols <- list()


for(p in pv) {
  f <- function(x) x*(1-x^p)
  advection <- function(x) f(x) - dD(x)
  G <- fvade(advection,D,xi,'r')
  
  sol <- PolicyIterationSingular(G,G1,k,uopt,do.minimize = FALSE)
  
  sols[[length(sols)+1]] <- list(p=p,sol=sol)
}


## -----------------------------------------------------------------------------
x11()
par(mfrow=c(2,1))
matplot(xc,cbind(sols[[1]]$sol$V,sols[[2]]$sol$V,sols[[3]]$sol$V),
        type="l",xlab="x",ylab="V")
matplot(xc,cbind(sols[[1]]$sol$u,sols[[2]]$sol$u,sols[[3]]$sol$u),
        type="l",xlab="x",ylab="u",ylim=c(0,20))