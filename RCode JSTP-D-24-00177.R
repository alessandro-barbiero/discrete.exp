# R Code for "A new discrete exponential distribution: properties and applications"
# Accompanies the manuscript currently under review (as of January, 2025)
# Author: Alessandro Barbiero

rm(list=ls())
library(bbmle) # for Maximum Likelihood estimation
# Proposed Discrete Exponential Distribution
# cdf
pdexp <- Vectorize(function(x, lambda) ifelse(x<0,0,1-exp(-lambda*x)/lambda*(1-exp(-lambda))))
# pmf
ddexp <- function(x, lambda) pdexp(x,lambda)-pdexp(x-1,lambda)
# quantile function
qdexp <- function(alpha, lambda) ceiling(-1/lambda*log(lambda*(1-alpha)/(1-exp(-lambda))))
# random generation
rdexp <- function(n, lambda) 
{
  u <- runif(n)
  qdexp(u, lambda)
}

# some graphs of the pmf
op<-par()
par(mai=c(0.8,0.8,0.1,0.1), mgp=c(2.5,1,0))
x <- 0:10
lambda <- c(2,1,0.5,0.2,0.1)
plot(x, ddexp(x, lambda[1]), pch=0, col=8,axes=FALSE, ylab="pmf")
for(j in 2:5)
{
points(x,ddexp(x,lambda[j]),pch=(13+j),col=j+7)
}
axis(1,x)
axis(2)
box()
legend_text <- as.expression(lapply(lambda, function(d) {
  bquote(lambda==.(d))
} ))
legend("topright",legend_text,pch=c(0,15:18),col=c(8:12))
par(op)

# failure rate function
frdexp <- function(x, lambda)
{
  ddexp(x, lambda)/(1-pdexp(x-1,lambda))
}
# moments
E1 <- function(lambda) 1/lambda
E2 <- function(lambda) (1+exp(-lambda))/lambda/(1-exp(-lambda))
E3 <- function(lambda) (4*exp(-lambda)+exp(-2*lambda)+1)/lambda/(1-exp(-lambda))^2
E4 <- function(lambda) (11*exp(-lambda)+11*exp(-2*lambda)+exp(-3*lambda)+1)/lambda/(1-exp(-lambda))^3
V  <- function(lambda) E2(lambda)-(E1(lambda))^2
Skewness <- function(lambda) (E3(lambda)-(E1(lambda))^3-3*E1(lambda)*V(lambda))/(V(lambda))^1.5
Kurtosis <- function(lambda) (E4(lambda)-4*E1(lambda)*E3(lambda)+3*E1(lambda)^2*(E1(lambda)^2+2*V(lambda)))/V(lambda)^2

# graphs of Skewness and Kurtosis
par(mfrow=c(1,2))
plot(function(lambda) Skewness(lambda), xlim=c(0.01,10), xlab=expression(lambda), ylab="Skewness")
plot(function(lambda) Kurtosis(lambda), xlim=c(0.01,10), xlab=expression(lambda), ylab="Kurtosis")
par(op)

# comparison with the geometric distribution
p <- 0.2
lambda <- -log(1-p)
lambda
x <- 0:10
pg <- dgeom(x,p)
pe <- ddexp(x,lambda)
prob<-cbind(pg,pe)
op <- par()
par(mai=c(0.8,0.8,0.2,0.2))
barplot(t(prob), pch=8, ylab=expression(p[x]), beside=TRUE, col=c("red","blue"),names.arg=0:10,density=c(-1,30),angle=45)
text=c("Geometric", "Discrete Exp.")
legend("topright",text,col=c("red","blue"),pch=c(15,12))
par(op)
lines(c(0,0),c(pg[1],pe[1]),lwd=1,col="red")
for(i in 1:10)
{
lines(c(i,i),c(pg[i+1],pe[i+1]),lwd=1,col="green")
}
# comparison with the geometric distribution (cdf)
par(mai=c(.7,.7,.1,.1),mgp=c(2.25,1,0),cex=1)
lambda <- 1/3
p <- 1-exp(-lambda)
h <- 0:10
plot(function(x) pexp(x, lambda), xlim=c(0,9.75),ylab="pdf and pmf",ylim=c(0,1))
abline(h=0)
points(c(-.5,h), pgeom(c(-.5,h), p), type="s",col="red")
points(c(-.5,h), pgeom(c(-.5,h), p), type="p",col="red",pch=17)
points(c(-.5,h), pdexp(c(-.5,h), lambda), type="s", col="blue")
points(c(-.5,h), pdexp(c(-.5,h), lambda), type="p", col="blue",pch=19)
abline(h=1, lty=2)
legend.text <- c("geometric","discrete exp")
legend("bottomright", legend=legend.text, pch=c(17,19), col=c("red","blue"),text.col=c("red","blue"), bty="n")
par(op)

# Dispersion Index (DI)
DI <- function(lambda) (exp(-lambda)*(lambda+1)+lambda-1)/lambda/(1-exp(-lambda))
plot(function(lambda) DI(lambda), xlim=c(0.2,10), ylab="DI", xlab=expression(lambda))
abline(h=1,lty=3)
lambda.d <- optimize(function(lambda) DI(lambda)-1, lower=2, upper=4)$minimum
axis(1,at=lambda.d,labels=expression(lambda[d]))

# Zero Modification Index (ZMI)
ZMI <- Vectorize(function(lambda) {1+log(ddexp(0,lambda))*lambda})
plot(function(lambda) ZMI(lambda), xlim=c(0.01,10), xlab=expression(lambda), ylab="ZMI",axes=FALSE)
axis(2)
abline(h=0)
lines(c(1,1),c(-1,0), lty=3)
uniroot(function(lambda) ZMI(lambda),lower=0.01,upper=10)
opt<-optimize(ZMI,lower=0,upper=10)
lines(c(opt$minimum,opt$minimum),c(-1,opt$objective), lty=3)
axis(1,at=c(0,1,opt$minimum,5,10),labels = c(0,1,expression(lambda[m]),5,10))
box()

# Stress-strength Reliability Parameter
RelSS <- function(lambda1,lambda2)
{
  (lambda1-1+exp(-lambda1))/lambda1 + (1-exp(-lambda1))^2*(1-exp(-lambda2))/(lambda1*lambda2)/(1-exp(-lambda1-lambda2))
}

# Right-Tail Index
RTI <- function(lambda)
{
  sqrt(lambda*(1-exp(-lambda)))/(1-exp(-lambda/2)) - 1
}
plot(function(lambda) RTI(lambda), xlim=c(0,5), xlab=expression(lambda), ylab="RTI")

# Shanno Entropy
ShE <- function(lambda)
{
  -ddexp(0,lambda)*log(ddexp(0,lambda)) + exp(-lambda) -(1-exp(-lambda))/lambda*log((1-exp(-lambda))^2/lambda)
}
plot(function(lambda) ShE(lambda), ylim=c(0,6),xlim =c(0.01,10), xlab=expression(lambda), ylab="Shannon entropy")
abline(h=0,lty=3)
plot(function(lambda) 1-log(lambda), col="red", lty=3, add=TRUE, ylim=c(0,6),xlim =c(0.01,10))

##################################################
# GENERALIZED TWO-PARAMETER DISCRETE EXPONENTIAL #
##################################################
pdgenexp <- function(x,lambda,alpha) pdexp(x,lambda)^alpha
ddgenexp <- function(x,lambda,alpha) pdexp(x,lambda)^alpha - pdexp(x-1,lambda)^alpha
qdgenexp <- function(u, lambda, alpha) ceiling(-1/lambda*log(lambda*(1-u^(1/alpha))/(1-exp(-lambda))))

# Plot of pmf
par(mfrow=c(3,4), mai=c(0.3,0.5,0.1,0.1),mgp=c(2,1,0),cex=.75)
lambda <- c(0.2,0.5,1)
alpha <- c(.5, .75, 1.25, 2)
ylim <- c(0.31,0.46,0.61)
x<-0:15
for(i in 1:3)
{
for (j in 1:4)
{
plot(x, ddgenexp(x, lambda[i], alpha[j]), type="h",ylab="pmf", ylim=c(0,ylim[i]))
leg <- c(as.expression(bquote(alpha==.(alpha[j]))),as.expression(bquote(lambda==.(lambda[i]))))
legend("topright", legend = leg, bty="n")
}
}

# RANDOM GENERATION
rdgenexp <- function(n, lambda,alpha) 
{
  u <- runif(n)
  qdgenexp(u, lambda,alpha)
}
# EXPECTED VALUE
E1gen <- function(lambda,alpha)
{
  if(lambda<=0 | alpha <=0) return(9999)
  k <- 5*qdgenexp(.999,lambda,alpha)
  if(k>1e5) return(integrate(function(t) 1-pexp(t,lambda)^alpha,lower=0,upper=Inf)$value)
  else
  {
    x <- 0:k
    return(sum(x*ddgenexp(x,lambda,alpha)))
  }
}
E1gen(1/4,2)
# SECOND RAW MOMENT
E2gen <- function(lambda,alpha)
{
  if(lambda<=0 | alpha <=0) return(9999)
  k <- 5*qdgenexp(.999,lambda,alpha)
  if(k>1e5) return(2*integrate(function(t) t*(1-pexp(t,lambda)^alpha),lower=0,upper=Inf)$value)
  else
  {
      x <- 0:k
      return(sum(x^2*ddgenexp(x,lambda,alpha)))
  }
}
E2gen(1/4,2)


# example
x <- rdgenexp(100, 1, 2)
hat.p0<-mean(x==0)
hat.p1<-mean(x==1)
#c <- (log(pdexp(1,lambda))-log(pdexp(0,lambda)))/log(pdexp(0,lambda))
c <- (log(hat.p0+hat.p1)-log(hat.p0))/log(hat.p0)
lambda.hat <- uniroot(function(t) t-exp(-t)*(1-exp(-t))-(t-(1-exp(-t)))^(1+c)/t^(c),lower=0.001,upper=10)$root
lambda.hat
alpha.hat <- log(hat.p0)/log(1-1/lambda.hat*(1-exp(-lambda.hat)))
alpha.hat

# LOSS FUNCTION for the method of moments
x <- rdexp(100,1/2)
fn.mom <- function(par)
{
  lambda <- par[1]
  alpha <-  par[2]
  (mean(x)-E1gen(lambda,alpha))^2 + (mean(x^2)-E2gen(lambda,alpha))^2
}
optim(c(1/mean(x),1), fn.mom)


# LOG-LIKELIHOOD function discrete exp
dexp.loglik <- function(lambda,x) 
{
  -sum(log(ddexp(x,lambda)))
}
# LOG-LIKELIHOOD function generalized discrete exp
dgenexp.loglik <- function(lambda,alpha,x) 
{
  -sum(log(ddgenexp(x,lambda,alpha)))
}
######################
# ESTIMATION ROUTINE #
####### START ########
######################
est.dgenexp <- function(x, method="ML")
{
  eps <- 0.0001
  hat.p0 <- mean(x==0)
  hat.p1 <- mean(x==1)
    if(method=="ML")
  {
    if(hat.p0==0 | hat.p1==0 | hat.p0+hat.p1==1) {return(list(c(NA,NA),matrix(NA,2,2)))}
    else
    {
    res <- mle2(dgenexp.loglik,start=list(lambda=1/mean(x),alpha=1),data=list(x=x))
    return(list(res@coef,confint(res)))
    }
  }
  if(method=="M")
  {
    res <- optim(c(1/mean(x), 1), fn.mom)
    return(res$par)
  }
  else
  {
      if(hat.p0==0 | hat.p1==0 | hat.p0+hat.p1==1)
      return(c(NA,NA)) else
    {
    c <- (log(hat.p0+hat.p1)-log(hat.p0))/log(hat.p0)
    lambda.hat <- uniroot(function(t) t-exp(-t)*(1-exp(-t))-(t-(1-exp(-t)))^(1+c)/t^(c),lower=1e-8,upper=10)$root
    alpha.hat <- log(hat.p0)/log(1-1/lambda.hat*(1-exp(-lambda.hat)))
    alpha.hat
    return(c(lambda.hat,alpha.hat))
    }
  }
}
######################
######## END #########
# ESTIMATION ROUTINE #
######################

##############################################
######### Monte Carlo simulation plan ########
################### START ####################
##############################################
set.seed(12345)
n      <- 100
lambda <- 1
alpha  <- 1
S      <- 100
v.par.P  <-  matrix(0,S,2)
v.par.M  <-  matrix(0,S,2)
v.par.ML  <- matrix(0,S,2)
IC.lambda<- matrix(0,S,2)
IC.alpha <- matrix(0,S,2)
for(i in 1:S)
{
cat("Simulation #",i,"\n")
    x <- rdgenexp(n, lambda, alpha)
  v.par.P[i,] <- est.dgenexp(x, "P")
  v.par.M[i,] <- est.dgenexp(x, "M")
  ML.res <- est.dgenexp(x, "ML")
  v.par.ML[i,]<- ML.res[[1]]
  IC.lambda[i,] <- ML.res[[2]][1,]
  IC.alpha[i,]  <- ML.res[[2]][2,]
}
filename <- paste("lambda",lambda,"alpha",alpha,"n",n)
filetxt <- paste(filename,"txt",sep=".")
filename <- paste(filename,"Rdata",sep=".")
sink(file=filetxt)

colnames(v.par.P)<-c("Plambda","Palpha")
colnames(v.par.M)<-c("Mlambda","Malpha")
colnames(v.par.ML)<-c("MLlambda","MLalpha")

cat("Sample means:","\n")
apply(v.par.P,2,mean,na.rm=TRUE)
cat("Proportion of non-feasible samples for MP","\n")
mean(is.na(v.par.P))  
apply(v.par.M,2,mean) 
apply(v.par.ML,2,mean,na.rm=TRUE)
cat("Sample RMSEs:","\n")
rmse.theta <- function(v, theta) sqrt(mean((v-theta)^2,na.rm=TRUE))
cat("P","\n")
rmse.theta(v.par.P[,1],lambda)
rmse.theta(v.par.P[,2],alpha)
cat("M","\n")
rmse.theta(v.par.M[,1],lambda)
rmse.theta(v.par.M[,2],alpha)
cat("ML","\n")
rmse.theta(v.par.ML[,1],lambda)
rmse.theta(v.par.ML[,2],alpha)

cat("CI coverage:","\n")
mean(IC.lambda[,1]<lambda & lambda<IC.lambda[,2],na.rm=TRUE)
mean(IC.alpha[,1]<alpha & alpha<IC.alpha[,2],na.rm=TRUE)
cat("CI length:","\n")
mean(IC.lambda[,2] - IC.lambda[,1], na.rm=TRUE)
mean(IC.alpha[,2]  - IC.alpha[,1], na.rm=TRUE)
save(list = ls(all.names = TRUE),file=filename)
sink()
##############################
par(mfrow=c(1,2), mai=c(0.4,0.4,0.4,0.4))
boxplot(v.par.P[,1],v.par.M[,1],v.par.ML[,1],names=c("MP","MM","ML"),main=expression(hat(lambda)))
abline(h=0.5, lty=3)
boxplot(v.par.P[,2],v.par.M[,2],v.par.ML[,2],names=c("MP","MM","ML"),main=expression(hat(alpha)))
abline(h=0.75, lty=3)
par(op)
############################################
####### Monte Carlo simulation plan ########
################### END ####################
############################################

# DATA ANALYSIS
freq <- c(58,34,14,8,6,2) # CONSUL TABLE 5.3 PAG.119 OK d-gen-exp
h    <- 0:5
data <- rep(h,freq)
x    <- data
# sample moments
m <- mean(x)
s <- mean((x-m)^2)
mean((x-m)^3)/s^1.5
mean((x-m)^4)/s^2
# fitting the two-parameter discrete exp via MLE
res.two <- mle2(dgenexp.loglik,start=list(lambda=1,alpha=1),data=list(x=x))
summary(res.two)
res.two@coef
ddgenexp(0:5,res.two@coef[1],res.two@coef[2])
#
plot(ecdf(x), verticals=TRUE, ylab="ecdf and cdf",main="")
lines(c(5,7),c(1,1))
points(0:7,pdgenexp(0:7,res.two@coef[1],res.two@coef[2]),col="red",type="s")
points(0:7,pdgenexp(0:7,res.two@coef[1],res.two@coef[2]),col="red")
legend("right",legend=c("ecdf","cdf"),pch=c(19,1),col=c("black","red"))
#
ntheo <- ddgenexp(0:5,res.two@coef[1],res.two@coef[2])*length(x)
ntheo[6] <- length(x)-sum(ntheo[1:5])
ntheo
ntheo[5]<- ntheo[5] + ntheo[6]
ntheo <- head(ntheo,-1)
nobs<-table(x)
nobs[6] <- nobs[5] + nobs[6]
nobs <- head(nobs,-1)
cbind(nobs,ntheo)
chi2 <- sum((ntheo-nobs)^2/ntheo)
1-pchisq(chi2,5-2-1)
# fitting the one-parameter discrete exp via MLE
res.one <- mle2(dexp.loglik,start=list(lambda=1),data=list(x=x))
summary(res.one)
res.one@coef
ntheo <- ddexp(0:5,res.one@coef)*length(x)
ntheo[6] <- length(x)-sum(ntheo[1:5])
ntheo
nobs<-table(x)
ntheo[5]<- ntheo[5] + ntheo[6]
ntheo <- head(ntheo,-1)
nobs[6] <- nobs[5] + nobs[6]
nobs <- head(nobs,-1)
chi2 <- sum((ntheo-nobs)^2/ntheo)
1-pchisq(chi2,5-1-1)
# fitting the negative binomial distribution
L.NB <- function(x, par)
{
  size<-par[1]
  prob<-par[2]
  return(-sum(log(dnbinom(x, size, prob))))
}
#
res.nb <- optim(par = c(2,0.5), L.NB, x = x, method = "L-BFGS-B", lower = c(1,0.001), upper = c(Inf,.999))
res.nb$value*2 + 2*2 # AIC
dnbinom(0:5, res.nb$par[1], res.nb$par[2])*length(x) 
# alternatively, log-likelihood function for NegBin
dnegbin.loglik <- function(n,p,x) 
{
  -sum(log(dnbinom(x,n,p)))
}
res.nb.alt <- mle2(dnegbin.loglik,start=list(n=2,p=0.5),data=list(x=x))
summary(res.nb.alt)

ntheo <- dnbinom(0:5, res.nb.alt@coef[1], res.nb.alt@coef[2])*length(x) 
ntheo[6] <- length(x)-sum(ntheo[1:5])
ntheo
ntheo[5]<- ntheo[5] + ntheo[6]
ntheo <- head(ntheo,-1)
nobs<-table(x)
nobs[6] <- nobs[5] + nobs[6]
nobs <- head(nobs,-1)
cbind(nobs,ntheo)
chi2 <- sum((ntheo-nobs)^2/ntheo)
1-pchisq(chi2,5-2-1)

#################################
########## APPLICATION ##########
#################################
### discretization and use of ###
## Panjer recursive formula #####
##### compound distribution #####
#### GEOMETRIC + EXPONENTIAL ####
# sum of N iid non-negative rvs X
# N ~ Geom(p), sitting on N={0,1,2,3,...}
# X ~ Exp(beta)

p      <- 0.2
beta   <- 1/4
nmax <- 500 # the maximum value of S we will consider (at least for the graphs)
# discretization of the exponential rv
f.C   <-  ddexp(0:nmax, beta)      # min. Cramer distance
f.S   <- dgeom(0:nmax,1-exp(-beta)) # preservation of the cdf at integer values
# approximation-by-discretization
g.C   <- numeric(nmax+1)
g.S   <- numeric(nmax+1)
#
g.C[1]   <- dgeom(0,p) + sum(dgeom(1:nmax,p)*f.C[1]^(1:nmax))
g.S[1]   <- dgeom(0,p) + sum(dgeom(1:nmax,p)*f.S[1]^(1:nmax))
# g_0 = dgeom(0,p) +... 
# approximate value of P(S=x) via
# Panjer recursive formula
for(i in 1:nmax)
{
  # for a Geometric(p)
  # a=1-p, b=0, f0=0
  g.C[i+1]   <- 1/(1-(1-p)*f.C[1])*sum((1-p)*f.C[2:(i+1)]*g.C[i:1])
  g.S[i+1]   <- 1/(1-(1-p)*f.S[1])*sum((1-p)*f.S[2:(i+1)]*g.S[i:1])
}
#####################################################
# approximate cumulative distribution function of S #
#####################################################
gg.C   <- cumsum(g.C)
gg.S   <- cumsum(g.S)
# the exact cdf of S is a mixed-type rv
G      <- numeric(nmax+1)
# with probability mass = p at zero
# and a continuous part at R^+ which is exponential with parameter beta*p
# so its expression for x>=0 is F_S(x) = p + (1-beta*p*exp(-beta*p*x))*(1-p)
G[1]   <- dgeom(0,p) # P(S<=0)=p
G[2:(nmax+1)] <- p + (1-exp(-beta*p*(1:nmax)))*(1-p)
##### ERRORS #####
par(mai=c(.6,.8,.1,.5),mgp=c(2,1,0),las=2)
plot(0:nmax, abs(gg.S-G), col="red",pch=19,ylab="",xlab="x",xlim=c(0,500),las=2)
points(0:nmax, abs(gg.C-G), col="blue",pch=0)
abline(h=0, lty=3)
legend.text <-c("geometric","new discr.exp","true cdf")
legend.col <- c("red","blue","green")
legend(legend.text,x="right",col=legend.col,pch=c(19,0,1),bty="n",inset=0.05)
mtext(text=expression(abs(F[S]-F[tilde(S)])),side=2,line=1)
sc <- 1/max(abs(gg.S-G))
#
points(0:nmax, G[1:(nmax+1)]/sc, col="green")
ax <- seq(0,0.05,0.01)
axis(4, at=c(0,1/2/sc,1/sc), labels=c(0,0.5,1), las=2, col.axis="green")
