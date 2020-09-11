library(mgcv)
library(maxLik)

n_obs<-1000
n_levels<-10
x1 = runif(n_obs,-2,2)
rfx<-rnorm(n_levels,4.5,3)
mu = rep(rfx,(n_obs/n_levels)) + 2.3*x1 - 0.5*x1^2 + 1.9*x1^3
y = rnorm(1000,mu,2)
plot(x1,y,col="gray")
practice.dat <- data.frame(x1=x1,
                           rfx=as.factor(rep(letters[1:n_levels],(n_obs/n_levels))),
                           y=y)
testgam <- gam(y~s(x1)+s(rfx,bs="re"),data=practice.dat)
## prediction ignoring random effect
pred <- predict.gam(testgam,
                    newdata = data.frame(x1=seq(-2,2,0.01),
                                         rfx="foo"), 
                    exclude = "s(rfx)")
lines(seq(-2,2,0.01),pred,col="red",lwd=5)

# coefficients
beta <- testgam$coefficients
# linear predictor
Xp <- predict.gam(testgam,type = "lpmatrix")

plot(x1,Xp%*%beta) ## no the rfx levels apparent here

## can I "re-fit" this spline by ML using the design matrix?
LogLik=function(pars,response,U){
  pars1 = pars[1:ncol(U)]; pars2=pars[-(1:ncol(U))]
  mu = U%*%pars1;  
  val = dnorm(x = response, 
             mean=mu,
             sd = exp(pars2),
             log=T) 
  return(val); 
}
practice.out <- maxLik(logLik=LogLik,
                       start=rnorm((ncol(Xp)+1),0,1), 
                       response=y,
                       U=Xp,
                       method="BHHH",
                       control=list(iterlim=5000,printLevel=2),
                       finalHessian=FALSE)

plot(beta,practice.out$estimate[1:ncol(Xp)])
abline(0,1)

plot(x1,y,col="lightgray")
points(x1,Xp%*%practice.out$estimate[1:ncol(Xp)],pch=".",col="blue")
lines(seq(-2,2,0.01),pred,col="red",lwd=6)
## I can pull out the average across rfx levels by setting those coefficients to zero
beta_mean <- practice.out$estimate
beta_mean[colnames(Xp) %in% paste0("s(rfx).",1:n_levels)] <- 0
points(x1,Xp%*%beta_mean[1:ncol(Xp)],col="blue")
# WORKS
