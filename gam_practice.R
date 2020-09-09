
  x1 = runif(1000,-2,2)
  x2 = runif(1000,4,15)
  mu = 4.5 + 2.3*x1 + 0.7*x1^2 - 0.09*x1^3
  y = rnorm(1000,mu,2)

plot(x1,y)
testgam <- gam(y~s(x1),data=data.frame(y=y,x1=x1))
pred <- predict.gam(testgam,newdata = data.frame(x1=seq(-2,2,0.01)))
lines(seq(-2,2,0.01),pred,col="red",lwd=3)

# coefficients
beta <- testgam$coefficients
# linear predictor
Xp <- predict.gam(testgam,type = "lpmatrix")

plot(x1,Xp%*%beta)

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
points(x1,Xp%*%practice.out$estimate[1:ncol(Xp)],pch=16,col="blue")
lines(seq(-2,2,0.01),pred,col="red",lwd=3)
