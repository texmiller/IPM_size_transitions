graphics.off(); 

library(mgcv);library(maxLik); library(fda); 

n_obs<-1000
n_levels<-10
x1 = runif(n_obs,-2,2)
rfx<-rnorm(n_levels,0,3)
mu = rep(rfx,(n_obs/n_levels)) + 2.3*x1 - 0.5*x1^2 + 1.9*x1^3
y = rnorm(n_obs,mu,5)
plot(x1,y,col="gray")
practice.dat <- data.frame(x1=x1,
                           rfx=as.factor(rep(letters[1:n_levels],(n_obs/n_levels))),
                           y=y)
testgam <- gam(y~s(x1)+s(rfx,bs="re"),data=practice.dat)

## prediction ignoring random effect
pred <- predict.gam(testgam,newdata = data.frame(x1=seq(-2,2,0.01),rfx="foo"), exclude = "s(rfx)")
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


############### Choose model complexity based on AIC 
AICs=rep(10^12,12); fits=list(12);
for(k in 4:12) {
    testgam <- gam(y~s(x1,k=k)+s(rfx,bs="re"),data=practice.dat)
    Xp <- predict.gam(testgam,type = "lpmatrix")
    practice.outk <- maxLik(logLik=LogLik,
                       start=rnorm((ncol(Xp)+1),0,1), 
                       response=y,
                       U=Xp,
                       method="BHHH",
                       control=list(iterlim=5000,printLevel=2),
                       finalHessian=FALSE)
    fits[[k]]= practice.outk;                   
    AICs[k]=2*length(practice.outk$estimate) - 2*practice.outk$maximum; 
cat(k,length(practice.outk$estimate),"\n");                    
}                       

### plot the AIC-selected choice of dimension 
k.best = which.min(AICs); 
testgam <- gam(y~s(x1,k=k.best)+s(rfx,bs="re"),data=practice.dat)
beta <- testgam$coefficients
Xp <- predict.gam(testgam,type = "lpmatrix")
beta_mean <- fits[[k.best]]$estimate
beta_mean[colnames(Xp) %in% paste0("s(rfx).",1:n_levels)] <- 0
points(x1,Xp%*%beta_mean[1:ncol(Xp)],col="forestgreen",pch=16)

### AICs for the ML fits, as a function of basis dimension 
dev.new(); plot(4:12,AICs[4:12]); 

########################################################### 
# Now do the same process, using fda for the spline basis 
###########################################################
dev.new(); 
plot(x1,y,col="gray")
fac.rfx = as.factor(rep(letters[1:n_levels],(n_obs/n_levels)))

practice.dat <- data.frame(y=y,x1=x1,rfx=fac.rfx)
U = model.matrix(~practice.dat$rfx-1);

AICs = rep(10^6,12); 
for(k in 4:12) {
    bspBasis = create.bspline.basis(rangeval=range(x1),nbasis=k); 
    B = eval.basis(x1,bspBasis); 
    Xb = cbind(B,U); 
    fitk = lm(y~Xb-1); # it's Gaussian, so we can do it the quick way. 
    AICs[k]=AIC(fitk); 
}    
    
k.best = which.min(AICs); 
bspBasis = create.bspline.basis(rangeval=range(x1),nbasis=k.best); 
B = eval.basis(x1,bspBasis); 
Xb = cbind(B,U); 
fitk = lm(y~Xb-1); 
pars = coef(fitk)[1:ncol(B)]
yhat = B%*%pars; 
points(x1,yhat,type="p",pch=16,col="forestgreen"); 








