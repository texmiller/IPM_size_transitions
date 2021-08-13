rm(list=ls(all=TRUE)); 

require(minqa); require(fda); require(viridisLite); 
require(moments); require(zoo); require(fields); 

setwd("c:/repos/IPM_size_transitions/gam practice"); 
source("JPfuns.R"); 
source("../Diagnostics.R"); 

notExp2 = function (x, d, b=1/d) exp(d * sin(x * b))  # from mgcv

graphics.off(); 

#######################################################
#  Display how NP skew and kurtosis depend on parameters
#######################################################

SkewMat = KurtMat = matrix(NA,150,151);
lambda = seq(-2,2,length=150)
tau = seq(-1.5,1,length=151);  
for(i in 1:150) {
for(j in 1:151){
    delta=exp(-tau[j]); epsilon = lambda[i]*delta; 
    SkewMat[i,j]=JP_NPskewness(epsilon,delta);
    KurtMat[i,j]=JP_NPkurtosis(epsilon,delta);
}}

dev.new(); 
image.plot(lambda,tau,SkewMat);  title(main="NP Skewness"); 
contour(lambda,tau,SkewMat,add=TRUE,labcex=1) 

dev.new(); 
image.plot(lambda,tau,KurtMat); title(main = "NP Kurtosis"); 
contour(lambda,tau,KurtMat,add=TRUE,labcex=1) 


################################################################
#  Display how ordinary skew and kurtosis depend on parameters
################################################################

SkewMat = KurtMat = matrix(NA,70,71);
lambda = seq(-2,2,length=70)
tau = seq(-1.5,1,length=71);  
for(i in 1:70) {
for(j in 1:71){
    delta=exp(-tau[j]); epsilon = lambda[i]*delta; 
    out = SJP_moments(epsilon,delta); 
    SkewMat[i,j]=out$skew;
    KurtMat[i,j]=out$excess.kurtosis;
}}

# graphics.off(); 
dev.new(); 
image.plot(lambda,tau,SkewMat);  title(main="Skewness"); 
contour(lambda,tau,SkewMat,add=TRUE,labcex=1) 

dev.new(); 
image.plot(lambda,tau,KurtMat); title(main = "Excess Kurtosis"); 
contour(lambda,tau,KurtMat,add=TRUE,labcex=1) 


######################################################################### 
## TEST DATA: can we recover known epsilon and delta?  
## Bayesian fit, using the (lambda, tau) parameterization with BayesianTools 
#########################################################################
epsilon.true = 1; delta.true = 0.5;  
lambda.true=epsilon.true/delta.true; tau.true = - log(delta.true); 
data = rSJP(500,epsilon = epsilon.true, delta=delta.true); 

ll = function(par) {
   lambda=par[1]; tau=par[2]; 
   delta=exp(-tau); epsilon=lambda*delta; 
   Lik = dSJP(data,epsilon,delta);
   if(sum(!is.finite(Lik)) > 0){
        val=0; 
        }else{
        val = prod(Lik)
    }    
    return(log(val));  
} 

bayesianSetup = createBayesianSetup(likelihood = ll, lower = rep(-10, 2), upper = rep(10, 2))
settings <- list(iterations = 100000, adapt = F, DRlevels = 1, nrChains = 3, gibbsProbabilities = c(0,1), temperingFunction = NULL, 
optimize = F,  message = TRUE, startValue=c(0,2))

### VERY SLOW compared to doing MH "by hand" below 
out <- runMCMC(bayesianSetup, sampler="Metropolis", settings = settings)


 
################################################################## 
## TEST DATA: can we recover known epsilon and delta?  
## Bayesian fit, using the (lambda, tau) parameterization 
##################################################################
epsilon.true = 1; delta.true = 0.5;  
lambda.true=epsilon.true/delta.true; tau.true = - log(delta.true); 
SJP_moments(epsilon.true,delta.true); 

y = rSJP(500,epsilon = epsilon.true, delta=delta.true); 
out = RSJP_ML(y); out; lambda.true; tau.true; 

# Likelihood times prior, pars=c(nu,tau)  
# These are closer than (epsilon,delta) to: nu measures skew, tau measures kurtosis 
Fun = function(par, data) {
   lambda=par[1]; tau=par[2]; 
   delta=exp(-tau); epsilon=lambda*delta; ### PAY CLOSE ATTENTION! 
   Lik = dSJP(data,epsilon,delta);
   if(sum(!is.finite(Lik)) > 0){
        val=0; 
        }else{
        val = prod(Lik)# *prior.eps(epsilon)*prior.del(delta)
    }    
    return(val);  
} 

### Do Metropolis-Hastings with bivariate Gaussian proposal 
### with pars = c(lambda,tau)) 
N = 10^6; scale=c(.2,0.2); 
pars = matrix(NA,N,2); pars[1,]=rnorm(2,0,0.1); oldFun = Fun(pars[1,],data=y); 
accept = 0; 
for(j in 2:N) {
    if(j%%10000==0) cat(j, "\n"); 
    newPars = pars[j-1,] + rnorm(2,0,scale); 
    newFun = Fun(newPars,data=y); 
    if(newFun >= oldFun) { 
        oldFun=newFun; pars[j,]=newPars; accept=accept+1; 
    }else{  
        r = newFun/oldFun; 
        if(runif(1)<r) {
            oldFun=newFun; pars[j,]=newPars; accept = accept + 1; 
         }else{
             pars[j,]=pars[j-1,]   
         }   
    }
}

keep = seq(N/4,N,by=50); 
pars=pars[keep,]; 
lambda=pars[,1]; tau=pars[,2]; ### PAY CLOSE ATTENTION! 
delta=exp(-tau); epsilon=delta*lambda; 

graphics.off(); 
dev.new(); par(mfrow=c(1,2)); plot(lambda); plot(tau); 

graphics.off(); 
dev.new(); plot(lambda,tau); 
points(lambda.true,tau.true,pch=1,col="red",lwd=2,cex=2); 

graphics.off(); 
dev.new(); plot(epsilon,delta); points(epsilon.true, delta.true,pch=1,col="red",lwd=2,cex=2); 

dev.new();
xx = 1.2*range(y); px = seq(xx[1],xx[2],length=1000); 
py1 = dSJP(px, epsilon=-1, delta=1); 
py2 = dSJP(px, epsilon=-5, delta=1.5); 
par(yaxs="i",bty="l",cex=1.2); 
matplot(px,cbind(py1,py2), col=c("black","red"),type="l",lwd=2, lty=c(1,2), xlab="Value",ylab="Scaled JP density function") 
legend("topleft",legend=c("epsilon = -1, delta = 1", "epsilon = -5, delta = 1.5"), col=c("black","red"),lty=c(1,2), lwd=2,pch=NULL,bty="n"); 



############################################################
#  Explore what diagnostic plots can reveal about residuals
############################################################    
    
######### Create covariate for residuals 
z = rbeta(500,2,2); z=sort(z); hist(z); 

########### Create artificial "residuals" with known parameters 
true.lambda = -z^2; true.tau = 0.25  
true.delta=exp(-true.tau); true.epsilon=true.delta*true.lambda; 

resids = rSJP(length(z), epsilon=true.epsilon, delta = true.delta); 
par(mfrow=c(2,1)); 
hist(resids); plot(z,resids); 

jarque.test(resids)$p.value # normality test
agostino.test(resids)$p.value # skew 
anscombe.test(resids)$p.value # kurtosis 

graphics.off(); dev.new(); 
out=rollMomentsNP(z,resids,windows=4);

###########################################################################
# Fit sJP to binned data
###########################################################################
require(tidyverse); 
X <- data.frame(init=z,resids=resids); 
X <- X %>% mutate(size_bin = cut_number(init,n=6))
bins = levels(X$size_bin);
nbins=length(bins); 

fits = list(nbins); 
fittedPars=matrix(NA,nbins,3); 
for(j in 1:nbins){
    Xj=subset(X,size_bin==bins[j]);
    fits[[j]] = SJP_ML(Xj$resids); 
    fittedPars[j,1]=median(Xj$init); 
    fittedPars[j,2:3] = fits[[j]]$estimate; 
}    
dev.new(); matplot(fittedPars[,1],fittedPars[,2:3], type="o", lty=1, pch=1); 






