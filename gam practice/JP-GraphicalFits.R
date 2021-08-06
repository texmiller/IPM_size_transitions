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
nu = seq(-3,3,length=150)
tau = seq(-1,1,length=151);  
for(i in 1:150) {
for(j in 1:151){
    delta=exp(-0.5*tau[j]); epsilon=delta*nu[i]; 
    SkewMat[i,j]=JP_NPskewness(epsilon,delta);
    KurtMat[i,j]=JP_NPkurtosis(epsilon,delta);
}}

graphics.off(); 
dev.new(); 
image.plot(nu,tau,SkewMat,col=plasma(64));  title(main="NP Skewness"); 
contour(nu,tau,SkewMat,add=TRUE) 

dev.new(); 
image.plot(nu,tau,KurtMat,col=plasma(64)); title(main = "NP Kurtosis"); 
contour(nu,tau,KurtMat,add=TRUE) 


################################################################
#  Display how ordinary skew and kurtosis depend on parameters
################################################################

SkewMat = KurtMat = matrix(NA,50,51);
nu = seq(-3,3,length=50)
tau = seq(-1,1,length=51);  
for(i in 1:50) {
for(j in 1:51){
    delta=exp(-0.5*tau[j]); epsilon=delta*nu[i]; 
    out = SJP_moments(epsilon,delta); 
    SkewMat[i,j]=out$skew;
    KurtMat[i,j]=out$excess.kurtosis;
}}

graphics.off(); 
dev.new(); 
image.plot(nu,tau,SkewMat,col=plasma(64));  title(main="Skewness"); 
contour(nu,tau,SkewMat,add=TRUE) 

dev.new(); 
image.plot(nu,tau,KurtMat,col=plasma(64)); title(main = "Excess Kurtosis"); 
contour(nu,tau,KurtMat,add=TRUE) 

   

######### TEST DATA: can we recover known epsilon and delta?  
prior.eps=function(epsilon) dnorm(epsilon,mean=0,sd=100); 
prior.del=approxfun(c(0,0.02,50,51),c(0,1,1,0),rule=2); 

nu.true = -0.5; tau.true = 1.5; 
delta.true=exp(-0.5*tau.true); epsilon.true=delta.true*nu.true;
SJP_moments(epsilon.true,delta.true); 


y = rSJP(500,epsilon = epsilon.true, delta=delta.true); 

# Likelihood times prior, pars=c(nu,tau)  
# These are closer than (epsilon,delta) to: nu measures skew, tau measures kurtosis 
Fun = function(par, data) {
   nu=par[1]; tau=par[2]; 
   delta=exp(-0.5*tau); epsilon=delta*nu; ### PAY CLOSE ATTENTION! 
   Lik = dSJP(data,epsilon,delta);
   if(sum(!is.finite(Lik)) > 0){
        val=0; 
        }else{
        val = prod(Lik)*prior.eps(epsilon)*prior.del(delta)
    }    
    return(val);  
} 

### Do Metropolis-Hastings with bivariate Gaussian proposal 
### with pars = c(nu,tau)) 
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
nu=pars[,1]; tau=pars[,2]; ### PAY CLOSE ATTENTION! 
delta=exp(-0.5*tau); epsilon=delta*nu; 

graphics.off(); 
dev.new(); par(mfrow=c(1,2)); plot(epsilon); plot(delta); 

graphics.off(); 
dev.new(); plot(nu,tau); 
points(nu.true,tau.true,pch=1,col="red",lwd=2,cex=2); 

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

z = -1+2*rbeta(500,2,2); z=sort(z); hist(z); 

########### Create artificial "residuals" with known parameters 
True.epsilon= rep(-0.8,length(z)); True.delta = exp(-0.5*z); 
resids = rSJP(length(z), epsilon=True.epsilon, delta = True.delta); 
par(mfrow=c(2,1)); 
hist(resids); plot(z,resids); 

jarque.test(resids)$p.value # normality test
agostino.test(resids)$p.value # skew 
anscombe.test(resids)$p.value # kurtosis 

dev.new(); 
out=rollMoments(z,resids,windows=5);

###########################################################################
# Fit sJP to binned data
###########################################################################
require(tidyverse); 
X <- data.frame(init=z,resids=resids); 
X <- X %>% mutate(size_bin = cut_number(init,n=5))
bins = levels(X$size_bin); 
X2=subset(X,size_bin==bins[2]); y = X2$resids; 





