####################################################################################
## Exploring the feasibility of quantile-regression diagnostics for scaled residuals
## Fake 'scaled residuals' are generated from JP distribution
####################################################################################


library(quantreg); library(fda); library(sn); 
library(moments); library(mgcv); 
graphics.off(); 

root = "c:/repos"; # edit as needed 
setwd(root); setwd("IPM_size_transitions/gam practice"); 

source("JPfuns.R");

### Generate data from RSJP distribution
##  JP with mean=0, var=1, in (lambda, tau) parameters. 
nx = 500; 
x = sort(2*rbeta(nx,3,3)); 
lambdas = (x-1)^2; lambdas=lambdas-mean(lambdas); 
y = numeric(nx); 
for (i in 1:nx) y[i]=rRSJP(1,lambdas[i],tau=0) 

## Quantiles of the distribution 
Y = matrix(NA,nx,3)
for(i in 1:nx) {
    out = qRSJP(c(0.1,0.5,0.9),lambda=lambdas[i], tau=0);
    Y[i,]=out; 
    if(i%%20==0) cat(i, "\n"); 
} 

# plot data and true quantiles 
plot(x,y); 
matpoints(x,Y, type="l",lty=2,lwd=2,col="red"); 


xydata = data.frame(x = x, y = y); 
S.10 = qgam(y~s(x),data=xydata,qu=0.1,argGam=list(gamma=1)); 
   q.10 = predict(S.10,newdata=xydata); 
S.50 = qgam(y~s(x),data=xydata,qu=0.5,argGam=list(gamma=1)); 
   q.50 = predict(S.50,newdata=xydata); 
S.90 = qgam(y~s(x),data=xydata,qu=0.9,argGam=list(gamma=1)); 
   q.90 = predict(S.90,newdata=xydata);
matpoints(x, cbind(q.10,q.50,q.90), type="l", lty=1, lwd=1, col="blue"); 

NPS_hat = matrix(NA,nx,10); 
NPS_hat[,1] = (q.10 + q.90 - 2*q.50)/(q.90 - q.10);

y0 = y; 
for(k in 2:10) {
  for (i in 1:nx) y[i]=rRSJP(1,lambdas[i],tau=0) 
  xydata = data.frame(x = x, y = y); 
  S.10 = qgam(y~s(x),data=xydata,qu=0.1,argGam=list(gamma=1)); 
    q.10 = predict(S.10,newdata=xydata); 
  S.50 = qgam(y~s(x),data=xydata,qu=0.5,argGam=list(gamma=1)); 
    q.50 = predict(S.50,newdata=xydata); 
  S.90 = qgam(y~s(x),data=xydata,qu=0.9,argGam=list(gamma=1)); 
    q.90 = predict(S.90,newdata=xydata);
  matpoints(x, cbind(q.10,q.50,q.90), type="l", lty=1, lwd=1, col="blue"); 
  NPS_hat[,k] = (q.10 + q.90 - 2*q.50)/(q.90 - q.10);
 
}
matpoints(x,Y, type="l",lty=2,lwd=2,col="red"); 

NPS = (Y[,1]+Y[,3]-2*Y[,2])/(Y[,3]-Y[,1]); 
dev.new(); 
matplot(x,cbind(NPS,NPS_hat),type="l", lty=c(1,rep(2,10)), col=c("black", rep("blue",10)),lwd=c(2,rep(1,10))); 




X = ns(x,df=4, intercept=TRUE, Boundary.knots=c(0,2)); 

boot.10=boot.rq(X,y,tau=0.1, R=501);
boot.50=boot.rq(X,y,tau=0.5,  R=501, U=boot.10$U);
boot.90=boot.rq(X,y,tau=0.9, R=501, U=boot.10$U);

bootq.10=X%*%t(boot.10$B);
bootq.50=X%*%t(boot.50$B);
bootq.90=X%*%t(boot.90$B);

matplot(x, cbind(bootq.10,bootq.50, bootq.90), type="l"); 

bootNPS = (bootq.10 + bootq.90 - 2*bootq.50)/(bootq.90 - bootq.10); 
LS = apply(bootNPS,1,function(x) quantile(x,0.05)); 
MS = apply(bootNPS,1,function(x) quantile(x,0.5));
US = apply(bootNPS,1,function(x) quantile(x,0.95)); 
matplot(x,cbind(LS,MS,US),type="l",col="red",lty=2,lwd=2); 
points(x,NPS,type="l",col="black",lwd=2); 
abline(0,0,col="blue",lty=2); 

b0 = apply(bootNPS,2,mean); 
quantile(b0,c(0.01, 0.05,0.95,0.99)); 

b1 = apply(bootNPS,2,function(s) lm(s~x)$coef[2]);
quantile(b1,c(0.01, 0.05,0.95,0.99)); 

b2 = apply(bootNPS,2,function(s) lm(s~x + I(x^2))$coef[3]);
quantile(b2,c(0.01, 0.05,0.95,0.99)); 

require(mgcv); 
gam_fit = gam(list(y~s(x),~s(x), ~s(x), ~s(x)), family="shash"); 
summary(gam_fit); 