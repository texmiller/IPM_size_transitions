####################################################################################
## Exploring the feasibility of quantile-regression diagnostics for scaled residuals
## Fake 'scaled residuals' are generated from JP distribution
####################################################################################

library(quantreg); library(fda); library(sn); 
library(moments); library(mgcv); library(qgam); 
graphics.off(); 

root = "c:/repos"; # edit as needed 
setwd(root); setwd("IPM_size_transitions/gam practice"); 

source("JPfuns.R");
source("rqAICfuns.R");
source("../Diagnostics.R");  

### Generate data from RSJP distribution
##  JP with mean=0, var=1, in (lambda, tau) parameters. 
nx = 500; x = sort(2*rbeta(nx,3,3)); 
lambdas = 1 + (x-1); lambdas=lambdas-min(lambdas);
lambdas = (x-1); 
y = matrix(NA,nx,10); 
for(i in 1:nx) y[i,]=rRSJP(10,lambdas[i],tau=0)     
xydata = data.frame(x = x, y = y); 

## Quantiles of the distribution 
Y = matrix(NA,nx,3)
for(i in 1:nx) {
    out = qRSJP(c(0.1,0.5,0.9),lambda=lambdas[i], tau=0);
    Y[i,]=out; 
    if(i%%20==0) cat(i, "\n"); 
} 

dev.new(height=7,width=10); 
par(mfcol=c(2,3),mar=c(4,4,2,1),mgp=c(2,1,0),bty="l"); 

# plot a typical data set, and true quantiles 
plot(x,y[,1],xlab="Fitted",ylab="Residuals"); 
matpoints(x,Y, type="l",lty=2,lwd=2,col="red"); 
title(main="qGAM"); 

NPS_hat = matrix(NA,nx,10); 
for(k in 1:10) {
    z=y[,k]; xzdata = data.frame(x=x,z=z); 
    S.10 = qgam(z~s(x,k=5), data=xzdata,qu=0.05,argGam=list(gamma=2)); 
        q.10 = predict(S.10,newdata=xzdata); 
    S.50 = qgam(z~s(x,k=5),data=xzdata,qu=0.5,argGam=list(gamma=2)); 
        q.50 = predict(S.50,newdata=xzdata); 
    S.90 = qgam(z~s(x,k=5),data=xzdata,qu=0.95,argGam=list(gamma=2)); 
        q.90 = predict(S.90,newdata=xydata);
    matpoints(x, cbind(q.10,q.50,q.90), type="l", lty=1, lwd=1, col="grey50"); 
    NPS_hat[,k] = (q.10 + q.90 - 2*q.50)/(q.90 - q.10);
}

matpoints(x,Y, type="l",lty=1,lwd=2,col="red"); 

NPS = (Y[,1]+Y[,3]-2*Y[,2])/(Y[,3]-Y[,1]); 
matplot(x,cbind(NPS,NPS_hat),type="l", lty=c(1,rep(2,10)), col=c("black", rep("blue",10)),lwd=c(2,rep(1,10))); 
err = scale(t(NPS_hat),center=NPS,scale=FALSE); RMSE = round(sqrt(mean(err^2)),digits=3); 
title(main=paste("RMSE=", RMSE)); 



################## REPEAT WITH SPLINES, df=1 to 3, by AIC  
plot(x,y[,1],xlab="Fitted",ylab="Residuals"); 
matpoints(x,Y, type="l",lty=2,lwd=2,col="red"); 
title(main="B-spline with AIC 1.4"); 

NPS_hat = matrix(NA,nx,10); 
for(k in 1:10) {
    z=y[,k]; 
    f.10 <- rqAIC3(x,z,tau=0.1, L=0, U=2); q.10 = fitted(f.10); 
    f.50 <- rqAIC3(x,z,tau=0.5, L=0, U=2); q.50 = fitted(f.50); 
    f.90 <- rqAIC3(x,z,tau=0.9, L=0, U=2); q.90 = fitted(f.90); 
    matpoints(x, cbind(q.10,q.50,q.90), type="l", lty=1, lwd=1, col="grey50"); 
    NPS_hat[,k] = (q.10 + q.90 - 2*q.50)/(q.90 - q.10);
}

matpoints(x,Y, type="l",lty=1,lwd=2,col="red"); 

NPS = (Y[,1]+Y[,3]-2*Y[,2])/(Y[,3]-Y[,1]); 
matplot(x,cbind(NPS,NPS_hat),type="l", lty=c(1,rep(2,10)), col=c("black", rep("blue",10)),lwd=c(2,rep(1,10))); 
err = scale(t(NPS_hat),center=NPS,scale=FALSE); RMSE = round(sqrt(mean(err^2)),digits=3); 
title(main=paste("RMSE=", RMSE)); 

################## REPEAT WITH SMOOTHED BINNED ESTIMATES  
plot(x,y[,1],xlab="Fitted",ylab="Residuals"); 
matpoints(x,Y, type="l",lty=2,lwd=2,col="red"); 
title(main="Smoothing Binned Estimates"); 

NPS_hat = matrix(NA,nx,10); 
binmeans = binskews = matrix(NA,10,10); 
for(k in 1:10) {
    z=y[,k]; 
    out = rollMomentsNP(x,z,windows=6,smooth=FALSE,scaled=TRUE,xlab=NULL)
    binmeans[,k]=out$rollxmedian; 
    binskews[,k]=out$rollskew; 
    gamdata = data.frame(px = out$rollx, py = out$rollskew); 
    fit=gam(py~s(px,k=6),gamma=2,data=gamdata);
    NPS_hat[,k] = predict(fit,newdata=data.frame(px=x), type="response"); 
}

NPS = (Y[,1]+Y[,3]-2*Y[,2])/(Y[,3]-Y[,1]); 
matplot(x,cbind(NPS,NPS_hat),type="l", lty=c(1,rep(2,10)), col=c("black", rep("blue",10)),lwd=c(2,rep(1,10))); 
matpoints(binmeans,binskews,type="p",col="red",pch=1); 
err = scale(t(NPS_hat),center=NPS,scale=FALSE); RMSE = round(sqrt(mean(err^2)),digits=3); 
title(main=paste("RMSE=", RMSE)); 

