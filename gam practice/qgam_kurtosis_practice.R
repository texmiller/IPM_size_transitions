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
taus = 0.75*(x-1); 
y = matrix(NA,nx,10); 
for(i in 1:nx) y[i,]=rRSJP(10,tau=taus[i])     
xydata = data.frame(x = x, y = y); 

#fit = kde2d(x,y[,1],n=20, h=c(0.2,1.5*bandwidth.nrd(y[,1])),lims=c(0,2,-4,4));
#z2 = fit$z; 
#for(i in 1:20){
#    z2[i,]= dRSJP(fit$y,tau=0.5*(fit$x[i]-1)) * dbeta(fit$x[i]/2,3,3); 
#}    
#graphics.off(); 
#dev.new(width=8,height=5); par(mfrow=c(1,2)); 
#contour(fit$x,fit$y,fit$z/max(fit$z));
#contour(fit$x,fit$y,z2/max(z2)); 

## Quantiles of the distribution 
Y = matrix(NA,nx,5)
for(i in 1:nx) {
    out = qRSJP(c(0.05,0.25,0.5,0.75,0.95),lambda=0, tau=taus[i]);
    Y[i,]=out; 
    if(i%%20==0) cat(i, "\n"); 
}

## NP kurtosis of a Gaussian 
qN = qnorm(c(0.05,0.25,0.75,0.95))
KG = (qN[4]-qN[1])/(qN[3]-qN[2]); 

## NP kurtosis of the JP distribution 
NPK = (Y[,5]- Y[,1])/(Y[,4]-Y[,2]); 
NPK = NPK/KG - 1; 

dev.new(height=7,width=10); 
par(mfcol=c(2,4),mar=c(4,4,2,1),mgp=c(2,1,0),bty="l"); 

# plot a typical data set, and true quantiles 
plot(x,y[,1],xlab="Fitted",ylab="Residuals"); 
matpoints(x,Y, type="l",lty=2,lwd=2,col="red"); 
title(main="qGAM"); 

NPK_hat = matrix(NA,nx,10); 
for(k in 1:10) {
    z=y[,k]; xzdata = data.frame(x=x,z=z); 
    S.05 = qgam(z~s(x,k=4),data=xzdata,qu=0.05,argGam=list(gamma=2)); 
        q.05 = predict(S.05,newdata=xzdata); 
    S.25 = qgam(z~s(x,k=4),data=xzdata,qu=0.25,argGam=list(gamma=2)); 
        q.25 = predict(S.25,newdata=xzdata);         
    S.75 = qgam(z~s(x,k=4),data=xzdata,qu=0.75,argGam=list(gamma=2)); 
        q.75 = predict(S.75,newdata=xzdata); 
    S.95 = qgam(z~s(x,k=4),data=xzdata,qu=0.95,argGam=list(gamma=2)); 
        q.95 = predict(S.95,newdata=xydata);
    matpoints(x, cbind(q.05,q.25,q.75,q.95), type="l", lty=1, lwd=1, col="grey50"); 
    NPK_hat[,k] = ((q.95-q.05)/(q.75-q.25))/KG - 1; 
}

matpoints(x,Y, type="l",lty=1,lwd=2,col="red"); 
matplot(x,cbind(NPK,NPK_hat),type="l", lty=c(1,rep(2,10)), col=c("black", rep("blue",10)),lwd=c(2,rep(1,10))); 
err = scale(t(NPK_hat),center=NPK,scale=FALSE); RMSE = round(sqrt(mean(err^2)),digits=3); 
rug(x); 
title(main=paste("RMSE=", RMSE)); 


################## REPEAT WITH SPLINES, df=1 to 3, by AIC  
plot(x,y[,1],xlab="Fitted",ylab="Residuals"); 
matpoints(x,Y, type="l",lty=2,lwd=2,col="red"); 
title(main="B-spline Q-regression, AIC 1.4"); 

NPK_hat = matrix(NA,nx,10); 
for(k in 1:10) {
    z=y[,k]; 
    f.05 <- rqAIC3(x,z,tau=0.05, L=0, U=2); q.05 = fitted(f.05); 
    f.25 <- rqAIC3(x,z,tau=0.25, L=0, U=2); q.25 = fitted(f.25); 
    f.75 <- rqAIC3(x,z,tau=0.75, L=0, U=2); q.75 = fitted(f.75); 
    f.95 <- rqAIC3(x,z,tau=0.95, L=0, U=2); q.95 = fitted(f.95); 
    matpoints(x, cbind(q.05,q.25,q.75,q.95), type="l", lty=1, lwd=1, col="grey50"); 
    NPK_hat[,k] = ((q.95-q.05)/(q.75-q.25))/KG - 1; 
}

matpoints(x,Y, type="l",lty=1,lwd=2,col="red"); 
matplot(x,cbind(NPK,NPK_hat),type="l", lty=c(1,rep(2,10)), col=c("black", rep("blue",10)),lwd=c(2,rep(1,10))); 
err = scale(t(NPK_hat),center=NPK,scale=FALSE); RMSE = round(sqrt(mean(err^2)),digits=3); 
rug(x); 
title(main=paste("RMSE=", RMSE))


################## REPEAT WITH SMOOTHED BINNED ESTIMATES  
plot(x,y[,1],xlab="Fitted",ylab="Residuals"); 
matpoints(x,Y, type="l",lty=2,lwd=2,col="red"); 
title(main="Smoothing Binned Estimates"); 

NPK_hat = matrix(NA,nx,10); 
for(k in 1:10) {
    z=y[,k]; 
    out = rollMomentsNP(x,z,windows=6,smooth=FALSE,scaled=TRUE,xlab=NULL)
    gamdata = data.frame(px = out$rollx, py = out$rollkurt); 
    fit=gam(py~s(px,k=6),gamma=2,data=gamdata);
    NPK_hat[,k] = predict(fit,newdata=data.frame(px=x), type="response"); 
}

matpoints(x,Y, type="l",lty=1,lwd=2,col="red"); 
matplot(x,cbind(NPK,NPK_hat),type="l", lty=c(1,rep(2,10)), col=c("black", rep("blue",10)),lwd=c(2,rep(1,10))); 
err = scale(t(NPK_hat),center=NPK,scale=FALSE); RMSE = round(sqrt(mean(err^2)),digits=3); 
rug(x); 
title(main=paste("RMSE=", RMSE))


################## REPEAT WITH SPLINE DENSITY ESTIMATE (gss package); 
require(gss); 
plot(x,y[,1],xlab="Fitted",ylab="Residuals"); 
matpoints(x,Y, type="l",lty=2,lwd=2,col="red"); 
title(main="GSS Spline Density Estimate"); 

NPK_hat = matrix(NA,nx,10); 
for(k in 1:10) {
    z=y[,k]; X = data.frame(x=x,y=z); 
    fit <- sscden1(~x + y + x*y, ~y,data=X)
    quan <- qsscden(fit,c(.05,.25,.5,.75,.95),data.frame(x=X$x)); 
    quan = t(quan); 
    matpoints(x,quan, type="l",lty=1, lwd=1, col="grey50"); 
    NPK_hat[,k] = ((quan[,5]-quan[,1])/(quan[,4]-quan[,2]))/KG - 1; 
}

matpoints(x,Y, type="l",lty=1,lwd=2,col="red"); 
matplot(x,cbind(NPK,NPK_hat),type="l", lty=c(1,rep(2,10)), col=c("black", rep("blue",10)),lwd=c(2,rep(1,10))); 
err = scale(t(NPK_hat),center=NPK,scale=FALSE); RMSE = round(sqrt(mean(err^2)),digits=3); 
rug(x); 
title(main=paste("RMSE=", RMSE))

