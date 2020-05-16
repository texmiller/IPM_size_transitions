##############################################################################
# A simple example of shrinkage based on Pseudoroegneria spicata, USSES Idaho. 
#
# Original: SPE May 2020 
#
##############################################################################

rm(list=ls(all=TRUE));
setwd("c:/repos/IPM_size_transitions/idaho"); 

require(lme4); require(tidyverse); require(maxLik); 
source("Utilities.R");

##############################################################
# 1. Generate artificial data for the example 
##############################################################

#############  Pull in real data 
allD<-read.csv(file="PSSP_growth_data.csv"); 
allD$year <- factor(allD$year); 

## drop out these arbitrarily small individuals
e1 = which(abs(allD$area.t0-0.25)<0.001); 
e2 = which(abs(allD$area.t1-0.25)<0.001); 
dropD = allD[-c(e1,e2),]; 

# Limit to historical control quadrats
e = which(dropD$Treatment=="Control");
dropD = droplevels(dropD[e,]);  
plot(dropD$logarea.t0, dropD$logarea.t1) 

################# Fit an SN1 distribution by maxLik 
# Note, sigma depends on fitted value. Sigma depending on initial
# size was much worse (Delta AIC of about 30)
LogLik=function(pars,new.size,U){
	pars1 = pars[1:ncol(U)]; pars2=pars[-(1:ncol(U))];
	mu = U%*%pars1;  
	dSN1(new.size, mu=mu, sigma=exp(pars2[1]+pars2[2]*mu),nu=pars2[3],log=TRUE)
}
U=model.matrix(~year + logarea.t0:year - 1, data=dropD);
start=c(runif(ncol(U)), rep(0,3))
out=maxLik(logLik=LogLik,start=start, new.size=dropD$logarea.t1,U=U,
			method="BHHH",control=list(iterlim=5000,printLevel=2),finalHessian=FALSE); 
out=maxLik(logLik=LogLik,start=out$estimate, new.size=dropD$logarea.t1,U=U,
			method="NM",control=list(iterlim=10000,printLevel=1),finalHessian=FALSE); 
out=maxLik(logLik=LogLik,start=out$estimate, new.size=dropD$logarea.t1,U=U,
			method="BHHH",control=list(iterlim=5000,printLevel=2),finalHessian=FALSE); 
pars=out$estimate; 

####################### Simulate a data set from the fitted model 
simData = data.frame(init.size = dropD$logarea.t0, year = factor(dropD$year));
pars1 = pars[1:ncol(U)]; pars2=pars[-(1:ncol(U))];
mu = U%*%pars1; nu=pars2[3]
simData$new.size = rSN1(nrow(simData),mu=mu,sigma=exp(pars2[1]+pars2[2]*mu),nu=pars2[3]);

plot(new.size~init.size,data=simData); 
truePars = pars; 
trueRanIntercept = truePars[1:22]-mean(truePars[1:22])
trueRanSlope = truePars[23:44]-mean(truePars[23:44])

##############################################################
# 2. Fit in a fixed-effects framework 
##############################################################
LogLik=function(pars,new.size,U){
	pars1 = pars[1:ncol(U)]; pars2=pars[-(1:ncol(U))];
	mu = U%*%pars1;  
	dSN1(new.size, mu=mu, sigma=exp(pars2[1]+pars2[2]*mu),nu=pars2[3],log=TRUE)
}

U=model.matrix(~year + init.size:year - 1, data=simData);
start=c(runif(ncol(U)), rep(0,3))
out=maxLik(logLik=LogLik,start=start, new.size=simData$new.size,U=U,
			method="BHHH",control=list(iterlim=5000,printLevel=2),finalHessian=FALSE); 
out=maxLik(logLik=LogLik,start=out$estimate, new.size=simData$new.size,U=U,
			method="NM",control=list(iterlim=10000,printLevel=2),finalHessian=FALSE); 
out=maxLik(logLik=LogLik,start=out$estimate, new.size=simData$new.size,U=U,
			method="BHHH",control=list(iterlim=5000,printLevel=2),finalHessian=TRUE); 

#####save results of ML fit  
coefs = out$estimate; names(coefs)<-colnames(U);
V = vcov(out); SEs = sqrt(diag(vcov(out))); 

##############################################################
# 3. Shrink the fixed-effects estimates 
##############################################################
T = length(unique(simData$year)); 
V1 = V[1:T,1:T]; V2 = V[(T+1):(2*T),(T+1):(2*T)]; 

fixed.fx = coefs[1:T]; fixed.fx = fixed.fx-mean(fixed.fx); 
var.hat = mean(fixed.fx^2) - mean(diag(V1)) + (sum(V1)-sum(diag(V1)))/(2*T*(T-1)); 
shrinkRanIntercept = fixed.fx*sqrt(var.hat/(var.hat + diag(V1)));

fixed.fx2 = coefs[(T+1):(2*T)]; fixed.fx2 = fixed.fx2-mean(fixed.fx2); 
var2.hat = mean(fixed.fx2^2) - mean(diag(V2)) + (sum(V2)-sum(diag(V2)))/(2*T*(T-1)); 
shrinkRanSlope = fixed.fx2*sqrt(var2.hat/(var2.hat + diag(V2))); 

graphics.off(); dev.new(height=4,width=8)
par(mfrow=c(1,2),bty="l",mar=c(4,4,1,1),mgp=c(2,1,0)); 
plot(trueRanIntercept,shrinkRanIntercept,type="p",xlab="True year-specific intercept",ylab="Estimated"); abline(0,1); 
plot(trueRanSlope,shrinkRanSlope,type="p",xlab="True year-specific slope",ylab="Estimated"); abline(0,1); 

save.image(file="SimpleShrinkageExample.Rdata"); 
