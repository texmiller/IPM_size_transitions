##############################################################################
# Growth kernel modeling for Pseudoroegneria spicata, USSES Idaho. 
# Compares f(initial size) vs. f(fitted value) specifications for SD, 
# skewness, and kurtosis parameters in JSU distribution 
#
# Original: SPE May 
#
##############################################################################

rm(list=ls(all=TRUE));
setwd("c:/repos/IPM_size_transitions/idaho"); 

require(car); require(lme4); require(zoo); require(moments); require(mgcv); 
require(gamlss); require(gamlss.tr); require(AICcmodavg); 
require(lmerTest); require(tidyverse); require(maxLik); 

source("Utilities.R");

##############################################################
# 1. Read in the data
##############################################################
allD<-read.csv(file="PSSP_growth_data.csv"); 
allD$year <- factor(allD$year); 

######################################################################
# Clean up the data
######################################################################
library(tidyverse)
## drop out these arbitrarily small individuals
e = which(abs(allD$area.t0-0.25)<0.001); 
dropD = allD[-e,]; 
e = which(abs(dropD$area.t1-0.25)<0.001); 
dropD = dropD[-e,]; 
plot(dropD$logarea.t0, dropD$logarea.t1) 

# condense treatments based on prior analysis: no difference between
# historical and modern control quadrats 
e = which(dropD$Treatment=="ControlModern");
dropD$Treatment[e] = "Control"; 
dropD = droplevels(dropD);  

################################################################################################################
# Fitting gamlss JSU to log size data, with (sigma,tau,nu) depending on the linear predictor for mu
################################################################################################################

# Model matrix for the fixed and random effects, specified so that each year gets its own coefficient
# rather than fitting contrasts against a baseline year (the default in R's regression functions) 
U=model.matrix(~year + logarea.t0:year + I(logarea.t0^2) +  
	W.ARTR + I(W.HECO + W.POSE) + W.PSSP+  I(W.allcov + W.allpts) + Group + Treatment - 1, data=dropD);

LogLik=function(pars,response,U){
	pars1 = pars[1:ncol(U)]; pars2=pars[-(1:ncol(U))];
	mu = U%*%pars1;  
	val = dJSU(response, mu=mu,
	sigma=exp(pars2[1]+pars2[2]*mu),
	nu = pars2[3]+pars2[4]*mu,
	tau = exp(pars2[5]+pars2[6]*mu), log=TRUE)
	return(val); 
}

############ Fit the data. Paranoid as usual about convergence
coefs = list(5); LL=numeric(5);  
for(j in 1:5) {
	out=maxLik(logLik=LogLik,start=c(runif(ncol(U)), rep(0,6)), response=dropD$logarea.t1,U=U,
		method="BHHH",control=list(iterlim=5000,printLevel=2),finalHessian=FALSE); 

	out=maxLik(logLik=LogLik,start=out$estimate,response=dropD$logarea.t1,U=U,
		method="NM",control=list(iterlim=5000,printLevel=1),finalHessian=FALSE); 

	out=maxLik(logLik=LogLik,start=out$estimate,response=dropD$logarea.t1,U=U,
		method="BHHH",control=list(iterlim=5000,printLevel=2),finalHessian=FALSE); 

	coefs[[j]] = out$estimate; LL[j] = out$maximum;
	cat(j, "#--------------------------------------#",out$maximum,"\n"); 
}

j = min(which(LL==max(LL)));
GLM_style_fit=maxLik(logLik=LogLik,start=coefs[[j]],response=dropD$logarea.t1,U=U,
		method="BHHH",control=list(iterlim=5000,printLevel=2),finalHessian=TRUE); 

################################################################################################################
# Fitting gamlss JSU to log size data, with (sigma,tau,nu) depending on initial size 
################################################################################################################
LogLik2=function(pars,response,init,U){
	pars1 = pars[1:ncol(U)]; pars2=pars[-(1:ncol(U))];
	mu = U%*%pars1;  
	val = dJSU(response, mu=mu,
	sigma=exp(pars2[1]+pars2[2]*init),	
	nu = pars2[3]+pars2[4]*init,
	tau = exp(pars2[5]+pars2[6]*init), log=TRUE)
	return(val); 
}

############ Fit the data. Paranoid as usual about convergence
coefs = list(5); LL=numeric(5);  
for(j in 1:5) {
	out=maxLik(logLik=LogLik2,start=c(runif(ncol(U)), rep(0,6)), response=dropD$logarea.t1, init=dropD$logarea.t0, U=U,
		method="BHHH",control=list(iterlim=5000,printLevel=2),finalHessian=FALSE); 

	out=maxLik(logLik=LogLik2,start=out$estimate,response=dropD$logarea.t1,init=dropD$logarea.t0, U=U,
		method="NM",control=list(iterlim=5000,printLevel=1),finalHessian=FALSE); 

	out=maxLik(logLik=LogLik2,start=out$estimate,response=dropD$logarea.t1,init=dropD$logarea.t0, U=U,
		method="BHHH",control=list(iterlim=5000,printLevel=2),finalHessian=FALSE); 

	coefs[[j]] = out$estimate; LL[j] = out$maximum;
	cat(j, "#--------------------------------------#",out$maximum,"\n"); 
}

j = min(which(LL==max(LL)));
lmer_style_fit=maxLik(logLik=LogLik2,start=coefs[[j]],response=dropD$logarea.t1,init=dropD$logarea.t0, U=U,
		method="BHHH",control=list(iterlim=5000,printLevel=2),finalHessian=TRUE); 

AIC(lmer_style_fit)
AIC(GLM_style_fit); 













