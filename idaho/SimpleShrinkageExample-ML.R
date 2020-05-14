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

set.seed(1774); 

##############################################################
# 1. Generate artificial data for the example 
##############################################################

#############  Pull in real data 
allD<-read.csv(file="PSSP_growth_data.csv"); 
allD$year <- factor(allD$year); 

## drop out these arbitrarily small individuals
e = which(abs(allD$area.t0-0.25)<0.001); 
dropD = allD[-e,]; 
e = which(abs(dropD$area.t1-0.25)<0.001); 
dropD = dropD[-e,]; 

# Limit to control quadrats, historical and modern
e = which(dropD$Treatment=="ControlModern");
dropD$Treatment[e] = "Control"; 
dropD = subset(dropD,Treatment=="Control"); 
dropD = droplevels(dropD);  
plot(dropD$logarea.t0, dropD$logarea.t1) 


## Fit with nonconstant variance
log_model <- lmer(logarea.t1 ~ logarea.t0 + (logarea.t0|year),
			control=lmerControl(optimizer="bobyqa"),data=dropD,REML=TRUE); 

varPars = function(pars) {
    return(-sum(dnorm(resids, mean=0, sd=exp(pars[1] + pars[2]*fitted_vals),log=TRUE)))
}	

err = 1; rep=0; while(err > 0.000001) {
	rep=rep+1; 
	fitted_vals = fitted(log_model); resids = residuals(log_model); 
	out=optim(c(sd(resids),0),varPars,method="BFGS",control=list(maxit=5000)); 
	pars=out$par; 
	new_sigma = exp(pars[1] + pars[2]*fitted_vals); new_weights = 1/((new_sigma)^2)
	new_weights = 0.5*(weights(log_model) + new_weights);  
	new_model <- update(log_model,weights=new_weights); 
	err = weights(log_model)-weights(new_model); err=sqrt(mean(err^2)); 
	cat(rep,err,"\n") # check on convergence of estimated weights 
	log_model <- new_model; 
}

e  = sample(1:nrow(dropD),3000,replace=FALSE)
simData = data.frame(init.size = dropD$logarea.t0[e], year = dropD$year[e]);
simData$year=factor(simData$year); 
mu = fitted(log_model)[e]; 
simData$new.size = rnorm(3000,mean=mu,sd=exp(pars[1] + pars[2]*mu));

plot(new.size~init.size,data=simData); 

##############################################################
# 2. Fit in a fixed-effects framework 
##############################################################

### Log Likelihood function as required by maxLik 
LogLik=function(pars,new.size,U){
	pars1 = pars[1:ncol(U)]; pars2=pars[-(1:ncol(U))];
	mu = U%*%pars1;  
	dnorm(new.size, mean=mu, sd=exp(pars2[1]+pars2[2]*mu),log=TRUE)
}

U=model.matrix(~year + init.size:year - 1, data=simData);

start=c(runif(ncol(U)), rep(0,2))
out=maxLik(logLik=LogLik,start=start, new.size=simData$new.size,U=U,
			method="BHHH",control=list(iterlim=5000,printLevel=2),finalHessian=FALSE); 
out=maxLik(logLik=LogLik,start=out$estimate, new.size=simData$new.size,U=U,
			method="NM",control=list(iterlim=10000,printLevel=2),finalHessian=FALSE); 
out=maxLik(logLik=LogLik,start=start, new.size=simData$new.size,U=U,
			method="BHHH",control=list(iterlim=5000,printLevel=2),finalHessian=TRUE); 

#####save results of ML fit  
coefs = out$estimate; SEs = sqrt(diag(vcov(out))); names(coefs)<-colnames(U);

##############################################################
# 3. Shrink the fixed-effects estimates 
##############################################################

fixed.fx = coefs[1:30]; fixed.fx = fixed.fx-mean(fixed.fx); 
fixed.se = SEs[1:30]; 
var.hat = mean(fixed.fx^2)-mean(fixed.se^2)
shrinkRanIntercept = fixed.fx*sqrt(var.hat/(var.hat + fixed.se^2));

fixed.fx2 = coefs[31:60]; fixed.fx2 = fixed.fx2-mean(fixed.fx2); 
fixed.se2 = SEs[31:60]; 
var2.hat = mean(fixed.fx2^2)-mean(fixed.se2^2)
shrinkRanSlope = fixed.fx2*sqrt(var2.hat/(var2.hat + fixed.se2^2)); 

graphics.off(); dev.new(height=4,width=8)
par(mfrow=c(1,2),bty="l",mar=c(4,4,1,1),mgp=c(2,1,0)); 
plot(ranef(log_model)$year[,1],shrinkRanIntercept,type="p",xlab="True year-specific intercept",ylab="Estimated"); abline(0,1); 
plot(ranef(log_model)$year[,2],shrinkRanSlope,type="p",xlab="True year-specific slope",ylab="Estimated"); abline(0,1); 

