##############################################################################
# Growth kernel modeling for Artemisia tridentata, USSES Idaho. 
#
# Historical and modern data from by Adler et al. removals experiments paper.
# Modern data are just the Control and Grass-removal treaments.  
#
# Original: SPE April 2020
#
##############################################################################

rm(list=ls(all=TRUE));
## Steve wd
setwd("c:/repos/IPM_size_transitions/idaho"); 
## Tom wd
setwd("./idaho"); 

require(car); require(lme4); require(zoo); require(moments); require(mgcv); 
require(gamlss); require(gamlss.tr); require(AICcmodavg); require(lmerTest); 

source("Utilities.R");

### Spline scatterplot smoothing function
### Adjustable dof cost gamma, for use with rollapply outputs 
spline.scatter.smooth=function(x,y,gamma=10,show.quadratic=FALSE,...) {
  fit=gam(y~s(x),gamma=gamma,method="REML")
  plot(x,y,type="p",...);
  out=predict(fit,type="response"); 
  points(x,out,type="l",lwd=1)
  if(show.quadratic){
    fit2 = lm(y~x+I(x^2));
    points(x,fit2$fitted,type="l",col="red",lwd=2,lty=2);
  }
}

##############################################################
# 1. Read in the data
##############################################################
allD<-read.csv(file="ARTR_growth_data.csv"); 
allD$year <- factor(allD$year); 

######################################################################
# Tom's work: exploring the idaho data
######################################################################
library(tidyverse)
## drop out these arbitrarily small individuals
arb_small <- unique(c(which(allD$area.t0==0.25),which(allD$area.t1==0.25)))
dropD <- allD[-arb_small,]
plot(dropD$logarea.t0, dropD$logarea.t1) 

# condense treatments based on prior analysis
e = which(dropD$Treatment=="ControlModern");
dropD$Treatment[e] = "Control"; 
dropD = droplevels(dropD);  
table(dropD$Treatment) ## n=1040 vs 41 -- is it worth even fitting this effect?
					   ## Well, it was the main point of our 2018 paper on this experiment! 
					   ## But for IPM-building, maybe it's not important. 	

########################################################################## 
## Pilot fits with log transformation, constant variance 
########################################################################## 
log_models <- list()
log_models[[1]] <- lmer(logarea.t1~ logarea.t0 + W.ARTR + W.HECO + W.POSE + W.PSSP+  W.allcov + W.allpts + Treatment + 
             Group+(logarea.t0|year),control=lmerControl(optimizer="bobyqa"),data=dropD,REML=F); 

log_models[[2]] <- lmer(logarea.t1~logarea.t0 + W.ARTR  + Treatment + 
             Group+(logarea.t0|year),control=lmerControl(optimizer="bobyqa"),data=dropD,REML = F)

log_models[[3]] <- lmer(logarea.t1~logarea.t0 + Treatment + 
             Group+(logarea.t0|year),control=lmerControl(optimizer="bobyqa"),data=dropD,REML = F)

log_models[[4]] <- lmer(logarea.t1~logarea.t0 + Treatment + 
             Group+(1|year),control=lmerControl(optimizer="bobyqa"),data=dropD,REML = F)
 			  
log_models[[5]] <- lmer(logarea.t1~logarea.t0 + Treatment + 
             Group+(0+logarea.t0|year),control=lmerControl(optimizer="bobyqa"),data=dropD,REML = F)

log_models[[6]] <- lmer(logarea.t1~logarea.t0 + Treatment + 
             (0+logarea.t0|year),control=lmerControl(optimizer="bobyqa"),data=dropD,REML = F)

aictab(log_models); 

########################################################################## 
## Use iterative re-weighting to fit with nonconstant variance,  
## then AIC for model selection.  
##########################################################################

## NegLogLik function to fit variance model for residuals 
varPars = function(pars) {
    return(-sum(dnorm(resids, mean=0, sd=pars[1]*exp(pars[2]*fitted_vals),log=TRUE)))
}	

for(mod in 1:6) {
err = 1; rep=0; 
while(err > 0.001) {
	rep=rep+1; 
	log_model = log_models[[mod]];
	fitted_vals = fitted(log_model);resids = residuals(log_model); 
	out=optim(c(sd(resids),0),varPars,control=list(maxit=25000)); 
	pars=out$par; 
	new_sigma = pars[1]*exp(pars[2]*fitted_vals); new_weights = 1/((new_sigma)^2)
	new_weights = 0.5*(weights(log_model) + new_weights); # cautious update 
	new_model <- update(log_model,weights=new_weights); 
	err = weights(log_model)-weights(new_model); err=sqrt(mean(err^2)); 
	cat(mod,rep,err,"\n") # check on convergence of estimated weights 
	log_models[[mod]]<-new_model; 
}}
aictab(log_models); #model 3 is the winner 


######### For a fair AIC comparison, fit all models with the same weights 
aics = unlist(lapply(log_models,AIC)); best_model=which(aics==min(aics)); 
best_weights=weights(log_models[[best_model]]); 
for(mod in 1:6) {log_models[[mod]] <- update(log_models[[mod]],weights=best_weights)}
aictab(log_models); # still model 3, so done 

######### Here's the best Gaussian model ########################################
aics = unlist(lapply(log_models,AIC)); best_model=which(aics==min(aics)); 
log_model = log_models[[best_model]]; 

######################################################################
# Interrogate the scaled residuals - Gaussian? NO. 
###################################################################### 
plot(fitted(log_model), sqrt(1/weights(log_model))); 
log_scaledResids = residuals(log_model)*sqrt(weights(log_model))
plot(fitted(log_model), log_scaledResids); 

qqPlot(log_scaledResids); # bad in just the lower tail, but REALLY bad 

jarque.test(log_scaledResids) # normality test: FAILS, P < 0.001 
anscombe.test(log_scaledResids) # kurtosis: FAILS, P < 0.001 
agostino.test(log_scaledResids) # skewness: FAILS, P<0.001 

########################################################################
## See what family gamlss thinks the scaled residuals conform to: ST5. 
########################################################################
log_gamlss_fit <- fitDist(log_scaledResids,type="realline")
log_gamlss_fit$fits
log_gamlss_fit$failed

########################################################################
## Rollapply diagnostics on the scaled residuals 
########################################################################
px = fitted(log_model); py=log_scaledResids; 
e = order(px); px=px[e]; py=py[e];  

rollx=rollapply(px,50,mean,by=25);
rollmean = rollapply(py,50,mean,by=25); 
rollsd=rollapply(py,50,sd,by=25); 
rollkurt=rollapply(py,50,kurtosis,by=25);
rollskew=rollapply(py,50,skewness,by=25);

par(mfrow=c(2,2),mar=c(5,5,2,1),cex.axis=1.3,cex.lab=1.3); 
spline.scatter.smooth(rollx,rollmean,gamma=2,xlab="Fitted values",ylab="Mean");
spline.scatter.smooth(rollx,rollsd,gamma=2,xlab="Fitted values",ylab="Std Dev"); 
spline.scatter.smooth(rollx,rollskew,gamma=2,xlab="Fitted values",ylab="Skew"); 
abline(h=0,col="blue",lty=2) 
spline.scatter.smooth(rollx,rollkurt/3-1,gamma=2,xlab="Fitted values",ylab="Excess kurtosis"); 
abline(h=0,col="blue",lty=2)

