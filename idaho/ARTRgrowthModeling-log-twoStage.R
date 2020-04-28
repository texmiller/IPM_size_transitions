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
setwd("c:/repos/IPM_size_transitions/idaho"); 

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
drop_allD <- allD[-arb_small,]
plot(drop_allD$logarea.t0, drop_allD$logarea.t1) 

# condense treatments based on prior analysis
e = which(drop_allD$Treatment=="ControlModern");
drop_allD$Treatment[e] = "Control"; 
drop_allD = droplevels(drop_allD);  

########################################################################## 
## Pilot fits with log transformation, constant variance 
########################################################################## 
log_models <- list()
vf <- varExp(form =~ logarea.t0)
log_models[[1]] <- lmer(logarea.t1~ logarea.t0 + W.ARTR + W.HECO + W.POSE + W.PSSP+  W.allcov + W.allpts + Treatment + 
             Group+(logarea.t0|year),control=lmerControl(optimizer="bobyqa"),data=drop_allD,REML=F); 

log_models[[2]] <- lmer(logarea.t1~logarea.t0 + W.ARTR  + Treatment + 
             Group+(logarea.t0|year),control=lmerControl(optimizer="bobyqa"),data=drop_allD,REML = F)

log_models[[3]] <- lmer(logarea.t1~logarea.t0 + Treatment + 
             Group+(logarea.t0|year),control=lmerControl(optimizer="bobyqa"),data=drop_allD,REML = F)

log_models[[4]] <- lmer(logarea.t1~logarea.t0 + Treatment + 
             Group+(1|year),control=lmerControl(optimizer="bobyqa"),data=drop_allD,REML = F)
 			  
log_models[[5]] <- lmer(logarea.t1~logarea.t0 + Treatment + 
             Group+(0+logarea.t0|year),control=lmerControl(optimizer="bobyqa"),data=drop_allD,REML = F)

log_models[[6]] <- lmer(logarea.t1~logarea.t0 + Treatment + 
             (0+logarea.t0|year),control=lmerControl(optimizer="bobyqa"),data=drop_allD,REML = F)

aictab(log_models); 

########################################################################## 
## Use iterative re-weighting to fit with nonconstant variance,  
## then AIC for model selection.  
##########################################################################
for(mod in 1:6) {
for(rep in 1:25) {
	log_model = log_models[[mod]];
	fitted_vals = fitted(log_model);
	resids = residuals(log_model); 
	resid_model = lm(log(abs(resids))~fitted_vals); 
	new_weights = 0.5*(weights(log_model) + exp(-2*fitted(resid_model))); # cautious update 
	new_model <- update(log_model,weights=new_weights); 
	err = weights(log_model)-weights(new_model); err=sqrt(mean(err^2)); 
	cat(mod,rep,err,"\n") # check on convergence of estimated weights 
	log_models[[mod]]<-new_model; 
}}
aictab(log_models); 

######### For a fair AIC comparison, fit all models with the same weights 
aics = unlist(lapply(log_models,AIC)); best_model=which(aics==min(aics)); 
best_weights=weights(log_models[[best_model]]); 
for(mod in 1:6) {log_models[[mod]] <- update(log_models[[mod]],weights=best_weights)}
aictab(log_models); 

########## Make sure that the model choice hasn't changed
aics = unlist(lapply(log_models,AIC)); best_model=which(aics==min(aics)); 
best_weights=weights(log_models[[best_model]]); 
for(mod in 1:6) {log_models[[mod]] <- update(log_models[[mod]],weights=best_weights)}
aictab(log_models); 

######### Here's the best Gaussian model ########################################
aics = unlist(lapply(log_models,AIC)); best_model=which(aics==min(aics)); 
log_model = log_models[[best_model]]; 

summary(log_model); ## no Group effect, so re-fit 
log_model <- lmer(logarea.t1~logarea.t0 + Treatment + 
             (logarea.t0|year),control=lmerControl(optimizer="bobyqa"),data=drop_allD,
			 weights=weights(log_model),REML = F)
summary(log_model); 

######################################################################
# Interrogate the scaled residuals - Gaussian? NO. 
###################################################################### 
plot(fitted(log_model), sqrt(1/weights(log_model))); 
log_scaledResids = residuals(log_model)*sqrt(weights(log_model))
plot(fitted(log_model), log_scaledResids); 

qqPlot(log_scaledResids); # bad in just the lower tail

jarque.test(log_scaledResids) # normality test: FAILS, P < 0.001 
anscombe.test(log_scaledResids) # kurtosis: FAILS, P < 0.001 
agostino.test(log_scaledResids) # skewness: FAILS, P<0.001 

########################################################################
## See what family gamlss thinks the scaled residuals conform to: JSU. 
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

################################################################################################################
## Try fitting gamlss JSU to log size data. Use the fixed effect structure corresponding to the best lmer fit,
## and believe the rollaply diagnostics saying that skew and kurtosis are pretty much constant.
#################################################################################################################

#### Extract random effects from the pilot Gaussian fit, and use them as an offset 
coefs=fixef(log_model)
fixed.effect = coefs[1] + coefs[2]*drop_allD$logarea.t0;
fixed.effect[drop_allD$Treatment=="No_grass"]  =  fixed.effect[drop_allD$Treatment=="No_grass"] + coefs[3]; 
random.effect = fitted(log_model)-fixed.effect; 
drop_allD$random.effect = random.effect; 

### Pilot fit: variance depends on initial logarea. 
theFamily = "ST5"
fit_LSS <- gamlss(logarea.t1 ~ logarea.t0 +  Treatment + offset(random.effect), 
                  data=drop_allD, family=theFamily, method=RS(250), 
                  sigma.formula = ~logarea.t0, 
				  nu.formula = ~1, tau.formula = ~1)

### Iterative re-fit, with variance depending on fitted values 
drop_allD$fitted=fitted(fit_LSS);
for(k in 1:5) {
	fit_vals= fitted(fit_LSS);
	drop_allD$fitted <- 0.5*(drop_allD$fitted + fit_vals); # cautious update 
	fit_LSS <- gamlss(logarea.t1 ~ logarea.t0 + Treatment + offset(random.effect),
                  data=drop_allD, family= theFamily, method=RS(250),start.from=fit_LSS,
                  sigma.formula = ~fitted, 
				  nu.formula = ~1, tau.formula = ~1)
	new_fit = fitted(fit_LSS); 
	err = (new_fit - fit_vals)^2; 
	cat(k, mean(err)^0.5,"\n"); 
}				  

###############################################################################
## Well, it fit. does it describe the data well?
## Simulate data from best model
###############################################################################
n_sim <- 250
idaho_sim<-matrix(NA,nrow=nrow(drop_allD),ncol=n_sim)
for(i in 1:n_sim){
  idaho_sim[,i] <- rST5(n = nrow(drop_allD), 
                        mu = predict(fit_LSS), sigma = fit_LSS$sigma.fv,
                        nu = fit_LSS$nu.fv, tau = fit_LSS$tau.fv)
}

################################################################
#  Rollapply diagnostics applied to the fitted model 
#  Compare mean, sd, skewness, excess kurtosis vs. fitted value
################################################################
px = fitted(fit_LSS); e = order(px); px=px[e]; 

## Real data 
py = drop_allD$logarea.t1; py=py[e]; 
rollx=rollapply(px,100,mean,by=50);
rollmean = rollapply(py,100,mean,by=50); 
rollsd=rollapply(py,100,sd,by=50); 
rollskew=rollapply(py,100,skewness,by=50);
rollkurt=rollapply(py,100,kurtosis,by=50)/3-1;

thePlot=function(x,Y,xlab,ylab,gamma=1.4){
    pY= matrix(0,nrow(Y),ncol(Y))
	for(j in 1:ncol(Y)) {
	  fit=gam(Y[,j]~s(x),gamma=gamma,method="REML")
      pY[,j]=predict(fit,type="response"); 
    } 
	matplot(x,pY,col="grey50",type="o",pch=1,xlab=xlab,ylab=ylab);	
	points(x,pY[,1],col="red",cex=1.5,pch=16); 
}

par(mfrow=c(2,2),mar=c(5,5,2,1),cex.axis=1.3,cex.lab=1.3); 

Y = matrix(NA,length(rollx),n_sim); 
for(j in 1:n_sim) Y[,j] = rollapply(idaho_sim[e,j],100,mean,by=50); 
thePlot(rollx,cbind(rollmean,Y),xlab="Fitted value",ylab="Mean");

for(j in 1:n_sim) Y[,j] = rollapply(idaho_sim[e,j],100,sd,by=50); 
thePlot(rollx,cbind(rollsd,Y),xlab="Fitted value",ylab="Std Dev");

for(j in 1:n_sim) Y[,j] = rollapply(idaho_sim[e,j],100,skewness,by=50); 
thePlot(rollx,cbind(rollskew,Y),xlab="Fitted value",ylab="Skewness");

for(j in 1:n_sim) Y[,j] = rollapply(idaho_sim[e,j],100,kurtosis,by=50)/3-1; 
thePlot(rollx,cbind(rollkurt,Y),xlab="Fitted value",ylab="Excess Kurtosis");
