##############################################################################
# Growth kernel modeling for Pseudoroegneria spicata, USSES Idaho. 
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
require(gamlss); require(gamlss.tr); require(AICcmodavg); 
require(lmerTest); 

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
allD<-read.csv(file="PSSP_growth_data.csv"); 
allD$year <- factor(allD$year); 

######################################################################
# Tom's work: exploring the idaho data
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

########################################################################## 
## Pilot fits with log transformation, constant variance 
########################################################################## 
log_models <- list()
log_models[[1]] <- lmer(logarea.t1~ logarea.t0 + W.ARTR + W.HECO + W.POSE + W.PSSP+  W.allcov + W.allpts + Treatment + 
             Group+(logarea.t0|year),control=lmerControl(optimizer="bobyqa"),data=dropD,REML=F); 

## coefficients on HECO and POSE are nearly identical, so group them: deltaAIC = 2, i.e. no change at all in logLik.   
log_models[[2]] <- lmer(logarea.t1~ logarea.t0 + W.ARTR + I(W.HECO + W.POSE) + W.PSSP+  W.allcov + W.allpts + Treatment + 
             Group+(logarea.t0|year),control=lmerControl(optimizer="bobyqa"),data=dropD,REML=F); 

## W.allpts is non-significant. Consider two options: drop, group with all other cover 
log_models[[3]] <- lmer(logarea.t1~ logarea.t0 + W.ARTR + I(W.HECO + W.POSE) + W.PSSP+  W.allcov + Treatment + 
             Group+(logarea.t0|year),control=lmerControl(optimizer="bobyqa"),data=dropD,REML=F);  			 

## Based on AIC, the winner by a hair is to group all heterospecific grasses, and group all cover besides the 'big 4'
log_models[[4]] <- lmer(logarea.t1~ logarea.t0 + W.ARTR + I(W.HECO + W.POSE) + W.PSSP+  I(W.allcov + W.allpts) + Treatment + 
             Group + (logarea.t0|year), control=lmerControl(optimizer="bobyqa"),data=dropD,REML=F); 	

########################################################################## 
## Use iterative re-weighting to fit with nonconstant variance,  
## then AIC for model selection.  
##########################################################################

## NegLogLik function to fit variance model for residuals 
varPars = function(pars) {
    return(-sum(dnorm(resids, mean=0, sd=pars[1]*exp(pars[2]*fitted_vals),log=TRUE)))
}	

for(mod in 1:4) {
err = 1; rep=0; 
while(err > 0.000001) {
	rep=rep+1; 
	log_model = log_models[[mod]];
	fitted_vals = fitted(log_model);resids = residuals(log_model); 
	out=optim(c(sd(resids),0),varPars,control=list(maxit=5000)); 
	pars=out$par; 
	new_sigma = pars[1]*exp(pars[2]*fitted_vals); new_weights = 1/((new_sigma)^2)
	new_weights = 0.5*(weights(log_model) + new_weights); # cautious update 
	new_model <- update(log_model,weights=new_weights); 
	err = weights(log_model)-weights(new_model); err=sqrt(mean(err^2)); 
	cat(mod,rep,err,"\n") # check on convergence of estimated weights 
	log_models[[mod]]<-new_model; 
}}
aictab(log_models); # model 4 is the winner, by a hair 

######### For a fair AIC comparison, fit all models with the same weights 
aics = unlist(lapply(log_models,AIC)); best_model=which(aics==min(aics)); 
best_weights=weights(log_models[[best_model]]); 
for(mod in 1:4) {log_models[[mod]] <- update(log_models[[mod]],weights=best_weights)}
aictab(log_models); # still model 4, by a hair 

######### Here's the best Gaussian model ########################################
aics = unlist(lapply(log_models,AIC)); best_model=which(aics==min(aics)); 
log_model = log_models[[best_model]]; 
summary(log_model); 

######################################################################
# Interrogate the scaled residuals - Gaussian? NO. 
###################################################################### 
plot(fitted(log_model), 1/sqrt(weights(log_model)),xlab="Fitted",ylab="Estimated residual Std Dev"); 
log_scaledResids = residuals(log_model)*sqrt(weights(log_model))
plot(fitted(log_model), log_scaledResids); 

qqPlot(log_scaledResids); # really bad in lower tail, not too bad in upper 

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
## import the random effects, and believe rollaply diagnostics saying skew and kurtosis are nearly linear.
#################################################################################################################

#### Extract random effects from the pilot Gaussian fit, and use them as an offset 
coefs=fixef(log_model)
fixed.effect = coefs[1] + coefs[2]*dropD$logarea.t0;
fixed.effect[dropD$Treatment=="No_shrub"]  =  fixed.effect[dropD$Treatment=="No_shrub"] + coefs[3]; 
random.effect = fitted(log_model)-fixed.effect; 
dropD$random.effect = random.effect; 

### Pilot fit: variance depends on initial logarea. 
theFamily = "JSU"
fit_LSS <- gamlss(logarea.t1 ~ logarea.t0 +  Treatment + offset(random.effect), 
                  data=dropD, family=theFamily, method=RS(250), 
                  sigma.formula = ~logarea.t0, 
				  nu.formula = ~logarea.t0, tau.formula = ~logarea.t0)

### Iterative re-fit, with variance depending on fitted values 
dropD$fitted=fitted(fit_LSS);
for(k in 1:15) {
	fit_vals= fitted(fit_LSS);
	dropD$fitted <- 0.5*(dropD$fitted + fit_vals); # cautious update 
	fit_LSS <- gamlss(logarea.t1 ~ logarea.t0 + I(logarea.t0)^2 + Treatment + offset(random.effect),
                  data=dropD, family= theFamily, method=RS(250),start.from=fit_LSS,
                  sigma.formula = ~ fitted, 
				  nu.formula = ~fitted, tau.formula = ~fitted)
	new_fit = fitted(fit_LSS); 
	err = (new_fit - fit_vals)^2; 
	cat(k, mean(err)^0.5,"\n"); 
}				  


###############################################################################
## Well, it fit. does it describe the data well?
## Simulate data from best model
###############################################################################
n_sim <- 250
idaho_sim<-matrix(NA,nrow=nrow(dropD),ncol=n_sim)
for(i in 1:n_sim){
  idaho_sim[,i] <- rJSU(n = nrow(dropD), 
                        mu = predict(fit_LSS), sigma = fit_LSS$sigma.fv,
                        nu = fit_LSS$nu.fv, tau = fit_LSS$tau.fv)
}
e = idaho_sim < log(0.5); 
idaho_sim[e]  = log(0.5); 

################################################################
#  Rollapply diagnostics applied to the fitted model 
#  Compare mean, sd, skewness, excess kurtosis vs. fitted value
################################################################
px = fitted(fit_LSS); e = order(px); px=px[e]; 

## Real data 
py = dropD$logarea.t1; py=py[e]; 
rollx=rollapply(px,100,mean,by=50);
rollmean = rollapply(py,100,mean,by=50); 
rollsd=rollapply(py,100,sd,by=50); 
rollskew=rollapply(py,100,skewness,by=50);
rollkurt=rollapply(py,100,kurtosis,by=50)/3-1;

thePlot=function(x,Y,xlab,ylab,gamma=1.4,...){
    pY= matrix(0,nrow(Y),ncol(Y))
	for(j in 1:ncol(Y)) {
	  fit=gam(Y[,j]~s(x),gamma=gamma,method="REML")
      pY[,j]=predict(fit,type="response"); 
    } 
	matplot(x,pY,col="grey50",type="o",pch=1,xlab=xlab,ylab=ylab,...);	
	points(x,pY[,1],col="red",cex=1.5,pch=16); 
}

par(mfrow=c(2,2),mar=c(5,5,2,1),cex.axis=1.3,cex.lab=1.3); 

Y = matrix(NA,length(rollx),n_sim); 
for(j in 1:n_sim) Y[,j] = rollapply(idaho_sim[e,j],100,mean,by=50); 
thePlot(rollx,cbind(rollmean,Y),xlab="Fitted value",ylab="Mean");

for(j in 1:n_sim) Y[,j] = rollapply(idaho_sim[e,j],100,sd,by=50); 
thePlot(rollx,cbind(rollsd,Y),xlab="Fitted value",ylab="Std Dev",ylim=c(0,max(cbind(rollsd,Y))));

for(j in 1:n_sim) Y[,j] = rollapply(idaho_sim[e,j],100,skewness,by=50); 
thePlot(rollx,cbind(rollskew,Y),xlab="Fitted value",ylab="Skewness");

for(j in 1:n_sim) Y[,j] = rollapply(idaho_sim[e,j],100,kurtosis,by=50)/3-1; 
thePlot(rollx,cbind(rollkurt,Y),xlab="Fitted value",ylab="Excess Kurtosis");

#######################################################################################
# Now do the Shrinkage approach, and compare to the lmer random effects
# Fit everything as fixed effects; used estimated s.e.'s to estimate the
# mixing sigma; then compute BLUPS for the unobserved true effects based 
# on estimated sigma and estimated s.e.'s 
#######################################################################################

### Pilot fit: variance depends on initial logarea. 
theFamily = "JSU"; XD=dropD; 

# The only way I know to avoid making fixed-effect year and logarea:year coefficients
# be contrasts with a baseline, so they align with the way lmer works. 
U=model.matrix(~year + logarea.t0:year - 1, data=XD); 

# fit using the model.matrix, without an intercept 
fit_LSS <- gamlss(logarea.t1 ~ U + I(logarea.t0^2) + Treatment -1, 
                  data=XD, family=theFamily, method=RS(250), 
                  sigma.formula = ~logarea.t0, 
				  nu.formula = ~logarea.t0, tau.formula = ~logarea.t0)

### Iterative re-fit, with variance depending on fitted values 
dropD$fitted=fitted(fit_LSS);
for(k in 1:15) {
	fit_vals= fitted(fit_LSS);
	XD$fitted <- 0.5*(dropD$fitted + fit_vals); # cautious update 
	fit_LSS <- gamlss(logarea.t1 ~  U + I(logarea.t0^2) + Treatment -1,
                  data=XD, family= theFamily, method=RS(250),start.from=fit_LSS,
                  sigma.formula = ~fitted, 
				  nu.formula = ~fitted, tau.formula = ~fitted)
	new_fit = fitted(fit_LSS); 
	err = (new_fit - fit_vals)^2; 
	cat(k, mean(err)^0.5,"\n"); 
}				  

S = summary(fit_LSS); 

par(mfrow=c(2,2),mar=c(4,4,1,1),mgp=c(2,1,0),bty="l"); 
# shrinkage random effects for (1|year) 
fixed.fx = S[1:30,1]; fixed.fx = fixed.fx-mean(fixed.fx); 
fixed.se = S[1:30,2]; 
sigma2.hat = mean(fixed.fx^2)-mean(fixed.se^2)
# BLUP based on sigma.hat and estimated s.e.
shrunk.fx = fixed.fx*(sigma2.hat/(sigma2.hat + fixed.se^2)); 

# lmer random effects for (1|year) 
ran.fx = ranef(log_model)[[1]]; ran.fx = ran.fx[,1]; 
sd(fixed.fx); sd(shrunk.fx); sd(ran.fx); 
plot(ran.fx,shrunk.fx,xlab="lmer year random effects",ylab="Shrunk year fixed effects",pch=1,lwd=2,cex=1.2); 
abline(0,1,col="blue",lty=2); 

# shrinkage random effects for (logarea.t0|year) 
fixed.fx2 = S[31:60,1]; fixed.fx2 = fixed.fx2-mean(fixed.fx2); 
fixed.se = S[31:60,2]; 
sigma2.hat = mean(fixed.fx^2)-mean(fixed.se^2)
# BLUP based on sigma.hat and estimated s.e.
shrunk.fx2 = fixed.fx2*(sigma2.hat/(sigma2.hat + fixed.se^2));

# lmer random effects for (1|year) 
ran.fx2 = ranef(log_model)[[1]]; ran.fx2 = ran.fx2[,2]; 
sd(fixed.fx2); sd(shrunk.fx2); sd(ran.fx2); 
plot(ran.fx2,shrunk.fx2,xlab="lmer size:year random effects",ylab="Shrunk size:year fixed effects",pch=1,lwd=2,cex=1.2); 
abline(0,1,col="blue",lty=2); 

matplot(cbind(ran.fx,shrunk.fx),cbind(ran.fx2,shrunk.fx2),col=c("blue","red"),pch=c(1,2),cex=1.4,lwd=2,xlab="Intercept Year-effect",
ylab="Slope year-effect"); 
legend("topright", legend=c("lmer","Shrinkage"),col=c("blue","red"),cex=1.4,pch=c(1,2),bty="n"); 


##################################################################
#  Rollapply diagnostics applied to the model fitted via Shrinkage
#  Compare mean, sd, skewness, excess kurtosis vs. fitted value
##################################################################
px = fitted(fit_LSS); e = order(px); px=px[e]; 

## Real data 
py = dropD$logarea.t1; py=py[e]; 
rollx=rollapply(px,100,mean,by=50);
rollmean = rollapply(py,100,mean,by=50); 
rollsd=rollapply(py,100,sd,by=50); 
rollskew=rollapply(py,100,skewness,by=50);
rollkurt=rollapply(py,100,kurtosis,by=50)/3-1;

thePlot=function(x,Y,xlab,ylab,gamma=1.4,...){
    pY= matrix(0,nrow(Y),ncol(Y))
	for(j in 1:ncol(Y)) {
	  fit=gam(Y[,j]~s(x),gamma=gamma,method="REML")
      pY[,j]=predict(fit,type="response"); 
    } 
	matplot(x,pY,col="grey50",type="o",pch=1,xlab=xlab,ylab=ylab,...);	
	points(x,pY[,1],col="red",cex=1.5,pch=16); 
}

par(mfrow=c(2,2),mar=c(5,5,2,1),cex.axis=1.3,cex.lab=1.3); 

Y = matrix(NA,length(rollx),n_sim); 
for(j in 1:n_sim) Y[,j] = rollapply(idaho_sim[e,j],100,mean,by=50); 
thePlot(rollx,cbind(rollmean,Y),xlab="Fitted value",ylab="Mean");

for(j in 1:n_sim) Y[,j] = rollapply(idaho_sim[e,j],100,sd,by=50); 
thePlot(rollx,cbind(rollsd,Y),xlab="Fitted value",ylab="Std Dev",ylim=c(0,max(cbind(rollsd,Y))));

for(j in 1:n_sim) Y[,j] = rollapply(idaho_sim[e,j],100,skewness,by=50); 
thePlot(rollx,cbind(rollskew,Y),xlab="Fitted value",ylab="Skewness");

for(j in 1:n_sim) Y[,j] = rollapply(idaho_sim[e,j],100,kurtosis,by=50)/3-1; 
thePlot(rollx,cbind(rollkurt,Y),xlab="Fitted value",ylab="Excess Kurtosis");
