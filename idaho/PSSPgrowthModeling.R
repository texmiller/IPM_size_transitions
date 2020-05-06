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
# Try fitting gamlss JSU to log size data. Use the fixed effect structure corresponding to the best lmer fit,
# and believe rollaply diagnostics saying skew is linear, and kurtosis is quadratic. 
# 
# The random effects terms in the lmer fit are fitted here as fixed effects, then adjusted by shrinkage.  
# Having fitted (1|year) and (0 + logarea.t0|yer) as fixed effects, we have year-specific coefficients
# and their estimated standard errors. We then compute BLUPS for the unobserved true effects based 
# on estimated mixing sigma and estimated s.e.'s.  
#################################################################################################################

theFamily = JSU(tau.link="identity"); XD=dropD; 

# This is the one way I know to avoid the fixed-effect year and logarea:year coefficients
# being reported as contrasts with a baseline. Doing it like this gives us one value for
# each year, corresponding to the ranef() values reported by lmer. 

U=model.matrix(~year + logarea.t0:year - 1, data=XD); 
# fit using the model.matrix, without an intercept 
fit_LSS <- gamlss(logarea.t1 ~ U + I(logarea.t0^2) + Treatment - 1, 
                  data=XD, family=theFamily, method=RS(250), 
                  sigma.formula = ~logarea.t0, 
				  nu.formula = ~logarea.t0, tau.formula = ~logarea.t0)

### Iterative re-fit, with variance depending on fitted values 
dropD$fitted=fitted(fit_LSS); err =1; k=0; 
while(err>0.0001) {
	fit_vals= fitted(fit_LSS);
	XD$fitted <- 0.5*(dropD$fitted + fit_vals); # cautious update 
	fit_LSS <- gamlss(logarea.t1 ~  U + I(logarea.t0^2) + Treatment -1,
                  data=XD, family= theFamily, method=RS(250),start.from=fit_LSS,
                  sigma.formula = ~ fitted, 
				  nu.formula = ~ fitted, tau.formula = ~ fitted )
	new_fit = fitted(fit_LSS); 
	err = (new_fit - fit_vals)^2; err = mean(err)^0.5; k=k+1; 
	cat(k, err, "\n"); 
}				  
dropD$fitted = fitted(fit_LSS); 

######################################################################
#  Compare lmer and Shrinkage random effects
######################################################################

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
fixed.se2 = S[31:60,2]; 
sigma2.hat = mean(fixed.fx2^2)-mean(fixed.se2^2)
# BLUP based on sigma.hat and estimated s.e.
shrunk.fx2 = fixed.fx2*(sigma2.hat/(sigma2.hat + fixed.se2^2));

# lmer random effects for (1|year) 
ran.fx2 = ranef(log_model)[[1]]; ran.fx2 = ran.fx2[,2]; 
sd(fixed.fx2); sd(shrunk.fx2); sd(ran.fx2); 
plot(ran.fx2,shrunk.fx2,xlab="lmer size:year random effects",ylab="Shrunk size:year fixed effects",pch=1,lwd=2,cex=1.2); 
abline(0,1,col="blue",lty=2); 

matplot(cbind(ran.fx,shrunk.fx),cbind(ran.fx2,shrunk.fx2),col=c("blue","red"),pch=c(1,2),cex=1.4,lwd=2,xlab="Intercept Year-effect",
ylab="Slope year-effect"); 
legend("topright", legend=c("lmer","Shrinkage"),col=c("blue","red"),cex=1.4,pch=c(1,2),bty="n"); 

save.image(file="PSSPgrowthModels.Rdata"); 

#############################################################################
#  Binned data diagnostics applied to the model fitted via Shrinkage. 
#  Compare mean, sd, skewness, excess kurtosis as a function of fitted value.
#  TO DO: functional programming instead of copy/paste/edit for each panel. 
#############################################################################

# Simulate data from best model
n_sim <- 500
idaho_sim<-matrix(NA,nrow=nrow(dropD),ncol=n_sim)
for(i in 1:n_sim){
  idaho_sim[,i] <- rJSU(n = nrow(dropD), 
                        mu = predict(fit_LSS), sigma = fit_LSS$sigma.fv,
                        nu = fit_LSS$nu.fv, tau = fit_LSS$tau.fv)
}

# Round the output to approximate the rounding in the data recording 

e = idaho_sim < log(0.6); idaho_sim[e]  = log(0.5); 
e = (idaho_sim > log(0.65))&(idaho_sim < log(0.85)); idaho_sim[e] <- log(0.75); 
e = (idaho_sim > log(0.9))&(idaho_sim < log(1.1)); idaho_sim[e] <- log(1); 

## moments of the real data by size bin
n_bins = 12
alpha_scale = 0.7
idaho_moments <- dropD %>% 
  arrange(fitted) %>% 
  mutate(size_bin = cut_number(fitted,n=n_bins)) %>% 
  group_by(size_bin) %>% 
  summarise(mean_t1 = mean(logarea.t1),
            sd_t1 = sd(logarea.t1),
            skew_t1 = skewness(logarea.t1),
            kurt_t1 = kurtosis(logarea.t1),
            bin_mean = mean(fitted),
            bin_n = n()) 


par(mfrow=c(2,2),mar=c(4,4,1,1),mgp=c(2,1,0)); 
sim_bin_means=sim_moment_means = matrix(NA,12,n_sim); 
for(i in 1:n_sim){
    sim_moments <- bind_cols(dropD,data.frame(sim=idaho_sim[,i])) %>% 
    arrange(fitted) %>% 
    mutate(size_bin = cut_number(fitted,n=n_bins)) %>% 
    group_by(size_bin) %>% 
    summarise(mean_t1 = mean(sim),
              bin_mean = mean(fitted))
	sim_bin_means[,i]=sim_moments$bin_mean; sim_moment_means[,i]=sim_moments$mean_t1; 		  
}
matplot(idaho_moments$bin_mean, sim_moment_means,col=alpha("gray",0.5),pch=16,xlab="Mean size t0",ylab="mean(Size t1)",cex=1.4); 
points(idaho_moments$bin_mean, idaho_moments$mean_t1,pch=1,lwd=2,col=alpha("red",alpha_scale),cex=1.4)
points(idaho_moments$bin_mean, apply(sim_moment_means,1,median),pch=1,lwd=2,col=alpha("black",alpha_scale),cex=1.4)
legend("topleft",legend=c("Model simulations","Median of simulations","Data"),
col=c(alpha("gray",0.5),alpha("black",alpha_scale), alpha("red",alpha_scale)),pch=1,lwd=2,cex=1.1,bty="n"); 


for(i in 1:n_sim){
    sim_moments <- bind_cols(dropD,data.frame(sim=idaho_sim[,i])) %>% 
    arrange(fitted) %>% 
    mutate(size_bin = cut_number(fitted,n=n_bins)) %>% 
    group_by(size_bin) %>% 
    summarise(mean_t1 = sd(sim),
              bin_mean = mean(fitted))
	sim_bin_means[,i]=sim_moments$bin_mean; sim_moment_means[,i]=sim_moments$mean_t1; 		  
}
matplot(idaho_moments$bin_mean, sim_moment_means,col=alpha("gray",0.5),pch=16,xlab="Mean size t0",ylab="SD(Size t1)",cex=1.4); 
points(idaho_moments$bin_mean, idaho_moments$sd_t1,pch=1,lwd=2,col=alpha("red",alpha_scale),cex=1.4)
points(idaho_moments$bin_mean, apply(sim_moment_means,1,median),pch=1,lwd=2,col=alpha("black",alpha_scale),cex=1.4)

for(i in 1:n_sim){
    sim_moments <- bind_cols(dropD,data.frame(sim=idaho_sim[,i])) %>% 
    arrange(fitted) %>% 
    mutate(size_bin = cut_number(fitted,n=n_bins)) %>% 
    group_by(size_bin) %>% 
    summarise(mean_t1 = skewness(sim),
              bin_mean = mean(fitted))
	sim_bin_means[,i]=sim_moments$bin_mean; sim_moment_means[,i]=sim_moments$mean_t1; 		  
}
matplot(idaho_moments$bin_mean, sim_moment_means,col=alpha("gray",0.5),pch=16,xlab="Mean size t0",ylab="Skew(Size t1)",cex=1.4); 
points(idaho_moments$bin_mean, idaho_moments$skew_t1,pch=1,lwd=2,col=alpha("red",alpha_scale),cex=1.4)
points(idaho_moments$bin_mean, apply(sim_moment_means,1,median),pch=1,lwd=2,col=alpha("black",alpha_scale),cex=1.4)

for(i in 1:n_sim){
    sim_moments <- bind_cols(dropD,data.frame(sim=idaho_sim[,i])) %>% 
    arrange(fitted) %>% 
    mutate(size_bin = cut_number(fitted,n=n_bins)) %>% 
    group_by(size_bin) %>% 
    summarise(mean_t1 = kurtosis(sim),
              bin_mean = mean(fitted))
	sim_bin_means[,i]=sim_moments$bin_mean; sim_moment_means[,i]=sim_moments$mean_t1; 		  
}
matplot(idaho_moments$bin_mean, sim_moment_means,col=alpha("gray",0.5),pch=16,xlab="Mean size t0",ylab="Kurtosis(Size t1)",cex=1.4); 
points(idaho_moments$bin_mean, idaho_moments$kurt_t1,pch=1,lwd=2,col=alpha("red",alpha_scale),cex=1.4)
points(idaho_moments$bin_mean, apply(sim_moment_means,1,median),pch=1,lwd=2,col=alpha("black",alpha_scale),cex=1.4)


#################################################################################
#   Simulation: how well do we recover known random effects? 
#################################################################################

U=model.matrix(~year + logarea.t0:year + I(logarea.t0^2) + Treatment - 1, data=dropD); 
LogLik=function(pars,response,U){
	pars1 = pars[1:ncol(U)]; pars2=pars[-(1:ncol(U))];
	mu = U%*%pars1;  
	val = dJSU(response, mu=mu,
	sigma=exp(pars2[1]+pars2[2]*mu),
	nu = pars2[3]+pars2[4]*mu,
	tau = exp(pars2[5]+pars2[6]*mu), log=TRUE)
	return(val); 
}

# try to fit the real data

out=maxLik(logLik=LogLik,start=runif(ncol(U)+6),response=dropD$logarea.t1,U=U,
		method="BFGS",control=list(iterlim=5000,printLevel=1),finalHessian=FALSE); 
out=maxLik(logLik=LogLik,start=out$estimate,response=dropD$logarea.t1,U=U,
		method="NM",control=list(iterlim=5000,printLevel=1),finalHessian=FALSE); 
out=maxLik(logLik=LogLik,start=out$estimate,response=dropD$logarea.t1,U=U,
		method="BFGS",control=list(iterlim=5000,printLevel=1),finalHessian=FALSE); 

simRanIntercept = simRanSlope = matrix(NA,30,50); 
for(j in 1:1) {
		U=model.matrix(~year + logarea.t0:year + I(logarea.t0^2) + Treatment - 1, data=dropD); 

		# fit to the simulated responses 
		out=maxLik(logLik=LogLik,start=runif(ncol(U)+6),response=idaho_sim[,j],U=U,
			method="BFGS",control=list(iterlim=5000,printLevel=1),finalHessian=FALSE); 
		
		out=maxLik(logLik=LogLik,start=out$estimate,response=idaho_sim[,j],,U=U,
			method="NM",control=list(iterlim=2500,printLevel=1),finalHessian=FALSE); 
			
		out=maxLik(logLik=LogLik,start=out$estimate,response=idaho_sim[,j],U=U,
			method="BFGS",control=list(iterlim=5000,printLevel=1),finalHessian=FALSE); 
		
		
		S = summary(fit_j); ### this won't work any more, need to use out$estimate. 

		# shrinkage random effects for (1|year) 
		fixed.fx = S[1:30,1]; fixed.fx = fixed.fx-mean(fixed.fx); 
		fixed.se = S[1:30,2]; 
		sigma2.hat = mean(fixed.fx^2)-mean(fixed.se^2)
		shrunk.fxj = fixed.fx*(sigma2.hat/(sigma2.hat + fixed.se^2)); 
		simRanIntercept[,j]= shrunk.fxj; 

		# shrinkage random effects for (logarea.t0|year) 
		fixed.fx2 = S[31:60,1]; fixed.fx2 = fixed.fx2-mean(fixed.fx2); 
		fixed.se2 = S[31:60,2]; 
		sigma2.hat = mean(fixed.fx2^2)-mean(fixed.se2^2)
		shrunk.fx2j = fixed.fx2*(sigma2.hat/(sigma2.hat + fixed.se2^2));
		simRanSlope[,j] = shrunk.fx2j; 
		
		cat("#####################################################","\n")
		cat("Done with simulated data set ",j,"\n"); 
}
	




