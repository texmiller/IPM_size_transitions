##############################################################################
# Growth kernel modeling for Pseudoroegneria spicata, USSES Idaho. 
#
# Historical and modern data from by Adler et al. removals experiments paper.
# Modern data are just the Control and Grass-removal treaments.  
#
# In this script, gam with family=gaulss is used for the pilot model 
#
# Original: SPE April 2020
# Last update, August 2020 
#
##############################################################################

rm(list=ls(all=TRUE));
setwd("c:/repos/IPM_size_transitions/idaho"); 

require(car); require(zoo); require(moments); require(mgcv); 
require(gamlss); require(AICcmodavg); 
require(tidyverse); require(maxLik); 

source("../Diagnostics.R"); 

##############################################################
# 1. Read in the data
##############################################################
allD<-read.csv(file="PSSP_growth_data.csv"); 
allD$year <- factor(allD$year); 

######################################################################
# Clean up the data
######################################################################
if(FALSE) {
    library(tidyverse)
## drop out these arbitrarily small individuals
    e = which(abs(allD$area.t0-0.25)<0.001); 
    dropD = allD[-e,]; 
    e = which(abs(dropD$area.t1-0.25)<0.001); 
    dropD = dropD[-e,]; 
    plot(dropD$logarea.t0, dropD$logarea.t1) 
} else {
    dropD = allD
}    

# condense treatments based on prior analysis: no difference between
# historical and modern control quadrats 
e = which(dropD$Treatment=="ControlModern");
dropD$Treatment[e] = "Control"; 
dropD = droplevels(dropD);  

########################################################################## 
## Pilot Gaussian fits with log transformation, nonconstant variance
## depending on initial size, fitted using gam(family=gaulss) 
##
## So that the models can be updated iteratively, a variable "sigma_covar" is
## added to the data frame holding the data. For these initial fits,
## sigma_covar equals initial size. It will eventually be updated iteratively
## so that it equals the fitted value (predicted mean of subsequent size)
########################################################################## 
log_models <- list()
dropD$sigma_covar = dropD$logarea.t0; 
dropD$clogarea.t0 = dropD$logarea.t0 - mean(dropD$logarea.t0); 
log_models[[1]] <- gam(list(logarea.t1~ s(logarea.t0) + W.ARTR + W.HECO + W.POSE + W.PSSP+  W.allcov + W.allpts + Treatment + 
             Group + s(year,bs="re") + s(clogarea.t0,year,bs="re"), ~s(sigma_covar,k=6)), family=gaulss,gamma=1.4,data=dropD);  

## coefficients on HECO and POSE are nearly identical, so group them: deltaAIC = 2, i.e. no change at all in logLik.   
log_models[[2]] <- gam(list(logarea.t1~ s(logarea.t0) + W.ARTR + I(W.HECO + W.POSE) + W.PSSP+  W.allcov + W.allpts + Treatment + 
             Group + s(year,bs="re") + s(clogarea.t0,year,bs="re"), ~s(sigma_covar,k=6)), family=gaulss,gamma=1.4,data=dropD);  

## W.allpts is non-significant. Consider two options: drop, group with all other cover 
log_models[[3]] <- 	gam(list(logarea.t1~ s(logarea.t0) + W.ARTR + I(W.HECO + W.POSE) + W.PSSP+  W.allpts + Treatment + 
             Group + s(year,bs="re") + s(clogarea.t0,year,bs="re"), ~s(sigma_covar,k=6)), family=gaulss,gamma=1.4,data=dropD);  	 

## Based on AIC, the winner by a hair is to group all heterospecific grasses, and group all cover besides the 'big 4'
log_models[[4]] <- gam(list(logarea.t1~ s(logarea.t0) + W.ARTR + I(W.HECO + W.POSE) + W.PSSP +  I(W.allcov + W.allpts) + Treatment + 
             Group + s(year,bs="re") + s(clogarea.t0,year,bs="re"), ~s(sigma_covar,k=6)), family=gaulss,gamma=1.4,data=dropD);  

                   
for(j in 1:4) cat(j, AIC(log_models[[j]]), "\n"); # model 4, by a hair 

### Save these models, for AIC comparison with models where sigma depends on the fitted value  
init_models <- log_models; 

########################################################################################## 
## iterate to fit models where sigma depends on fitted value
##########################################################################################
for(mod in 1:4) {
  fitGAU = log_models[[mod]]
  fitted_all = predict(fitGAU,type="response",data=dropD);                  
  fitted_vals = new_fitted_vals = fitted_all[,1]; 
  weights = fitted_all[,2]; # these are 1/sigma values; see ?gaulss for details. 

  err=100; k=0; 
  while(err>10^(-6)) {
    dropD$sigma_covar = new_fitted_vals; 
    fitGAU <- update(fitGAU); 
    fitted_all = predict(fitGAU,type="response",data=dropD);   
    new_fitted_vals = fitted_all[,1]; new_weights = fitted_all[,2];
    err = weights - new_weights; err=sqrt(mean(err^2)); 
    weights = new_weights; 
    k=k+1; cat(k,err,"\n"); 
  }   
  log_models[[mod]] =  fitGAU;
}

for(j in 1:4) cat(j, AIC(log_models[[j]]), "\n"); 
for(j in 1:4) cat(j, AIC(init_models[[j]]), "\n"); 
## Models based on fitted value are strongly preferred (Delta AIC \approx 50)

######### Here's the best Gaussian model ########################################
aics = unlist(lapply(log_models,AIC)); best_model=which(aics==min(aics)); 
log_model = log_models[[best_model]]; 
summary(log_model); 

######################################################################
# Interrogate the scaled residuals - Gaussian? NO!  
###################################################################### 
fitted_all = predict(log_model,type="response",data=dropD); 
plot(fitted_all[,1], 1/fitted_all[,2],xlab="Fitted",ylab="Estimated residual Std Dev"); 

scaledResids = residuals(log_model,type="response")*fitted_all[,2]
plot(fitted_all[,1], scaledResids); 

qqPlot(scaledResids); # really bad in lower tail, not too bad in upper 

jarque.test(scaledResids) # normality test: FAILS, P < 0.001 
anscombe.test(scaledResids) # kurtosis: FAILS, P < 0.001 
agostino.test(scaledResids) # skewness: FAILS, P<0.001 

########################################################################
## Rollapply diagnostics on the scaled residuals 
########################################################################
px = fitted_all[,1]; py=scaledResids; 

graphics.off(); dev.new(width=8,height=6); 
par(mfrow=c(2,2),bty="l",mar=c(4,4,2,1),mgp=c(2.2,1,0),cex.axis=1.4,cex.lab=1.4);   
z = rollMomentsNP(px,py,windows=20,smooth=TRUE,scaled=TRUE) 
# mean and SD look OK; skew changes sign, kurtosis is always positive but close to 0 

dev.copy2pdf(file="../manuscript/figures/RollingNPMomentsPSSP.pdf") 

###########################################################################
# Fit suitable distributions to binned data
# NOTE: we bin with respect to fitted values of the pilot Gaussian model
# because that is going to be the form of the final model. 
###########################################################################
logResids <- data.frame(init=fitted_all[,1],resids=scaledResids); 
logResids <- logResids %>% mutate(size_bin = cut_number(init,n=25))

source("../fitChosenDists.R"); 

tryDists=c("SHASH","SHASHo2","SEP1","SEP3","SEP4"); 
# JSU and ST distributions removed because they can't do exactly zero kurtosis. 
# EGB2 was initially included, but removed as it was slow and failed to converge for many bins; the  
# gamlss book notes that it is only good for mildly leptokurtosis and low skewness. 

bins = levels(logResids$size_bin); maxVals = matrix(NA,length(bins),length(tryDists)); 
for(j in 1:length(bins)){
for(k in 1:length(tryDists)) {
	Xj=subset(logResids,size_bin==bins[j])
	fitj = gamlssMaxlik(y=Xj$resids,DIST=tryDists[k]); 
	maxVals[j,k] = fitj$out[[1]]$maximum;
	cat("Finished ", tryDists[k]," bin ",j, fitj$maximum,"\n"); 
}
}

## ranking top 3 for each bin 
for(j in 1:length(bins)){
	e = order(-maxVals[j,]); 
	cat(j, tryDists[e][1:3],"\n"); 
}	

# Overall ranking: SHASH is the winner 
e = order(-colSums(maxVals)); 
rbind(tryDists[e],round(colSums(maxVals)[e],digits=3)); 

Results with 20 bins 
# [1,] "SHASH"     "SEP4"      "SEP1"      "SHASHo2"   "SEP3"     
# [2,] "-6489.268" "-6491.345" "-6493.011" "-6494.701" "-6497.499"

Results with 25 bins       
#[1,] "SHASH"     "SEP1"      "SEP4"      "SEP3"      "SHASHo2"  
#[2,] "-6477.462" "-6477.806" "-6478.314" "-6482.432" "-6482.462"
 

########################################################################
## Binned-data SHASH parameter estimates on subsequent size 
########################################################################
source("../Diagnostics.R"); 
source("../fitChosenDists.R"); 
binPars = binnedPars(y=dropD$logarea.t1,sortVar = fitted_all[,1],nBins=20,DIST="SHASH",rolling=FALSE) 

graphics.off(); dev.new(width=8,height=6); 
par(mfrow=c(2,2),bty="l",mar=c(4,4,2,1),mgp=c(2.2,1,0),cex.axis=1.4,cex.lab=1.4);   

with(binPars,{
spline.scatter.smooth(bin_means,mus,gamma=1.4,xlab="Fitted value",ylab=expression(paste("Location parameter  ", mu ))); #OK
add_panel_label("a"); 
spline.scatter.smooth(mus,sigmas,gamma=1.4,xlab=expression(paste("Location parameter  ", mu )),
			ylab=expression(paste("log Scale parameter  ", sigma))); #linear 
add_panel_label("b"); 
spline.scatter.smooth(mus,nus,gamma=1.4,xlab=expression(paste("Location parameter  ", mu )),
			ylab=expression(paste("log Shape parameter  ", nu ))); #quadratic
add_panel_label("c"); 
spline.scatter.smooth(mus,taus,gamma=1.4,xlab=expression(paste("Location parameter  ", mu )),
			ylab=expression(paste("log Shape parameter  ", tau))); #linear 
add_panel_label("d"); 
})
  
dev.copy2pdf(file="../manuscript/figures/RollingSHASHparsPSSP.pdf") 
savePlot(file="../manuscript/figures/RollingSHASHparsPSSP.png",type="png"); 

save.image(file="PSSPmodeling.Rdata");   

##############################################################################################
# We want to fit an SHASH growth model, based on the pilot Gaussian model. 
# We can't use splines in this model, so we interrogate the fitted gam model to choose
# functional forms that can accomplish the same shapes. 
##############################################################################################

plot(log_model,scale=0); # remind ourselves what the splines look like. No wagging tails! 

### sigma as a function of fitted values 
fitted_all = predict(log_model,type="response",data=dropD); 
zvals = dropD$logarea.t0; yhat=fitted_all[,1]; sigma_hat = 1/fitted_all[,2]; 
sd_fit2 = lm(log(sigma_hat) ~ yhat + I(yhat^2)); # major trend in residuals 
sd_fit3 = lm(log(sigma_hat) ~ yhat + I(yhat^2) + I(yhat^3)); # looks pretty good, R^2 = 0.99 

### nonlinear effect of initial size 
fitted_terms = predict(log_model,type="terms",data=dropD); 
mu_hat = fitted_terms[,"s(logarea.t0)"]; 
mu_fit1 = lm(mu_hat ~ zvals); # major trend in residuals 
mu_fit2 = lm(mu_hat ~ zvals + I(zvals^2)); # looks good - use this one.  
mu_fit3 = lm(mu_hat ~ zvals + I(zvals^2) + I(zvals^3)); # nearly identical to quadratic. 

# sigma_hat2 = fitted_terms[,"s.1(sigma_covar)"]
# plot(sigma_hat,exp( -0.3564609 + sigma_hat2)+0.01);  # Note, the fitted_terms column for the smooth omits the intercept. 


##################################################################################################
# Final model fitting is done by maximum likelihood using maxLik 
# 
# The random effects terms are fitted here as fixed effects, then adjusted by shrinkage.  
# After fitting (1|year) and (0 + logarea.t0|year) as fixed effects by ML, we have year-specific 
# coefficients and their estimated standard errors, which are all we need to do shrinkage. 
#####################################################################################################

# Model matrix for the fixed and random effects, specified so that each year gets its own coefficient
# rather than fitting contrasts against a baseline year (the default in R's regression functions) 
# Note that logarea.t0:year below takes care of the linear term in logarea.t0 

U=model.matrix(~year + logarea.t0:year + I(logarea.t0^2) +  
	W.ARTR + I(W.HECO + W.POSE) + W.PSSP+  I(W.allcov + W.allpts) + Group + Treatment - 1, data=dropD);

## funtional forms of the distribution parameters vs. linear predictor 
f.sigma=function(p,x) exp(p[1] + p[2]*x +p[3]*x^2 + p[4]*x^3)  # from lm fits
f.nu = function(p,x) exp(p[5] + p[6]*x + p[7]*x^2)  # from the rolling pars plot 
f.tau = function(p,x) exp(p[8] + p[9]*x)  # from the rolling pars plot 


## Log-likelihood for the model. 
LogLik=function(pars,response,U){
	pars1 = pars[1:ncol(U)]; pars2=pars[-(1:ncol(U))];
	mu = U%*%pars1;  # "fitted values" -- actually, the linear predictor. 
	val = dSHASH(response, mu=mu, sigma = f.sigma(pars2,mu), 
                nu = f.nu(pars2,mu), tau = f.tau(pars2,mu), log=TRUE)
	return(val); 
}

############ Fit the data. Paranoid as usual about convergence. 
coefs = list(5); LL=numeric(5);  

# To get starting values for the linear predictor, do a pilot fit
# using the spline for sigma from the gam fit.   
pfit = lm(dropD$logarea.t1 ~ U-1, weights=1/sigma_hat^2); 
fixed_start = as.numeric(coef(pfit)); 

# shape coefficients from the rolling parameters plot 
fit_nu = with(binPars,lm(nus ~ mus + I(mus^2)));
fit_tau = with(binPars,lm(taus~ mus)); 
p0=c(fixed_start, coef(sd_fit3), coef(fit_nu), coef(fit_tau));

for(j in 1:5) {
	out=maxLik(logLik=LogLik,start=p0*exp(0.1*rnorm(length(p0))), response=dropD$logarea.t1,U=U,
		method="BHHH",control=list(iterlim=5000,printLevel=2),finalHessian=FALSE); 

	out=maxLik(logLik=LogLik,start=out$estimate,response=dropD$logarea.t1,U=U,
		method="NM",control=list(iterlim=5000,printLevel=1),finalHessian=FALSE); 

	out=maxLik(logLik=LogLik,start=out$estimate,response=dropD$logarea.t1,U=U,
		method="BHHH",control=list(iterlim=5000,printLevel=2),finalHessian=FALSE); 

	coefs[[j]] = out$estimate; LL[j] = out$maximum;
	cat(j, "#--------------------------------------#",out$maximum,"\n"); 
}

j = min(which(LL==max(LL)));
MLfit=maxLik(logLik=LogLik,start=coefs[[j]],response=dropD$logarea.t1,U=U,
		method="BHHH",control=list(iterlim=5000,printLevel=2),finalHessian=TRUE); 

######### save results of SHASH ML fit. 
names(MLfit$estimate)<-c(colnames(U),c("sigma.0","sigma.1","sigma.2","sigma.3","nu.0","nu.1","nu.2","tau.0","tau.1"));  
coefs=MLfit$estimate; SEs = sqrt(diag(vcov(MLfit))); 

save.image(file="PSSPmodeling.Rdata");   

#######################################################################################
#  Plot the fits 
#######################################################################################
graphics.off(); dev.new(width=8,height=6); 
par(mfrow=c(2,2),bty="l",mar=c(4,4,2,1),mgp=c(2.2,1,0),cex.axis=1.4,cex.lab=1.4);   

pars1 = coefs[1:ncol(U)]; pars2=coefs[-(1:ncol(U))];
mu = U%*%pars1; 
x = seq(min(mu),max(mu),length=100); 
plot(x,f.sigma(pars2,x), xlab="Fitted value", ylab="Scale parameter sigma",type="l"); 
plot(x, f.nu(pars2,x), xlab="Fitted value",ylab="Shape parameter nu",type="l");
plot(x,	f.tau(pars2,x), xlab="Fitted value",ylab="Shape parameter tau",type="l");


#############################################################################
#  Binned quantiles diagnostic applied to the model fitted by ML 
#############################################################################
# Simulate data from fitted SHASH model
pars1 = coefs[1:ncol(U)]; 
pars2 = coefs[-(1:ncol(U))];
MLmu = U%*%pars1;  

n_sim <- 500
idaho_sim<-matrix(NA,nrow=nrow(dropD),ncol=n_sim)
for(i in 1:n_sim){
  idaho_sim[,i] <- rSHASH(n = nrow(dropD), mu = MLmu, 
                    sigma = f.sigma(pars2,MLmu),
                    nu = f.nu(pars2,MLmu),
                    tau = f.tau(pars2,MLmu))
   if(i%%10==0) cat(i,"\n");                  
}

# Round the output to approximate the rounding in the data recording 
e = idaho_sim < log(0.6); idaho_sim[e]  = log(0.5); 
e = (idaho_sim > log(0.65))&(idaho_sim < log(0.85)); idaho_sim[e] <- log(0.75); 
e = (idaho_sim > log(0.9))&(idaho_sim < log(1.1)); idaho_sim[e] <- log(1); 

dev.new(width=7,height=9); 
out=quantileComparePlot(dropD$logarea.t0,dropD$logarea.t1,idaho_sim,12);

dev.copy2pdf(file="../manuscript/figures/QuantileComparePlotPSSP.pdf")
save.image(file="PSSPmodeling.Rdata"); 

#############################################################################
#  Do the same for the pilot Gaussian fit 
#############################################################################
n_sim <- 500
idaho_Gsim<-matrix(NA,nrow=nrow(dropD),ncol=n_sim)
fitted_all = predict(log_model,type="response"); # you remember log_model... 
yhat = fitted_all[,1]; sigma_hat = 1/fitted_all[,2]; 
for(i in 1:n_sim){
  idaho_Gsim[,i] <- rnorm(n = nrow(dropD), mean=yhat, sd=sigma_hat); 
  if(i%%50==0) cat(i,"\n");                  
}

# Round the output to approximate the rounding in the data recording 
e = idaho_Gsim < log(0.6); idaho_Gsim[e]  = log(0.5); 
e = (idaho_Gsim > log(0.65))&(idaho_Gsim < log(0.85)); idaho_Gsim[e] <- log(0.75); 
e = (idaho_Gsim > log(0.9))&(idaho_Gsim < log(1.1)); idaho_Gsim[e] <- log(1); 

dev.new(width=7,height=9); 
out=quantileComparePlot(dropD$logarea.t0,dropD$logarea.t1,idaho_Gsim,nBins=12);
title(main="Gaussian pilot model"); 

dev.copy2pdf(file="../manuscript/figures/QuantileComparePlotPSSP-pilot.pdf")

#############################################################################
#  Binned moments comparison plot 
#############################################################################
out=momentsComparePlot(dropD$logarea.t0,dropD$logarea.t1,idaho_sim, idaho_Gsim,nBins=12);
dev.copy2pdf(file="../manuscript/figures/PSSPMomentsComparePlot.pdf")

out=momentsComparePlot(dropD$logarea.t0,dropD$logarea.t1,idaho_sim,idaho_Gsim,fittedMean=fitted_all[,1],nBins=12);
dev.copy2pdf(file="../manuscript/figures/PSSPMomentsComparePlot2.pdf")

########################################################################################
#   Simulation: how well do we recover known random effects? 
#
#   Use simulated replicates as "data" and compare shrinkage estimates with the "truth". 
#   Key question: how well does shrinkage do at getting the variance right? 
#########################################################################################
U=model.matrix(~year + logarea.t0:year + I(logarea.t0^2) +  
	W.ARTR + I(W.HECO + W.POSE) + W.PSSP+  I(W.allcov + W.allpts) + Group + Treatment - 1, data=dropD);

##### Re-do simulations, so there is no rounding of small values 
pars1 = coefs[1:ncol(U)]; pars2 = coefs[-(1:ncol(U))]; MLmu = U%*%pars1;  
n_sim <- 250 
idaho_sim2<-matrix(NA,nrow=nrow(dropD),ncol=n_sim)
for(i in 1:n_sim){
  idaho_sim2[,i] <- rSHASH(n = nrow(dropD), 
                    mu = MLmu, sigma = f.sigma(pars2,MLmu),
                    nu = f.nu(pars2,MLmu),tau = f.tau(pars2,MLmu))
   if(i%%10==0) cat(i,"\n");                  
}

LogLik=function(pars,response,U){
	pars1 = pars[1:ncol(U)]; pars2=pars[-(1:ncol(U))];
	mu = U%*%pars1;  
	val = dSHASH(response, mu=mu,
	sigma=f.sigma(pars2,mu), nu = f.nu(pars2,mu),
	tau = f.tau(pars2,mu), log=TRUE)
	return(val); 
}

#### Estimate the random effects in various ways, compare the results 
shrinkRanIntercept = shrinkRanSlope = matrix(NA,30,250); 
shrinkRanIntercept2 = shrinkRanSlope2 = matrix(NA,30,250); 
fixRanIntercept = fixRanSlope = matrix(NA,30,250); 
lmerRanIntercept = lmerRanSlope = matrix(NA,30,250); 

# Use the same weights for lmer fits 
MLweights =  fitted_all[,2];

T = length(unique(dropD$year)); 

for(j in 1:250) {

		# fit by ML to the simulated responses 
		outj=maxLik(logLik=LogLik,start=coefs,response=idaho_sim2[,j],U=U,
			method="BHHH",control=list(iterlim=5000,printLevel=2),finalHessian=FALSE); 
		
		outj=maxLik(logLik=LogLik,start=outj$estimate,response=idaho_sim2[,j],,U=U,
			method="NM",control=list(iterlim=5000,printLevel=1),finalHessian=FALSE); 
			
		outj=maxLik(logLik=LogLik,start=outj$estimate,response=idaho_sim2[,j],U=U,
			method="BHHH",control=list(iterlim=5000,printLevel=1),finalHessian=TRUE); 
		
		coefj=outj$estimate; SEj = sqrt(diag(vcov(outj))); V = vcov(outj); 
		V1 = V[1:T,1:T]; V2 = V[(T+1):(2*T),(T+1):(2*T)]; 
		
		# shrinkage random effects for (1|year) 
		fixed.fxj = coefj[1:30]; fixed.fxj = fixed.fxj-mean(fixed.fxj); 
		fixed.sej = SEj[1:30]; 
		sigma2.hat = mean(fixed.fxj^2)-mean(fixed.sej^2)+ (sum(V1)-sum(diag(V1)))/(2*T*(T-1)); 
	    shrinkRanIntercept[,j] = fixed.fxj*(sigma2.hat/(sigma2.hat + fixed.sej^2)); 
		shrinkRanIntercept2[,j] = fixed.fxj*sqrt(sigma2.hat/(sigma2.hat + fixed.sej^2));
		fixRanIntercept[,j] = fixed.fxj; 


		# shrinkage random effects for (logarea.t0|year) 
		fixed.fx2j = coefj[42:71]; fixed.fx2j = fixed.fx2j-mean(fixed.fx2j); 
		fixed.se2j = SEj[42:71]; 
		sigma2.hat = mean(fixed.fx2j^2)-mean(fixed.se2j^2)+ (sum(V2)-sum(diag(V2)))/(2*T*(T-1)); 
		shrinkRanSlope[,j] = fixed.fx2j*(sigma2.hat/(sigma2.hat + fixed.se2j^2));
		shrinkRanSlope2[,j] = fixed.fx2j*sqrt(sigma2.hat/(sigma2.hat + fixed.se2j^2));
		fixRanSlope[,j] = fixed.fx2j; 

		# lmer random effects, using the true variance function  
		lmerFit <- lmer(idaho_sim2[,j]~ logarea.t0 + I(logarea.t0^2) + W.ARTR + I(W.HECO + W.POSE) + W.PSSP+  I(W.allcov + W.allpts) + Treatment + 
             Group + (logarea.t0|year), control=lmerControl(optimizer="bobyqa"),weights=MLweights,data=dropD,REML=TRUE); 	
		ran.fx = ranef(lmerFit)[[1]]; 
		lmerRanIntercept[,j] = ran.fx[,1]; 
		lmerRanSlope[,j] = ran.fx[,2]; 


		cat("#####################################################","\n")
		cat("Done with simulated data set ",j,"\n"); 
}
	
par(mfrow=c(2,2),bty="l",mgp=c(2,1,0),mar=c(4,4,1,1),cex.axis=1.3,cex.lab=1.3);
trueRanIntercept = coefs[1:30]-mean(coefs[1:30]);
matplot(trueRanIntercept,fixRanIntercept,type="p",pch=1,col="black");abline(0,1,col="blue"); 
add_panel_label("a"); 
matplot(trueRanIntercept,shrinkRanIntercept2,type="p",pch=1,col="black");abline(0,1,col="blue"); 
add_panel_label("b"); 
matplot(trueRanIntercept,shrinkRanIntercept,type="p",pch=1,col="black");abline(0,1,col="blue"); 
add_panel_label("c"); 
matplot(trueRanIntercept,lmerRanIntercept,type="p",pch=1,col="black");abline(0,1,col="blue"); 
add_panel_label("d"); 
dev.copy2pdf(file="../manuscript/figures/PSSPShrinkageTest.pdf"); 


# compare estimates of the between-year mixing sigma, averaging over reps 
sd(trueRanIntercept); 
mean(apply(fixRanIntercept,2,sd)); 
mean(apply(shrinkRanIntercept,2,sd)); 
mean(apply(shrinkRanIntercept2,2,sd)); # winner
mean(apply(lmerRanIntercept,2,sd)); 

# compare root-mean-square error
RMSE = function(target, ests) {
	out=apply(ests, 2, function(x) sum((x-target)^2) )
	sqrt(mean(out)/length(target)) 
}
RMSE(trueRanIntercept,fixRanIntercept); 
RMSE(trueRanIntercept,shrinkRanIntercept); 
RMSE(trueRanIntercept,shrinkRanIntercept2); # winner 
RMSE(trueRanIntercept,lmerRanIntercept); 

par(mfrow=c(2,2),bty="l",mgp=c(2,1,0),mar=c(4,4,1,1),cex.axis=1.3,cex.lab=1.3);
trueRanSlope = coefs[42:71]-mean(coefs[42:71]);
matplot(trueRanSlope,fixRanSlope,type="p",pch=1,col="black");abline(0,1,col="blue"); 
matplot(trueRanSlope,shrinkRanSlope,type="p",pch=1,col="black");abline(0,1,col="blue"); 
matplot(trueRanSlope,shrinkRanSlope2,type="p",pch=1,col="black");abline(0,1,col="blue"); 
matplot(trueRanSlope,lmerRanSlope,type="p",pch=1,col="black");abline(0,1,col="blue"); 

# compare estimates of the between-year mixing sigma, averaging over reps 
sd(trueRanSlope); 
mean(apply(fixRanSlope,2,sd)); 
mean(apply(shrinkRanSlope,2,sd)); 
mean(apply(shrinkRanSlope2,2,sd)); # winner 
mean(apply(lmerRanSlope,2,sd)); 

# compare sum of squared errors
RMSE(trueRanSlope,fixRanSlope); 
RMSE(trueRanSlope,shrinkRanSlope); 
RMSE(trueRanSlope,shrinkRanSlope2); # winner 
RMSE(trueRanSlope,lmerRanSlope); 

save(file="ShrinkageTest.Rdata",list=c("idaho_sim2","trueRanIntercept","fixRanIntercept","shrinkRanIntercept","shrinkRanIntercept2","lmerRanIntercept",
"trueRanSlope","fixRanSlope","shrinkRanSlope","shrinkRanSlope2","lmerRanSlope","dropD")); 

save.image(file="PSSPmodeling.Rdata"); 


