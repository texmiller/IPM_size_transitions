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
require(gamlss); require(gamlss.tr); require(AICcmodavg); 
require(tidyverse); require(maxLik); 

source("Utilities.R");
source("../Diagnostics.R"); 

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

########################################################################## 
## Pilot Gaussian fits with log transformation, nonconstant variance
## depending on initial size 
########################################################################## 
log_models <- list()
dropD$sigma_covar = dropD$logarea.t0; # for iterative fitting, later 
log_models[[1]] <- gam(list(logarea.t1~ s(logarea.t0) + W.ARTR + W.HECO + W.POSE + W.PSSP+  W.allcov + W.allpts + Treatment + 
             Group + s(year,bs="re") + s(logarea.t0,year,bs="re"), ~s(sigma_covar)), family=gaulss,gamma=1.4,data=dropD);  

## coefficients on HECO and POSE are nearly identical, so group them: deltaAIC = 2, i.e. no change at all in logLik.   
log_models[[2]] <- gam(list(logarea.t1~ s(logarea.t0) + W.ARTR + I(W.HECO + W.POSE) + W.PSSP+  W.allcov + W.allpts + Treatment + 
             Group + s(year,bs="re") + s(logarea.t0,year,bs="re"), ~s(sigma_covar)), family=gaulss,gamma=1.4,data=dropD);  

## W.allpts is non-significant. Consider two options: drop, group with all other cover 
log_models[[3]] <- 	gam(list(logarea.t1~ s(logarea.t0) + W.ARTR + I(W.HECO + W.POSE) + W.PSSP+  W.allpts + Treatment + 
             Group + s(year,bs="re") + s(logarea.t0,year,bs="re"), ~s(sigma_covar)), family=gaulss,gamma=1.4,data=dropD);  	 

## Based on AIC, the winner by a hair is to group all heterospecific grasses, and group all cover besides the 'big 4'
log_models[[4]] <- gam(list(logarea.t1~ s(logarea.t0) + W.ARTR + I(W.HECO + W.POSE) + W.PSSP+  I(W.allcov + W.allpts) + Treatment + 
             Group + s(year,bs="re") + s(logarea.t0,year,bs="re"), ~s(sigma_covar)), family=gaulss,gamma=1.4,data=dropD);  

for(j in 1:4) cat(j, AIC(log_models[[j]]), "\n"); # model 4, by a hair 
init_models <- log_models; 

fitted_all = predict(init_models[[4]],type="response",data=dropD); 
plot(dropD$logarea.t0,fitted_all[,2]); 

########################################################################################## 
## iterate to fit a model where sigma depends on fitted value, not on initial size 
##########################################################################################
for(mod in 1:4) {
  fitGAU = log_models[[mod]]
  fitted_all = predict(fitGAU,type="response",data=dropD);                  
  fitted_vals = new_fitted_vals = fitted_all[,1]; 
  weights = fitted_all[,2]; # what I call "weights" here are 1/sigma values; see ?gaulss for details. 

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
## so we will go with those 

######### Here's the best Gaussian model ########################################
aics = unlist(lapply(log_models,AIC)); best_model=which(aics==min(aics)); 
log_model = log_models[[best_model]]; 
summary(log_model); 

######################################################################
# Interrogate the scaled residuals - Gaussian? NO. 
###################################################################### 
fitted_all = predict(log_model,type="response",data=dropD); 
plot(fitted_all[,1], 1/fitted_all[,2],xlab="Fitted",ylab="Estimated residual Std Dev"); 

scaledResids = residuals(log_model,type="response")*fitted_all[,2]
plot(fitted_all[,1], scaledResids); 

qqPlot(log_scaledResids); # really bad in lower tail, not too bad in upper 

jarque.test(scaledResids) # normality test: FAILS, P < 0.001 
anscombe.test(scaledResids) # kurtosis: FAILS, P < 0.001 
agostino.test(scaledResids) # skewness: FAILS, P<0.001 

########################################################################
## Rollapply diagnostics on the scaled residuals 
########################################################################
px = fitted_all[,1]; py=scaledResids; 

graphics.off(); dev.new(width=8,height=6); 
par(mfrow=c(2,2),bty="l",mar=c(4,4,2,1),mgp=c(2.2,1,0),cex.axis=1.4,cex.lab=1.4);   
z = rollMomentsNP(px,py,windows=10,smooth=TRUE,scaled=TRUE) 
# mean and SD look OK; skew changes sign, kurtosis is always positive but close to 0 


dev.copy2pdf(file="../manuscript/figures/RollingNPMomentsPSSP.pdf") 

####################################
#  OKAY DOWN TO HERE 
####################################

###########################################################################
# Fit suitable distributions to binned data 
###########################################################################
logResids <- data.frame(init=dropD$logarea.t0,resids=scaledResids); 
logResids <- logResids %>% mutate(size_bin = cut_number(init,n=10))

source("../fitChosenDists.R"); 

tryDists=c("EGB2","GT","JSU", "SHASHo","SEP1","SEP3","SEP4"); 

bins = levels(logResids$size_bin); maxVals = matrix(NA,length(bins),length(tryDists)); 
for(j in 1:length(bins)){
for(k in 1:length(tryDists)) {
	Xj=subset(logResids,size_bin==bins[j])
	fitj = gamlssMaxlik(y=Xj$resids,DIST=tryDists[k]); 
	maxVals[j,k] = fitj$maximum;
	cat("Finished ", tryDists[k]," ",j,k, fitj$maximum,"\n"); 
}
}

## best two for each bin 
for(j in 1:length(bins)){
	e = order(-maxVals[j,]); 
	cat(j, tryDists[e][1:8],"\n"); 
}	

# overall ranking 
e = order(-colSums(maxVals)); 
rbind(tryDists[e],round(colSums(maxVals)[e],digits=3)); 










###########################################################################
# Try fitDist on binned data, various ways (only the last version remains)
###########################################################################
logResids<- data.frame(fitted=fitted(log_model),resids=dropD$logarea.t1); 
logResids <- logResids %>% mutate(size_bin = cut_number(fitted,n=12))

bins = levels(logResids$size_bin); dists=list(length(bins)); 
for(j in 1:length(bins)){
	Xj=subset(logResids,size_bin==bins[j])
	dists[[j]]=fitDist(resids,data=Xj,type="realline"); 
	cat(j,"\n"); 
}
### Results are not consistent, lots of convergence errors, 
### and best for a bin often has skew but no kurtosis  

########################################################################
## Binned data SHASH parameter estimates on subsequent size 
########################################################################
dropD$fitted <- fitted(log_model); 
dropD <- dropD %>% mutate(size_bin = cut_number(fitted,n=12))
bins = levels(dropD$size_bin); 
mus = sigmas = nus = taus = bin_means = numeric(length(bins)); 
for(j in 1:length(bins)){
	Xj=subset(dropD,size_bin==bins[j])
	fitj = gamlssML(Xj$logarea.t1 ~ 1,family="SHASHo") # lazy, but OK for now 
	mus[j]=fitj$mu; sigmas[j]=fitj$sigma; nus[j]=fitj$nu; taus[j]=fitj$tau; 
	bin_means[j]=mean(Xj$fitted); 
}	

graphics.off(); dev.new(width=8,height=6); 
par(mfrow=c(2,2),bty="l",mar=c(4,4,2,1),mgp=c(2.2,1,0),cex.axis=1.4,cex.lab=1.4);   
spline.scatter.smooth(bin_means,mus,xlab="Fitted value",ylab=expression(paste("Location parameter  ", mu ))); #OK
add_panel_label("a"); 
spline.scatter.smooth(mus,log(sigmas),xlab=expression(paste("Location parameter  ", mu )),
			ylab=expression(paste("log Scale parameter  ", sigma))); #linear 
add_panel_label("b"); 
spline.scatter.smooth(mus,nus,xlab=expression(paste("Location parameter  ", mu )),
			ylab=expression(paste("Skewness parameter  ", nu ))); #quadratic
add_panel_label("c"); 
spline.scatter.smooth(mus,log(taus),xlab=expression(paste("Location parameter  ", mu )),
			ylab=expression(paste("log Kurtosis parameter  ", tau))); #linear 
add_panel_label("d"); 
  
dev.copy2pdf(file="../manuscript/figures/RollingSHASHparsPSSP.pdf") 
savePlot(file="../manuscript/figures/RollingSHASHparsPSSP.png",type="png"); 


 
################################################################################################################
# Try fitting gamlss SHASHo to log size data. Use the fixed effect structure corresponding to the best lmer fit,
# and believe rollaply diagnostics 
#
# Fitting is done by maximum likelihood using maxLik (easier than mle) 
# 
# The random effects terms in the lmer fit are fitted here as fixed effects, then adjusted by shrinkage.  
# Having fitted (1|year) and (0 + logarea.t0|year) as fixed effects, we have year-specific coefficients
# and their estimated standard errors, which are all we need to do shrinkage. 
#################################################################################################################

# Model matrix for the fixed and random effects, specified so that each year gets its own coefficient
# rather than fitting contrasts against a baseline year (the default in R's regression functions) 
U=model.matrix(~year + logarea.t0:year + I(logarea.t0^2) +  
	W.ARTR + I(W.HECO + W.POSE) + W.PSSP+  I(W.allcov + W.allpts) + Group + Treatment - 1, data=dropD);

LogLik=function(pars,response,U){
	pars1 = pars[1:ncol(U)]; pars2=pars[-(1:ncol(U))];
	mu = U%*%pars1;  
	val = dSHASHo(response, mu=mu,
	sigma=exp(pars2[1]+pars2[2]*mu),
	nu = pars2[3]+pars2[4]*mu + pars2[5]*mu^2,
	tau = exp(pars2[6]+pars2[7]*mu), log=TRUE)
	return(val); 
}

############ Fit the data. Paranoid as usual about convergence
coefs = list(5); LL=numeric(5);  

# Starting values from the pilot model are jittered to do multi-start optimization). 
# Using good starting values really speeds up convergence in the ML fits  

# Linear predictor coefficients extracted from the lmer model 
lmer_coef = coef(log_model)$year; 
fixed_start = unlist(c(lmer_coef[1:30,1], lmer_coef[1,-(1:2)], lmer_coef[1:30,2])); 

# Shape and scale coefficients from the rollaply diagnostic plots 
fit_sigma = lm(log(sigmas)~mus);
fit_nu = lm(nus~mus + I(mus^2));
fit_tau = lm(log(taus)~mus); 
p0=c(fixed_start, coef(fit_sigma), coef(fit_nu),coef(fit_tau));


for(j in 1:5) {
	out=maxLik(logLik=LogLik,start=p0*exp(0.2*rnorm(length(p0))), response=dropD$logarea.t1,U=U,
		method="BHHH",control=list(iterlim=5000,printLevel=2),finalHessian=FALSE); 

	out=maxLik(logLik=LogLik,start=out$estimate,response=dropD$logarea.t1,U=U,
		method="NM",control=list(iterlim=5000,printLevel=1),finalHessian=FALSE); 

	out=maxLik(logLik=LogLik,start=out$estimate,response=dropD$logarea.t1,U=U,
		method="BHHH",control=list(iterlim=5000,printLevel=2),finalHessian=FALSE); 

	coefs[[j]] = out$estimate; LL[j] = out$maximum;
	cat(j, "#--------------------------------------#",out$maximum,"\n"); 
}

j = min(which(LL==max(LL)));
out=maxLik(logLik=LogLik,start=coefs[[j]],response=dropD$logarea.t1,U=U,
		method="BHHH",control=list(iterlim=5000,printLevel=2),finalHessian=TRUE); 

######### save results of ML fit.  
names(out$estimate)<-colnames(U); coefs=out$estimate; SEs = sqrt(diag(vcov(out))); 

###########################################################################
#  Compare lmer and Shrinkage random effects
#  With SHASH this is not really meaningful, as mu is not the mean response
###########################################################################

par(mfrow=c(2,2),mar=c(4,4,1,1),mgp=c(2,1,0),bty="l"); 
# shrinkage random effects for (1|year) 
fixed.fx = coefs[1:30]; fixed.fx = fixed.fx-mean(fixed.fx); 
fixed.se = SEs[1:30]; 
sigma2.hat = mean(fixed.fx^2)-mean(fixed.se^2)
shrunk.fx = fixed.fx*sqrt(sigma2.hat/(sigma2.hat + fixed.se^2)); 

# lmer random effects for (1|year) 
ran.fx = ranef(log_model)[[1]]; ran.fx = ran.fx[,1]; 
sd(fixed.fx); sd(shrunk.fx); sd(ran.fx); 
plot(ran.fx,shrunk.fx,xlab="lmer year random effects",ylab="Shrunk year fixed effects",pch=1,lwd=2,cex=1.2); 
abline(0,1,col="blue",lty=2); 

# shrinkage random effects for (logarea.t0|year) 
fixed.fx2 = coefs[42:71]; fixed.fx2 = fixed.fx2-mean(fixed.fx2); 
fixed.se2 = SEs[42:71]; 
sigma2.hat = mean(fixed.fx2^2)-mean(fixed.se2^2)
shrunk.fx2 = fixed.fx2*sqrt(sigma2.hat/(sigma2.hat + fixed.se2^2));

# lmer random effects for (logarea.t0|year) 
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
# Simulate data from fitted SHASH model
pars1 = coefs[1:ncol(U)]; 
pars2 = coefs[-(1:ncol(U))];
MLmu = U%*%pars1;  

n_sim <- 500
idaho_sim<-matrix(NA,nrow=nrow(dropD),ncol=n_sim)
for(i in 1:n_sim){
  idaho_sim[,i] <- rSHASHo(n = nrow(dropD), 
                    mu = MLmu, 
					sigma = exp(pars2[1]+pars2[2]*MLmu), 
					nu=pars2[3]+pars2[4]*MLmu + pars2[5]*MLmu^2,
					tau=exp(pars2[6]+pars2[7]*MLmu))
}

# Round the output to approximate the rounding in the data recording 
e = idaho_sim < log(0.6); idaho_sim[e]  = log(0.5); 
e = (idaho_sim > log(0.65))&(idaho_sim < log(0.85)); idaho_sim[e] <- log(0.75); 
e = (idaho_sim > log(0.9))&(idaho_sim < log(1.1)); idaho_sim[e] <- log(1); 

## moments of the real data by size bin
dropD$fitted  = dropD$logarea.t0; ### NOTE! these plots are really structured by initial size

Lkurtosis=function(x) log(kurtosis(x)); 

n_bins = 10
alpha_scale = 0.7
idaho_moments <- dropD %>% 
  arrange(fitted) %>% 
  mutate(size_bin = cut_number(fitted,n=n_bins)) %>% 
  group_by(size_bin) %>% 
  summarise(mean_t1 = mean(logarea.t1),
            sd_t1 = sd(logarea.t1),
            skew_t1 = skewness(logarea.t1),
            kurt_t1 = Lkurtosis(logarea.t1),
            bin_mean = mean(fitted),
            bin_n = n()) 


par(mfrow=c(2,2),mar=c(4,4,2,1),cex.axis=1.3,cex.lab=1.3,mgp=c(2,1,0),bty="l"); 
sim_bin_means=sim_moment_means = matrix(NA,10,n_sim); 
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
add_panel_label("a")

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
add_panel_label("b")

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
add_panel_label("c")

for(i in 1:n_sim){
    sim_moments <- bind_cols(dropD,data.frame(sim=idaho_sim[,i])) %>% 
    arrange(fitted) %>% 
    mutate(size_bin = cut_number(fitted,n=n_bins)) %>% 
    group_by(size_bin) %>% 
    summarise(mean_t1 = Lkurtosis(sim),
              bin_mean = mean(fitted))
	sim_bin_means[,i]=sim_moments$bin_mean; sim_moment_means[,i]=sim_moments$mean_t1; 		  
}
matplot(idaho_moments$bin_mean, sim_moment_means,col=alpha("gray",0.5),pch=16,xlab="Mean size t0",ylab="Log kurtosis(Size t1)",cex=1.4); 
points(idaho_moments$bin_mean, idaho_moments$kurt_t1,pch=1,lwd=2,col=alpha("red",alpha_scale),cex=1.4)
points(idaho_moments$bin_mean, apply(sim_moment_means,1,median),pch=1,lwd=2,col=alpha("black",alpha_scale),cex=1.4)
add_panel_label("d")

dev.copy2pdf(file="../manuscript/figures/BinnedConditionalMoments.pdf")

#############################################################################
# Ditto, for quantiles 
#############################################################################
source("quantileComparePlot.R"); 

dev.new(width=7,height=9); 
out=quantileComparePlot(dropD$logarea.t0,dropD$logarea.t1,idaho_sim,10);

dev.copy2pdf(file="../manuscript/figures/BinnedConditionalQuantiles.pdf")

########################################################################################
#   Simulation: how well do we recover known random effects? 
#
#   Use simulated replicates as "data" and compare shrinkage estimates with the "truth". 
#   Key question: how well does shrinkage do at getting the variance right? 
#########################################################################################

##### Re-do simulations, so there is no rounding of small values 
pars1 = coefs[1:ncol(U)]; pars2 = coefs[-(1:ncol(U))]; MLmu = U%*%pars1;  
n_sim <- 250 
idaho_sim2<-matrix(NA,nrow=nrow(dropD),ncol=n_sim)
for(i in 1:n_sim){
  idaho_sim2[,i] <- rSHASHo(n = nrow(dropD), 
                    mu = MLmu, 
					sigma = exp(pars2[1]+pars2[2]*MLmu), 
					nu=pars2[3]+pars2[4]*MLmu + pars2[5]*MLmu^2,
					tau=exp(pars2[6]+pars2[7]*MLmu))
}

LogLik=function(pars,response,U){
	pars1 = pars[1:ncol(U)]; pars2=pars[-(1:ncol(U))];
	mu = U%*%pars1;  
	val = dSHASHo(response, mu=mu,
	sigma=exp(pars2[1]+pars2[2]*mu),
	nu = pars2[3]+pars2[4]*mu + pars2[5]*mu^2,
	tau = exp(pars2[6]+pars2[7]*mu), log=TRUE)
	return(val); 
}

#### Estimate the random effects in various ways, compare the results 
shrinkRanIntercept = shrinkRanSlope = matrix(NA,30,250); 
shrinkRanIntercept2 = shrinkRanSlope2 = matrix(NA,30,250); 
fixRanIntercept = fixRanSlope = matrix(NA,30,250); 
lmerRanIntercept = lmerRanSlope = matrix(NA,30,250); 

# Model matrix for ML fits 
U=model.matrix(~year + logarea.t0:year + I(logarea.t0^2) +  
	W.ARTR + I(W.HECO + W.POSE) + W.PSSP+  I(W.allcov + W.allpts) + Group + Treatment - 1, data=dropD);

# Use the correct weights for lmer fits 
pars1 = coefs[1:ncol(U)]; pars2 = coefs[-(1:ncol(U))];
MLmu = U%*%pars1; MLsigma = exp(pars2[1]+pars2[2]*MLmu); 
MLweights = 1/MLsigma^2

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
dev.copy2pdf(file="../manuscript/figures/ShrinkageTest.pdf"); 


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






