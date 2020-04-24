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
require(gamlss); require(gamlss.tr); require(AICcmodavg); 

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
# Tom's work starting here: exploring the idaho data
######################################################################
library(tidyverse)
## drop out these arbitrarily small individuals
arb_small <- unique(c(which(allD$area.t0==0.25),which(allD$area.t1==0.25)))
drop_allD <- allD[-arb_small,]
plot(drop_allD$logarea.t0, drop_allD$logarea.t1) 

## model selection with log transformation
log_models <- list()
log_models[[1]] <- lmer(logarea.t1~logarea.t0+W.ARTR + W.HECO + W.POSE + W.PSSP+  W.allcov + W.allpts + Treatment + 
             (1|Group)+(logarea.t0|year),control=lmerControl(optimizer="bobyqa"),data=drop_allD,REML = F) 

log_models[[2]] <- lmer(logarea.t1~logarea.t0 + W.ARTR  + Treatment + 
             (1|Group)+(logarea.t0|year),control=lmerControl(optimizer="bobyqa"),data=drop_allD,REML = F)

log_models[[3]] <- lmer(logarea.t1~logarea.t0 + Treatment + 
             (1|Group)+(logarea.t0|year),control=lmerControl(optimizer="bobyqa"),data=drop_allD,REML = F)

log_models[[4]] <- lmer(logarea.t1~logarea.t0 + Treatment + 
             (1|Group)+(1|year),control=lmerControl(optimizer="bobyqa"),data=drop_allD,REML = F)
 			  
log_models[[5]] <- lmer(logarea.t1~logarea.t0 + Treatment + 
             (1|Group)+(0+logarea.t0|year),control=lmerControl(optimizer="bobyqa"),data=drop_allD,REML = F)

log_models[[6]] <- lmer(logarea.t1~logarea.t0 + Treatment + 
             (0+logarea.t0|year),control=lmerControl(optimizer="bobyqa"),data=drop_allD,REML = F)
			 
for(mod in 1:6) {
cat(mod,"\n"); 
for(rep in 1:5) {
	log_model = log_models[[mod]];
	fitted_vals = fitted(log_model);
	resids = residuals(log_model); 
	resid_model = lm(log(abs(resids))~fitted_vals); 
	log_models[[mod]] <- update(log_model,weights=exp(-fitted(resid_model))); 
}}
aictab(log_models); 

aics = unlist(lapply(log_models,AIC)); best_model=which(aics==min(aics)); 
best_weights=weights(log_models[[best_model]]); 
for(mod in 1:6) {log_models[[mod]] <- update(log_models[[mod]],weights=best_weights)}
aictab(log_models); 

aics = unlist(lapply(log_models,AIC)); best_model=which(aics==min(aics)); 
log_model = log_models[[best_model]]; 

plot(fitted(log_model), 1/weights(log_model)); 
log_scaledResids = residuals(log_model)*weights(log_model)
plot(fitted(log_model), log_scaledResids); 

qqPlot(log_scaledResids); # bad in just the lower tail

jarque.test(log_scaledResids) # normality test: FAILS, P < 0.001 
anscombe.test(log_scaledResids) # kurtosis: FAILS, P < 0.001 
agostino.test(log_scaledResids) # skewness: FAILS, P<0.001 

## See what family gamlss thinks the scaled residuals conform to. 
log_gamlss_fit <- fitDist(log_scaledResids,type="realline")
log_gamlss_fit$fits
log_gamlss_fit$failed

################################################################################################################
## try fitting gamlss ST2 to log size data. use the fixed effect structure corresponding to log_models[[3]]
## SPE: ST2 is suggested by using fitDist on the scaled residuals, and it gives a better (lower deviance) fit
#################################################################################################################
fit_ST2 <- gamlss(logarea.t1 ~ as.factor(year) + logarea.t0:as.factor(year) + Treatment + Group, 
                  data=drop_allD, family="ST2", method=RS(250),
                  sigma.formula = ~logarea.t0, nu.formula = ~1, tau.formula = ~1)

## well, it fit. does it describe the data well?
## simulate data from best model
n_sim <- 500
idaho_sim<-matrix(NA,nrow=nrow(drop_allD),ncol=n_sim)
for(i in 1:n_sim){
  idaho_sim[,i] <- rST2(n = nrow(drop_allD), 
                        mu = predict(fit_ST2), 
                        sigma = exp(fit_ST2$sigma.coefficients[1] + fit_ST2$sigma.coefficients[2] * drop_allD$logarea.t0),
                        nu = exp(fit_ST2$nu.coefficients),
                        tau = exp(fit_ST2$tau.coefficients))
}


## moments of the real data by size bin
n_bins = 10
alpha_scale = 0.7
idaho_moments <- drop_allD %>% 
  arrange(logarea.t0) %>% 
  mutate(size_bin = cut_number(logarea.t0,n=n_bins)) %>% 
  group_by(size_bin) %>% 
  summarise(mean_t1 = mean(logarea.t1),
            sd_t1 = sd(logarea.t1),
            skew_t1 = skewness(logarea.t1),
            kurt_t1 = kurtosis(logarea.t1),
            bin_mean = mean(logarea.t0),
            bin_n = n()) 

par(mfrow=c(2,2))
plot(idaho_moments$bin_mean,idaho_moments$mean_t1,type="n",xlab="Mean size t0",ylab="mean(Size t1)",ylim=c(1,8))
for(i in 1:n_sim){
  sim_moments <- bind_cols(drop_allD,data.frame(sim=idaho_sim[,i])) %>% 
    arrange(logarea.t0) %>% 
    mutate(size_bin = cut_number(logarea.t0,n=n_bins)) %>% 
    group_by(size_bin) %>% 
    summarise(mean_t1 = mean(sim),
              bin_mean = mean(logarea.t0))
  points(sim_moments$bin_mean,sim_moments$mean_t1,col=alpha("gray",0.5))
}
points(idaho_moments$bin_mean,idaho_moments$mean_t1,cex=2,pch=16,col=alpha("red",alpha_scale))

plot(idaho_moments$bin_mean,idaho_moments$sd_t1,type="n",xlab="Mean size t0",ylab="sd(Size t1)",ylim=c(0,4))
for(i in 1:n_sim){
  sim_moments <- bind_cols(drop_allD,data.frame(sim=idaho_sim[,i])) %>% 
    arrange(logarea.t0) %>% 
    mutate(size_bin = cut_number(logarea.t0,n=n_bins)) %>% 
    group_by(size_bin) %>% 
    summarise(sd_t1 = sd(sim),
              bin_mean = mean(logarea.t0))
  points(sim_moments$bin_mean,sim_moments$sd_t1,col=alpha("gray",0.5))
}
points(idaho_moments$bin_mean,idaho_moments$sd_t1,cex=2,pch=16,col=alpha("red",alpha_scale))

plot(idaho_moments$bin_mean,idaho_moments$skew_t1,type="n",xlab="Mean size t0",ylab="skewness(Size t1)",ylim=c(-10,1))
abline(h=0,col="blue",lty=2) 
for(i in 1:n_sim){
  sim_moments <- bind_cols(drop_allD,data.frame(sim=idaho_sim[,i])) %>% 
    arrange(logarea.t0) %>% 
    mutate(size_bin = cut_number(logarea.t0,n=n_bins)) %>% 
    group_by(size_bin) %>% 
    summarise(skew_t1 = skewness(sim),
              bin_mean = mean(logarea.t0))
  points(sim_moments$bin_mean,sim_moments$skew_t1,col=alpha("gray",0.5))
}
points(idaho_moments$bin_mean,idaho_moments$skew_t1,cex=2,pch=16,col=alpha("red",alpha_scale))

plot(idaho_moments$bin_mean,idaho_moments$kurt_t1,type="n",xlab="Mean size t0",ylab="kurtosis(Size t1)",ylim=c(0,25))
abline(h=0,col="blue",lty=2) 
for(i in 1:n_sim){
  sim_moments <- bind_cols(drop_allD,data.frame(sim=idaho_sim[,i])) %>% 
    arrange(logarea.t0) %>% 
    mutate(size_bin = cut_number(logarea.t0,n=n_bins)) %>% 
    group_by(size_bin) %>% 
    summarise(kurt_t1 = kurtosis(sim),
              bin_mean = mean(logarea.t0))
  points(sim_moments$bin_mean,sim_moments$kurt_t1,col=alpha("gray",0.5))
}
points(idaho_moments$bin_mean,idaho_moments$kurt_t1,cex=2,pch=16,col=alpha("red",alpha_scale))



