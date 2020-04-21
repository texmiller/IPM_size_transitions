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
require(gamlss); require(gamlss.tr); 

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

######################################################################
# Tom's work starting here: exploring the idaho data
######################################################################
library(tidyverse)
plot(allD$area.t0, allD$area.t1) 
hist(log(allD$area.t1)) 
hist(sqrt(allD$area.t1)) 
plot(allD$logarea.t0, allD$logarea.t1) 

## looks like there are arbitrarily small values
table(allD$area.t0);
table(allD$area.t1) ##32 and 28 0.25's in t0 and t1 data, all other values unique
## drop out these arbitrarily small individuals
arb_small <- unique(c(which(allD$area.t0==0.25),which(allD$area.t1==0.25)))
drop_allD <- allD[-arb_small,]
plot(drop_allD$logarea.t0, drop_allD$logarea.t1) 

## for an IPM, one could calculate the probability of graduating out of "small":
smalls <- allD[arb_small,]
length(which(smalls$area.t0==0.25 & smalls$area.t1!=0.25))/nrow(smalls)
## and the size distribution given that you graduate out of small
hist(log(smalls$area.t1[smalls$area.t0==0.25 & smalls$area.t1!=0.25]))

## have another look at rolling moments with smalls dropped
e = order(drop_allD$area.t0); 
## first column, log transform
rollmean=rollapply(drop_allD$logarea.t0[e],50,mean,by=25);
rollvar=rollapply(drop_allD$logarea.t1[e],50,sd,by=25); max(rollvar)/min(rollvar); 
rollkurt=rollapply(drop_allD$logarea.t1[e],50,kurtosis,by=25);
rollskew=rollapply(drop_allD$logarea.t1[e],50,skewness,by=25);
par(mfcol=c(3,2),mar=c(5,5,2,1),cex.axis=1.3,cex.lab=1.3); 
spline.scatter.smooth(rollmean,rollvar,gamma=2,xlab="",ylab="Std Dev",ylim=c(0,max(rollvar))); 
title(main="Log-transformed size"); 
spline.scatter.smooth(rollmean,rollskew,gamma=2,xlab="",ylab="Skewness"); 
abline(h=0,col="blue",lty=2) 
spline.scatter.smooth(rollmean,rollkurt/3-1,gamma=2,xlab="Sqrt area t0",ylab="Excess kurtosis"); 
abline(h=0,col="blue",lty=2)
## second column, sqrt transform
rollmean=rollapply(drop_allD$area.t0[e]^0.5,50,mean,by=25);
rollvar=rollapply(drop_allD$area.t1[e]^0.5,50,sd,by=25); max(rollvar)/min(rollvar);
rollkurt=rollapply(drop_allD$area.t1[e]^0.5,50,kurtosis,by=25);
rollskew=rollapply(drop_allD$area.t1[e]^0.5,50,skewness,by=25);
spline.scatter.smooth(rollmean,rollvar,xlab="",ylab="Std Dev",ylim=c(0,max(rollvar))); 
title(main="Sqrt-transformed size"); 
spline.scatter.smooth(rollmean,rollskew,xlab="",ylab="Skewness"); 
abline(h=0,col="blue",lty=2) 
spline.scatter.smooth(rollmean,rollkurt/3-1,xlab="Sqrt area t0",ylab="Excess kurtosis"); 
abline(h=0,col="blue",lty=2)
## The sqrt transformation *does* look better in terms of skewness and kurtosis.

## model selection with log transformation
log_models <- list()
log_models[[1]] <- lmer(logarea.t1~logarea.t0+W.ARTR + W.HECO + W.POSE + W.PSSP+  W.allcov + W.allpts + Treatment + 
             (1|Group)+(logarea.t0|year),control=lmerControl(optimizer="bobyqa"),data=drop_allD,REML = F) 

log_models[[2]] <- lmer(logarea.t1~logarea.t0 + W.ARTR + W.POSE + W.PSSP + Treatment + 
             (1|Group)+(logarea.t0|year),control=lmerControl(optimizer="bobyqa"),data=drop_allD,REML = F) 

log_models[[3]] <- lmer(logarea.t1~logarea.t0 + Treatment + 
             (1|Group)+(logarea.t0|year),control=lmerControl(optimizer="bobyqa"),data=drop_allD,REML = F) 

log_models[[4]] <- lmer(logarea.t1~logarea.t0 + Treatment + 
             (1|Group)+(1|year),control=lmerControl(optimizer="bobyqa"),data=drop_allD,REML = F) 			  

log_models[[5]] <- lmer(logarea.t1~logarea.t0 + Treatment + 
             (1|Group)+(0+logarea.t0|year),control=lmerControl(optimizer="bobyqa"),data=drop_allD,REML = F)

log_models[[6]] <- lmer(logarea.t1~logarea.t0 + Treatment + 
             (0+logarea.t0|year),control=lmerControl(optimizer="bobyqa"),data=drop_allD,REML = F)
AICtab(log_models) 
# diagnostics of best model
spline.scatter.smooth(drop_allD$logarea.t0,abs(residuals(log_models[[3]]))); 
fit = lm(abs(residuals(log_models[[3]]))~logarea.t0,data=drop_allD); 
abline(fit,col="blue");

scaledResids = residuals(log_models[[3]])/fit$fitted; 
qqPlot(scaledResids); # now bad in just the lower tail

jarque.test(scaledResids) # normality test: FAILS, P < 0.001 
anscombe.test(scaledResids) # kurtosis: FAILS, P < 0.001 
agostino.test(scaledResids) # skewness: FAILS, P<0.001 

## model selection with sqrt transformation -- I am using log area for the random effects like Steve did above
sqrt_models <- list()
sqrt_models[[1]] <- lmer(sqrt(area.t1)~sqrt(area.t0)+W.ARTR + W.HECO + W.POSE + W.PSSP+  W.allcov + W.allpts + Treatment + 
                          (1|Group)+(logarea.t0|year),control=lmerControl(optimizer="bobyqa"),data=drop_allD,REML = F) 

sqrt_models[[2]] <- lmer(sqrt(area.t1)~sqrt(area.t0) + W.ARTR + W.POSE + W.PSSP + Treatment + 
                          (1|Group)+(logarea.t0|year),control=lmerControl(optimizer="bobyqa"),data=drop_allD,REML = F) 

sqrt_models[[3]] <- lmer(sqrt(area.t1)~sqrt(area.t0) + Treatment + 
                          (1|Group)+(logarea.t0|year),control=lmerControl(optimizer="bobyqa"),data=drop_allD,REML = F) 

sqrt_models[[4]] <- lmer(sqrt(area.t1)~sqrt(area.t0) + Treatment + 
                          (1|Group)+(1|year),control=lmerControl(optimizer="bobyqa"),data=drop_allD,REML = F) 			  

sqrt_models[[5]] <- lmer(sqrt(area.t1)~sqrt(area.t0) + Treatment + 
                          (1|Group)+(0+logarea.t0|year),control=lmerControl(optimizer="bobyqa"),data=drop_allD,REML = F)

sqrt_models[[6]] <- lmer(sqrt(area.t1)~sqrt(area.t0) + Treatment + 
                          (0+logarea.t0|year),control=lmerControl(optimizer="bobyqa"),data=drop_allD,REML = F)
AICtab(sqrt_models)
# diagnostics of best model
spline.scatter.smooth(drop_allD$area.t0^0.5,abs(residuals(sqrt_models[[3]]))); 
fit = lm(abs(residuals(sqrt_models[[3]]))~sqrt(area.t0),data=drop_allD); 
abline(fit,col="blue");

scaledResids = residuals(sqrt_models[[3]])/fit$fitted; 
qqPlot(scaledResids); # bad in both tails

jarque.test(scaledResids) # normality test: FAILS, P < 0.001 
anscombe.test(scaledResids) # kurtosis: FAILS, P < 0.001 
agostino.test(scaledResids) # skewness: FAILS, P<0.001 
## both log and sqrt transforms support model 3 (Steve's m2), which is different from Steve's result (support for his m5 = my model6)
## I made several changes including REML=F, dropping arbitrarily small plants, and using same size variable in fixed and random effects.
## The latter change made the biggest difference as far as I can tell. 

## Even with my adjustments we arrive at the same place -- size residuals are not gaussian
## Can I get gamlss to work with the log transform, dropped smalls?

## first see what family gamlss thinks these data conform to, just based on the marginal distribution of size_t1
sqrt_gamlss_fit <- fitDist(drop_allD$area.t1^0.5,type="realplus")
sqrt_gamlss_fit$fits
sqrt_gamlss_fit$failed

log_gamlss_fit <- fitDist(drop_allD$logarea.t1,type="realAll")
log_gamlss_fit$fits
log_gamlss_fit$failed

## try fitting gamlss ST4 to log size data. use the fixed effect structure corresponding to log_models[[3]]
fit_ST4 <- gamlss(logarea.t1 ~ as.factor(year) + logarea.t0:as.factor(year) + Treatment + Group,
                  data=drop_allD, family="ST4",
                  sigma.formula = ~logarea.t0, nu.formula = ~1, tau.formula = ~1)

## well, it fit. does it describe the data well?
## simulate data from best model
n_sim <- 500
idaho_sim<-matrix(NA,nrow=nrow(drop_allD),ncol=n_sim)
for(i in 1:n_sim){
  idaho_sim[,i] <- rST4(n = nrow(drop_allD), 
                        mu = predict(fit_ST4), 
                        sigma = exp(fit_ST4$sigma.coefficients[1] + fit_ST4$sigma.coefficients[2] * drop_allD$logarea.t0),
                        nu = exp(fit_ST4$nu.coefficients),
                        tau = exp(fit_ST4$tau.coefficients))
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



