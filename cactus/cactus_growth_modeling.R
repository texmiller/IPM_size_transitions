##############################################################################
# Growth kernel modeling for Cylindropuntia imbricata - modeled off of Steve's PSSPgrowthModeling-ML.R script. 
#
# Data collection described in Miller et al. 2009, Ohm and Miller 2014, Czachura and Miller 2020, and elsewhere  
# This data set includes 2018 data, which has not been previously published. 
#
# Last update: 38 May, 2020
#
##############################################################################

rm(list=ls(all=TRUE));
setwd("./cactus"); 

require(car); require(lme4); require(zoo); require(moments); require(mgcv); 
require(gamlss); require(gamlss.tr); require(AICcmodavg); 
require(lmerTest); require(tidyverse); require(maxLik); 

# Steve's diagnostics functions
source("../Diagnostics.R") 
# function for converting cactus size measurements to volume
volume <- function(h, w, p){
  (1/3)*pi*h*(((w + p)/2)/2)^2
}


##############################################################
# 1. Read in the data and pare down data set to just size transitions
##############################################################

cactus<-read_csv("cholla_demography_20042018_EDI.csv") %>% 
  ## filter out transplants and drop seed addition plots (which start with letter H)
  ## also drop Year_t==2018 because I don't have 2019 size data (not entered yet). 2018 data still included in 2017-2018 transition
    filter(Transplant == 0,
         str_sub(Plot,1,1)!="H",
         Year_t!=2018) %>% 
  ## convert height, max width, perp width to volume of cone, take natural log
  mutate(vol_t = volume(Height_t,Width_t,Perp_t),
         vol_t1 = volume(Height_t1,Width_t1,Perp_t1),
         plot = as.factor(Plot),
         year_t = as.factor(Year_t)) %>%
  select(year_t,plot,vol_t,vol_t1) %>% 
  filter() %>% 
  ## sort by initial size
  arrange(vol_t) %>% 
  ## drop rows with NAs
  drop_na()

## how much coverage do I have across years and plots?
table(cactus$plot,cactus$year_t) ## first four years has fewer plants from fewer plots

########################################################################## 
## 2. Gaussian fits: log transformation, constant variance, test quadratic term for initial size
## There are no other fixed effects screaming out to be included in the models, but we could test for climate drivers
## and ant defense effects, as we have done elsewhere
########################################################################## 
cactus_lmer_models <- list()
## random effects are intercept-only - convergence troubles otherwise
cactus_lmer_models[[1]] <- lmer(log(vol_t1) ~ log(vol_t) + (1|year_t) + (1|plot), data=cactus,REML=F)
cactus_lmer_models[[2]] <- lmer(log(vol_t1) ~ log(vol_t) + I(log(vol_t)^2) + (1|year_t) + (1|plot), data=cactus,REML=F)

########################################################################## 
## 3. Use iterative re-weighting to fit with nonconstant variance,  
## then AIC for model selection.  
##########################################################################

## NegLogLik function to fit variance model for residuals 
varPars = function(pars) {
  return(-sum(dnorm(resids, mean=0, sd=exp(pars[1] + pars[2]*fitted_vals),log=TRUE)))
}	

for(mod in 1:length(cactus_lmer_models)) {
  err = 1; rep=0; 
  while(err > 0.000001) {
    rep=rep+1; model = cactus_lmer_models[[mod]];
    fitted_vals = fitted(model);resids = residuals(model); 
    out=optim(c(sd(resids),0),varPars,control=list(maxit=5000)); 
    pars=out$par; 
    new_sigma = exp(pars[1] + pars[2]*fitted_vals); new_weights = 1/((new_sigma)^2)
    new_weights = 0.5*(weights(model) + new_weights); # cautious update 
    new_model <- update(model,weights=new_weights); 
    err = weights(model)-weights(new_model); err=sqrt(mean(err^2)); 
    cat(mod,rep,err,"\n") # check on convergence of estimated weights 
    cactus_lmer_models[[mod]]<-new_model; 
  }}
aictab(cactus_lmer_models) ## no strong support for quadratic term

######### For a fair AIC comparison, fit all models with the same weights 
aics = unlist(lapply(cactus_lmer_models,AIC)); best_model=which(aics==min(aics)); 
best_weights=weights(cactus_lmer_models[[best_model]]); 
for(mod in 1:length(cactus_lmer_models)) {cactus_lmer_models[[mod]] <- update(cactus_lmer_models[[mod]],weights=best_weights)}
aictab(cactus_lmer_models); # still support for model 1

######### Here's the best Gaussian model ########################################
aics = unlist(lapply(cactus_lmer_models,AIC)); best_model=which(aics==min(aics)); 
cactus_lmer_best = cactus_lmer_models[[best_model]]; 
summary(cactus_lmer_best); 

### refit with REML, as recommended for estimating random effects 
best_weights = weights(cactus_lmer_best)
cactus_lmer_best <- lmer(log(vol_t1) ~ log(vol_t) + (1|year_t) + (1|plot), data=cactus,weights=best_weights,REML=TRUE) 

######################################################################
# Interrogate the scaled residuals - Gaussian? NO. 
###################################################################### 
plot(fitted(cactus_lmer_best), 1/sqrt(weights(cactus_lmer_best)),xlab="Fitted",ylab="Estimated residual Std Dev"); 
log_scaledResids = residuals(cactus_lmer_best)*sqrt(weights(cactus_lmer_best))
plot(fitted(cactus_lmer_best), log_scaledResids); 

qqPlot(log_scaledResids); # really bad in lower tail, not too bad in upper 

jarque.test(log_scaledResids) # normality test: FAILS, P < 0.001 
anscombe.test(log_scaledResids) # kurtosis: FAILS, P < 0.001 
agostino.test(log_scaledResids) # skewness: FAILS, P<0.001 

