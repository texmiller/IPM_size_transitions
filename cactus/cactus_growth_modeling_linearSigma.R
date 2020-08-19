rm(list=ls(all=TRUE));

setwd("c:/repos/IPM_size_transitions/cactus"); #Steve
# setwd("C:/Users/tm9/Desktop/git local/IPM_size_transitions/cactus"); #Tom

require(car); require(lme4); require(zoo); require(moments); require(mgcv); 
require(gamlss); require(gamlss.tr); require(AICcmodavg); 
require(lmerTest); require(tidyverse); require(maxLik); 
require(actuar); require(lattice); require(grid); require(scales);
require(sgt); require(formatR); require(popbio); require(bbmle)

# function for converting cactus size measurements to volume
volume <- function(h, w, p){
  (1/3)*pi*h*(((w + p)/2)/2)^2
}

invlogit <- function(x){exp(x)/(1+exp(x))}

# log kurtosis function for diagnostics
Lkurtosis=function(x) log(kurtosis(x)); 

# Steve's diagnostics functions
source("../Diagnostics.R")
source("../fitChosenDists.R"); 

# read in data
CYIM_full<-read_csv("cholla_demography_20042018_EDI.csv") %>% 
  ## filter out transplants and drop seed addition plots (which start with letter H)
  ## also drop Year_t==2018 because I don't have 2019 size data (not entered yet). 2018 data still included in 2017-2018 transition.
  ## For now I am dropping the early years from plots T1-T3 because these cause problems with shrinkage parameter estimates.
    filter(Transplant == 0,
         str_sub(Plot,1,1)!="H",
         Year_t > 2007,
         Year_t!=2018) %>% 
  ## convert height, max width, perp width to volume of cone, take natural log
  mutate(vol_t = volume(Height_t,Width_t,Perp_t),
         vol_t1 = volume(Height_t1,Width_t1,Perp_t1),
         plot = as.factor(Plot),
         year_t = as.factor(Year_t),
         ID = interaction(TagID,plot)) %>%
  select(ID,year_t,plot,vol_t,vol_t1,Survival_t1,Goodbuds_t1) %>% 
  ## sort by initial size
  arrange(vol_t) 

## pull out and na.omit size transitions for growth modeling
CYIM_full %>% 
  select(ID,year_t,plot,vol_t,vol_t1) %>% 
  ## drop rows with NAs
  drop_na() -> CYIM

table(CYIM$plot,CYIM$year_t) ## note that plots 7 and 8 started in 2011
## how many unique individuals?
CYIM_full %>% select(ID) %>% unique() %>% summarise(n())
## how many observation-years?
CYIM_full %>% nrow()

###################################################################################################
# Gaussian fits with sigma depending on initial size ---------------------------------------------
###################################################################################################
CYIM_init_models <- list() 
CYIM_init_models[[1]] <- lmer(log(vol_t1) ~ log(vol_t) + (1|year_t) + (1|plot), data=CYIM,REML=F,control=lmerControl(optimizer="bobyqa"))
CYIM_init_models[[2]] <- lmer(log(vol_t1) ~ log(vol_t) + I(log(vol_t)^2) + (1|year_t) + (1|plot), data=CYIM,REML=F,control=lmerControl(optimizer="bobyqa"))

## Now use the iterative re-weighting approach to re-fit these with non-constant variance:
## NegLogLik function to fit variance model for residuals 
## SPE: based on the residual diagnostic plot for a linear model of log(sigma), a quadratic term is added. 
varPars1 = function(pars) {
  return(-sum(dnorm(resids, mean=0, sd=exp(pars[1] + pars[2]*log(vol_t)),log=TRUE)))
}	

pars<-list(); vol_t = CYIM$vol_t; 
for(mod in 1:length(CYIM_init_models)) {
  err = 1; rep=0; 
  while(err > 0.000001) {
    rep=rep+1; model = CYIM_init_models[[mod]];
    resids = residuals(model); 
    out=optim(c(sd(resids),0),varPars1,control=list(maxit=25000)); 
    pars[[mod]]=out$par; 
    new_sigma = exp(pars[[mod]][1] + pars[[mod]][2]*log(vol_t)); new_weights = 1/((new_sigma)^2)
    new_weights = 0.5*(weights(model) + new_weights); # cautious update 
    new_model <- update(model,weights=new_weights); 
    err = weights(model)-weights(new_model); err=sqrt(mean(err^2)); 
    cat(mod,rep,err,"\n") # check on convergence of estimated weights 
    CYIM_init_models[[mod]]<-new_model; 
  }}


######### For a fair AIC comparison, fit all models with the same weights 
aics = unlist(lapply(CYIM_init_models,AIC)); best_model=which(aics==min(aics)); 
best_weights=weights(CYIM_init_models[[best_model]]); 
for(mod in 1:length(CYIM_init_models)) {
  CYIM_init_models[[mod]] <- update(CYIM_init_models[[mod]],weights=best_weights)
  }
AIC(CYIM_init_models[[1]],CYIM_init_models[[2]]) # quadratic model wins 


# Gaussian fits with sigma depending on fitted value ----------------------------------------------------- 
CYIM_lmer_models <- list() 
CYIM_lmer_models[[1]] <- lmer(log(vol_t1) ~ log(vol_t) + (1|year_t) + (1|plot), data=CYIM,REML=F,control=lmerControl(optimizer="bobyqa"))
CYIM_lmer_models[[2]] <- lmer(log(vol_t1) ~ log(vol_t) + I(log(vol_t)^2) + (1|year_t) + (1|plot), data=CYIM,REML=F,control=lmerControl(optimizer="bobyqa"))

## Now use the iterative re-weighting approach to re-fit these with non-constant variance:
## NegLogLik function to fit variance model for residuals 
## SPE: based on the residual diagnostic plot for a linear model of log(sigma), a quadratic term is added. 
varPars = function(pars) {
  return(-sum(dnorm(resids, mean=0, sd=exp(pars[1] + pars[2]*fitted_vals),log=TRUE)))
}	

pars<-list()
for(mod in 1:length(CYIM_lmer_models)) {
  err = 1; rep=0; 
  while(err > 0.00001) {
    rep=rep+1; model = CYIM_lmer_models[[mod]];
    fitted_vals = fitted(model); resids = residuals(model); 
    out=optim(c(sd(resids),0),varPars,control=list(maxit=5000)); 
    pars[[mod]]=out$par; 
    new_sigma = exp(pars[[mod]][1] + pars[[mod]][2]*fitted_vals); new_weights = 1/((new_sigma)^2)
    new_weights = 0.75*(weights(model)) + 0.25*new_weights; # cautious update 
    new_model <- update(model,weights=new_weights); 
    err = weights(model)-weights(new_model); err=sqrt(mean(err^2)); 
    cat(mod,rep,err,"\n") # check on convergence of estimated weights 
    CYIM_lmer_models[[mod]]<-new_model; 
  }}

######### For a fair AIC comparison, fit all models with the same weights 
aics = unlist(lapply(CYIM_lmer_models,AIC)); best_model=which(aics==min(aics)); 
best_weights=weights(CYIM_lmer_models[[best_model]]); 
for(mod in 1:length(CYIM_lmer_models)) {
  CYIM_lmer_models[[mod]] <- update(CYIM_lmer_models[[mod]],weights=best_weights)
  }
AIC(CYIM_lmer_models[[1]],CYIM_lmer_models[[2]]) # The quadratic model is favored.

####### Should sigma depend on initial size or fitted value? 
AIC(CYIM_init_models[[2]],CYIM_lmer_models[[2]]) # Initial size! Delta AIC \approx 32. 

## swap in the initial-size models;
CYIM_lmer_models = CYIM_init_models; 

## Last step is to re-fit with REML=T.
aics = unlist(lapply(CYIM_lmer_models,AIC)); best_model=which(aics==min(aics));
CYIM_lmer_best = CYIM_lmer_models[[best_model]] 
best_weights = weights(CYIM_lmer_best)
CYIM_lmer_best <- lmer(log(vol_t1) ~ log(vol_t) + I(log(vol_t)^2) + (1|year_t) + (1|plot), data=CYIM,weights=best_weights,REML=TRUE) 

##refit the residuals as a function of initial size 
resids = residuals(CYIM_lmer_best)
best_pars <- optim(c(sd(resids),0,0),varPars1,control=list(maxit=5000))

scaledResids = residuals(CYIM_lmer_best)*sqrt(weights(CYIM_lmer_best))
par(mfrow=c(1,2))
plot(fitted(CYIM_lmer_best), scaledResids) 
qqPlot(scaledResids) # really bad in both tails
jarque.test(scaledResids) # normality test: FAILS, P < 0.001 
anscombe.test(scaledResids) # kurtosis: FAILS, P < 0.001 
agostino.test(scaledResids) # skewness: FAILS, P<0.001 

# One last look at the standardized residuals. They are roughly mean zero and unit variance -- 
# so that checks out. But there is negative skew and excess kurtosis, especially at large sizes. 
px = fitted(CYIM_lmer_best); py=scaledResids; 

# print rolling moments figure
par(mfrow=c(2,2),bty="l",mar=c(4,4,2,1),mgp=c(2.2,1,0),cex.axis=1.4,cex.lab=1.4);   
z = rollMomentsNP(px,py,windows=12,smooth=TRUE,scaled=TRUE) 
## Very clear quadratic trend in sigma! The model for sigma needs to be quadratic. 
