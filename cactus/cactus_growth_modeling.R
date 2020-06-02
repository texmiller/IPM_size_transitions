##############################################################################
# Growth kernel modeling for Cylindropuntia imbricata - modeled off of Steve's PSSPgrowthModeling-ML.R script. 
#
# Data collection described in Miller et al. 2009, Ohm and Miller 2014, Czachura and Miller 2020, and elsewhere  
# This data set includes 2018 data, which has not been previously published. 
#
# Last update: 2 June, 2020
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

CYIM<-read_csv("cholla_demography_20042018_EDI.csv") %>% 
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
table(CYIM$plot,CYIM$year_t) ## first four years has fewer plants from fewer plots

########################################################################## 
## 2. Gaussian fits: log transformation, constant variance, test quadratic term for initial size
## There are no other fixed effects screaming out to be included in the models, but we could test for climate drivers
## and ant defense effects, as we have done elsewhere
########################################################################## 
CYIM_lmer_models <- list()
## random effects are intercept-only - convergence troubles otherwise
CYIM_lmer_models[[1]] <- lmer(log(vol_t1) ~ log(vol_t) + (1|year_t) + (1|plot), data=CYIM,REML=F,control=lmerControl(optimizer="bobyqa"))
CYIM_lmer_models[[2]] <- lmer(log(vol_t1) ~ log(vol_t) + I(log(vol_t)^2) + (1|year_t) + (1|plot), data=CYIM,REML=F,control=lmerControl(optimizer="bobyqa"))

########################################################################## 
## 3. Use iterative re-weighting to fit with nonconstant variance,  
## then AIC for model selection.  
##########################################################################

## NegLogLik function to fit variance model for residuals 
varPars = function(pars) {
  return(-sum(dnorm(resids, mean=0, sd=exp(pars[1] + pars[2]*fitted_vals),log=TRUE)))
}	

for(mod in 1:length(CYIM_lmer_models)) {
  err = 1; rep=0; 
  while(err > 0.000001) {
    rep=rep+1; model = CYIM_lmer_models[[mod]];
    fitted_vals = fitted(model);resids = residuals(model); 
    out=optim(c(sd(resids),0),varPars,control=list(maxit=5000)); 
    pars=out$par; 
    new_sigma = exp(pars[1] + pars[2]*fitted_vals); new_weights = 1/((new_sigma)^2)
    new_weights = 0.5*(weights(model) + new_weights); # cautious update 
    new_model <- update(model,weights=new_weights); 
    err = weights(model)-weights(new_model); err=sqrt(mean(err^2)); 
    cat(mod,rep,err,"\n") # check on convergence of estimated weights 
    CYIM_lmer_models[[mod]]<-new_model; 
  }}
aictab(CYIM_lmer_models) ## no strong support for quadratic term

######### For a fair AIC comparison, fit all models with the same weights 
aics = unlist(lapply(CYIM_lmer_models,AIC)); best_model=which(aics==min(aics)); 
best_weights=weights(CYIM_lmer_models[[best_model]]); 
for(mod in 1:length(CYIM_lmer_models)) {CYIM_lmer_models[[mod]] <- update(CYIM_lmer_models[[mod]],weights=best_weights)}
aictab(CYIM_lmer_models); # still support for model 1

######### Here's the best Gaussian model ########################################
CYIM_lmer_best = CYIM_lmer_models[[best_model]]; 
summary(CYIM_lmer_best); 

### refit with REML, as recommended for estimating random effects 
best_weights = weights(CYIM_lmer_best)
CYIM_lmer_best <- lmer(log(vol_t1) ~ log(vol_t) + (1|year_t) + (1|plot), data=CYIM,weights=best_weights,REML=TRUE) 
plot(log(vol_t1) ~ log(vol_t),data=CYIM)
abline(fixef(CYIM_lmer_best))

######################################################################
# Interrogate the scaled residuals - Gaussian? NO. 
###################################################################### 
plot(fitted(CYIM_lmer_best), 1/sqrt(weights(CYIM_lmer_best)),xlab="Fitted",ylab="Estimated residual Std Dev"); 
scaledResids = residuals(CYIM_lmer_best)*sqrt(weights(CYIM_lmer_best))
plot(fitted(CYIM_lmer_best), scaledResids); 

qqPlot(scaledResids); # really bad in both tails

jarque.test(scaledResids) # normality test: FAILS, P < 0.001 
anscombe.test(scaledResids) # kurtosis: FAILS, P < 0.001 
agostino.test(scaledResids) # skewness: FAILS, P<0.001 


########################################################################
## Rollapply diagnostics on the scaled residuals 
########################################################################
px = fitted(CYIM_lmer_best); py=scaledResids; 

graphics.off(); dev.new(width=8,height=6); 
par(mfrow=c(2,2),bty="l",mar=c(4,4,2,1),mgp=c(2.2,1,0),cex.axis=1.4,cex.lab=1.4);   
z = rollMoments(px,py,windows=10,smooth=TRUE,scaled=TRUE) 
dev.copy2pdf(file="../manuscript/figures/RollingMomentsCYIM.pdf") 
## scaled residuals look about mean zero and unit variance (but maybe some weird trends with fitted values)
## clear negative skew and excess kurtosis

## what if we did not do the interative re-weighting and just examined the residuals of the original gaussian model?
orig_mod <- lmer(log(vol_t1) ~ log(vol_t) + (1|year_t) + (1|plot), data=CYIM,REML=F,control=lmerControl(optimizer="bobyqa"))
graphics.off(); dev.new(width=8,height=6); 
par(mfrow=c(2,2),bty="l",mar=c(4,4,2,1),mgp=c(2.2,1,0),cex.axis=1.4,cex.lab=1.4);
rollMoments(fitted(orig_mod),residuals(orig_mod),windows=10,smooth=TRUE,scaled=TRUE) 
dev.copy2pdf(file="../manuscript/figures/RollingMomentsCYIM_rawresids.pdf") 

########################################################################
## what is the likely distribution of the residuals?
## because there are trends wrt to fitted value, I will do this in slices. 
########################################################################
## I should re-write this code as a loop because it runs fitDist 3 times more than necessary
n_bins <- 8
select_dist <- tibble(fit_best = fitted(CYIM_lmer_best),
                      scale_resid = residuals(CYIM_lmer_best)*sqrt(weights(CYIM_lmer_best))) %>% 
  mutate(bin = as.integer(cut_number(fit_best,n_bins))) %>% 
  group_by(bin) %>% 
  mutate(best_dist = names(fitDist(scale_resid,type="realline")$fits[1]),
         secondbest_dist = names(fitDist(scale_resid,type="realline")$fits[2]),
         aic_margin = fitDist(scale_resid,type="realline")$fits[2] - fitDist(scale_resid,type="realline")$fits[1]) %>% 
  summarise(n_bin = n(),
            best_dist = unique(best_dist),
            secondbest_dist = unique(secondbest_dist),
            aic_margin = unique(aic_margin))
# a little bit of everything but some version of the skewed t is commonly favored, though the SHASH is strongly favored for the smaller bins

## again again, what if I used the raw residuals? -- still support for the skewed t
select_dist_orig <- tibble(fit_best = fitted(orig_mod),
                      scale_resid = residuals(orig_mod)*sqrt(weights(orig_mod))) %>% 
  mutate(bin = as.integer(cut_number(fit_best,n_bins))) %>% 
  group_by(bin) %>% 
  mutate(best_dist = names(fitDist(scale_resid,type="realline")$fits[1]),
         secondbest_dist = names(fitDist(scale_resid,type="realline")$fits[2]),
         aic_margin = fitDist(scale_resid,type="realline")$fits[2] - fitDist(scale_resid,type="realline")$fits[1]) %>% 
  summarise(n_bin = n(),
            best_dist = unique(best_dist),
            secondbest_dist = unique(secondbest_dist),
            aic_margin = unique(aic_margin))

########################################################################
## Fit parameters of the skewed t to binned size data -- to get a visual sense of relationships between mu and other params
########################################################################
CYIM_bin_fit <-CYIM %>% 
  mutate(fitted = fitted(CYIM_lmer_best),
         bin = as.integer(cut_number(fitted,n_bins))) %>% 
  mutate(mu=NA, sigma=NA,nu=NA,tau=NA)
for(b in 1:n_bins){
  ## I get the most stable tau estimates with ST3
  bin_fit <- gamlssML(log(CYIM_bin_fit$vol_t1[CYIM_bin_fit$bin==b]) ~ 1,family="ST3")
  CYIM_bin_fit$mu[CYIM_bin_fit$bin==b] <- bin_fit$mu
  CYIM_bin_fit$sigma[CYIM_bin_fit$bin==b] <- bin_fit$sigma
  CYIM_bin_fit$nu[CYIM_bin_fit$bin==b] <- bin_fit$nu
  CYIM_bin_fit$tau[CYIM_bin_fit$bin==b] <- bin_fit$tau
}
CYIM_bin_fit %>% 
  group_by(bin) %>% 
  summarise(N = n(),
            mean_fitted = mean(fitted),
            mu = unique(mu),
            sigma=unique(sigma),
            nu=unique(nu),
            tau=unique(tau)) -> CYIM_bin_fit

par(mfrow=c(2,2),bty="l",mar=c(4,4,2,1),mgp=c(2.2,1,0),cex.axis=1.4,cex.lab=1.4);
## Steve's spline.scatter.smooth() function not working for me and I did not both trying to figure out why
plot(CYIM_bin_fit$mean_fitted,CYIM_bin_fit$mu,xlab="Fitted value",ylab=expression(paste("Location parameter  ", mu )),type="b")
plot(CYIM_bin_fit$mu,CYIM_bin_fit$sigma,xlab=expression(paste("Location parameter  ", mu )),
                      ylab=expression(paste("Scale parameter  ", sigma)),type="b")
plot(CYIM_bin_fit$mu,CYIM_bin_fit$nu,xlab=expression(paste("Location parameter  ", mu )),
                      ylab=expression(paste("Skewness parameter  ", nu )),type="b")
plot(CYIM_bin_fit$mu,CYIM_bin_fit$tau,xlab=expression(paste("Location parameter  ", mu )),
                      ylab=expression(paste("Kurtosis parameter  ", tau)),type="b") ## garbage tau estimate for group 1
