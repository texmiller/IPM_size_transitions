rm(list=ls(all=TRUE));

setwd("c:/repos/IPM_size_transitions/creosote"); #Steve
# setwd("C:/Users/tm9/Desktop/git local/IPM_size_transitions/creosote"); #Tom

require(car); require(lme4); require(zoo); require(moments); require(mgcv); 
require(gamlss); require(gamlss.tr); require(AICcmodavg); 
require(lmerTest); require(tidyverse); require(maxLik); 
require(actuar); require(lattice); require(grid); require(scales);
require(sgt); require(formatR); require(popbio); require(bbmle)

# misc functions
volume <- function(h, w, p){
  (1/3)*pi*h*(((w + p)/2)/2)^2
}

invlogit <- function(x){exp(x)/(1+exp(x))}
 
# log kurtosis function for diagnostics
Lkurtosis=function(x) log(kurtosis(x)); 

# Steve's diagnostics functions
source("../Diagnostics.R")

## read in data for Larrea tridentata (LATR)
LATR <- read.csv("creosote_growth_density.csv") %>% 
  #calculate volume
  mutate(vol_t = volume(max.ht_t,max.w_t,perp.w_t),
         vol_t1 = volume(max.ht_t1,max.w_t1,perp.w_t1),
         #standardize weighted density to mean zero
         d.stand = (weighted.dens - mean(weighted.dens, na.rm = TRUE)) / sd(weighted.dens, na.rm = TRUE),
         #create unique transect as interaction of transect and site
         unique.transect = interaction(transect, site)) %>% 
  drop_na(vol_t,vol_t1)

# closer look at the outliers
LATR %>% mutate(change=(vol_t1)-(vol_t)) %>% 
  filter(change > quantile(change,0.99,na.rm=T) | change < quantile(change,0.01,na.rm=T)) %>% 
  select(X,site,transect,designated.window,plant,year_t,max.ht_t,max.w_t,perp.w_t,max.ht_t1,max.w_t1,perp.w_t1)
## FPS-3-500 and MOD-3-200 both look like errors. These are the suspect row numbers
outliers <- c(617,684,686,688,882)
LATR %>% filter(!(X %in% outliers)) -> LATR

# first look at size transitions
plot(log(LATR$vol_t),log(LATR$vol_t1))

############################################################################
# Gaussian fits and model selection using mgcv 
############################################################################
LATR_gam_models=list()
## Pilot fits, where sigma depends on initial size only
LATR$fitted_vals = log(LATR$vol_t); 
LATR_gam_models[[1]] <- gam(list(log(vol_t1) ~s(log(vol_t)) + s(unique.transect,bs="re"), ~s(fitted_vals)), 
                            data=LATR, gamma=1.4, family=gaulss())
LATR_gam_models[[2]] <- gam(list(log(vol_t1) ~s(log(vol_t)) + s(d.stand) + s(unique.transect,bs="re"), ~s(fitted_vals)), 
                            data=LATR, gamma=1.4, family=gaulss())                
LATR_gam_models[[3]] <- gam(list(log(vol_t1) ~s(log(vol_t)) + s(d.stand) + d.stand:log(vol_t) + s(unique.transect,bs="re"), ~s(fitted_vals)), 
                            data=LATR, gamma=1.4, family=gaulss())   
## these models will be iterated to fit sigma as f(fitted value)
LATR_gam_models[[4]] <- gam(list(log(vol_t1) ~s(log(vol_t)) + s(unique.transect,bs="re"), ~s(fitted_vals)), 
                            data=LATR, gamma=1.4, family=gaulss())
LATR_gam_models[[5]] <- gam(list(log(vol_t1) ~s(log(vol_t)) + s(d.stand) + s(unique.transect,bs="re"), ~s(fitted_vals)), 
                            data=LATR, gamma=1.4, family=gaulss())                
LATR_gam_models[[6]] <- gam(list(log(vol_t1) ~s(log(vol_t)) + s(d.stand) + d.stand:log(vol_t) + s(unique.transect,bs="re"), ~s(fitted_vals)), 
                            data=LATR, gamma=1.4, family=gaulss())  
for(mod in 4:6) {
  fitGAU = LATR_gam_models[[mod]]
  fitted_all = predict(fitGAU,type="response",data=LATR);                  
  fitted_vals = new_fitted_vals = fitted_all[,1]; 
  weights = fitted_all[,2]; # what I call "weights" here are 1/sigma values; see ?gaulss for details.

  err=100; k=0; 
  while(err>10^(-6)) {
    LATR$fitted_vals = new_fitted_vals; 
    fitGAU <- update(fitGAU); 
    fitted_all = predict(fitGAU,type="response",data=LATR);   
    new_fitted_vals = fitted_all[,1]; new_weights = fitted_all[,2];
    err = weights - new_weights; err=sqrt(mean(err^2)); 
    weights = new_weights; 
    k=k+1; cat(k,err,"\n"); 
  }   
  LATR_gam_models[[mod]] =  fitGAU;
}

AIC(LATR_gam_models[[1]]); AIC(LATR_gam_models[[2]]); AIC(LATR_gam_models[[3]]); 
AIC(LATR_gam_models[[4]]); AIC(LATR_gam_models[[5]]); AIC(LATR_gam_models[[6]]); 
## Models 2 and 3 are top and basically equivalent. I''ll proceed with 2, the simpler one
## It surprises me that initial size was the better covariated than initial value, because I expected the
## additional effect of density to affect sigma. But maybe size and density are correlated.
plot(LATR$d.stand,log(LATR$vol_t)) ## yes, they are - onward
LATR_gam_model <- LATR_gam_models[[2]]; 
LATR$fitted_vals = new_fitted_vals; 


##################################################################  
# Extract values of the fitted splines to explore their properties 
##################################################################
fitted_terms = predict(LATR_gam_model,type="terms");   

##### effect of initial size on mean of final size 
plot(log(LATR$vol_t), fitted_terms[,1]); 
mean_fit1 = lm(fitted_terms[,1]~log(LATR$vol_t)); ## linear initial size is a perfect fit

##### effect of d.stand on mean of final size 
plot(LATR$d.stand, fitted_terms[,2]); ## complicated 
d.stand_fit1 = lm(fitted_terms[,2]~LATR$d.stand); 
d.stand_fit2 = lm(fitted_terms[,2]~LATR$d.stand + I(LATR$d.stand^2)); # R^2 = 0.88 
d.stand_fit3 = lm(fitted_terms[,2]~LATR$d.stand + I(LATR$d.stand^2) + I(LATR$d.stand^3)); # R^2 = 0.99
points(LATR$d.stand,fitted(d.stand_fit3), col="red"); 

##### sigma versus fitted values 
fitted_all = predict(LATR_gam_model,type="response");   
sigma.hat = 1/fitted_all[,2]; 
plot(LATR$fitted_vals, log(sigma.hat)) ; 
sigma_fit1 = lm(log(sigma.hat)~LATR$fitted_vals); ## linear in fitted value is a perfect fit 

##################################################################  
# Inspect scaled residuals to evaluate the pilot model: FAILS 
##################################################################
scaledResids = residuals(LATR_gam_model,type="response")/sigma.hat;  # note the 'type' argument is needed
par(mfrow=c(1,2))
plot(fitted_all[,1], scaledResids,xlab="Fitted values", ylab="Scaled residuals") 

qqPlot(scaledResids) # really bad in both tails
jarque.test(scaledResids) # normality test: FAILS, P < 0.001 
anscombe.test(scaledResids) # kurtosis: FAILS, P < 0.001 
agostino.test(scaledResids) # skewness: FAILS, P<0.001 

######## rolling NP moments diagnostic: skew is small but variable, tails are fat. 
px = LATR$fitted_vals; py=scaledResids; 
z = rollMomentsNP(px,py,windows=8,smooth=TRUE,scaled=TRUE) 

################################################################################
## Finding a better distribution. Must at least allow positive excess kurtosis 
#################################################################################
source("../fitChosenDists.R")

cand_dist=c("GT","JSU","LO","SEP3","SEP4", "SHASHo", "TF") ##NET not behaving, SEP1 takes too long and appears unstable
n_bins <- 6
select_dist <- tibble(init_size = log(LATR$vol_t),
                      scale_resid = scaledResids) %>% 
  mutate(bin = as.integer(cut_number(init_size,n_bins)),
         best_dist = NA,
         secondbest_dist = NA,
         aic_margin = NA) 
for(b in 1:n_bins){
  bin_fit <- gamlssMaxlik(y=select_dist$scale_resid[select_dist$bin==b],DIST=cand_dist)
  select_dist$best_dist[select_dist$bin==b] <- cand_dist[which(bin_fit$aics==sort(bin_fit$aics)[1])]#names(bin_fit$fits[1])
  select_dist$secondbest_dist[select_dist$bin==b] <- cand_dist[which(bin_fit$aics==sort(bin_fit$aics)[2])]#names(bin_fit$fits[2])
  select_dist$aic_margin[select_dist$bin==b] <- sort(bin_fit$aics)[2] - sort(bin_fit$aics)[1]
}
select_dist %>% 
  group_by(bin) %>% 
  summarise(n_bin = n(),
            best_dist = unique(best_dist),
            secondbest_dist = unique(secondbest_dist),
            aic_margin = unique(aic_margin))
## t family is best or second best for almost every bin

## visualize TF parameters in relation to initial size bin
LATR_bin_fit <-LATR %>% 
  mutate(init_size = log(LATR$vol_t),
         bin = as.integer(cut_number(init_size,n_bins))) %>% 
  mutate(mu=NA, sigma=NA,nu=NA)
for(b in 1:n_bins){
  bin_fit <- gamlssMaxlik(y=log(LATR_bin_fit$vol_t1[LATR_bin_fit$bin==b]),DIST="TF")
  LATR_bin_fit$mu[LATR_bin_fit$bin==b] <- bin_fit$out[[1]]$estimate["eta.mu"] ## identity link
  LATR_bin_fit$sigma[LATR_bin_fit$bin==b] <- exp(bin_fit$out[[1]]$estimate["eta.sigma"]) ## log link
  LATR_bin_fit$nu[LATR_bin_fit$bin==b] <- exp(bin_fit$out[[1]]$estimate["eta.nu"]) ## log link
}
LATR_bin_fit %>% 
  group_by(bin) %>% 
  summarise(N = n(),
            mean_size = mean(init_size),
            mu = unique(mu),
            sigma=unique(sigma),
            nu=unique(nu)) -> LATR_bin_fit

pdf("../manuscript/figures/creosote_binned_ST.pdf",height = 4,width = 10,useDingbats = F)
par(mfrow=c(1,3),bty="l",mar=c(4,4,2,1),mgp=c(2.2,1,0),cex.axis=1.4,cex.lab=1.4);
plot(LATR_bin_fit$mean_size,LATR_bin_fit$mu,xlab="Initial size",ylab=expression(paste("Skewed t parameter  ", mu )),type="b",pch=16)
plot(LATR_bin_fit$mean_size,LATR_bin_fit$sigma,xlab="Initial size",
     ylab=expression(paste("Skewed t  parameter  ", sigma)),type="b",pch=16)
plot(LATR_bin_fit$mean_size,LATR_bin_fit$nu,xlab="Initial size",
     ylab=expression(paste("Skewed t  parameter  ", nu )),type="b",pch=16)
dev.off()
