### move to the right local directory 
tom = "C:/Users/tm9/Dropbox/github/IPM_size_transitions"
steve = "c:/repos/IPM_size_transitions" 
home = ifelse(Sys.info()["user"] == "Ellner", steve, tom)
setwd(home); setwd("creosote"); 

library(lme4)
library(mgcv)
library(tidyverse)
library(quantreg)
library(gamlss.dist)
library(maxLik)
library(bbmle)
library(popbio)
library(moments)
library(minpack.lm)
library(sqldf)
library(SuppDists)
library(Rage)

## functions
Q.mean<-function(q.25,q.50,q.75){(q.25+q.50+q.75)/3}
Q.sd<-function(q.25,q.75){(q.75-q.25)/1.35}
Q.skewness<-function(q.10,q.50,q.90){(q.10 + q.90 - 2*q.50)/(q.90 - q.10)}
Q.kurtosis<-function(q.05,q.25,q.75,q.95){
  qN = qnorm(c(0.05,0.25,0.75,0.95))
  KG = (qN[4]-qN[1])/(qN[3]-qN[2])
  return(((q.95-q.05)/(q.75-q.25))/KG - 1)
}
invlogit<-function(x){exp(x)/(1+exp(x))}

## grab the creosote demography data from github
source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/04_CDataPrep.R")
LATR_full <- CData %>% 
  mutate(unique.transect = interaction(transect, site)) %>%
  mutate(log_volume_t = log(volume_t),
         log_volume_t1 = log(volume_t1),
         dens_scaled = weighted.dens/100)

## color pallete for growth figure -- thanks, colorbrewer
dens_pallete<-c("#bdd7e7","#6baed6","#2171b5")

# Prepare a data subset for growth that drops rows missing either t or t1 size data
## update to the growth data -- dropping a few unbelievable outliers following additional QA/QC
outliers<-c("MOD.2.50.3.2016","MOD.3.200.1.2015","MOD.3.200.1.2014","PDC.2.0.5.2014")
LATR_grow <- LATR_full %>% 
  ## this ID will help us drop outliers below
  mutate(ID=interaction(site,transect,actual.window,plant,year_t)) %>% 
  drop_na(volume_t, volume_t1) %>% 
  ##need to scale weighted density because 1st and 2nd order variables were hugely different in range
  filter(!ID%in%outliers) %>% 
  ##bin density variation to make a nice plot
  mutate(dens_bin=cut_interval(dens_scaled,n=length(dens_pallete),labels=F),
         dens_col=dens_pallete[dens_bin],
         size_bin=cut_interval(log_volume_t,n=length(dens_pallete),labels=F),
         size_col=dens_pallete[size_bin])

e = order(LATR_grow$log_volume_t); 
LATR_grow = LATR_grow[e,]; 

# fit candidate gaussian growth models
LATR_GAU<-list()
LATR_GAU[[1]] <- lmer(log_volume_t1 ~ log_volume_t + (1|unique.transect), data=LATR_grow, REML=F)
LATR_GAU[[2]] <- lmer(log_volume_t1 ~ log_volume_t + dens_scaled + (1|unique.transect), data=LATR_grow, REML=F)
LATR_GAU[[3]] <- lmer(log_volume_t1 ~ log_volume_t*dens_scaled + (1|unique.transect), data=LATR_grow, REML=F)
LATR_GAU[[4]] <- lmer(log_volume_t1 ~ log_volume_t + dens_scaled + I(dens_scaled^2) + (1|unique.transect), data=LATR_grow, REML=F)
LATR_GAU[[5]] <- lmer(log_volume_t1 ~ log_volume_t*dens_scaled + log_volume_t*I(dens_scaled^2) + (1|unique.transect), data=LATR_grow, REML=F)
AICctab(LATR_GAU,sort=F)

## now use iterative re-weighting to fit sd as function of expected value
sdloglik = function(pars) {
  dnorm(resids, mean=0, sd=exp(pars[1]+pars[2]*fitted_vals),log=TRUE)
}	
pars<-list()
for(mod in 1:length(LATR_GAU)) {
  err = 1; rep=0; 
  while(err > 0.000001) {
    rep=rep+1
    model = LATR_GAU[[mod]]
    fitted_vals = fitted(model)
    resids = residuals(model) 
    out=maxLik(logLik=sdloglik,start=c(exp(sd(resids)),0))
    pars[[mod]]=out$estimate 
    new_sigma = exp(pars[[mod]][1] + pars[[mod]][2]*fitted_vals)
    new_weights = 1/((new_sigma)^2)
    new_weights = 0.5*(weights(model) + new_weights) # cautious update 
    new_model <- update(model,weights=new_weights) 
    err = weights(model)-weights(new_model)
    err=sqrt(mean(err^2))
    cat(mod,rep,err,"\n") # check on convergence of estimated weights 
    LATR_GAU[[mod]]<-new_model 
  }}

##re-do model selection using the best weights for all models
best_weights<-weights(LATR_GAU[[which.min(AICctab(LATR_GAU,sort=F)$dAICc)]])
for(mod in 1:length(LATR_GAU)) {
  LATR_GAU[[mod]]<-update(LATR_GAU[[mod]],weights=best_weights)
}

AICctab(LATR_GAU,sort=F); ## Model 4 wins by a hair 

## finally, re-fit with REML=T and best weights
LATR_GAU_best<-update(LATR_GAU[[4]],weights=best_weights,REML=T)
best_weights<-weights(LATR_GAU_best)
LATR_grow$GAU_fitted <- fitted(LATR_GAU_best)
LATR_grow$GAU_resids <- residuals(LATR_GAU_best)
LATR_grow$GAU_scaled_resids <- LATR_grow$GAU_resids*sqrt(best_weights) ##sqrt(weights)=1/sd

## should be mean zero unit variance
mean(LATR_grow$GAU_scaled_resids);
sd(LATR_grow$GAU_scaled_resids);

###################################################################################
############### Residual diagnostics on the lmer fit with nonconstant variance.  
###################################################################################
library(car); 

gvals = seq(5,20,by=1);
pvals = numeric(length(gvals)); 

for(j in 1:length(pvals)) {
    e = order(LATR_grow$GAU_fitted); 
    sresids = LATR_grow$GAU_scaled_resids[e]; 
    indexx = c(1:length(sresids)); 
    u = gvals[j]*indexx/(1 + max(indexx)); 
    u = floor(u); u = factor(u); 
    out = leveneTest(sresids,u,center=median)
    pvals[j] = as.numeric(out["Pr(>F)"][[1]][1]) 
   
}    

###################################################################################
############### Residual diagnostics on gam fit, family = gaulss.  
###################################################################################

fig_gaulss <- gam(list(log_volume_t1~log_volume_t + s(dens_scaled) + s(unique.transect,bs="re"),~s(log_volume_t)), 
    family="gaulss", data=LATR_grow, method="REML",gamma=1) 

## Now use iterative re-fitting to fit gam model with SD=f(fitted)
  fitted_all = predict(fit_gaulss,type="response",data=LATR_grow);                  
  new_fitted_vals = fitted_all[,1];
  LATR_grow$fitted_vals = new_fitted_vals; 
  weights = fitted_all[,2]; # what I call "weights" here are 1/sigma values; see ?gaulss for details.
  fit_gaulss = gam(list(log_volume_t1~log_volume_t + s(dens_scaled) + s(unique.transect,bs="re"),~s(fitted_vals)), 
    family="gaulss", data=LATR_grow, method="REML",gamma=1.2) 
  
  err=100; k=0; 
  while(err>10^(-8)) {
    LATR_grow$fitted_vals = new_fitted_vals; 
    fit_gaulss <- gam(list(log_volume_t1~log_volume_t + s(dens_scaled) + s(unique.transect,bs="re"),~s(fitted_vals)), 
    family="gaulss", data=LATR_grow, method="REML",gamma=1.2)  
    fitted_all = predict(fit_gaulss,type="response",data=LATR_grow);   
    new_fitted_vals = fitted_all[,1]; new_weights = fitted_all[,2];
    err = weights - new_weights; err=sqrt(mean(err^2)); 
    weights = new_weights; 
    k=k+1; cat(k,err,"\n"); 
  }   
out = predict(fit_gaulss,type="response"); 
plot(out[,1],1/out[,2]); 

pvals2 = numeric(length(gvals)); 
for(j in 1:length(pvals)) {
    x = predict(fit_gaulss,type="response",data=LATR_grow)[,1];  
    sresids = residuals(fit_gaulss,type="pearson") ## equivalent to our scaling 
    e = order(x); sresids = sresids[e]; 
    indexx = c(1:length(sresids)); 
    u = gvals[j]*indexx/(1 + max(indexx)); 
    u = floor(u); u = factor(u); 
    out = leveneTest(sresids,u,center=median); 
    pvals2[j] = as.numeric(out["Pr(>F)"][[1]][1]) 
}  
matplot(gvals,cbind(pvals,pvals2),type="o",lty=1,col=c("black","red")); 
abline(h=0.05,lty=2, col="blue"); 




