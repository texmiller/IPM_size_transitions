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

### SPE: how good is that log-linear model for the standard deviation? 
LATR_GAU[[6]] <- gam(list(log_volume_t1~log_volume_t + s(dens_scaled) + s(unique.transect,bs="re"),~s(log_volume_t)), 
    family="gaulss", data=LATR_grow, method="ML",gamma=1.4) 
## TM: model 6 has smooth terms for both density and SD(size). 
## Try a model with only smooth term for SD, same linear model for mean
LATR_GAU[[7]] <- gam(list(log_volume_t1~log_volume_t + dens_scaled + I(dens_scaled^2) + s(unique.transect,bs="re"),~s(log_volume_t)), 
                     family="gaulss", data=LATR_grow, method="ML",gamma=1.4) 

### SPE: not very good. The gam winds by 33 AIC units. 
AICctab(LATR_GAU,sort=F); 

### Now use iterative re-weighting to fit gam model with SD=f(fitted)
for(i in 6:7){
  fit_gaulss = LATR_GAU[[i]]
  fitted_all = predict(fit_gaulss,type="response",data=LATR_grow);                  
  new_fitted_vals = fitted_all[,1];
  LATR_grow$fitted_vals = new_fitted_vals; 
  weights = fitted_all[,2]; # what I call "weights" here are 1/sigma values; see ?gaulss for details.
  if(i==6){fit_gaulss = gam(list(log_volume_t1~log_volume_t + s(dens_scaled) + s(unique.transect,bs="re"),~s(fitted_vals)), 
    family="gaulss", data=LATR_grow, method="ML",gamma=1.4) }
  if(i==7){fit_gaulss = gam(list(log_volume_t1~log_volume_t + dens_scaled + I(dens_scaled^2) + s(unique.transect,bs="re"),~s(fitted_vals)), 
                            family="gaulss", data=LATR_grow, method="ML",gamma=1.4) }  
  err=100; k=0; 
  while(err>10^(-6)) {
    LATR_grow$fitted_vals = new_fitted_vals; 
    #fit_gaulss <- gam(list(log_volume_t1~log_volume_t + s(dens_scaled) + s(unique.transect,bs="re"),~s(fitted_vals)), 
    #family="gaulss", data=LATR_grow, method="ML",gamma=1.4)  
    fit_gaulss<-update(fit_gaulss)
    fitted_all = predict(fit_gaulss,type="response",data=LATR_grow);   
    new_fitted_vals = fitted_all[,1]; new_weights = fitted_all[,2];
    err = weights - new_weights; err=sqrt(mean(err^2)); 
    weights = new_weights; 
    k=k+1; cat(k,err,"\n"); 
  }   
LATR_GAU[[i]] = fit_gaulss; 
}

AICctab(LATR_GAU,sort=F); 
## the gam wins by 30.7 AIC units 
## TM: yes, looks like the smooth for SD(size) is much better
## I'm curious what the density smooth looks like
plot(LATR_grow$dens_scaled,predict.gam(LATR_GAU[[6]],type="terms",terms="s(dens_scaled)"))


### SPE: let's see what the gam says about standard deviation vs. fitted 
out = predict(LATR_GAU[[6]],type="response"); 
plot(out[,1],1/out[,2]); 

### SPE: back to using the log-linear model for SD. 
LATR_GAU[[6]] = LATR_GAU[[1]]   # kludge, so the gam isn't selected as the best model below 
## TM: not sure what Steve is doing here

## finally, re-fit with REML=T and best weights
LATR_GAU_best<-update(LATR_GAU[[which.min(AICctab(LATR_GAU,sort=F)$dAICc)]],weights=best_weights,REML=T)
best_weights<-weights(LATR_GAU_best)
LATR_grow$GAU_fitted <- fitted(LATR_GAU_best)
LATR_grow$GAU_resids <- residuals(LATR_GAU_best)
LATR_grow$GAU_scaled_resids <- LATR_grow$GAU_resids*sqrt(best_weights) ##sqrt(weights)=1/sd

## should be mean zero unit variance
mean(LATR_grow$GAU_scaled_resids);sd(LATR_grow$GAU_scaled_resids)

## get parameters for sd as f(fitted)
## not using these (SPE) -- so don't compute themm and see what fails. 
GAU_sd_coef<-maxLik(logLik=sdloglik,start=c(exp(sd(LATR_grow$GAU_resids)),0))  

##are the standardized residuals gaussian? -- no
jarque.test(LATR_grow$GAU_scaled_resids) # normality test: FAILS, P < 0.001 
anscombe.test(LATR_grow$GAU_scaled_resids) # kurtosis: FAILS, P < 0.001 
agostino.test(LATR_grow$GAU_scaled_resids) # skewness: FAILS, P<0.001 

## fit quantile regression with quantreg
q.fit<-predict(rq(GAU_scaled_resids~GAU_fitted, data=LATR_grow,tau=c(0.05,0.1,0.25,0.5,0.75,0.9,0.95)))
q.05<-q.fit[,1]
q.10<-q.fit[,2] 
q.25<-q.fit[,3] 
q.50<-q.fit[,4]
q.75<-q.fit[,5]
q.90<-q.fit[,6] 
q.95<-q.fit[,7]

size_means<-LATR_grow %>% group_by(size_bin) %>% summarise(mean=mean(log_volume_t),col=unique(size_col))
dens_means<-LATR_grow %>% group_by(dens_bin) %>% summarise(mean=mean(dens_scaled),col=unique(dens_col))
## make a function for plotting
predict_size<-function(size_t,dens){
  fixef(LATR_GAU_best)[1]+fixef(LATR_GAU_best)[2]*LATR_grow$log_volume_t+
    fixef(LATR_GAU_best)[3]*size_means$mean_size[1]+fixef(LATR_GAU_best)[4]*(size_means$mean_size[1])^2
}

plot(LATR_grow$log_volume_t,LATR_grow$log_volume_t1,
     col=alpha(LATR_grow$dens_col,0.75),pch=16)
points(LATR_grow$log_volume_t,
      fixef(LATR_GAU_best)[1]+fixef(LATR_GAU_best)[2]*LATR_grow$log_volume_t+
        fixef(LATR_GAU_best)[3]*dens_means$mean[1]+fixef(LATR_GAU_best)[4]*(dens_means$mean[1])^2,
      pch=".")
points(LATR_grow$log_volume_t,
       fixef(LATR_GAU_best)[1]+fixef(LATR_GAU_best)[2]*LATR_grow$log_volume_t+
         fixef(LATR_GAU_best)[3]*dens_means$mean[2]+fixef(LATR_GAU_best)[4]*(dens_means$mean[2])^2,
       pch=".")
points(LATR_grow$log_volume_t,
       fixef(LATR_GAU_best)[1]+fixef(LATR_GAU_best)[2]*LATR_grow$log_volume_t+
         fixef(LATR_GAU_best)[3]*dens_means$mean[3]+fixef(LATR_GAU_best)[4]*(dens_means$mean[3])^2,
       pch=".")

dens_dummy<-seq(min(LATR_grow$dens_scaled),
                max(LATR_grow$dens_scaled),0.1)

pdf("./manuscript/figures/creosote_diagnostics.pdf",height = 3, width = 9,useDingbats = F)
par(mar = c(5, 4, 2, 2), oma=c(0,0,0,2), mfrow=c(1,3)) 
plot(LATR_grow$dens_scaled*100,LATR_grow$log_volume_t1,
     col=alpha(LATR_grow$size_col,0.5),pch=16,
     xlab="Weighted density",ylab="Size at t+1")
lines(dens_dummy*100,
      fixef(LATR_GAU_best)[1]+fixef(LATR_GAU_best)[2]*size_means$mean[1]+
        fixef(LATR_GAU_best)[3]*dens_dummy+fixef(LATR_GAU_best)[4]*dens_dummy^2,
      pch=".",col=dens_means$col[1],lwd=2)
lines(dens_dummy*100,
      fixef(LATR_GAU_best)[1]+fixef(LATR_GAU_best)[2]*size_means$mean[2]+
        fixef(LATR_GAU_best)[3]*dens_dummy+fixef(LATR_GAU_best)[4]*dens_dummy^2,
      pch=".",col=dens_means$col[2],lwd=2)
lines(dens_dummy*100,
      fixef(LATR_GAU_best)[1]+fixef(LATR_GAU_best)[2]*size_means$mean[3]+
        fixef(LATR_GAU_best)[3]*dens_dummy+fixef(LATR_GAU_best)[4]*dens_dummy^2,
      pch=".",col=dens_means$col[3],lwd=2)
legend("bottomright",title="Size at t",legend=round(size_means$mean,2),
       bty="n",cex=0.8,pch=16,col=dens_pallete)
title("A",adj=0,font=3)

plot(sort(LATR_grow$GAU_fitted),
     exp(GAU_sd_coef$estimate[1]+GAU_sd_coef$estimate[2]*sort(LATR_grow$GAU_fitted)),
     xlab="Expected size at t+1",ylab="SD of size at t+1",type="l",lwd=3)
title("B",adj=0,font=3)

plot(LATR_grow$GAU_fitted,LATR_grow$GAU_scaled_resids,col=alpha("black",0.25),
     xlab="Expected size at t+1",ylab="Scaled residuals of size at t+1")
points(LATR_grow$GAU_fitted,q.05,col="black",pch=".")
points(LATR_grow$GAU_fitted,q.10,col="black",pch=".")
points(LATR_grow$GAU_fitted,q.25,col="black",pch=".")
points(LATR_grow$GAU_fitted,q.50,col="black",pch=".")
points(LATR_grow$GAU_fitted,q.75,col="black",pch=".")
points(LATR_grow$GAU_fitted,q.90,col="black",pch=".")
points(LATR_grow$GAU_fitted,q.95,col="black",pch=".")
par(new = TRUE)                           
plot(c(LATR_grow$GAU_fitted,LATR_grow$GAU_fitted),
     c(Q.skewness(q.10,q.50,q.90),Q.kurtosis(q.05,q.25,q.75,q.95)),
     col=c(rep(alpha("blue",0.25),nrow(LATR_grow)),rep(alpha("red",0.25),nrow(LATR_grow))),
     pch=16,cex=.5,axes = FALSE, xlab = "", ylab = "")
abline(h=0,col="lightgray",lty=3)
axis(side = 4,cex.axis=0.8,at = pretty(range(c(Q.skewness(q.10,q.50,q.90),Q.kurtosis(q.05,q.25,q.75,q.95)))))
mtext("Skewness", side = 4, line = 2,col="blue",cex=0.7)
mtext("Excess Kurtosis", side = 4, line =3,col="red",cex=0.7)
title("C",adj=0,font=3)
dev.off()

##################################################################################################
### SPE: let's make the diagnostics figure again, using the spline function for SD 
###      but sticking to the lmer() model
##################################################################################################

best_weights = out[,2]^2; 

## finally, re-fit with REML=T and best weights
LATR_GAU_best<-update(LATR_GAU[[which.min(AICctab(LATR_GAU,sort=F)$dAICc)]],weights=best_weights,REML=T)
LATR_grow$GAU_fitted <- fitted(LATR_GAU_best)
LATR_grow$GAU_resids <- residuals(LATR_GAU_best)
LATR_grow$GAU_scaled_resids <- LATR_grow$GAU_resids*sqrt(best_weights) ##sqrt(weights)=1/sd
## should be mean zero unit variance
mean(LATR_grow$GAU_scaled_resids);sd(LATR_grow$GAU_scaled_resids)

## fit quantile regression with quantreg
q.fit<-predict(rq(GAU_scaled_resids~GAU_fitted, data=LATR_grow,tau=c(0.05,0.1,0.25,0.5,0.75,0.9,0.95)))
q.05<-q.fit[,1]
q.10<-q.fit[,2] 
q.25<-q.fit[,3] 
q.50<-q.fit[,4]
q.75<-q.fit[,5]
q.90<-q.fit[,6] 
q.95<-q.fit[,7]

size_means<-LATR_grow %>% group_by(size_bin) %>% summarise(mean=mean(log_volume_t),col=unique(size_col))
dens_means<-LATR_grow %>% group_by(dens_bin) %>% summarise(mean=mean(dens_scaled),col=unique(dens_col))
## make a function for plotting
predict_size<-function(size_t,dens){
  fixef(LATR_GAU_best)[1]+fixef(LATR_GAU_best)[2]*LATR_grow$log_volume_t+
    fixef(LATR_GAU_best)[3]*size_means$mean_size[1]+fixef(LATR_GAU_best)[4]*(size_means$mean_size[1])^2
}

plot(LATR_grow$log_volume_t,LATR_grow$log_volume_t1,
     col=alpha(LATR_grow$dens_col,0.75),pch=16)
points(LATR_grow$log_volume_t,
      fixef(LATR_GAU_best)[1]+fixef(LATR_GAU_best)[2]*LATR_grow$log_volume_t+
        fixef(LATR_GAU_best)[3]*dens_means$mean[1]+fixef(LATR_GAU_best)[4]*(dens_means$mean[1])^2,
      pch=".")
points(LATR_grow$log_volume_t,
       fixef(LATR_GAU_best)[1]+fixef(LATR_GAU_best)[2]*LATR_grow$log_volume_t+
         fixef(LATR_GAU_best)[3]*dens_means$mean[2]+fixef(LATR_GAU_best)[4]*(dens_means$mean[2])^2,
       pch=".")
points(LATR_grow$log_volume_t,
       fixef(LATR_GAU_best)[1]+fixef(LATR_GAU_best)[2]*LATR_grow$log_volume_t+
         fixef(LATR_GAU_best)[3]*dens_means$mean[3]+fixef(LATR_GAU_best)[4]*(dens_means$mean[3])^2,
       pch=".")

dens_dummy<-seq(min(LATR_grow$dens_scaled),
                max(LATR_grow$dens_scaled),0.1)

pdf("./manuscript/figures/creosote_diagnostics_gamSD.pdf",height = 3, width = 9,useDingbats = F)
par(mar = c(5, 4, 2, 2), oma=c(0,0,0,2), mfrow=c(1,3)) 
plot(LATR_grow$dens_scaled*100,LATR_grow$log_volume_t1,
     col=alpha(LATR_grow$size_col,0.5),pch=16,
     xlab="Weighted density",ylab="Size at t+1")
lines(dens_dummy*100,
      fixef(LATR_GAU_best)[1]+fixef(LATR_GAU_best)[2]*size_means$mean[1]+
        fixef(LATR_GAU_best)[3]*dens_dummy+fixef(LATR_GAU_best)[4]*dens_dummy^2,
      pch=".",col=dens_means$col[1],lwd=2)
lines(dens_dummy*100,
      fixef(LATR_GAU_best)[1]+fixef(LATR_GAU_best)[2]*size_means$mean[2]+
        fixef(LATR_GAU_best)[3]*dens_dummy+fixef(LATR_GAU_best)[4]*dens_dummy^2,
      pch=".",col=dens_means$col[2],lwd=2)
lines(dens_dummy*100,
      fixef(LATR_GAU_best)[1]+fixef(LATR_GAU_best)[2]*size_means$mean[3]+
        fixef(LATR_GAU_best)[3]*dens_dummy+fixef(LATR_GAU_best)[4]*dens_dummy^2,
      pch=".",col=dens_means$col[3],lwd=2)
legend("bottomright",title="Size at t",legend=round(size_means$mean,2),
       bty="n",cex=0.8,pch=16,col=dens_pallete)
title("A",adj=0,font=3)

e = order(out[,1]); 
plot(out[e,1], 1/out[e,2],
     xlab="Expected size at t+1",ylab="SD of size at t+1",type="l",lwd=3)
title("B",adj=0,font=3)

plot(LATR_grow$GAU_fitted,LATR_grow$GAU_scaled_resids,col=alpha("black",0.25),
     xlab="Expected size at t+1",ylab="Scaled residuals of size at t+1")
points(LATR_grow$GAU_fitted,q.05,col="black",pch=".")
points(LATR_grow$GAU_fitted,q.10,col="black",pch=".")
points(LATR_grow$GAU_fitted,q.25,col="black",pch=".")
points(LATR_grow$GAU_fitted,q.50,col="black",pch=".")
points(LATR_grow$GAU_fitted,q.75,col="black",pch=".")
points(LATR_grow$GAU_fitted,q.90,col="black",pch=".")
points(LATR_grow$GAU_fitted,q.95,col="black",pch=".")
par(new = TRUE)                           
plot(c(LATR_grow$GAU_fitted,LATR_grow$GAU_fitted),
     c(Q.skewness(q.10,q.50,q.90),Q.kurtosis(q.05,q.25,q.75,q.95)),
     col=c(rep(alpha("blue",0.25),nrow(LATR_grow)),rep(alpha("red",0.25),nrow(LATR_grow))),
     pch=16,cex=.5,axes = FALSE, xlab = "", ylab = "")
abline(h=0,col="lightgray",lty=3)
axis(side = 4,cex.axis=0.8,at = pretty(range(c(Q.skewness(q.10,q.50,q.90),Q.kurtosis(q.05,q.25,q.75,q.95)))))
mtext("Skewness", side = 4, line = 2,col="blue",cex=0.7)
mtext("Excess Kurtosis", side = 4, line =3,col="red",cex=0.7)
title("C",adj=0,font=3)
dev.off()
########################## END of re-making the diagnostics graph 

##################################################################################################
### SPE: what about **studentized** residuals rather than scaled residuals? 
### The punchline is: ALMOST NO DIFFERENCE, with this sample size no points 
###      have especially high leverage.  
###################################################################################################
best_weights = out[,2]^2; 

## finally, re-fit with REML=T and best weights
LATR_GAU_best<-update(LATR_GAU[[which.min(AICctab(LATR_GAU,sort=F)$dAICc)]],weights=best_weights,REML=T)

LATR_grow$GAU_resids <- residuals(LATR_GAU_best)

LATR_grow$GAU_scaled_resids1 <- LATR_grow$GAU_resids*sqrt(best_weights) ##sqrt(weights)=1/sd
LATR_grow$GAU_scaled_resids <- rstudent(LATR_GAU_best);

plot(LATR_grow$GAU_scaled_resids, LATR_grow$GAU_scaled_resids1); abline(0,1); 
########################## END of studentized residuals 


##################################################################################################
### SPE: let's make the last panel again, using the spline model throughout.  
### Prereq: run lines 1 - 130. Loads data and creates the spline model, fit_gaulss. 
### This has SD as a function of fitted, through iterative refitting. 
### Be careful to get the right type of residuals! See ?residuals.gam  
###################################################################################################
out = predict(fit_gaulss,type="response"); 
e = order(out[,1]); plot(out[e,1], 1/out[e,2],type="p"); 
best_weights = out[,2]^2; 

LATR_grow$GAU_fitted <- out[,1]; 
LATR_grow$GAU_resids <- residuals(fit_gaulss,type="response")
LATR_grow$GAU_scaled_resids <- LATR_grow$GAU_resids*sqrt(best_weights) 
## should be mean zero, unit variance
mean(LATR_grow$GAU_scaled_resids);sd(LATR_grow$GAU_scaled_resids)

######### gam() will do the scaling for us! 
pearson_resids = residuals(fit_gaulss,type="pearson");
range(LATR_grow$GAU_scaled_resids - pearson_resids); 

require(qgam); 
S.05<-qgam(GAU_scaled_resids~s(GAU_fitted,k=6), data=LATR_grow,qu=0.05)  
S.10<-qgam(GAU_scaled_resids~s(GAU_fitted,k=6), data=LATR_grow,qu=0.1)   
S.25<-qgam(GAU_scaled_resids~s(GAU_fitted,k=6), data=LATR_grow,qu=0.25) 
S.50<-qgam(GAU_scaled_resids~s(GAU_fitted,k=6), data=LATR_grow,qu=0.5)   
S.75<-qgam(GAU_scaled_resids~s(GAU_fitted,k=6), data=LATR_grow,qu=0.75) 
S.90<-qgam(GAU_scaled_resids~s(GAU_fitted,k=6), data=LATR_grow,qu=0.9)   
S.95<-qgam(GAU_scaled_resids~s(GAU_fitted,k=6), data=LATR_grow,qu=0.95)  

## NP skewness
q.10<-predict(S.10);q.50<-predict(S.50);q.90<-predict(S.90)
NPS_hat = (q.10 + q.90 - 2*q.50)/(q.90 - q.10)

## NP kurtosis (relative to Gaussian)
q.05<-predict(S.05);q.25<-predict(S.25);q.75<-predict(S.75);q.95<-predict(S.95)
qN = qnorm(c(0.05,0.25,0.75,0.95))
KG = (qN[4]-qN[1])/(qN[3]-qN[2])
NPK_hat = ((q.95-q.05)/(q.75-q.25))/KG - 1

pdf("./manuscript/figures/creosote_diagnostics_gaulss_qgam.pdf",height = 6, width = 8,useDingbats = F)
par(mar = c(5, 4, 2, 3), oma=c(0,0,0,4)) 
plot(LATR_grow$GAU_fitted,LATR_grow$GAU_scaled_resids,col=alpha("black",0.25),
     xlab="Expected size at t+1",ylab="Scaled residuals of size at t+1")
points(LATR_grow$GAU_fitted,q.05,col="black",pch=".")
points(LATR_grow$GAU_fitted,q.10,col="black",pch=".")
points(LATR_grow$GAU_fitted,q.25,col="black",pch=".")
points(LATR_grow$GAU_fitted,q.50,col="black",pch=".")
points(LATR_grow$GAU_fitted,q.75,col="black",pch=".")
points(LATR_grow$GAU_fitted,q.90,col="black",pch=".")
points(LATR_grow$GAU_fitted,q.95,col="black",pch=".")
par(new = TRUE)                           
plot(c(LATR_grow$GAU_fitted,LATR_grow$GAU_fitted),
     c(Q.skewness(q.10,q.50,q.90),Q.kurtosis(q.05,q.25,q.75,q.95)),
     col=c(rep(alpha("blue",0.25),nrow(LATR_grow)),rep(alpha("red",0.25),nrow(LATR_grow))),
     pch=16,cex=.5,axes = FALSE, xlab = "", ylab = "")
abline(h=0,col="lightgray",lty=3)
axis(side = 4,cex.axis=0.8,at = pretty(range(c(Q.skewness(q.10,q.50,q.90),Q.kurtosis(q.05,q.25,q.75,q.95)))))
mtext("Skewness", side = 4, line = 2,col="blue")
mtext("Excess Kurtosis", side = 4, line =3,col="red")
dev.off()


######### END of re-making the diagnostics graph with gam and qgam. 



















## based on these results I will fit a JSU distribution to the residuals
## will need to fit variance, skew, and kurtosis as functions of the mean
## makes random effects tricker -- will fit as fixed effects and use Steve's "shrinkage" methods

## try a max-lik fit using the fitted values from lmer as mu
## and other params as functions of mu
JSULogLik=function(pars){
  dJSU(LATR_grow$log_volume_t1, 
          mu=LATR_grow$GAU_fitted,
          sigma=exp(GAU_sd_coef$estimate[1]+GAU_sd_coef$estimate[2]*LATR_grow$GAU_fitted),
          nu = pars[3]+pars[4]*LATR_grow$GAU_fitted,
          tau = exp(pars[5]+pars[6]*LATR_grow$GAU_fitted), log=TRUE)
}
## starting parameters
p0<-c(0,0,0,0,0,0)

## pass this through several ML algorithms
JSUout=maxLik(logLik=JSULogLik,start=p0*exp(0.2*rnorm(length(p0))),method="BHHH",control=list(iterlim=5000,printLevel=2),finalHessian=FALSE) 
JSUout=maxLik(logLik=JSULogLik,start=JSUout$estimate,method="NM",control=list(iterlim=5000,printLevel=1),finalHessian=FALSE)
JSUout=maxLik(logLik=JSULogLik,start=JSUout$estimate,method="BHHH",control=list(iterlim=5000,printLevel=2),finalHessian=FALSE)

## simulate from fitted model
n_sim<-100
JSUsim_mean<-JSUsim_sd<-JSUsim_skew<-JSUsim_kurt<-matrix(NA,nrow=nrow(LATR_grow),ncol=n_sim)
GAUsim_mean<-GAUsim_sd<-GAUsim_skew<-GAUsim_kurt<-matrix(NA,nrow=nrow(LATR_grow),ncol=n_sim)
for(i in 1:n_sim){
  cat(" --------------- Simulation number", i, "\n"); 
  LATR_grow$log_volume_t1.simJSU <- rJSU(n=nrow(LATR_grow),
                                      mu=LATR_grow$GAU_fitted,
                                      sigma=exp(GAU_sd_coef$estimate[1]+GAU_sd_coef$estimate[2]*LATR_grow$GAU_fitted),
                                      nu=JSUout$estimate[3]+JSUout$estimate[4]*LATR_grow$GAU_fitted,
                                      tau=exp(JSUout$estimate[5]+JSUout$estimate[6]*LATR_grow$GAU_fitted))
  q.fit<-predict(rq(log_volume_t1.simJSU~log_volume_t, data=LATR_grow,tau=c(0.05,0.1,0.25,0.5,0.75,0.9,0.95)))
  JSUsim_mean[,i]<-Q.mean(q.fit[,3],q.fit[,4],q.fit[,5])
  JSUsim_sd[,i]<-Q.sd(q.fit[,3],q.fit[,5])
  JSUsim_skew[,i]<-Q.skewness(q.fit[,2],q.fit[,4],q.fit[,6])
  JSUsim_kurt[,i]<-Q.kurtosis(q.fit[,1],q.fit[,3],q.fit[,5],q.fit[,7])
  
  LATR_grow$log_volume_t1.simGAU <- rnorm(n=nrow(LATR_grow),
                                         mean=LATR_grow$GAU_fitted,
                                         sd=exp(GAU_sd_coef$estimate[1]+GAU_sd_coef$estimate[2]*LATR_grow$GAU_fitted))
  q.fit<-predict(rq(log_volume_t1.simGAU~log_volume_t, data=LATR_grow,tau=c(0.05,0.1,0.25,0.5,0.75,0.9,0.95)))
  GAUsim_mean[,i]<-Q.mean(q.fit[,3],q.fit[,4],q.fit[,5])
  GAUsim_sd[,i]<-Q.sd(q.fit[,3],q.fit[,5])
  GAUsim_skew[,i]<-Q.skewness(q.fit[,2],q.fit[,4],q.fit[,6])
  GAUsim_kurt[,i]<-Q.kurtosis(q.fit[,1],q.fit[,3],q.fit[,5],q.fit[,7])
}

## and now the real data
q.fit<-predict(rq(log_volume_t1~log_volume_t, data=LATR_grow,tau=c(0.05,0.1,0.25,0.5,0.75,0.9,0.95)))
q.05<-q.fit[,1]
q.10<-q.fit[,2] 
q.25<-q.fit[,3] 
q.50<-q.fit[,4]
q.75<-q.fit[,5]
q.90<-q.fit[,6] 
q.95<-q.fit[,7]

pdf("./manuscript/figures/creosote_JSU_fit.pdf",height = 6, width = 6,useDingbats = F)
par(mfrow=c(2,2),mar=c(4,4,1,1),mgp=c(2.1,1,0),cex.lab=1.2); 
plot(LATR_grow$log_volume_t,Q.mean(q.25,q.50,q.75),type="n",
     xlab="Size(t)",ylab="NP mean size(t+1)",
     ylim=c(min(c(GAUsim_mean,JSUsim_mean)),1 + max(c(GAUsim_mean,JSUsim_mean))))
     matpoints(LATR_grow$log_volume_t,GAUsim_mean,col=alpha("tomato",0.25),pch=".",cex=2)
     matpoints(LATR_grow$log_volume_t,1 + JSUsim_mean,col=alpha("cornflowerblue",0.25),pch=".",cex=2)
points(LATR_grow$log_volume_t,Q.mean(q.25,q.50,q.75),col="black",pch=".",cex=2)
points(LATR_grow$log_volume_t,1 + Q.mean(q.25,q.50,q.75),col="black",pch=".",cex=2)
legend("topleft",legend=c("Real data","Gaussian simulation","JSU simulation + offset"),
       lty=1,col=c("black","tomato","cornflowerblue"),cex=0.8,bty="n")

plot(LATR_grow$log_volume_t,Q.sd(q.25,q.75),type="n",
     xlab="Size(t)",ylab="NP SD Size(t+1)",
     ylim=c(min(c(GAUsim_sd,JSUsim_sd)),1 + max(c(GAUsim_sd,JSUsim_sd))))
     matpoints(LATR_grow$log_volume_t,GAUsim_sd,col=alpha("tomato",0.25),pch=".",cex=2)
     matpoints(LATR_grow$log_volume_t,1 + JSUsim_sd,col=alpha("cornflowerblue",0.25),pch=".",cex=2)
points(LATR_grow$log_volume_t,Q.sd(q.25,q.75),col="black",pch=".",cex=2)
points(LATR_grow$log_volume_t,1 + Q.sd(q.25,q.75),col="black",pch=".",cex=2)

plot(LATR_grow$log_volume_t,Q.skewness(q.10,q.50,q.90),type="n",
     xlab="Size(t)",ylab="NP skewness size(t+1)",
     ylim=c(min(c(GAUsim_skew,JSUsim_skew)),1 + max(c(GAUsim_skew,JSUsim_skew))))
     matpoints(LATR_grow$log_volume_t,GAUsim_skew,col=alpha("tomato",0.25),pch=".",cex=2)
     matpoints(LATR_grow$log_volume_t,1+ JSUsim_skew,col=alpha("cornflowerblue",0.25),pch=".",cex=2)
    points(LATR_grow$log_volume_t,Q.skewness(q.10,q.50,q.90),col="black",pch=".",cex=2)
    points(LATR_grow$log_volume_t,1 + Q.skewness(q.10,q.50,q.90),col="black",pch=".",cex=2)

plot(LATR_grow$log_volume_t,Q.kurtosis(q.05,q.25,q.75,q.95),type="n",
     xlab="Size(t)",ylab="NP kurtosis size(t+1)",
     ylim=c(min(GAUsim_kurt,JSUsim_kurt),1 + max(GAUsim_kurt,JSUsim_kurt)))
     matpoints(LATR_grow$log_volume_t,GAUsim_kurt,col=alpha("tomato",0.25),pch=".",cex=2)
     matpoints(LATR_grow$log_volume_t,1 + JSUsim_kurt,col=alpha("cornflowerblue",0.25),pch=".",cex=2)
points(LATR_grow$log_volume_t,Q.kurtosis(q.05,q.25,q.75,q.95),col="black",pch=".",cex=2)
points(LATR_grow$log_volume_t,1 + Q.kurtosis(q.05,q.25,q.75,q.95),col="black",pch=".",cex=2)

dev.off()


# compare IPM results between Gaussian and JSU growth kernel --------------
## fit other vital rates -- this is all done with mgcv, because 
## I already had these models at my fingertips
## set the gamma argument of gam() -- gamma>1 generates smoother fits, less likely to be overfit
gamma = 1.8
##### Flowering probability model -------------------------------------------------------------------------

# Populate year t of 2017-2018 transition year
# There are no 2018 data but this way we get all four years in the reproduction models
# Do this by creating the 2017-18 data as a stand-alone df then bind rows
LATR_dat_201718 <- LATR_full[LATR_full$year_t == 2016 & LATR_full$survival_t1 == 1, ]
# These are the 2017 survivors; make their year t demography last year's data
LATR_dat_201718$year_t <- 2017
LATR_dat_201718$year_t1 <- 2018
LATR_dat_201718$log_volume_t <- LATR_dat_201718$log_volume_t1
LATR_dat_201718$flowers_t <- LATR_dat_201718$flowers_t1
LATR_dat_201718$fruits_t <- LATR_dat_201718$fruits_t1
LATR_dat_201718$reproductive_fraction_t <- LATR_dat_201718$reproductive_fraction_t1
LATR_dat_201718$total.reproduction_t <- LATR_dat_201718$total.reproduction_t1
# Now set all the t1 data to NA
LATR_dat_201718$log_volume_t1 <- NA
LATR_dat_201718$flowers_t1 <- NA
LATR_dat_201718$fruits_t1 <- NA
LATR_dat_201718$reproductive_fraction_t1 <- NA
LATR_dat_201718$total.reproduction_t1 <- NA

# Bind rows and create log_vol as new variables (easier for GAMs)
LATR_flow_dat <- bind_rows(LATR_full,LATR_dat_201718) %>% 
  dplyr::select(unique.transect,log_volume_t,total.reproduction_t,dens_scaled) %>% drop_na()

# Create empty list to populate with model results
LATR_flower <- list()
# Three candidate models for the mean: size only, size + density, or size, density, and size:density
LATR_flower[[1]] <- gam(total.reproduction_t > 0 ~ s(log_volume_t) + s(unique.transect, bs = "re"),
                        data = LATR_flow_dat, gamma = gamma, family = "binomial")
LATR_flower[[2]] <- gam(total.reproduction_t > 0 ~ s(log_volume_t) + s(dens_scaled) + s(unique.transect, bs = "re"),
                        data = LATR_flow_dat, gamma = gamma, family = "binomial")
LATR_flower[[3]] <- gam(total.reproduction_t > 0 ~ s(log_volume_t) + s(dens_scaled) + ti(log_volume_t,dens_scaled) + s(unique.transect, bs = "re"),
                        data = LATR_flow_dat, gamma = gamma, family = "binomial")
# Collect model AICs into a single table
flower_aic<-AICtab(LATR_flower, base = TRUE, sort = FALSE)
# Set top model as "best"
LATR_flower_best <- LATR_flower[[which.min(flower_aic$AIC)]]

##### Fruit production model ------------------------------------------------------------------------------

# Create new df with plants that have produced at least one reproductive structure
LATR_fruits_dat <- subset(LATR_flow_dat, total.reproduction_t > 0)

# Create empty list to populate with model results
LATR_fruits <- list()
# Three candidate models for the mean: size only, size + density, or size, density, and size:density
LATR_fruits[[1]] <- gam(total.reproduction_t ~ s(log_volume_t) + s(unique.transect, bs = "re"),
                        data = LATR_fruits_dat, gamma = gamma, family = "nb")
LATR_fruits[[2]] <- gam(total.reproduction_t ~ s(log_volume_t) + s(dens_scaled) + s(unique.transect, bs = "re"),
                        data = LATR_fruits_dat, gamma = gamma, family = "nb")
LATR_fruits[[3]] <- gam(total.reproduction_t ~ s(log_volume_t) + s(dens_scaled) + ti(log_volume_t,dens_scaled) + s(unique.transect, bs = "re"),
                        data = LATR_fruits_dat, gamma = gamma, family = "nb")
# Collect model AICs into a single table
fruits_aic <- AICtab(LATR_fruits, base = TRUE, sort = FALSE)
# Set top model as "best"
LATR_fruits_best <- LATR_fruits[[which.min(fruits_aic$AIC)]]

##### Survival model --------------------------------------------------------------------------------------

# Combine transplants with large shrubs; keep only location info, survival, volume, and density
CData.Transplants %>% 
  dplyr::select("site", "transect", "actual.window", 
                "spring_survival_t1", "volume_t", "weighted.dens", "transplant") %>% 
  rename("survival_t1" = "spring_survival_t1") %>% 
  mutate(unique.transect = interaction(transect, site)) %>% 
  rbind(dplyr::select(LATR_full, "site", "transect", "actual.window", 
                      "survival_t1", "volume_t", "weighted.dens", "transplant","unique.transect")) %>% 
  mutate(log_volume_t = log(volume_t),
         dens_scaled = weighted.dens/100) %>% 
  drop_na() -> LATR_surv_dat

# Create empty list to populate with model results
LATR_surv <- list()
# Three candidate models for the mean: size only, size + density, or size, density, and size:density
LATR_surv[[1]] <- gam(survival_t1 ~ s(log_volume_t,by=as.factor(transplant)) + s(unique.transect, bs = "re"),
                      data = LATR_surv_dat, gamma = gamma, family = "binomial")
LATR_surv[[2]] <- gam(survival_t1 ~ s(log_volume_t,by=as.factor(transplant)) + s(dens_scaled,by=as.factor(transplant))  + s(unique.transect, bs = "re"),
                      data = LATR_surv_dat, gamma = gamma, family = "binomial")
LATR_surv[[3]] <- gam(survival_t1 ~ s(log_volume_t,by=as.factor(transplant)) + s(dens_scaled,by=as.factor(transplant)) + ti(log_volume_t,dens_scaled,by=as.factor(transplant)) + s(unique.transect, bs = "re"),
                      data = LATR_surv_dat, gamma = gamma, family = "binomial")
# Collect model AICs into a single table
surv_aic <- AICtab(LATR_surv, base = TRUE, sort = FALSE)
# Set top model as "best"
LATR_surv_best <- LATR_surv[[which.min(surv_aic$AIC)]]

##### Per-seed recruitment probability model --------------------------------------------------------------

## number of seeds per fruit
seeds_per_fruit <- 5
# Create subset df of recruits
LATR_recruits <- LATR_full %>% 
  mutate(unique.transect = interaction(transect, site)) %>% 
  group_by(year_t1, unique.transect, actual.window) %>% 
  filter(seedling_t1 == 1)
suppressMessages(summarise(LATR_recruits, recruits = n())) %>% 
  rename(window = actual.window) -> LATR_recruits

# Estimate total seeds produced in each window
# This is computed using the known plant sizes and the fitted flowering and fruiting models
LATR_transects <- Cdata.Transects.Windows %>% 
  mutate(unique.transect = interaction(transect, site),
         log_volume_t = log(volume),
         dens_scaled = weighted.dens/100)
LATR_transects$seeds = ceiling(invlogit(predict.gam(LATR_flower_best, newdata = LATR_transects))* 
                                 seeds_per_fruit*exp(predict.gam(LATR_fruits_best, newdata = LATR_transects)))
LATR_transects <- group_by(LATR_transects, unique.transect, window)
suppressMessages(summarise(LATR_transects, total_seeds = sum(seeds),
                           dens_scaled = unique(dens_scaled))) -> LATR_transects

# Take three copies of this df, assigning each one to a different year and assigning recruits to zero (for now)
LATR_recruitment <- bind_rows(LATR_transects %>% filter(unique.transect == "1.FPS" | unique.transect == "2.FPS" | unique.transect == "3.FPS") %>% 
                                mutate(year_t1 = 2014, recruits = 0), ## only FPS for 2013-2014
                              LATR_transects %>% mutate(year_t1 = 2015, recruits = 0),
                              LATR_transects %>% mutate(year_t1 = 2016, recruits = 0),
                              LATR_transects %>% mutate(year_t1 = 2017, recruits = 0)) %>% 
  left_join(., LATR_recruits, by = c("year_t1", "unique.transect", "window")) %>% 
  mutate(recruits.y = replace_na(recruits.y, 0),
         recruits = pmax(recruits.x, recruits.y, na.rm = T)) %>% 
  drop_na()

# Create empty list to populate with model results
LATR_recruit <- list()
# Two candidate models for the mean: no effect, or weighted density only
LATR_recruit[[1]] <- gam(cbind(recruits,total_seeds - recruits) ~ s(unique.transect, bs = "re"),
                         data = LATR_recruitment, gamma = gamma, family = "binomial")
LATR_recruit[[2]] <- gam(cbind(recruits,total_seeds - recruits) ~ s(dens_scaled) + s(unique.transect, bs = "re"),
                         data = LATR_recruitment, gamma = gamma, family = "binomial")

# Collect model AICs into a single table
recruit_aic <- AICtab(LATR_recruit, base = TRUE, sort = FALSE)
# Set top model as "best"
LATR_recruit_best <- LATR_recruit[[which.min(recruit_aic$AIC)]]
## annoying but necessary index wrangling
recruit_size_sd_index <- which(as.factor(names(coef(LATR_recruitsize_best)))=="(Intercept).1") ## this is where the sd coefficients start
recruit_size_coef_length <- length(coef(LATR_recruitsize_best))

##### Recruit sizes and integration limits (size bounds) --------------------------------------------------
  # Filter out seedlings and get their sizes
  LATR_recruit_size <- LATR_full %>% 
    filter(seedling_t1 == 1) %>% 
    mutate(log_volume = log(volume_t1)) %>% 
    arrange(dens_scaled)

# test for density dependence in recruit size
LATR_recruit_size_mod <- list()
LATR_recruit_size_mod[[1]] <- gam(list(log_volume ~ 1 + s(unique.transect,bs = "re"),~1), 
                                  data = LATR_recruit_size, gamma = gamma, family = gaulss())
LATR_recruit_size_mod[[2]] <- gam(list(log_volume ~ s(dens_scaled) + s(unique.transect,bs = "re"),~1), 
                                  data = LATR_recruit_size, gamma = gamma, family = gaulss())
LATR_recruit_size_mod[[3]] <- gam(list(log_volume ~ 1 + s(unique.transect,bs = "re"), ~s(dens_scaled)), 
                                  data = LATR_recruit_size, gamma = gamma, family = gaulss())
LATR_recruit_size_mod[[4]] <- gam(list(log_volume ~ s(dens_scaled) + s(unique.transect,bs = "re"), ~s(dens_scaled)), 
                                  data = LATR_recruit_size, gamma = gamma, family = gaulss())
# Collect model AICs into a single table
recruitsize_aic <- AICtab(LATR_recruit_size_mod, base = TRUE, sort = FALSE)
# Set top model as "best"
LATR_recruitsize_best <- LATR_recruit_size_mod[[which.min(recruitsize_aic$AIC)]]

# Create maximum and minimum size bounds for the IPM
LATR_size_bounds <- data.frame(min_size = log(min(LATR_full$volume_t, LATR_full$volume_t1[LATR_full$transplant == FALSE], na.rm = TRUE)),
                               max_size = log(max(LATR_full$volume_t, LATR_full$volume_t1[LATR_full$transplant == FALSE], na.rm = TRUE)))

##### SIPM functions ---------------------------------------------------------------------------------------
# Eviction extensions for upper and lower size limits
lower.extension <- -8
upper.extension <- 2

## load source functions and some data elements for wind dispersal
## will take a hot minute and will throw warnings
source("creosote/creosote_SIPM_source_fns.R")

## explore effect of matrix dimensions
dims<-c(100,150,200,250,300,350,400)
lambda0.GAU<-lambda0.JSU<-c()
meanlife.GAU<-meanlife.JSU<-c()
cstar_GAU<-cstar_JSU<-c()
for(i in 1:length(dims)){
  hold.GAU<-ApproxMatrix(dens=0,dist="GAU",mat.size=dims[i])
  hold.JSU<-ApproxMatrix(dens=0,dist="JSU",mat.size=dims[i])
  lambda0.GAU[i]<-lambda(hold.GAU$IPMmat)
  meanlife.GAU[i]<-life_expect_mean(hold.GAU$Pmat,start=1)
  lambda0.JSU[i]<-lambda(hold.JSU$IPMmat)
  meanlife.JSU[i]<-life_expect_mean(hold.JSU$Pmat,start=1)
  
  params_GAU <- WALD_par(mat=hold.GAU)
  params_JSU <- WALD_par(mat=hold.JSU)
  #sample dispersal events for empirical MGF - generates a N*mat.size matrix
  D.samples.GAU <- WALD_samples(N=10000,seed=ran.seeds,params=params_GAU) 
  D.samples.JSU <- WALD_samples(N=10000,seed=ran.seeds,params=params_JSU) 
  cstar_GAU[i] <- optimize(cs,lower=0.05,upper=4,emp=F,mat=hold.GAU,
                           params=params_GAU,D.samples=D.samples.GAU)$objective
  cstar_JSU[i] <- optimize(cs,lower=0.05,upper=4,emp=F,mat=hold.JSU,
                           params=params_JSU,D.samples=D.samples.JSU)$objective
  
  print(i)
}
plot(dims,lambda0.GAU,col="tomato",ylim=range(c(lambda0.GAU,lambda0.JSU)))
points(dims,lambda0.JSU,col="cornflowerblue")
plot(dims,meanlife.GAU,col="tomato",ylim=range(c(meanlife.GAU,meanlife.JSU)))
points(dims,meanlife.JSU,col="cornflowerblue")
plot(dims,cstar_GAU,col="tomato",ylim=range(c(cstar_GAU,cstar_JSU)))
points(dims,cstar_JSU,col="cornflowerblue")
## looks like JSU needs higher dimensional matrix
matdim <- 400

## look at lambda as a function of density
dens=seq(min(LATR_full$dens_scaled),max(LATR_full$dens_scaled),0.1)
lambda.GAU<-lambda.JSU<-c()
for(i in 1:length(dens)){
  lambda.GAU[i]<-lambda(ApproxMatrix(dens=dens[i],dist="GAU")$IPMmat)
  lambda.JSU[i]<-lambda(ApproxMatrix(dens=dens[i],dist="JSU")$IPMmat)
  print(i)
}

# Construct transition matrix for minimum weighted density (zero)
mat_GAU <- ApproxMatrix(dens = 0, dist="GAU")
mat_JSU <- ApproxMatrix(dens = 0, dist="JSU")
colSums(mat_JSU$Pmat)

params_GAU <- WALD_par(mat=mat_GAU)
params_JSU <- WALD_par(mat=mat_JSU)
#sample dispersal events for empirical MGF - generates a N*mat.size matrix
D.samples.GAU <- WALD_samples(N=10000,seed=ran.seeds,params=params_GAU) 
D.samples.JSU <- WALD_samples(N=10000,seed=ran.seeds,params=params_JSU) 

# Find the asymptotic wave speed c*(s) 
cstar_GAU <- optimize(cs,lower=0.05,upper=4,emp=F,mat=mat_GAU,
                      params=params_GAU,D.samples=D.samples.GAU)$objective
cstar_JSU <- optimize(cs,lower=0.05,upper=4,emp=F,mat=mat_JSU,
                      params=params_JSU,D.samples=D.samples.JSU)$objective


pdf("./manuscript/figures/creosote_DD_lambda.pdf",height = 4, width = 4,useDingbats = F)
par(mar=c(4,4,1,1))
plot(dens*100,lambda.JSU,type="b",col=alpha("cornflowerblue",0.75),pch=16,cex=1.2,
     xlab="Weighted density",ylab=expression(paste(lambda)))
lines(dens*100,lambda.GAU,type="b",col=alpha("tomato",0.75),pch=16,cex=1.2)
legend("topright",bty="n",legend=c("Gaussian","JSU"),
       pch=16,col=c(alpha("tomato",1),alpha("cornflowerblue",1)))
text(140,1.025,paste("c*=",round(cstar_GAU,5),"m/yr"),col=alpha("tomato",1))
text(140,1.02,paste("c*=",round(cstar_JSU,5),"m/yr"),col=alpha("cornflowerblue",1))
dev.off()


life_expect_mean(mat_GAU$Pmat,start=1);life_expect_var(mat_GAU$Pmat,start=1)
life_expect_mean(mat_JSU$Pmat,start=1);life_expect_var(mat_JSU$Pmat,start=1)

## remaining life expectancy of median-sized shrub
median_index<-which.min(abs(mat_GAU$meshpts-median(LATR_full$log_volume_t,na.rm=T)))
life_expect_mean(mat_JSU$Pmat,start=median_index)

plot(flower(mat_JSU$meshpts,d=0))
plot(survival(mat_JSU$meshpts,d=0))

colSums(mat_GAU$Pmat)
colSums(mat_JSU$Pmat)
