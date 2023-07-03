library(lme4)
library(tidyverse)
library(quantreg)
library(gamlss.dist)
library(maxLik)
library(bbmle)
library(popbio)
library(moments)

## Tom's local directory
setwd("C:/Users/tm9/Dropbox/github/IPM_size_transitions")

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
  mutate(unique.transect = interaction(transect, site))

## color pallete for growth figure -- thanks, colorbrewer
dens_pallete<-c("#bdd7e7","#6baed6","#2171b5")

# Prepare a data subset for growth that drops rows missing either t or t1 size data
## update to the growth data -- dropping a few unbelievable outliers following additional QA/QC
outliers<-c("MOD.2.50.3.2016","MOD.3.200.1.2015","MOD.3.200.1.2014","PDC.2.0.5.2014")
LATR_grow <- LATR_full %>% 
  ## this ID will help us drop outliers below
  mutate(ID=interaction(site,transect,actual.window,plant,year_t)) %>% 
  drop_na(volume_t, volume_t1) %>%
  mutate(log_volume_t = log(volume_t),
         log_volume_t1 = log(volume_t1),
         dens_scaled = weighted.dens/100) %>% 
  ##need to scale weighted density because 1st and 2nd order variables were hugely different in range
  filter(!ID%in%outliers) %>% 
  ##bin density variation to make a nice plot
  mutate(dens_bin=cut_interval(dens_scaled,n=length(dens_pallete),labels=F),
         dens_col=dens_pallete[dens_bin],
         size_bin=cut_interval(log_volume_t,n=length(dens_pallete),labels=F),
         size_col=dens_pallete[size_bin])

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
AICctab(LATR_GAU,sort=F)

## finally, re-fit with REML=T and best weights
LATR_GAU_best<-update(LATR_GAU[[which.min(AICctab(LATR_GAU,sort=F)$dAICc)]],weights=best_weights,REML=T)
best_weights<-weights(LATR_GAU_best)
LATR_grow$GAU_fitted <- fitted(LATR_GAU_best)
LATR_grow$GAU_resids <- residuals(LATR_GAU_best)
LATR_grow$GAU_scaled_resids <- LATR_grow$GAU_resids*sqrt(best_weights) ##sqrt(weights)=1/sd
## should be mean zero unit variance
mean(LATR_grow$GAU_scaled_resids);sd(LATR_grow$GAU_scaled_resids)
## get parameters for sd as f(fitted)
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
par(mfrow=c(2,2),mar=c(4,4,1,1))
plot(LATR_grow$log_volume_t,Q.mean(q.25,q.50,q.75),type="n",
     xlab="size t",ylab="mean size t1",
     ylim=c(min(c(GAUsim_mean,JSUsim_mean)),max(c(GAUsim_mean,JSUsim_mean))))
for(i in 1:n_sim){
  points(LATR_grow$log_volume_t,GAUsim_mean[,i],col=alpha("tomato",0.25),pch=".")
  points(LATR_grow$log_volume_t,JSUsim_mean[,i],col=alpha("cornflowerblue",0.25),pch=".")
}
points(LATR_grow$log_volume_t,Q.mean(q.25,q.50,q.75),col="black",pch=".",cex=2)
legend("topleft",legend=c("Real data","Gaussian simulation","JSU simulation"),
       lty=1,col=c("black","tomato","cornflowerblue"),cex=0.8,bty="n")

plot(LATR_grow$log_volume_t,Q.sd(q.25,q.75),type="n",
     xlab="size t",ylab="sd size t1",
     ylim=c(min(c(GAUsim_sd,JSUsim_sd)),max(c(GAUsim_sd,JSUsim_sd))))
for(i in 1:n_sim){
  points(LATR_grow$log_volume_t,GAUsim_sd[,i],col=alpha("tomato",0.25),pch=".")
  points(LATR_grow$log_volume_t,JSUsim_sd[,i],col=alpha("cornflowerblue",0.25),pch=".")
}
points(LATR_grow$log_volume_t,Q.sd(q.25,q.75),col="black",pch=".",cex=2)

plot(LATR_grow$log_volume_t,Q.skewness(q.10,q.50,q.90),type="n",
     xlab="size t",ylab="skewness size t1",
     ylim=c(min(c(GAUsim_skew,JSUsim_skew)),max(c(GAUsim_skew,JSUsim_skew))))
for(i in 1:n_sim){
  points(LATR_grow$log_volume_t,GAUsim_skew[,i],col=alpha("tomato",0.25),pch=".")
  points(LATR_grow$log_volume_t,JSUsim_skew[,i],col=alpha("cornflowerblue",0.25),pch=".")
}
points(LATR_grow$log_volume_t,Q.skewness(q.10,q.50,q.90),col="black",pch=".",cex=2)

plot(LATR_grow$log_volume_t,Q.kurtosis(q.05,q.25,q.75,q.95),type="n",
     xlab="size t",ylab="kurtosis size t1",
     ylim=c(min(GAUsim_kurt,JSUsim_kurt),max(GAUsim_kurt,JSUsim_kurt)))
for(i in 1:n_sim){
  points(LATR_grow$log_volume_t,GAUsim_kurt[,i],col=alpha("tomato",0.25),pch=".")
  points(LATR_grow$log_volume_t,JSUsim_kurt[,i],col=alpha("cornflowerblue",0.25),pch=".")
}
points(LATR_grow$log_volume_t,Q.kurtosis(q.05,q.25,q.75,q.95),col="black",pch=".",cex=2)
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
LATR_dat_201718$max.ht_t <- LATR_dat_201718$max.ht_t1
LATR_dat_201718$max.w_t <- LATR_dat_201718$max.w_t1
LATR_dat_201718$volume_t <- LATR_dat_201718$volume_t1
LATR_dat_201718$perp.w_t <- LATR_dat_201718$perp.w_t1
LATR_dat_201718$flowers_t <- LATR_dat_201718$flowers_t1
LATR_dat_201718$fruits_t <- LATR_dat_201718$fruits_t1
LATR_dat_201718$reproductive_fraction_t <- LATR_dat_201718$reproductive_fraction_t1
LATR_dat_201718$total.reproduction_t <- LATR_dat_201718$total.reproduction_t1

# Now set all the t1 data to NA
LATR_dat_201718$max.ht_t1 <- NA
LATR_dat_201718$max.w_t1 <- NA
LATR_dat_201718$volume_t1 <- NA
LATR_dat_201718$perp.w_t1 <- NA
LATR_dat_201718$flowers_t1 <- NA
LATR_dat_201718$fruits_t1 <- NA
LATR_dat_201718$reproductive_fraction_t1 <- NA
LATR_dat_201718$total.reproduction_t1 <- NA

# Bind rows and create log_vol as new variables (easier for GAMs)
LATR_flow_dat <- bind_rows(LATR_full,LATR_dat_201718) %>% 
  dplyr::select(unique.transect,volume_t,total.reproduction_t,weighted.dens) %>% drop_na()
LATR_flow_dat$log_volume_t <- log(LATR_flow_dat$volume_t)

# Create empty list to populate with model results
LATR_flower <- list()

# Three candidate models for the mean: size only, size + density, or size, density, and size:density
LATR_flower[[1]] <- gam(total.reproduction_t > 0 ~ s(log_volume_t) + s(unique.transect, bs = "re"),
                        data = LATR_flow_dat, gamma = gamma, family = "binomial")
LATR_flower[[2]] <- gam(total.reproduction_t > 0 ~ s(log_volume_t) + s(weighted.dens) + s(unique.transect, bs = "re"),
                        data = LATR_flow_dat, gamma = gamma, family = "binomial")
LATR_flower[[3]] <- gam(total.reproduction_t > 0 ~ s(log_volume_t) + s(weighted.dens) + ti(log_volume_t,weighted.dens) + s(unique.transect, bs = "re"),
                        data = LATR_flow_dat, gamma = gamma, family = "binomial")

# Collect model AICs into a single table
flower_aic<-AICtab(LATR_flower, base = TRUE, sort = FALSE)

# Set top model as "best"
LATR_flower_best <- LATR_flower[[which.min(flower_aic$AIC)]]
LATR_flower_fitted_terms <- predict(LATR_flower_best, type = "terms") 
LATR_flow_dat$pred <- predict.gam(LATR_flower_best, newdata = LATR_flow_dat, exclude = "s(unique.transect)")

##### Fruit production model ------------------------------------------------------------------------------

# Create new df with plants that have produced at least one reproductive structure
LATR_fruits_dat <- subset(LATR_flow_dat, total.reproduction_t > 0)

# Create empty list to populate with model results
LATR_fruits <- list()

# Three candidate models for the mean: size only, size + density, or size, density, and size:density
LATR_fruits[[1]] <- gam(total.reproduction_t ~ s(log_volume_t) + s(unique.transect, bs = "re"),
                        data = LATR_fruits_dat, gamma = gamma, family = "nb")
LATR_fruits[[2]] <- gam(total.reproduction_t ~ s(log_volume_t) + s(weighted.dens) + s(unique.transect, bs = "re"),
                        data = LATR_fruits_dat, gamma = gamma, family = "nb")
LATR_fruits[[3]] <- gam(total.reproduction_t ~ s(log_volume_t) + s(weighted.dens) + ti(log_volume_t,weighted.dens) + s(unique.transect, bs = "re"),
                        data = LATR_fruits_dat, gamma = gamma, family = "nb")
# Collect model AICs into a single table
fruits_aic <- AICtab(LATR_fruits, base = TRUE, sort = FALSE)

# Set top model as "best"
LATR_fruits_best <- LATR_fruits[[which.min(fruits_aic$AIC)]]
LATR_fruits_fitted_terms <- predict(LATR_fruits_best, type = "terms") 
LATR_fruits_dat$pred <- predict.gam(LATR_fruits_best, newdata = LATR_fruits_dat, exclude = "s(unique.transect)")

# Plot effect of size on fruits
# plot(LATR_fruits_dat$log_volume_t, LATR_fruits_fitted_terms[, "s(log_volume_t)"]) 

# Plot effect of density on fruits 
# plot(LATR_fruits_dat$weighted.dens, LATR_fruits_fitted_terms[, "s(weighted.dens)"]) 

##### Survival model --------------------------------------------------------------------------------------

# Combine transplants with large shrubs; keep only location info, survival, volume, and density
CData.Transplants %>% 
  dplyr::select("site", "transect", "actual.window", 
                "spring_survival_t1", "volume_t", "weighted.dens", "transplant") %>% 
  rename("survival_t1" = "spring_survival_t1") %>% 
  mutate(unique.transect = interaction(transect, site)) %>% 
  rbind(dplyr::select(LATR_full, "site", "transect", "actual.window", 
                      "survival_t1", "volume_t", "weighted.dens", "transplant","unique.transect")) %>% 
  mutate(log_volume_t = log(volume_t)) %>% 
  drop_na() -> LATR_surv_dat

# Investigate size overlap between transplant experiment and observational census
#hist(log(LATR_surv_dat$volume_t[LATR_surv_dat$transplant == FALSE]), breaks = 25)
#hist(log(LATR_surv_dat$volume_t[LATR_surv_dat$transplant == TRUE]), breaks = 10, add = TRUE, col = alpha("gray", 0.5))

# Plot survival against volume, grouped by transplant status
#plot(log(LATR_surv_dat$volume_t[LATR_surv_dat$transplant == FALSE]),
#     LATR_surv_dat$survival_t1[LATR_surv_dat$transplant == FALSE])
#points(log(LATR_surv_dat$volume_t[LATR_surv_dat$transplant == TRUE]),
#       LATR_surv_dat$survival_t1[LATR_surv_dat$transplant == TRUE] - 0.025, pch = 2)

# Create empty list to populate with model results
LATR_surv <- list()
# Three candidate models for the mean: size only, size + density, or size, density, and size:density
LATR_surv[[1]] <- gam(survival_t1 ~ s(log_volume_t,by=as.factor(transplant)) + s(unique.transect, bs = "re"),
                      data = LATR_surv_dat, gamma = gamma, family = "binomial")
LATR_surv[[2]] <- gam(survival_t1 ~ s(log_volume_t,by=as.factor(transplant)) + s(weighted.dens,by=as.factor(transplant))  + s(unique.transect, bs = "re"),
                      data = LATR_surv_dat, gamma = gamma, family = "binomial")
LATR_surv[[3]] <- gam(survival_t1 ~ s(log_volume_t,by=as.factor(transplant)) + s(weighted.dens,by=as.factor(transplant)) + ti(log_volume_t,weighted.dens,by=as.factor(transplant)) + s(unique.transect, bs = "re"),
                      data = LATR_surv_dat, gamma = gamma, family = "binomial")

# Collect model AICs into a single table
surv_aic <- AICtab(LATR_surv, base = TRUE, sort = FALSE)

# Set top model as "best"
LATR_surv_best <- LATR_surv[[which.min(surv_aic$AIC)]]
LATR_surv_fitted_terms <- predict(LATR_surv_best, type = "terms") 
LATR_surv_dat$pred <- predict.gam(LATR_surv_best, newdata = LATR_surv_dat, exclude = "s(unique.transect)")

# Plot effect of size on pr(survival)
# plot(LATR_surv_dat$log_volume_t, LATR_surv_fitted_terms[, "s(log_volume_t)"]) 

# Plot effect of density on pr(survival)
# plot(LATR_surv_dat$weighted.dens, LATR_surv_fitted_terms[, "s(weighted.dens)"]) 

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
# Note: we assume 6 seeds per fruit -- updated 18 Apr 2021 to 5
LATR_transects <- Cdata.Transects.Windows %>% 
  mutate(unique.transect = interaction(transect, site),
         log_volume_t = log(volume))
LATR_transects$seeds = ceiling(invlogit(predict.gam(LATR_flower_best, newdata = LATR_transects))* 
                                 seeds_per_fruit*exp(predict.gam(LATR_fruits_best, newdata = LATR_transects)))
LATR_transects <- group_by(LATR_transects, unique.transect, window)
suppressMessages(summarise(LATR_transects, total_seeds = sum(seeds),
                           weighted.dens = unique(weighted.dens))) -> LATR_transects

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
LATR_recruit[[2]] <- gam(cbind(recruits,total_seeds - recruits) ~ s(weighted.dens) + s(unique.transect, bs = "re"),
                         data = LATR_recruitment, gamma = gamma, family = "binomial")

# Collect model AICs into a single table
recruit_aic <- AICtab(LATR_recruit, base = TRUE, sort = FALSE)

# Set top model as "best"
LATR_recruit_best <- LATR_recruit[[which.min(recruit_aic$AIC)]]

# Plot effect of density on pr(seedling recruitment); negative density dependence
# LATR_recruit_fitted_terms <- predict(LATR_recruit[[2]], type = "terms") 
# plot(LATR_recruitment$weighted.dens,LATR_recruit_fitted_terms[, "s(weighted.dens)"])

##### Recruit sizes and integration limits (size bounds) --------------------------------------------------


  
  # Filter out seedlings and get their sizes
  LATR_recruit_size <- LATR_full %>% 
    filter(seedling_t1 == 1) %>% 
    mutate(log_volume = log(volume_t1)) %>% 
    arrange(weighted.dens)

# Plot distribution of recruit sizes using hist(LATR_recruit_size$log_volume)
#hist(LATR_recruit_size$log_volume)
# test for density dependence in recruit size
LATR_recruit_size_mod <- list()
LATR_recruit_size_mod[[1]] <- gam(list(log_volume ~ 1 + s(unique.transect,bs = "re"),~1), 
                                  data = LATR_recruit_size, gamma = gamma, family = gaulss())
LATR_recruit_size_mod[[2]] <- gam(list(log_volume ~ s(weighted.dens) + s(unique.transect,bs = "re"),~1), 
                                  data = LATR_recruit_size, gamma = gamma, family = gaulss())
LATR_recruit_size_mod[[3]] <- gam(list(log_volume ~ 1 + s(unique.transect,bs = "re"), ~s(weighted.dens)), 
                                  data = LATR_recruit_size, gamma = gamma, family = gaulss())
LATR_recruit_size_mod[[4]] <- gam(list(log_volume ~ s(weighted.dens) + s(unique.transect,bs = "re"), ~s(weighted.dens)), 
                                  data = LATR_recruit_size, gamma = gamma, family = gaulss())
# Collect model AICs into a single table
recruitsize_aic <- AICtab(LATR_recruit_size_mod, base = TRUE, sort = FALSE)

# Set top model as "best"
LATR_recruitsize_best <- LATR_recruit_size_mod[[which.min(recruitsize_aic$AIC)]]
LATR_recruitsize_fitted_terms <- predict(LATR_recruitsize_best, type = "terms") 
LATR_recruit_size$pred <- predict.gam(LATR_recruitsize_best, newdata = LATR_recruit_size, exclude = "s(unique.transect)")
## annoying but necessary index wrangling
recruit_size_sd_index <- which(as.factor(names(coef(LATR_recruitsize_best)))=="(Intercept).1") ## this is where the sd coefficients start
recruit_size_coef_length <- length(coef(LATR_recruitsize_best))

# Create maximum and minimum size bounds for the IPM
LATR_size_bounds <- data.frame(min_size = log(min(LATR_full$volume_t, LATR_full$volume_t1[LATR_full$transplant == FALSE], na.rm = TRUE)),
                               max_size = log(max(LATR_full$volume_t, LATR_full$volume_t1[LATR_full$transplant == FALSE], na.rm = TRUE)))



##### IPM functions ---------------------------------------------------------------------------------------
# Give dimension of the size vector (# bins of the approximating matrix)
TM.matdim <- 200

# Eviction extensions for upper and lower size limits
TM.lower.extension <- -8
TM.upper.extension <- 2

# Growth from size x to y at density d, using best GAM -- GAUSSIAN
TM.growth.GAU <- function(x, y, d){
  xb = pmin(pmax(x, LATR_size_bounds$min_size), LATR_size_bounds$max_size)
  mu = predict(LATR_m1,
          newdata = data.frame(weighted.dens = d, log_volume_t = xb, unique.transect = "1.FPS"),
          re.form = NA)
  sigma = sdfit$estimate[1]*exp(sdfit$estimate[2]*mu)
  return(dnorm(y,mu,sigma))
}

TM.growth.JSU <- function(x, y, d){
  xb = pmin(pmax(x, LATR_size_bounds$min_size), LATR_size_bounds$max_size)
  mu = predict(LATR_m1,
               newdata = data.frame(weighted.dens = d, log_volume_t = xb, unique.transect = "1.FPS"),
               re.form = NA)
  return(dJSU(y,
              mu=mu,
              sigma=exp(out1$estimate[1]+out1$estimate[2]*mu),
              nu=out1$estimate[3]+out1$estimate[4]*mu+out1$estimate[7]*mu^2,
              #nu=ifelse(LATR_grow$fitted<=out1$estimate[7],out1$estimate[3]+out1$estimate[4]*LATR_grow$fitted,out1$estimate[8]),
              tau=exp(out1$estimate[5]+out1$estimate[6]*mu)))
  }


# Survival of size x at density d using best GAM
# For nnaturally occuring plants (transplant = FALSE)
TM.survival <- function(x, d, elas, sens){
  xb = pmin(pmax(x, LATR_size_bounds$min_size), LATR_size_bounds$max_size)
  lpmat <- predict.gam(LATR_surv_best,
                       newdata = data.frame(weighted.dens = d, log_volume_t = xb, transplant = FALSE,
                                            unique.transect = "1.FPS"),
                       type = "lpmatrix",
                       exclude = "s(unique.transect)")
  pred <- invlogit(lpmat %*% coef(LATR_surv_best))
  if(elas=="survival"){pred<-pred*(1+pert)}
  if(sens=="survival"){pred<-pred+pert}
  return(pred)
}

# Combined growth and survival at density d
TM.growsurv <- function(x, y, d, elas, sens, dist){
  result <- TM.survival(x, d, elas, sens) * do.call(paste0("TM.growth.",dist),list(x, y, d))
  return(result)
  }

# Flowering at size x and density d using best GAM
TM.flower <- function(x, d, elas, sens){
  xb = pmin(pmax(x, LATR_size_bounds$min_size), LATR_size_bounds$max_size)
  lpmat <- predict.gam(LATR_flower_best,
                       newdata = data.frame(weighted.dens = d, log_volume_t = xb, unique.transect = "1.FPS"),
                       type = "lpmatrix",
                       exclude = "s(unique.transect)")
  pred <- invlogit(lpmat %*% coef(LATR_flower_best))
  if(elas=="flower"){pred<-pred*(1+pert)}
  if(sens=="flower"){pred<-pred+pert}
  return(pred)}

# Seed production (fruits * seeds/fruit) at size x and density d using best GAM
# Note: we assume 6 seeds per fruit, and best GAM is actually not density dependent
TM.seeds <- function(x, d, elas, sens, seeds.per.fruit = 5){
  xb = pmin(pmax(x, LATR_size_bounds$min_size), LATR_size_bounds$max_size)
  lpmat <- predict.gam(LATR_fruits_best,
                       newdata = data.frame(weighted.dens = d, log_volume_t = xb, unique.transect = "1.FPS"),
                       type = "lpmatrix",
                       exclude = "s(unique.transect)")
  pred <- exp(lpmat %*% coef(LATR_fruits_best))
  if(elas=="fertility"){pred<-pred*(1+pert)}
  if(sens=="fertility"){pred<-pred+pert}
  return(pred*seeds.per.fruit)}

# Seed-to-Seedling recruitment probability at density d
TM.recruitment <- function(d, elas, sens){
  lpmat <- predict.gam(LATR_recruit_best,
                       newdata = data.frame(weighted.dens = d, unique.transect = "1.FPS"),
                       type = "lpmatrix",
                       exclude = "s(unique.transect)")
  pred <- lpmat%*% coef(LATR_recruit_best)
  pred <- invlogit(pred[[1]])
  if(elas=="recruitment"){pred<-pred*(1-pert)} ## note substraction here bc pred is negative
  if(sens=="recruitment"){pred<-pred+pert} 
  return(pred)}

# Recruit size distribution at size y
TM.recruitsize <- function(y,d,elas,sens){
  lpmat <- predict.gam(LATR_recruitsize_best,
                       newdata = data.frame(weighted.dens = d, unique.transect = "1.FPS"),
                       type = "lpmatrix",
                       exclude = "s(unique.transect)")
  recruitsize_mu <- lpmat[, 1:(recruit_size_sd_index-1)] %*% coef(LATR_recruitsize_best)[1:(recruit_size_sd_index-1)]
  if(elas=="recruitsize.mean"){recruitsize_mu<-recruitsize_mu*(1+pert)}
  if(sens=="recruitsize.mean"){recruitsize_mu<-recruitsize_mu+pert}
  recruitsize_sigma <- exp(lpmat[, recruit_size_sd_index:recruit_size_coef_length] %*% coef(LATR_recruitsize_best)[recruit_size_sd_index:recruit_size_coef_length])
  if(elas=="recruitsize.sd"){recruitsize_sigma<-recruitsize_sigma*(1+pert)}
  if(sens=="recruitsize.sd"){recruitsize_sigma<-recruitsize_sigma+pert}
  return(dnorm(x = y, mean = recruitsize_mu, sd = recruitsize_sigma))
}

# Combined flowering, fertility, and recruitment
TM.fertrecruit <- function(x, y, d, elas, sens){
  TM.flower(x, d, elas, sens) * TM.seeds(x, d, elas, sens) * TM.recruitment(d,elas,sens) * TM.recruitsize(y,d,elas,sens)}

# Put it all together; projection matrix is a function of weighted density (dens)
# We need a large lower extension because growth variance (gaussian) is greater for smaller plants
TransMatrix <- function(dens, ext.lower = TM.lower.extension, ext.upper = TM.upper.extension,
                        min.size = LATR_size_bounds$min_size, max.size = LATR_size_bounds$max_size,
                        mat.size = TM.matdim,elas="none",sens="none",dist){
  
  # Matrix size and size extensions (upper and lower integration limits)
  n <- mat.size
  L <- min.size + ext.lower
  U <- max.size + ext.upper
  
  # Bin size for n bins
  h <- (U - L)/n
  
  # Lower boundaries of bins 
  b <- L + c(0:n)*h
  
  # Bin midpoints
  y <- 0.5*(b[1:n] + b[2:(n + 1)])
  
  # Growth/Survival matrix
  Pmat <- t(outer(y, y, TM.growsurv, d = dens, elas=elas, sens=sens, dist=dist)) * h 
  
  # Fertility/Recruiment matrix
  Fmat <- t(outer(y, y, TM.fertrecruit, d = dens, elas=elas, sens=sens)) * h 
  
  # Put it all together
  IPMmat <- Pmat + Fmat
  
  #and transition matrix
  return(list(IPMmat = IPMmat, Fmat = Fmat, Pmat = Pmat, meshpts = y))}

## look at lambda as a function of density
dens=seq(min(LATR_full$weighted.dens),max(LATR_full$weighted.dens),10)
lambda.GAU<-lambda.JSU<-c()
for(i in 1:length(dens)){
  lambda.GAU[i]<-lambda(TransMatrix(dens=dens[i],dist="GAU")$IPMmat)
  lambda.JSU[i]<-lambda(TransMatrix(dens=dens[i],dist="JSU")$IPMmat)
}

pdf("./manuscript/figures/creosote_DD_lambda.pdf",height = 6, width = 6,useDingbats = F)
plot(dens,lambda.JSU,type="b",col=alpha("red",0.5),pch=16,cex=1.2,
     xlab="Weighted density",ylab=expression(paste(lambda)))
lines(dens,lambda.GAU,type="b",col=alpha("blue",0.5),pch=16,cex=1.2)
legend("topright",bty="n",legend=c("Gaussian","JSU"),
       pch=16,col=c(alpha("blue",0.5),alpha("red",0.5)))
dev.off()

## wavespeed
boot.on = FALSE
boot.num <- 1
boot.switch <- TRUE
seeds <- sample.int(100000,size=boot.num)
# "01_SeedVelocities"
# Calculate terminal velocities of seeds and fit to distributions
source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/01_SeedVelocities.R")
# "02_WindSpeeds"
# Load in wind speeds and fit to distributions
source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/02_WindSpeeds.R")
# "03_Dispersal"
# Construct dispersal kernel functions for seeds
source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/03_Dispersal.R")
# "06_BootRes"
# Run resampling subroutine for wind speeds, terminal velocities, and demography
source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/06_BootRes.R")

# Construct transition matrix for minimum weighted density (zero)
TM_GAU <- TransMatrix(dens = 0, mat.size = 100,dist="GAU")
TM_JSU <- TransMatrix(dens = 0, mat.size = 100,dist="JSU")
#get WALD parameters from this data bootstrap - generates a list of length mat.size containing
# plant heights and correpsonding WALD parameters
## Function to compute the WALD mgf
WALDmgf <- function(s,nu,lambda) {
  t1 <- (lambda/nu) 
  t2 <- 2*(nu^2)*s/lambda 
  mgf <- exp(t1*(1-sqrt(1-t2)))
  return(mgf)
}    

## Function to compute the marginalize WALD mgf
margWALDmgf <- function(s,nu,lambda) {
  (1/pi)*integrate(function(q) WALDmgf(s*cos(q),nu,lambda),0,pi)$value
}
WALD_par <- function(h=0.15,elas="none",sens="none",TM){
  
  # Fit equation to convert volume to height for dispersal kernel use
  LATR_full %>%
    dplyr::select(max.ht_t, volume_t) %>% 
    drop_na(max.ht_t, volume_t) %>% 
    rename("h" = max.ht_t, "v" = volume_t) %>% 
    arrange(v) %>% 
    nlsLM(h ~ A*v^(1/3),
          start = list(A = 0), data = .) %>% 
    coef() %>% 
    as.numeric() -> A
  
  # Function converting volume to height (embedded here bc LATR_full will change with bootstrap iterations)
  # returns height in centimeters
  vol.to.height <- function(v){A*v^(1/3)}
  
  # size vector (log(volume))
  zvals <- TM$meshpts
  ## eviction problem!! if the zval is below true min or above true max, set to true min and true max
  zvals[zvals<LATR_size_bounds$min_size]=LATR_size_bounds$min_size
  zvals[zvals>=LATR_size_bounds$max_size]=LATR_size_bounds$max_size
  
  # Vector of heights across which dispersal kernel will be evaluated
  heights <- sapply(exp(zvals), vol.to.height)/100
  WALD.par <- vector("list",length(heights))
  WALD.par[heights>=h] <- lapply(heights[heights>=h],WALD.b.tom,elas=elas,sens=sens)
  
  return(list(heights=heights,WALD.par=WALD.par))
}
WALD_samples<-function(N,h=0.15,elas="none",sens="none",seed=NULL,params){
  r=matrix(0,nrow=N,ncol=length(params$heights))
  r[,params$heights>h]=sapply(params$heights[params$heights>h],WALD.f.e.h.tom,n=N,elas=elas,sens=sens,seed=seed)
  alpha <- matrix(runif(N*length(params$heights),0,2*pi),dim(r))
  X=r*cos(alpha)
  return(X)
}
## Function to compute wave speeds c(s)
cs <- function(s,h=0.15,emp=F,params,D.samples,TM) {
  # survival-growth matrix
  P <- TM$Pmat 
  # fertility matrix
  Fs <- TM$Fmat 
  for(j in 1:length(params$heights)){
    if(params$heights[j]<h){next}
    if(emp==F){Fs[,j] <- Fs[,j]*margWALDmgf(s,nu=params$WALD.par[[j]]$nu,lambda=params$WALD.par[[j]]$lambda)}
    if(emp==T){Fs[,j] <- Fs[,j]*empiricalWALDmgf(s,D.samples[,j])}
  }
  Hs <- P+Fs 
  L1 = abs(eigen(Hs)$values[1]); 
  return((1/s)*log(L1))  #
}

params_GAU <- WALD_par(TM=TM_GAU)
params_JSU <- WALD_par(TM=TM_JSU)
#sample dispersal events for empirical MGF - generates a N*mat.size matrix
D.samples.GAU <- WALD_samples(N=10000,seed=seeds,params=params_GAU) 
D.samples.JSU <- WALD_samples(N=10000,seed=seeds,params=params_JSU) 

# Find the asymptotic wave speed c*(s) 
cstar_GAU <- optimize(cs,lower=0.05,upper=4,emp=F,TM=TM_GAU,
                      params=params_GAU,D.samples=D.samples.GAU)$objective
cstar_JSU <- optimize(cs,lower=0.05,upper=4,emp=F,TM=TM_JSU,
                      params=params_JSU,D.samples=D.samples.JSU)$objective

pdf("./manuscript/figures/creosote_DD_lambda.pdf",height = 6, width = 6,useDingbats = F)
plot(dens,lambda.JSU,type="b",col=alpha("red",0.5),pch=16,cex=1.2,
     xlab="Weighted density",ylab=expression(paste(lambda)))
lines(dens,lambda.GAU,type="b",col=alpha("blue",0.5),pch=16,cex=1.2)
legend("topright",bty="n",legend=c("Gaussian","JSU"),
       pch=16,col=c(alpha("blue",0.5),alpha("red",0.5)))
text(100,1.025,"c*_GAU = 0.033m/yr",col=alpha("blue",0.5))
text(100,1.02,"c*_JSU = 0.026m/yr",col=alpha("red",0.5))
dev.off()
