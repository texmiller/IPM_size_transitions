### move to the right local directory 
tom = "C:/Users/tm9/Dropbox/github/IPM_size_transitions"
steve = "c:/repos/IPM_size_transitions" 
home = ifelse(Sys.info()["user"] == "Ellner", steve, tom)
setwd(home); setwd("creosote"); 

library(lme4)
library(mgcv)
library(qgam)
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

# fit pilot gaussian growth model with mgcv
## define k parameter (basis number) for global use
k_param<-6
# specify sd as constant then re-fit with iterative re-weighting
LATR_GAU <- gam(list(log_volume_t1~s(log_volume_t,k=k_param) + s(dens_scaled,k=k_param) + s(unique.transect,bs="re"),~1), 
                     family="gaulss", data=LATR_grow, method="ML",gamma=1.4) 
fitted_GAU = predict(LATR_GAU,type="response",data=LATR_grow)                  
new_fitted_vals = fitted_GAU[,1]
LATR_grow$fitted_vals = new_fitted_vals
weights = fitted_GAU[,2] #1/sigma_i
#re-fit with sd~f(fitted)
LATR_GAU <- gam(list(log_volume_t1~s(log_volume_t,k=k_param) + s(dens_scaled,k=k_param) + s(unique.transect,bs="re"),~s(fitted_vals,k=k_param)), 
                family="gaulss", data=LATR_grow, method="ML",gamma=1.4) 
err=100; k=0; 
while(err>10^(-6)) {
  LATR_grow$fitted_vals = new_fitted_vals; 
  LATR_GAU<-update(LATR_GAU)
  fitted_all = predict(LATR_GAU,type="response",data=LATR_grow)   
  new_fitted_vals = fitted_all[,1]; new_weights = fitted_all[,2]
  err = weights - new_weights; err=sqrt(mean(err^2)) 
  weights = new_weights; 
  k=k+1; cat(k,err,"\n") 
}
out = predict(LATR_GAU,type="response")
plot(out[,1],1/out[,2]); 

LATR_GAU_best<-LATR_GAU
LATR_grow$GAU_mean<- fitted(LATR_GAU_best)[,1]
LATR_grow$GAU_sd<- 1/fitted(LATR_GAU_best)[,2]
LATR_grow$GAU_resids <- residuals(LATR_GAU_best,type="response")
LATR_grow$GAU_scaled_resids <- residuals(LATR_GAU_best,type="pearson") 
## should be mean zero unit variance
mean(LATR_grow$GAU_scaled_resids);sd(LATR_grow$GAU_scaled_resids)

saveRDS(LATR_grow,file="creosote_Gau_pilot.rds"); 

##are the standardized residuals gaussian? -- no
jarque.test(LATR_grow$GAU_scaled_resids) # normality test: FAILS, P < 0.001 
anscombe.test(LATR_grow$GAU_scaled_resids) # kurtosis: FAILS, P < 0.001 
agostino.test(LATR_grow$GAU_scaled_resids) # skewness: FAILS, P<0.001 


################################################################
## Steve inserts: tests for non-constant residual variance 
## Results: small p-value on one test, but minuscule trend
#################################################################
source("../code/variance_diagnostics.R"); 

stopCluster(c1); 
c1<- makeCluster(3); 
registerDoParallel(c1);
R = 2000; 

out_bartlett = multiple_bartlett_test(LATR_grow$log_volume_t,LATR_grow$GAU_scaled_resids, 3, 10, R) 
out_bartlett$p_value; ##  0.20


out_bs = multiple_bs_test(LATR_grow$log_volume_t,LATR_grow$GAU_scaled_resids, 4, 10, R) ## p = 0.47; 
out_bs$p_value; ## 0.34


out_bartlett = multiple_bartlett_test(LATR_grow$GAU_mean,LATR_grow$GAU_scaled_resids, 3, 10, R) 
out_bartlett$p_value; ## 0.02

out_bs = multiple_bs_test(LATR_grow$GAU_mean,LATR_grow$GAU_scaled_resids, 4, 10, R)  
out_bs$p_value; ## p = 0.47;

out_bartlett = multiple_bartlett_test(LATR_grow$dens_scaled,LATR_grow$GAU_scaled_resids, 3, 10, R);
out_bartlett$p_value;  # 0.371

out_bs = multiple_bs_test(LATR_grow$dens_scaled,LATR_grow$GAU_scaled_resids, 4, 10, R)  
out_bs$p_value; # 0.5415

stopCluster(c1);

require(mgcv); 
fit = gam(I(GAU_scaled_resids^2)~s(GAU_mean),data=LATR_grow)  ## not significant
fit = gam(abs(GAU_scaled_resids)~s(GAU_mean),data=LATR_grow)  ## p = 0.06, but <1% of variance
fit1 = lm(abs(GAU_scaled_resids)~(GAU_mean),data=LATR_grow) ## p = 0.06, but <1% of variance 

########### Back to Tom.... 

## quantile regressions on stand resids
S.05<-qgam(GAU_scaled_resids~s(GAU_mean,k=k_param), data=LATR_grow,qu=0.05)
S.10<-qgam(GAU_scaled_resids~s(GAU_mean,k=k_param), data=LATR_grow,qu=0.1)
S.25<-qgam(GAU_scaled_resids~s(GAU_mean,k=k_param), data=LATR_grow,qu=0.25)
S.50<-qgam(GAU_scaled_resids~s(GAU_mean,k=k_param), data=LATR_grow,qu=0.5) 
S.75<-qgam(GAU_scaled_resids~s(GAU_mean,k=k_param), data=LATR_grow,qu=0.75)
S.90<-qgam(GAU_scaled_resids~s(GAU_mean,k=k_param), data=LATR_grow,qu=0.9) 
S.95<-qgam(GAU_scaled_resids~s(GAU_mean,k=k_param), data=LATR_grow,qu=0.95)

## NP skewness
NPS_hat = Q.skewness(q.10=predict(S.10),
                     q.50=predict(S.50),
                     q.90=predict(S.90))

## NP kurtosis (relative to Gaussian)
NPK_hat = Q.kurtosis(q.05=predict(S.05),
                     q.25=predict(S.25),
                     q.75=predict(S.75),
                     q.95=predict(S.95))

## figure of sizet1/density conditional on sizet, sd(fitted), and scaled resids
dens_dummy<-seq(min(LATR_grow$dens_scaled),
                max(LATR_grow$dens_scaled),0.1)
size_means<-LATR_grow %>% group_by(size_bin) %>% summarise(mean=mean(log_volume_t),col=unique(size_col))
size1_pred<-predict.gam(LATR_GAU_best,newdata=data.frame(dens_scaled=dens_dummy,log_volume_t=size_means$mean[1],fitted_vals=10,unique.transect="1.FPS"),
                        type = "response",
                        exclude = "s(unique.transect)")[,1]
size2_pred<-predict.gam(LATR_GAU_best,newdata=data.frame(dens_scaled=dens_dummy,log_volume_t=size_means$mean[2],fitted_vals=10,unique.transect="1.FPS"),
                        type = "response",
                        exclude = "s(unique.transect)")[,1]
size3_pred<-predict.gam(LATR_GAU_best,newdata=data.frame(dens_scaled=dens_dummy,log_volume_t=size_means$mean[3],fitted_vals=10,unique.transect="1.FPS"),
                        type = "response",
                        exclude = "s(unique.transect)")[,1]

pdf("../manuscript/figures/creosote_diagnostics.pdf",height = 3, width = 9,useDingbats = F)
par(mar = c(5, 4, 2, 2), oma=c(0,0,0,2), mfrow=c(1,3)) 
plot(LATR_grow$dens_scaled*100,LATR_grow$log_volume_t1,
     col=alpha(LATR_grow$size_col,0.5),pch=16,
     xlab="Weighted density",ylab="Size at t+1")
lines(dens_dummy*100,size1_pred,col=dens_means$col[1],lwd=2)
lines(dens_dummy*100,size2_pred,col=dens_means$col[2],lwd=2)
lines(dens_dummy*100,size3_pred,col=dens_means$col[3],lwd=2)
legend("bottomright",title="Size at t",legend=round(size_means$mean,2),
       bty="n",cex=0.8,pch=16,col=dens_pallete)
title("A",adj=0,font=3)

plot(LATR_grow$GAU_mean,LATR_grow$GAU_sd,
     xlab="Expected size at t+1",ylab="SD of size at t+1",type="l",lwd=3)
title("B",adj=0,font=3)

plot(LATR_grow$GAU_mean,LATR_grow$GAU_scaled_resids,col=alpha("black",0.25),
     xlab="Expected size at t+1",ylab="Scaled residuals of size at t+1")
points(LATR_grow$GAU_mean,predict(S.05),col="black",pch=".")
points(LATR_grow$GAU_mean,predict(S.10),col="black",pch=".")
points(LATR_grow$GAU_mean,predict(S.25),col="black",pch=".")
points(LATR_grow$GAU_mean,predict(S.50),col="black",pch=".")
points(LATR_grow$GAU_mean,predict(S.75),col="black",pch=".")
points(LATR_grow$GAU_mean,predict(S.90),col="black",pch=".")
points(LATR_grow$GAU_mean,predict(S.95),col="black",pch=".")
par(new = TRUE)                           
plot(c(LATR_grow$GAU_mean,LATR_grow$GAU_mean),
     c(NPS_hat,NPK_hat),
     col=c(rep(alpha("blue",0.25),nrow(LATR_grow)),rep(alpha("red",0.25),nrow(LATR_grow))),
     pch=16,cex=.5,axes = FALSE, xlab = "", ylab = "")
abline(h=0,col="lightgray",lty=3)
axis(side = 4,cex.axis=0.8,at = pretty(range(c(NPS_hat,NPK_hat))))
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
          mu=LATR_grow$GAU_mean,
          sigma=LATR_grow$GAU_sd,
          nu = pars[1]+pars[2]*LATR_grow$GAU_mean,
          tau = exp(pars[3]+pars[4]*LATR_grow$GAU_mean), log=TRUE)
}
## starting parameters
p0<-c(0,0,0,0)

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
                                      mu=LATR_grow$GAU_mean,
                                      sigma=LATR_grow$GAU_sd,
                                      nu=JSUout$estimate[1]+JSUout$estimate[2]*LATR_grow$GAU_mean,
                                      tau=exp(JSUout$estimate[3]+JSUout$estimate[4]*LATR_grow$GAU_mean))
  S.05<-qgam(log_volume_t1.simJSU~s(log_volume_t,k=k_param), data=LATR_grow,qu=0.05)
  S.10<-qgam(log_volume_t1.simJSU~s(log_volume_t,k=k_param), data=LATR_grow,qu=0.1)
  S.25<-qgam(log_volume_t1.simJSU~s(log_volume_t,k=k_param), data=LATR_grow,qu=0.25)
  S.50<-qgam(log_volume_t1.simJSU~s(log_volume_t,k=k_param), data=LATR_grow,qu=0.5) 
  S.75<-qgam(log_volume_t1.simJSU~s(log_volume_t,k=k_param), data=LATR_grow,qu=0.75)
  S.90<-qgam(log_volume_t1.simJSU~s(log_volume_t,k=k_param), data=LATR_grow,qu=0.9) 
  S.95<-qgam(log_volume_t1.simJSU~s(log_volume_t,k=k_param), data=LATR_grow,qu=0.95)
  JSUsim_mean[,i]<-Q.mean(predict(S.25),predict(S.50),predict(S.75))
  JSUsim_sd[,i]<-Q.sd(predict(S.25),predict(S.75))
  JSUsim_skew[,i]<-Q.skewness(predict(S.10),predict(S.50),predict(S.90))
  JSUsim_kurt[,i]<-Q.kurtosis(predict(S.05),predict(S.25),predict(S.75),predict(S.95))
  
  LATR_grow$log_volume_t1.simGAU <- rnorm(n=nrow(LATR_grow),
                                         mean=LATR_grow$GAU_mean,
                                         sd=LATR_grow$GAU_sd)
  S.05<-qgam(log_volume_t1.simGAU~s(log_volume_t,k=k_param), data=LATR_grow,qu=0.05)
  S.10<-qgam(log_volume_t1.simGAU~s(log_volume_t,k=k_param), data=LATR_grow,qu=0.1)
  S.25<-qgam(log_volume_t1.simGAU~s(log_volume_t,k=k_param), data=LATR_grow,qu=0.25)
  S.50<-qgam(log_volume_t1.simGAU~s(log_volume_t,k=k_param), data=LATR_grow,qu=0.5) 
  S.75<-qgam(log_volume_t1.simGAU~s(log_volume_t,k=k_param), data=LATR_grow,qu=0.75)
  S.90<-qgam(log_volume_t1.simGAU~s(log_volume_t,k=k_param), data=LATR_grow,qu=0.9) 
  S.95<-qgam(log_volume_t1.simGAU~s(log_volume_t,k=k_param), data=LATR_grow,qu=0.95)
  GAUsim_mean[,i]<-Q.mean(predict(S.25),predict(S.50),predict(S.75))
  GAUsim_sd[,i]<-Q.sd(predict(S.25),predict(S.75))
  GAUsim_skew[,i]<-Q.skewness(predict(S.10),predict(S.50),predict(S.90))
  GAUsim_kurt[,i]<-Q.kurtosis(predict(S.05),predict(S.25),predict(S.75),predict(S.95))
}

## and now the real data
q.05<-predict(qgam(log_volume_t1~s(log_volume_t,k=k_param), data=LATR_grow,qu=0.05))
q.10<-predict(qgam(log_volume_t1~s(log_volume_t,k=k_param), data=LATR_grow,qu=0.1))
q.25<-predict(qgam(log_volume_t1~s(log_volume_t,k=k_param), data=LATR_grow,qu=0.25))
q.50<-predict(qgam(log_volume_t1~s(log_volume_t,k=k_param), data=LATR_grow,qu=0.5) )
q.75<-predict(qgam(log_volume_t1~s(log_volume_t,k=k_param), data=LATR_grow,qu=0.75))
q.90<-predict(qgam(log_volume_t1~s(log_volume_t,k=k_param), data=LATR_grow,qu=0.9) )
q.95<-predict(qgam(log_volume_t1~s(log_volume_t,k=k_param), data=LATR_grow,qu=0.95))


pdf("../manuscript/figures/creosote_JSU_fit.pdf",height = 6, width = 6,useDingbats = F)
par(mfrow=c(2,2),mar=c(4,4,1,1),mgp=c(2.1,1,0),cex.lab=1.2); 
plot(LATR_grow$log_volume_t,Q.mean(q.25,q.50,q.75),type="n",
     xlab="Size(t)",ylab="NP mean size(t+1)",
     ylim=c(min(c(GAUsim_mean,JSUsim_mean)),1 + max(c(GAUsim_mean,JSUsim_mean))))
     matpoints(LATR_grow$log_volume_t,GAUsim_mean,col=alpha("tomato",0.25),type="l",lty=1)
     matpoints(LATR_grow$log_volume_t,1 + JSUsim_mean,col=alpha("cornflowerblue",0.25),type="l",lty=1)
lines(LATR_grow$log_volume_t,Q.mean(q.25,q.50,q.75),col="black",lwd=1.5)
lines(LATR_grow$log_volume_t,1 + Q.mean(q.25,q.50,q.75),col="black",lwd=1.5)
legend("topleft",legend=c("Real data","Gaussian simulation","JSU simulation + offset"),
       lty=1,col=c("black","tomato","cornflowerblue"),cex=0.8,bty="n")

plot(LATR_grow$log_volume_t,Q.sd(q.25,q.75),type="n",
     xlab="Size(t)",ylab="NP SD Size(t+1)",
     ylim=c(min(c(GAUsim_sd,JSUsim_sd)),1 + max(c(GAUsim_sd,JSUsim_sd))))
     matpoints(LATR_grow$log_volume_t,GAUsim_sd,col=alpha("tomato",0.25),type="l",lty=1)
     matpoints(LATR_grow$log_volume_t,1 + JSUsim_sd,col=alpha("cornflowerblue",0.25),type="l",lty=1)
lines(LATR_grow$log_volume_t,Q.sd(q.25,q.75),col="black",lwd=1.5)
lines(LATR_grow$log_volume_t,1 + Q.sd(q.25,q.75),col="black",lwd=1.5)

plot(LATR_grow$log_volume_t,Q.skewness(q.10,q.50,q.90),type="n",
     xlab="Size(t)",ylab="NP skewness size(t+1)",
     ylim=c(min(c(GAUsim_skew,JSUsim_skew)),1 + max(c(GAUsim_skew,JSUsim_skew))))
     matpoints(LATR_grow$log_volume_t,GAUsim_skew,col=alpha("tomato",0.25),type="l",lty=1)
     matpoints(LATR_grow$log_volume_t,1+ JSUsim_skew,col=alpha("cornflowerblue",0.25),type="l",lty=1)
    lines(LATR_grow$log_volume_t,Q.skewness(q.10,q.50,q.90),col="black",lwd=1.5)
    lines(LATR_grow$log_volume_t,1 + Q.skewness(q.10,q.50,q.90),col="black",lwd=1.5)

plot(LATR_grow$log_volume_t,Q.kurtosis(q.05,q.25,q.75,q.95),type="n",
     xlab="Size(t)",ylab="NP kurtosis size(t+1)",
     ylim=c(min(GAUsim_kurt,JSUsim_kurt),1 + max(GAUsim_kurt,JSUsim_kurt)))
     matpoints(LATR_grow$log_volume_t,GAUsim_kurt,col=alpha("tomato",0.25),type="l",lty=1)
     matpoints(LATR_grow$log_volume_t,1 + JSUsim_kurt,col=alpha("cornflowerblue",0.25),type="l",lty=1)
lines(LATR_grow$log_volume_t,Q.kurtosis(q.05,q.25,q.75,q.95),col="black",lwd=1.5)
lines(LATR_grow$log_volume_t,1 + Q.kurtosis(q.05,q.25,q.75,q.95),col="black",lwd=1.5)
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
## annoying but necessary index wrangling
recruit_size_sd_index <- which(as.factor(names(coef(LATR_recruitsize_best)))=="(Intercept).1") ## this is where the sd coefficients start
recruit_size_coef_length <- length(coef(LATR_recruitsize_best))

# Create maximum and minimum size bounds for the IPM
LATR_size_bounds <- data.frame(min_size = log(min(LATR_full$volume_t, LATR_full$volume_t1[LATR_full$transplant == FALSE], na.rm = TRUE)),
                               max_size = log(max(LATR_full$volume_t, LATR_full$volume_t1[LATR_full$transplant == FALSE], na.rm = TRUE)))

# IPM functions -----------------------------------------------------------
# using the fitted model objects above

# Growth from size x to y at density d, using best GAUSSIAN
growth.GAU <- function(x, y, d){
  xb = pmin(pmax(x, LATR_size_bounds$min_size), LATR_size_bounds$max_size)
  ## need expected value before I can use predict.gam because fitted_val is an input
  pred_mat = predict.gam(LATR_GAU_best,
               newdata = data.frame(dens_scaled = d, log_volume_t = xb, fitted_vals = 0, unique.transect = "1.FPS"),
               type = "lpmatrix",
               exclude = "s(unique.transect)")
  mu = pred_mat[,1:19]%*%coef(LATR_GAU_best)[1:19]
  pred = predict.gam(LATR_GAU_best,
                     newdata = data.frame(dens_scaled = d, log_volume_t = xb, fitted_vals = mu, unique.transect = "1.FPS"),
                     type = "response",
                     exclude = "s(unique.transect)")
  return(dnorm(y,mean=pred[,1],1/pred[,2]))
}

# Growth from size x to y at density d, using JSU
growth.JSU <- function(x, y, d){
  xb = pmin(pmax(x, LATR_size_bounds$min_size), LATR_size_bounds$max_size)
  pred_mat = predict.gam(LATR_GAU_best,
                         newdata = data.frame(dens_scaled = d, log_volume_t = xb, fitted_vals = 0, unique.transect = "1.FPS"),
                         type = "lpmatrix",
                         exclude = "s(unique.transect)")
  mu = pred_mat[,1:19]%*%coef(LATR_GAU_best)[1:19]
  pred = predict.gam(LATR_GAU_best,
                     newdata = data.frame(dens_scaled = d, log_volume_t = xb, fitted_vals = mu, unique.transect = "1.FPS"),
                     type = "response",
                     exclude = "s(unique.transect)")
  return(dJSU(y,
              mu=pred[,1],
              sigma=1/pred[,2],
              nu=JSUout$estimate[1]+JSUout$estimate[2]*pred[,1],
              tau=exp(JSUout$estimate[3]+JSUout$estimate[4]*pred[,1])))
}

# Survival of size x at density d using best GAM
# For nnaturally occuring plants (transplant = FALSE)
survival <- function(x, d){
  xb = pmin(pmax(x, LATR_size_bounds$min_size), LATR_size_bounds$max_size)
  pred <- predict.gam(LATR_surv_best,
                      newdata = data.frame(dens_scaled = d, log_volume_t = xb, transplant = FALSE,
                                           unique.transect = "1.FPS"),
                      type = "response",
                      exclude = "s(unique.transect)")
  return(pred)
}

# Combined growth and survival at density d
growsurv <- function(x, y, d, dist){
  survival(x, d) * do.call(paste0("growth.",dist),list(x, y, d))
}

# Flowering at size x and density d using best GAM
flower <- function(x, d){
  xb = pmin(pmax(x, LATR_size_bounds$min_size), LATR_size_bounds$max_size)
  pred <- predict.gam(LATR_flower_best,
                      newdata = data.frame(dens_scaled = d, log_volume_t = xb, unique.transect = "1.FPS"),
                      type = "response",
                      exclude = "s(unique.transect)")
  return(pred)
}

# Seed production (fruits * seeds/fruit) at size x and density d using best GAM
seeds <- function(x, d, seeds.per.fruit = 5){
  xb = pmin(pmax(x, LATR_size_bounds$min_size), LATR_size_bounds$max_size)
  pred <- predict.gam(LATR_fruits_best,
                      newdata = data.frame(dens_scaled = d, log_volume_t = xb, unique.transect = "1.FPS"),
                      type = "response",
                      exclude = "s(unique.transect)")
  return(pred*seeds.per.fruit)
}

# Seed-to-Seedling recruitment probability at density d
recruitment <- function(d){
  pred <- predict.gam(LATR_recruit_best,
                      newdata = data.frame(dens_scaled = d, unique.transect = "1.FPS"),
                      type = "response",
                      exclude = "s(unique.transect)")
  return(pred[1])
}

# Recruit size distribution at size y
recruitsize <- function(y,d){
  lpmat <- predict.gam(LATR_recruitsize_best,
                       newdata = data.frame(dens_scaled = d, unique.transect = "1.FPS"),
                       type = "lpmatrix",
                       exclude = "s(unique.transect)")
  recruitsize_mu <- lpmat[, 1:(recruit_size_sd_index-1)] %*% coef(LATR_recruitsize_best)[1:(recruit_size_sd_index-1)]
  recruitsize_sigma <- exp(lpmat[, recruit_size_sd_index:recruit_size_coef_length] %*% coef(LATR_recruitsize_best)[recruit_size_sd_index:recruit_size_coef_length])
  return(dnorm(x = y, mean = recruitsize_mu, sd = recruitsize_sigma))
}

# Combined flowering, fertility, and recruitment
fertrecruit <- function(x, y, d){
  flower(x,d)*seeds(x,d)*recruitment(d)*recruitsize(y,d)
}

# Put it all together; projection matrix is a function of weighted density (dens)
# We need a large lower extension because growth variance is greater for smaller plants
ApproxMatrix <- function(dens,ext.lower,ext.upper,
                         min.size=LATR_size_bounds$min_size,max.size=LATR_size_bounds$max_size,
                         mat.size,dist){
  
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
  Pmat <- t(outer(y, y, growsurv, d = dens, dist=dist)) * h 
  
  # Fertility/Recruiment matrix
  Fmat <- t(outer(y, y, fertrecruit, d = dens)) * h 
  
  # Put it all together
  IPMmat <- Pmat + Fmat
  
  #and transition matrix
  return(list(IPMmat = IPMmat, Fmat = Fmat, Pmat = Pmat, meshpts = y))
}

## load SIPM source functions and some data elements for wind dispersal
## will take a hot minute and will throw warnings
source("creosote_SIPM_source_fns.R")

# Eviction extensions for upper and lower size limits
lower.extension <- -8
upper.extension <- 2
## explore effect of matrix dimensions -- or skip to matdim
dims<-c(100,150,200,250,300,350,400)
lambda0.GAU<-lambda0.JSU<-c()
meanlife.GAU<-meanlife.JSU<-c()
cstar_GAU<-cstar_JSU<-c()
for(i in 1:length(dims)){
  hold.GAU<-ApproxMatrix(dens=0,dist="GAU",mat.size=dims[i],ext.lower=lower.extension,ext.upper=upper.extension)
  hold.JSU<-ApproxMatrix(dens=0,dist="JSU",mat.size=dims[i],ext.lower=lower.extension,ext.upper=upper.extension)
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
  lambda.GAU[i]<-lambda(ApproxMatrix(dens=dens[i],dist="GAU",mat.size=matdim,ext.lower=lower.extension,ext.upper=upper.extension)$IPMmat)
  lambda.JSU[i]<-lambda(ApproxMatrix(dens=dens[i],dist="JSU",mat.size=matdim,ext.lower=lower.extension,ext.upper=upper.extension)$IPMmat)
  print(i)
}

# Construct transition matrix for minimum weighted density (zero)
mat_GAU <- ApproxMatrix(dens = 0, dist="GAU",mat.size=matdim,ext.lower=lower.extension,ext.upper=upper.extension)
mat_JSU <- ApproxMatrix(dens = 0, dist="JSU",mat.size=matdim,ext.lower=lower.extension,ext.upper=upper.extension)

## check eviction
plot(mat_GAU$meshpts,colSums(mat_GAU$Pmat))
lines(mat_GAU$meshpts,survival(mat_GAU$meshpts,d=0),col="red",lwd=2)##good
plot(mat_JSU$meshpts,colSums(mat_JSU$Pmat))
lines(mat_JSU$meshpts,survival(mat_JSU$meshpts,d=0),col="red",lwd=2)##good

## write out for cross-spp analysis
write_rds(list(mat_GAU=mat_GAU,mat_JSU=mat_JSU),file="creosote_out.rds")

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


pdf("../manuscript/figures/creosote_DD_lambda.pdf",height = 4, width = 4,useDingbats = F)
par(mar=c(4,4,1,1))
plot(dens*100,lambda.JSU,type="b",col=alpha("cornflowerblue",0.75),pch=16,cex=1.2,
     xlab="Weighted density",ylab=expression(paste(lambda)))
lines(dens*100,lambda.GAU,type="b",col=alpha("tomato",0.75),pch=16,cex=1.2)
legend("topright",bty="n",legend=c("Gaussian","JSU"),
       pch=16,col=c(alpha("tomato",1),alpha("cornflowerblue",1)))
text(140,1.025,paste("c*=",round(cstar_GAU,4),"m/yr"),col=alpha("tomato",1))
text(140,1.02,paste("c*=",round(cstar_JSU,4),"m/yr"),col=alpha("cornflowerblue",1))
dev.off()


life_expect_mean(mat_GAU$Pmat,start=1);life_expect_var(mat_GAU$Pmat,start=1)
life_expect_mean(mat_JSU$Pmat,start=1);life_expect_var(mat_JSU$Pmat,start=1)

## remaining life expectancy of median-sized shrub
median_index<-which.min(abs(mat_GAU$meshpts-median(LATR_full$log_volume_t,na.rm=T)))
life_expect_mean(mat_GAU$Pmat,start=median_index)
life_expect_mean(mat_JSU$Pmat,start=median_index)

plot(flower(mat_JSU$meshpts,d=0))
plot(survival(mat_JSU$meshpts,d=0))

colSums(mat_GAU$Pmat)
colSums(mat_JSU$Pmat)
