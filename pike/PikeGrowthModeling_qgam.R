## Purpose: analyze pike growth data and diagnose appropriate growth kernel
## Here using gam and qgam analysis of growth data and standardized residuals
## Author: Tom Miller (mash-up with code by Steve Ellner)

## set repo as wd
tom<-"C:/Users/tm9/Dropbox/github/IPM_size_transitions"
steve<-NULL
setwd(tom)

## load libraries
library(tidyverse)
library(lme4)
library(mgcv)
library(scales)
library(qgam)
library(gamlss.dist)
library(popbio)
library(moments)
library(maxLik)

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

## read in data
pike_dat <- read.csv("pike/data/PikeGrowthData1944_1995.csv")
## this will need some data manipulation, starting with wide format
pike_wide <- pike_dat %>% 
  arrange(Year,Ind) %>% dplyr::select(-RowID) %>% 
  pivot_wider(names_from = Year,values_from = Length)
## make sure that the column names are ordered years with no missing year
as.numeric(names(pike_wide)[-1]) == 1944:1995 #good
##now stack transition years
pike_trans_year <- pike_wide[,(1:3)]
pike_trans_year$year <- as.numeric(names(pike_wide)[2])
names(pike_trans_year)[2:3]<-c("t0","t1")
for(i in 3:(ncol(pike_wide)-1)){
  hold <- pike_wide[,c(1,i,(i+1))]
  hold$year <- as.numeric(names(pike_wide)[i])
  names(hold)[2:3]<-c("t0","t1")
  pike_trans_year <- bind_rows(pike_trans_year,hold)
}
## drop rows without t0 and t1 size data
pike_trans_year %>% filter(!is.na(t0) & !is.na(t1)) %>% 
  mutate(log_t0 = log(t0),
         log_t1 = log(t1),
         growth = log_t1 - log_t0) %>% 
  arrange(log_t0,log_t1)-> pike_final


plot(pike_final$t0,pike_final$t1)
plot(pike_final$log_t0,pike_final$log_t1)
abline(0,1)

## pilot gaussian with mgcv
pike_gau<-gam(list(log_t1 ~ s(log_t0,k=4), ~s(t0,k=4)), data=pike_final, family=gaulss())
pike_gau_pred<-predict(pike_gau,type="response")
## scale residuals by fitted sd
fitted_sd<-1/pike_gau_pred[,2]
pike_final$scaledResids=residuals(pike_gau,type="response")/fitted_sd

## residuals should have mean zero, unit variance...
## ...in aggregate
mean(pike_final$scaledResids);sd(pike_final$scaledResids)
## ...and also in size bins
pike_final %>% 
  mutate(size_bin=cut_number(log_t0,n=5)) %>% 
  group_by(size_bin) %>% 
  summarise(mean_stdresid = mean(scaledResids),
            sd_stdresid = sd(scaledResids))
## no too bad
  
## fit quantile gams to scaled resids
q.05<-predict(qgam(scaledResids~s(log_t0,k=4), data=pike_final,qu=0.05))
q.10<-predict(qgam(scaledResids~s(log_t0,k=4), data=pike_final,qu=0.1)) 
q.25<-predict(qgam(scaledResids~s(log_t0,k=4), data=pike_final,qu=0.25))
q.50<-predict(qgam(scaledResids~s(log_t0,k=4), data=pike_final,qu=0.5))
q.75<-predict(qgam(scaledResids~s(log_t0,k=4), data=pike_final,qu=0.75))
q.90<-predict(qgam(scaledResids~s(log_t0,k=4), data=pike_final,qu=0.9)) 
q.95<-predict(qgam(scaledResids~s(log_t0,k=4), data=pike_final,qu=0.95)) 
NPS_hat<-Q.skewness(q.10,q.50,q.90)
NPK_hat<-Q.kurtosis(q.05,q.25,q.75,q.95)

## mean and sd of size transition
plot(pike_final$log_t0,pike_final$log_t1,pch=".")
lines(pike_final$log_t0,pike_gau_pred[,1],col="red")
par(new = TRUE)                           
plot(pike_final$log_t0,1/pike_gau_pred[,2],col="blue",pch=16,cex=.5,
     axes = FALSE, xlab = "", ylab = "",type="l")
axis(side = 4, at = pretty(range(1/pike_gau_pred[,2])))

## scaled residuals plot
par(mar = c(5, 5, 2, 3), oma=c(0,0,0,2)) 
plot(pike_final$log_t0,pike_final$scaledResids,col=alpha("black",0.25),pch=".",
     xlab="log length, time t",ylab="Scaled residuals of log length at t+1")
lines(pike_final$log_t0,q.05,col="black",pch=".")
lines(pike_final$log_t0,q.10,col="black",pch=".")
lines(pike_final$log_t0,q.25,col="black",pch=".")
lines(pike_final$log_t0,q.50,col="black",pch=".")
lines(pike_final$log_t0,q.75,col="black",pch=".")
lines(pike_final$log_t0,q.90,col="black",pch=".")
lines(pike_final$log_t0,q.95,col="black",pch=".")
par(new = TRUE)                           
plot(c(pike_final$log_t0,pike_final$log_t0),c(NPS_hat,NPK_hat),
     col=c(rep(alpha("blue",0.25),nrow(pike_final)),rep(alpha("red",0.25),nrow(pike_final))),
     type="l",cex=.5, axes = FALSE, xlab = "", ylab = "")
abline(h=0,col="lightgray",lty=3)
axis(side = 4, at = pretty(range(c(NPS_hat,NPK_hat))),cex.axis=0.8)
mtext("Skewness", side = 4, line = 2,col="blue")
mtext("Excess Kurtosis", side = 4, line =3,col="red")

## gonna try a gam shash (bc that's what's easy in mgcv)
pike_gam_shash <- gam(list(log_t1 ~ s(log_t0,k=4), # <- model for location 
                           ~ s(log_t0,k=4),   # <- model for log-scale
                           ~ s(log_t0,k=4),   # <- model for skewness
                           ~ s(log_t0,k=4)), # <- model for log-kurtosis
                      data = pike_final, 
                      family = shash,  
                      optimizer = "efs")
pike_shash_pred <- predict(pike_gam_shash,type="response")

## also, out of curiosity, gonna fit a gamma to growth data (they can only increase)
pike_gam_gamma <- gam(list(growth ~ s(log_t0),
                           ~ s(log_t0)),
                      data=pike_final,family=gammals)
## I don't know why this is not working, but F it

## simulate data from fitted models
n_sim<-10
## four sets of simulated data!
GAUsim_mean<-GAUsim_sd<-GAUsim_skew<-GAUsim_kurt<-matrix(NA,nrow=nrow(pike_final),ncol=n_sim)
SHASHsim_mean<-SHASHsim_sd<-SHASHsim_skew<-SHASHsim_kurt<-matrix(NA,nrow=nrow(pike_final),ncol=n_sim)

for(i in 1:n_sim){
  cat("################# Simulation number ", i, "\n") 
  ## add this iteration of sim data to real df
  pike_final$log_t1.sim.GAU <- rnorm(n=nrow(pike_final),
                                      mean=pike_gau_pred[,1],
                                      sd=1/pike_gau_pred[,2])
  pike_final$log_t1.sim.SHASH <- rSHASHo2(n=nrow(pike_final),
                                             mu=pike_shash_pred[,1],
                                             sigma=exp(pike_shash_pred[,2]),
                                             nu=pike_shash_pred[,3],
                                             tau=exp(pike_shash_pred[,4]))
  ## Qreg on sim data
  q.05.GAU<-predict(qgam(log_t1.sim.GAU~s(log_t0,k=4),data=pike_final,qu=0.05)) 
  q.10.GAU<-predict(qgam(log_t1.sim.GAU~s(log_t0,k=4), data=pike_final,qu=0.10)) 
  q.25.GAU<-predict(qgam(log_t1.sim.GAU~s(log_t0,k=4), data=pike_final,qu=0.25))
  q.50.GAU<-predict(qgam(log_t1.sim.GAU~s(log_t0,k=4), data=pike_final,qu=0.5))
  q.75.GAU<-predict(qgam(log_t1.sim.GAU~s(log_t0,k=4), data=pike_final,qu=0.75)) 
  q.90.GAU<-predict(qgam(log_t1.sim.GAU~s(log_t0,k=4), data=pike_final,qu=0.90))
  q.95.GAU<-predict(qgam(log_t1.sim.GAU~s(log_t0,k=4), data=pike_final,qu=0.95))
  GAUsim_mean[,i]<-Q.mean(q.25.GAU,q.50.GAU,q.75.GAU)
  GAUsim_sd[,i]<-Q.sd(q.25.GAU,q.75.GAU)
  GAUsim_skew[,i]<-Q.skewness(q.10.GAU,q.50.GAU,q.90.GAU)
  GAUsim_kurt[,i]<-Q.kurtosis(q.05.GAU,q.25.GAU,q.75.GAU,q.95.GAU)
  
  q.05.SHASH<-predict(qgam(log_t1.sim.SHASH~s(log_t0,k=4),data=pike_final,qu=0.05)) 
  q.10.SHASH<-predict(qgam(log_t1.sim.SHASH~s(log_t0,k=4), data=pike_final,qu=0.10)) 
  q.25.SHASH<-predict(qgam(log_t1.sim.SHASH~s(log_t0,k=4), data=pike_final,qu=0.25))
  q.50.SHASH<-predict(qgam(log_t1.sim.SHASH~s(log_t0,k=4), data=pike_final,qu=0.5))
  q.75.SHASH<-predict(qgam(log_t1.sim.SHASH~s(log_t0,k=4), data=pike_final,qu=0.75)) 
  q.90.SHASH<-predict(qgam(log_t1.sim.SHASH~s(log_t0,k=4), data=pike_final,qu=0.90))
  q.95.SHASH<-predict(qgam(log_t1.sim.SHASH~s(log_t0,k=4), data=pike_final,qu=0.95))
  SHASHsim_mean[,i]<-Q.mean(q.25.SHASH,q.50.SHASH,q.75.SHASH)
  SHASHsim_sd[,i]<-Q.sd(q.25.SHASH,q.75.SHASH)
  SHASHsim_skew[,i]<-Q.skewness(q.10.SHASH,q.50.SHASH,q.90.SHASH)
  SHASHsim_kurt[,i]<-Q.kurtosis(q.05.SHASH,q.25.SHASH,q.75.SHASH,q.95.SHASH)
}

## and now the real data
q.05<-predict(qgam(log_t1~s(log_t0,k=4),data=pike_final,qu=0.05)) 
q.10<-predict(qgam(log_t1~s(log_t0,k=4), data=pike_final,qu=0.10)) 
q.25<-predict(qgam(log_t1~s(log_t0,k=4), data=pike_final,qu=0.25)) 
q.50<-predict(qgam(log_t1~s(log_t0,k=4), data=pike_final,qu=0.5))
q.75<-predict(qgam(log_t1~s(log_t0,k=4), data=pike_final,qu=0.75))
q.90<-predict(qgam(log_t1~s(log_t0,k=4), data=pike_final,qu=0.90))
q.95<-predict(qgam(log_t1~s(log_t0,k=4), data=pike_final,qu=0.95))

par(mfrow=c(2,2),mar=c(4,4,1,1))
plot(pike_final$log_t0,Q.mean(q.25,q.50,q.75),type="n",
     xlab="size t",ylab="mean size t1",ylim=c(min(GAUsim_mean),max(GAUsim_mean)))
for(i in 1:n_sim){
  lines(pike_final$log_t0,GAUsim_mean[,i],col=alpha("tomato",0.25),pch=".")
  lines(pike_final$log_t0,SHASHsim_mean[,i],col=alpha("cornflowerblue",0.25),pch=".")
}
lines(pike_final$log_t0,Q.mean(q.25,q.50,q.75),col="black",pch=".",cex=2)
legend("topleft",legend=c("Real data","Simulated from \nfitted NO gam",
                          "Simulated from \nfitted SHASH2"),
       lty=1,col=c("black","tomato","cornflowerblue"),cex=0.8,bty="n")

plot(pike_final$log_t0,Q.sd(q.25,q.75),type="n",
     xlab="size t",ylab="sd size t1",ylim=c(min(GAUsim_sd),max(GAUsim_sd)))
for(i in 1:n_sim){
  lines(pike_final$log_t0,GAUsim_sd[,i],col=alpha("tomato",0.25),pch=".")
  lines(pike_final$log_t0,SHASHsim_sd[,i],col=alpha("cornflowerblue",0.25),pch=".")
}
lines(pike_final$log_t0,Q.sd(q.25,q.75),col="black",pch=".",cex=2)

plot(pike_final$log_t0,Q.skewness(q.10,q.50,q.90),type="n",
     xlab="size t",ylab="skewness size t1",ylim=c(min(Q.skewness(q.10,q.50,q.90)),max(Q.skewness(q.10,q.50,q.90))))
for(i in 1:n_sim){
  lines(pike_final$log_t0,GAUsim_skew[,i],col=alpha("tomato",0.25),pch=".")
  lines(pike_final$log_t0,SHASHsim_skew[,i],col=alpha("cornflowerblue",0.25),pch=".")
}
lines(pike_final$log_t0,Q.skewness(q.10,q.50,q.90),col="black",pch=".",cex=2)

plot(pike_final$log_t0,Q.kurtosis(q.05,q.25,q.75,q.95),type="n",
     xlab="size t",ylab="kurtosis size t1",ylim=c(min(Q.kurtosis(q.05,q.25,q.75,q.95)),max(Q.kurtosis(q.05,q.25,q.75,q.95))))
for(i in 1:n_sim){
  lines(pike_final$log_t0,GAUsim_kurt[,i],col=alpha("tomato",0.25),pch=".")
  lines(pike_final$log_t0,SHASHsim_kurt[,i],col=alpha("cornflowerblue",0.25),pch=".")
}
lines(pike_final$log_t0,Q.kurtosis(q.05,q.25,q.75,q.95),col="black",pch=".",cex=2)

