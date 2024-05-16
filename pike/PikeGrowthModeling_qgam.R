## Purpose: analyze pike growth data and diagnose appropriate growth kernel
## Here using gam and qgam analysis of growth data and standardized residuals
## Author: Tom Miller (mash-up with code by Steve Ellner)


### move to the right local directory 
tom = "C:/Users/tm9/Dropbox/github/IPM_size_transitions"
steve = "c:/repos/IPM_size_transitions" 
home = ifelse(Sys.info()["user"] == "Ellner", steve, tom)

setwd(home);

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
library(exactLTRE)

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
         growth = t1-t0,
         log_growth = log_t1 - log_t0) %>% 
  arrange(log_t0,log_t1)-> pike_final


plot(pike_final$t0,pike_final$t1)
plot(pike_final$log_t0,pike_final$log_t1)
abline(0,1)

## pilot gaussian with mgcv
## I am using k=6 in the final SHASH model (see below) so keeping the same for pilot Gaussian
pike_gau<-gam(list(log_t1 ~ s(log_t0,k=5), ~s(log_t0,k=5)), data=pike_final, family=gaulss())
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
## not too bad
  

##############################################################
# Steve interposes: test for constant variance 
##############################################################
source("code/variance_diagnostics.R"); 

stopCluster(c1); 
c1<- makeCluster(3); 
registerDoParallel(c1);
R = 2000; 

out_bartlett = multiple_bartlett_test(pike_gau_pred,pike_final$scaledResids, 3, 10, R) 
out_bartlett$p_value; #0 

out_bs = multiple_bs_test(pike_gau_pred,pike_final$scaledResids, 3, 10, R) 
out_bs$p_value; #0 

require(mgcv); 
fit = gam(I(pike_final$scaledResids^2)~s(pike_gau_pred))  ## p<0.001, under 1\% of variance
fit = gam(abs(pike_final$scaledResids)~s(pike_gau_pred))  ## p<0.001, under 1\% of variance

fit = gam(I(scaledResids^2)~s(log_t0),data=pike_final)  ## p<0.001, under 1\% of variance
fit = gam(abs(scaledResids)~s(log_t0),data = pike_final)  ## p<0.001, under 2\% of variance
stopCluster(c1); 
 
##############################################################
# Back to Tom 
##############################################################
 
 
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

pdf("manuscript/figures/pike_resid_diagnostics.pdf",height = 5, width = 11,useDingbats = F)
## mean and sd of size transition
par(mfrow=c(1,2),mar = c(5, 5, 2, 3), oma=c(0,0,0,2)) 
plot(pike_final$log_t0,pike_final$log_t1,pch=".",col=alpha("black",0.25),
     xlab="log(length), time t",ylab="log(length), time t+1")
lines(pike_final$log_t0,pike_gau_pred[,1],col="red",lwd=2)
abline(0,1)
par(new = TRUE)                           
plot(pike_final$log_t0,1/pike_gau_pred[,2],col="blue",pch=16,cex=.5,
     axes = FALSE, xlab = "", ylab = "",type="l",lwd=2)
axis(side = 4, at = pretty(range(1/pike_gau_pred[,2])))
mtext("Standard deviation", side = 4, line = 2)
abline(0,1)
legend("bottom",legend=c("Fitted mean","Fitted sd"),bg="white",lwd=2,col=c("red","blue"),cex=0.8)
title("A",font=3,adj=0)

## scaled residuals plot
plot(pike_final$log_t0,pike_final$scaledResids,col=alpha("black",0.25),pch=".",
     xlab="log(length), time t",ylab="Scaled residuals of log(length) at t+1")
lines(pike_final$log_t0,q.05,col="black")
lines(pike_final$log_t0,q.10,col="black")
lines(pike_final$log_t0,q.25,col="black")
lines(pike_final$log_t0,q.50,col="black")
lines(pike_final$log_t0,q.75,col="black")
lines(pike_final$log_t0,q.90,col="black")
lines(pike_final$log_t0,q.95,col="black")
par(new = TRUE)                           
matplot(cbind(pike_final$log_t0,pike_final$log_t0),
        cbind(NPS_hat,NPK_hat),
        type="l",lwd=2,
        col=c("blue3","red3"), lty=1, axes = FALSE, xlab = "", ylab = "",cex.axis=0.8)
axis(side = 4, at = pretty(range(c(NPS_hat,NPK_hat))),cex.axis=0.8)
mtext("NP skewness or kurtosis", side = 4, line = 2)
legend("topright",legend=c("Skewness","Excess kurtosis"),bg="white",lwd=1.5,col=c("blue3","red3"),cex=0.8)
title("B",font=3,adj=0)
dev.off()

##do we ever get a non-zero probability of shrinkage?
size_dummy <- seq(min(pike_final$log_t0),max(pike_final$log_t1),0.01)
plot(size_dummy,dnorm(size_dummy,
      mean=pike_gau_pred[1,1],
      sd=1/pike_gau_pred[1,2]),type="l")
abline(v=pike_final$log_t0[1],col="red")
lines(size_dummy,dnorm(size_dummy,
           mean=pike_gau_pred[nrow(pike_final),1],
           sd=1/pike_gau_pred[nrow(pike_final),2]),type="l")
abline(v=pike_final$log_t0[nrow(pike_final)],col="blue")

## gonna try a gam shash (bc that's what's easy in mgcv)
## k=6 gave better correspondence between real and sim data than k=4
pike_gam_shash <- gam(list(log_t1 ~ s(log_t0,k=5), # <- model for location 
                           ~ s(log_t0,k=5),   # <- model for log-scale
                           ~ s(log_t0,k=5),   # <- model for skewness
                           ~ s(log_t0,k=5)), # <- model for log-kurtosis
                      data = pike_final, 
                      family = shash,  
                      optimizer = "efs")
pike_shash_pred <- predict(pike_gam_shash,type="response")

## also, out of curiosity, gonna fit a gamma to growth data (they can only increase)
#pike_gam_gamma <- gam(list(log_growth+0.0001 ~ s(log_t0,k=6), ~ s(t0,k=6)),
#                      data=pike_final,family=gammals)
#plot(pike_final$t0,pike_final$log_growth+0.0001)
#points(pike_final$t0,predict(pike_gam_gamma,type="response")[,1],col="red")

## simulate data from fitted models
n_sim<-10
## four sets of simulated data!
GAUsim_mean<-GAUsim_sd<-GAUsim_skew<-GAUsim_kurt<-matrix(NA,nrow=nrow(pike_final),ncol=n_sim)
SHASHsim_mean<-SHASHsim_sd<-SHASHsim_skew<-SHASHsim_kurt<-matrix(NA,nrow=nrow(pike_final),ncol=n_sim)
#GAMMAsim_mean<-GAMMAsim_sd<-GAMMAsim_skew<-GAMMAsim_kurt<-matrix(NA,nrow=nrow(pike_final),ncol=n_sim)

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
  
  #pike_final$log_t1.sim.GAMMA <- pike_final$log_t0 + 
  #  rgamma(n=nrow(pike_final),shape=1/exp(fitted(pike_gam_gamma)[,2]),scale=fitted(pike_gam_gamma)[,1]*exp(fitted(pike_gam_gamma)[,2]))
  
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

  #q.05.GAMMA<-predict(qgam(log_t1.sim.GAMMA~s(log_t0,k=4),data=pike_final,qu=0.05)) 
  #q.10.GAMMA<-predict(qgam(log_t1.sim.GAMMA~s(log_t0,k=4), data=pike_final,qu=0.10)) 
  #q.25.GAMMA<-predict(qgam(log_t1.sim.GAMMA~s(log_t0,k=4), data=pike_final,qu=0.25))
  #q.50.GAMMA<-predict(qgam(log_t1.sim.GAMMA~s(log_t0,k=4), data=pike_final,qu=0.5))
  #q.75.GAMMA<-predict(qgam(log_t1.sim.GAMMA~s(log_t0,k=4), data=pike_final,qu=0.75)) 
  #q.90.GAMMA<-predict(qgam(log_t1.sim.GAMMA~s(log_t0,k=4), data=pike_final,qu=0.90))
  #q.95.GAMMA<-predict(qgam(log_t1.sim.GAMMA~s(log_t0,k=4), data=pike_final,qu=0.95))
  #GAMMAsim_mean[,i]<-Q.mean(q.25.GAMMA,q.50.GAMMA,q.75.GAMMA)
  #GAMMAsim_sd[,i]<-Q.sd(q.25.GAMMA,q.75.GAMMA)
  #GAMMAsim_skew[,i]<-Q.skewness(q.10.GAMMA,q.50.GAMMA,q.90.GAMMA)
  #GAMMAsim_kurt[,i]<-Q.kurtosis(q.05.GAMMA,q.25.GAMMA,q.75.GAMMA,q.95.GAMMA)
}

## and now the real data
q.05<-predict(qgam(log_t1~s(log_t0,k=4),data=pike_final,qu=0.05)) 
q.10<-predict(qgam(log_t1~s(log_t0,k=4), data=pike_final,qu=0.10)) 
q.25<-predict(qgam(log_t1~s(log_t0,k=4), data=pike_final,qu=0.25)) 
q.50<-predict(qgam(log_t1~s(log_t0,k=4), data=pike_final,qu=0.5))
q.75<-predict(qgam(log_t1~s(log_t0,k=4), data=pike_final,qu=0.75))
q.90<-predict(qgam(log_t1~s(log_t0,k=4), data=pike_final,qu=0.90))
q.95<-predict(qgam(log_t1~s(log_t0,k=4), data=pike_final,qu=0.95))



##########################################################################
## Plotting: comparison of simulation results with real data
##########################################################################
par(mfrow=c(2,2),mar=c(4,4,1,1),cex.axis=1.3,cex.lab=1.3,mgp=c(2.2,1,0))

plot(pike_final$log_t0,Q.mean(q.25,q.50,q.75),type="n",
     xlab="Size(t)",ylab="NP mean size(t+1)",ylim=c(min(GAUsim_mean),max(GAUsim_mean)))
abline(0,1,lty=2,col="blue"); 
for(i in 1:n_sim){
  lines(pike_final$log_t0,GAUsim_mean[,i],col=alpha("tomato",0.5),pch=".")
  lines(pike_final$log_t0,SHASHsim_mean[,i],col=alpha("cornflowerblue",0.5),pch=".")
  #lines(pike_final$log_t0,GAMMAsim_mean[,i],col=alpha("green",0.25),pch=".")
}
lines(pike_final$log_t0,Q.mean(q.25,q.50,q.75),col="black",pch=".",cex=2)
legend("topleft",legend=c("Real data","Gaussian simulation",
                          "SHASH simulation"),
       lty=1,col=c("black","tomato","cornflowerblue"),cex=1,bty="n")

plot(pike_final$log_t0,Q.sd(q.25,q.75),type="n",
     xlab="Size(t)",ylab="NP SD size(t+1)",ylim=c(min(GAUsim_sd),max(GAUsim_sd)))
for(i in 1:n_sim){
  lines(pike_final$log_t0,GAUsim_sd[,i],col=alpha("tomato",0.5),pch=".")
  lines(pike_final$log_t0,SHASHsim_sd[,i],col=alpha("cornflowerblue",0.5),pch=".")
  #lines(pike_final$log_t0,GAMMAsim_sd[,i],col=alpha("green",0.25),pch=".")
}
lines(pike_final$log_t0,Q.sd(q.25,q.75),col="black",pch=".",cex=2)

plot(pike_final$log_t0,Q.skewness(q.10,q.50,q.90),type="n",
     xlab="Size(t)",ylab="NP skewness size(t+1)",ylim=c(min(Q.skewness(q.10,q.50,q.90)),max(Q.skewness(q.10,q.50,q.90))))
for(i in 1:n_sim){
  lines(pike_final$log_t0,GAUsim_skew[,i],col=alpha("tomato",0.5),pch=".")
  lines(pike_final$log_t0,SHASHsim_skew[,i],col=alpha("cornflowerblue",0.5),pch=".")
  #lines(pike_final$log_t0,GAMMAsim_skew[,i],col=alpha("green",0.25),pch=".")
}
lines(pike_final$log_t0,Q.skewness(q.10,q.50,q.90),col="black",pch=".",cex=2)

plot(pike_final$log_t0,Q.kurtosis(q.05,q.25,q.75,q.95),type="n",
     xlab="Size(t)",ylab="NP kurtosis size(t+1)",ylim=c(min(Q.kurtosis(q.05,q.25,q.75,q.95)),max(Q.kurtosis(q.05,q.25,q.75,q.95))))
for(i in 1:n_sim){
  lines(pike_final$log_t0,GAUsim_kurt[,i],col=alpha("tomato",0.2),pch=".")
  lines(pike_final$log_t0,SHASHsim_kurt[,i],col=alpha("cornflowerblue",0.5),pch=".")
  #lines(pike_final$log_t0,GAMMAsim_kurt[,i],col=alpha("green",0.25),pch=".")
}
lines(pike_final$log_t0,Q.kurtosis(q.05,q.25,q.75,q.95),col="black",pch=".",cex=2)
## looks like it's going to be a SHASH IPM

dev.copy2pdf(file="manuscript/figures/pike_growth_sims.pdf"); 
save.image(file="pike_growth_sims.Rdata"); 


####################################################################################################
## Make the Gaussian and SHASH IPMs 
####################################################################################################

## fit recruit size distribution (assuming fish that did not have a previous size are recruits)
pike_trans_year %>% filter(is.na(t0) & !is.na(t1)) %>% 
  mutate(log_t1 = log(t1)) %>% 
  arrange(log_t1)-> pike_recruits

pike_recruitsize_shash <- gam(list(log_t1 ~ 1, # <- model for location 
                           ~ 1,   # <- model for log-scale
                           ~ 1,   # <- model for skewness
                           ~ 1), # <- model for log-kurtosis
                      data = pike_recruits, 
                      family = shash,  
                      optimizer = "efs")
pike_recruitsize_gau <- gam(list(log_t1 ~ 1, # <- model for mean 
                                   ~ 1), # <- model for sd
                              data = pike_recruits, 
                              family = gaulss(),  
                              optimizer = "efs")
AIC(pike_recruitsize_shash); AIC(pike_recruitsize_gau)
## SHASH recruit size is better

## survival data and gam fit
pike_surv <- read_csv("pike/data/PikeSurvivalData1953_1990.csv") %>% 
  mutate(log_t0 = log(Length)) %>% 
  arrange(log_t0)

pike_surv_gam <- gam(Survival ~ s(log_t0,k=5),data = pike_surv,family = binomial)
plot(pike_surv$log_t0,pike_surv$Survival)
lines(pike_surv$log_t0,fitted(pike_surv_gam))

## fertility data and gam fit
pike_fert <- read_csv("pike/data/FecundityData1963_2002.csv") %>% 
  mutate(log_t0 = log(Length),
         log_eggs = log(Eggs)) %>% 
  arrange(log_t0)
pike_fert_gam <- gam(Eggs ~ s(log_t0,k=5),data = pike_fert,family = nb)
plot(pike_fert$log_t0,pike_fert$Eggs)
lines(pike_fert$log_t0,fitted(pike_fert_gam))


# IPM ---------------------------------------------------------------------
## GROWTH - SHASH
gxy_SHASH<-function(x,y){
  xb=pmin(pmax(x,min(pike_final$log_t0)),max(pike_final$log_t1)) 
  pred=predict(pike_gam_shash,
               newdata = data.frame(log_t0=xb))
  return(dSHASHo2(x=y, 
                  mu=pred[,1],
                  sigma = exp(pred[,2]), 
                  nu = pred[,3], 
                  tau = exp(pred[,4])))
}

## GROWTH - Gaussian
gxy_GAU<-function(x,y){
  xb=pmin(pmax(x,min(pike_final$log_t0)),max(pike_final$log_t1)) 
  pred = predict(pike_gau,newdata = data.frame(log_t0=xb))
  return(dnorm(y,mean=pred[,1],sd=exp(pred[,2])))
}

## SURVIVAL
sx<-function(x){
  xb=pmin(pmax(x,min(pike_final$log_t0)),max(pike_final$log_t1)) 
  pred = predict(pike_surv_gam,newdata = data.frame(log_t0=xb),type="response")
  return(pred)
}

## COMBINED GROWTH_SURVIVAL
pxy <- function(x,y,dist){
  result <- sx(x)*do.call(paste0("gxy_",dist),list(x,y))
  return(result)
}

## FERTILITY
## these parameter values come from Table 2
r<-0.9 #fertilization probability
q<-0.5 #fraction female
S1<-0.000623 #survival from fertilized egg to 1yo

fx<-function(x){
  xb=pmin(pmax(x,min(pike_final$log_t0)),max(pike_final$log_t1)) 
  eggs=predict(pike_fert_gam,newdata = data.frame(log_t0=xb),type="response")
  return(eggs*r*q*S1)  
}

#SIZE DISTRIBUTION OF RECRUITS
#recruit.size<-function(y){
#  dnorm(x=y,
#        mean=pike_recruitsize_gau$coefficients[1],
#        sd=exp(pike_recruitsize_gau$coefficients[2]))
#}

recruit.size<-function(y){
  dSHASHo2(y,
        pike_recruitsize_shash$coefficients[1],
        exp(pike_recruitsize_shash$coefficients[2]),
		pike_recruitsize_shash$coefficients[3],
		exp(pike_recruitsize_shash$coefficients[4]),
		)
}

##COMBINED FERTILITY/RECRUITMENT
fxy<-function(x,y){
  return(fx(x)*recruit.size(y))
}

##PUT IT ALL TOGETHER
## defaults here come from experimentation in the basement
ApproxMatrix <- function(ext.lower=0,ext.upper=0.1,mat.size=800,dist){
  # Matrix size and size extensions (upper and lower integration limits)
  n <- mat.size
  L <- min(pike_final$log_t0) - ext.lower
  U <- max(pike_final$log_t0) + ext.upper
  # Bin size for n bins
  h <- (U - L)/n
  # Lower boundaries of bins 
  b <- L + c(0:n)*h
  # Bin midpoints
  y <- 0.5*(b[1:n] + b[2:(n + 1)])
  
  # Growth/Survival matrix
  Pmat <- t(outer(y, y, pxy, dist=dist)) * h 
  # Fertility/Recruitment matrix
  Fmat <- t(outer(y, y, fxy)) * h 
  # Put it all together
  IPMmat <- Pmat + Fmat
  
  return(list(IPMmat = IPMmat, Fmat = Fmat, Pmat = Pmat, meshpts = y))
}

## check effect of matrix dimension
dims<-seq(100,1000,100)
lambda_GAU<-lambda_SHASH<-c()
for(i in 1:length(dims)){
  lambda_GAU[i]<-lambda(ApproxMatrix(dist="GAU",mat.size=dims[i])$IPMmat)
  lambda_SHASH[i]<-lambda(ApproxMatrix(dist="SHASH",mat.size=dims[i])$IPMmat)
}
plot(dims,lambda_SHASH)
points(dims,lambda_GAU,col="red")
## it looks very stable for lambda, but eviction remains a problem
## unless we go to very high dimension, because the growth kernel becomes very spiky at the upper size limit

GAU_IPM<-ApproxMatrix(dist="GAU",mat.size=1000)
SHASH_IPM<-ApproxMatrix(dist="SHASH",mat.size=1000)
plot(GAU_IPM$meshpts,colSums(GAU_IPM$Pmat) / sx(GAU_IPM$meshpts))
plot(SHASH_IPM$meshpts,colSums(SHASH_IPM$Pmat) / sx(SHASH_IPM$meshpts))
## It takes a matrix dimension of about 1000 to chop up the growth kernel appropriately

## write out for cross-spp analysis
write_rds(list(GAU_IPM=GAU_IPM,SHASH_IPM=SHASH_IPM),file="pike/pike_out.rds")

## source Steve/Chrissy's life history functions
source("code/metaluck_fns_CMH.R")
X = matrix(NA, 2,13)

### Gaussian
matU = GAU_IPM$Pmat; matF = GAU_IPM$Fmat; c0 = rep(0,nrow(matU)); c0[1]=1; 
X[1,] = c(lambda(GAU_IPM$IPMmat),
          mean_lifespan(matU, mixdist=c0),
          var_lifespan(matU, mixdist=c0)^0.5,
          skew_lifespan(matU, mixdist=c0),
          mean_LRO(matU,matF,mixdist=c0), 
          var_LRO_mcr(matU,matF,mixdist=c0)^0.5,
          skew_LRO(matU,matF,mixdist=c0), 
          prob_repro(matU,matF)[1], 
          mean_age_repro(matU,matF,mixdist=c0), 
          lifespan_reproducers(matU,matF,mixdist=c0), 
          gen_time_Ta(matU,matF), 
          gen_time_mu1_v(matU,matF), 
          gen_time_R0(matU,matF))  


### SHASH
matU = SHASH_IPM$Pmat; matF = SHASH_IPM$Fmat; c0 = rep(0,nrow(matU)); c0[1]=1; 

X[2,] = c(
  lambda(SHASH_IPM$IPMmat),
  mean_lifespan(matU, mixdist=c0),
  var_lifespan(matU, mixdist=c0)^0.5,
  skew_lifespan(matU, mixdist=c0),
  mean_LRO(matU,matF,mixdist=c0), 
  var_LRO_mcr(matU,matF,mixdist=c0)^0.5,
  skew_LRO(matU,matF,mixdist=c0), 
  prob_repro(matU,matF)[1], 
  mean_age_repro(matU,matF,mixdist=c0), 
  lifespan_reproducers(matU,matF,mixdist=c0), 
  gen_time_Ta(matU,matF), 
  gen_time_mu1_v(matU,matF), 
  gen_time_R0(matU,matF))  

X = data.frame(round(X,digits=4)); 
names(X) = c("lambda","Mean lifespan", "SD lifespan", "Skew lifespan", "Mean LRO", "SD LRO", "Skew LRO", "Prob repro", "Mean age repro", "Conditional lifespan", 
             "Gen time Ta", "Gen time mu1(v)", "Gen time R0"); 
row.names(X) = c("Gaussian", "SHASH"); 

View(X); 

# basement ----------------------------------------------------------------

## play with eviction
n <- 600
L <- min(pike_final$log_t0)
U <- max(pike_final$log_t0)+0.1
# Bin size for n bins
h <- (U - L)/n
# Lower boundaries of bins 
b <- L + c(0:n)*h
# Bin midpoints
y <- 0.5*(b[1:n] + b[2:(n + 1)])

plot(y,gxy_GAU(x=min(pike_final$log_t0),y=y))
plot(y,gxy_SHASH(x=min(pike_final$log_t0),y=y))

plot(y,gxy_GAU(x=mean(pike_final$log_t0),y=y))
plot(y,gxy_SHASH(x=mean(pike_final$log_t0),y=y))

plot(y,gxy_GAU(x=max(pike_final$log_t0),y=y))
plot(y,gxy_SHASH(x=max(pike_final$log_t0),y=y))

colSums(t(outer(y, y, gxy_GAU)) * h)
colSums(t(outer(y, y, gxy_SHASH)) * h)
