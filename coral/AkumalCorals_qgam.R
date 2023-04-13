#############################################################
# The corals case study, based on Bruno et al. (2011)
#
# Original code by Steve Ellner (using Tom's diagnostic plots)
#
# Last modified by Steve Ellner August 20, 2020 
#############################################################

rm(list=ls(all=TRUE))
setwd("c:/repos/IPM_size_transitions/coral"); 

require(car); require(zoo); require(moments); require(mgcv); 
require(gamlss); require(AICcmodavg); 
require(tidyverse); require(maxLik); require(qgam)
library(gamlss.dist)

source("../Diagnostics.R"); 
source("../fitChosenDists.R"); 
source("AkumalCoralsSetup.R"); # load the data frame on healthy corals 

## quartile-based estimates of mean and sd
## see https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/1471-2288-14-135
Q.mean<-function(q.25,q.50,q.75){(q.25+q.50+q.75)/3}
Q.sd<-function(q.25,q.75){(q.75-q.25)/1.35}
## Steve's functions for NP skew and kurtosis
Q.skewness<-function(q.10,q.50,q.90){(q.10 + q.90 - 2*q.50)/(q.90 - q.10)}
Q.kurtosis<-function(q.05,q.25,q.75,q.95){
  qN = qnorm(c(0.05,0.25,0.75,0.95))
  KG = (qN[4]-qN[1])/(qN[3]-qN[2])
  return(((q.95-q.05)/(q.75-q.25))/KG - 1)
}

######################################################################### 
# Steve first looked at cubic root transformation -- see AkumalCorals.R
# First step is to fit a pilot Gaussian model. 
#########################################################################
fitGAU <- gam(list(logarea.t1~s(logarea.t0),~s(logarea.t0)), data=XH, gamma=1.4,family=gaulss())
summary(fitGAU); plot(fitGAU); 

## the mean looks almost linear; is there evidence against this? 
fitGAU0 <- gam(list(logarea.t1~logarea.t0,~s(logarea.t0)), data=XH, gamma=1.4, family=gaulss())
AIC(fitGAU); AIC(fitGAU0); # yes, Delta AIC of about 9 in favor of the spline 

## the log(sigma) fit looks almost linear; is there evidence against this? 
fitGAU00 <- gam(list(logarea.t1~s(logarea.t0),~logarea.t0), data=XH, gamma=1.4,family=gaulss())
AIC(fitGAU); AIC(fitGAU00);  # Delta AIC < 2, so perhaps weak evidence 

## proceeding with fitGAU as "best" Gaussian model
fitted_sd<-1/predict(fitGAU,type="response")[,2]
XH$scaledResids=residuals(fitGAU,type="response")/fitted_sd

## quantile regressions on stand resids
S.05<-qgam(scaledResids~s(logarea.t0,k=4), data=XH,qu=0.05)
S.10<-qgam(scaledResids~s(logarea.t0,k=4), data=XH,qu=0.1)
S.25<-qgam(scaledResids~s(logarea.t0,k=4), data=XH,qu=0.25)
S.50<-qgam(scaledResids~s(logarea.t0,k=4), data=XH,qu=0.5) 
S.75<-qgam(scaledResids~s(logarea.t0,k=4), data=XH,qu=0.75)
S.90<-qgam(scaledResids~s(logarea.t0,k=4), data=XH,qu=0.9) 
S.95<-qgam(scaledResids~s(logarea.t0,k=4), data=XH,qu=0.95)

## NP skewness
NPS_hat = Q.skewness(q.10=predict(S.10),
                     q.50=predict(S.50),
                     q.90=predict(S.90))

## NP kurtosis (relative to Gaussian)
NPK_hat = Q.kurtosis(q.05=predict(S.05),
                     q.25=predict(S.25),
                     q.75=predict(S.75),
                     q.95=predict(S.95))

## view diagnostics of scaled residuals
pdf("../manuscript/figures/coral_qgam_diagnostics.pdf",height = 5, width = 11,useDingbats = F)
par(mfrow=c(1,2),mar = c(5, 5, 2, 3), oma=c(0,0,0,2)) 
plot(XH$logarea.t0,XH$logarea.t1,pch=1,col=alpha("black",0.25),cex.axis=0.8,
     xlab="log area, time t",ylab="log area, time t+1")
points(XH$logarea.t0,predict(fitGAU,type="response")[,1],col=alpha("red",0.25),pch=16,cex=.5)
par(new = TRUE)                           
plot(XH$logarea.t0,fitted_sd,col=alpha("blue",0.25),pch=16,cex=.5,
     axes = FALSE, xlab = "", ylab = "",cex.axis=0.8)
axis(side = 4, at = pretty(range(fitted_sd)),cex.axis=0.8)
mtext("std dev", side = 4, line = 2)
legend("topleft",legend=c("Fitted mean","Fitted sd"),bg="white",pch=1,col=c("red","blue"),cex=0.8)
title("A",font=3,adj=0)

plot(XH$logarea.t0,XH$scaledResids,col=alpha("black",0.25),cex.axis=0.8,
     xlab="log area, time t",ylab="Scaled residuals of size at time t+1")
points(XH$logarea.t0,predict(S.05),col="black",pch=".")
points(XH$logarea.t0,predict(S.10),col="black",pch=".")
points(XH$logarea.t0,predict(S.25),col="black",pch=".")
points(XH$logarea.t0,predict(S.50),col="black",pch=".")
points(XH$logarea.t0,predict(S.75),col="black",pch=".")
points(XH$logarea.t0,predict(S.90),col="black",pch=".")
points(XH$logarea.t0,predict(S.95),col="black",pch=".")
par(new = TRUE)                           
plot(c(XH$logarea.t0,XH$logarea.t0),c(NPS_hat,NPK_hat),
     col=c(rep(alpha("blue",0.25),nrow(XH)),rep(alpha("red",0.25),nrow(XH))),
     pch=16,cex=.5, axes = FALSE, xlab = "", ylab = "",cex.axis=0.8)
abline(h=0,col="lightgray",lty=3)
axis(side = 4, at = pretty(range(c(NPS_hat,NPK_hat))),cex.axis=0.8)
mtext("Skewness and kurtosis", side = 4, line = 2)
legend("topleft",legend=c("Quantiles","NP skewness","NP excess kurtosis"),bg="white",pch=c(20,1,1),col=c("black","red","blue"),cex=0.8)
title("B",font=3,adj=0)
dev.off()

## compare this to the "old" way
z = rollMomentsNP(XH$logarea.t0,XH$scaledResids,windows=8,smooth=TRUE,scaled=TRUE,xlab="Initial log area") 

## there is clearly some skew and kurtosis that is not accounted for in the gaussian model
## but could the gaussian model still be a reasonable approximation for size transitions?
## simulate data from Gaussian model
n_sim<-100
gau_mean<-gau_sd<-gau_skew<-gau_kurt<-matrix(NA,nrow=nrow(XH),ncol=n_sim)
for(i in 1:n_sim){
  ## add this iteration of sim data to real df
  XH$logarea.sim <- rnorm(n=nrow(XH),
                          mean=predict(fitGAU,type="response")[,1],
                          sd=(1/predict(fitGAU,type="response")[,2]))
  ## Qreg on sim data
  q.05<-predict(qgam(logarea.sim~s(logarea.t0,k=4),data=XH,qu=0.05)) 
  q.10<-predict(qgam(logarea.sim~s(logarea.t0,k=4),data=XH,qu=0.10)) 
  q.25<-predict(qgam(logarea.sim~s(logarea.t0,k=4),data=XH,qu=0.25))
  q.50<-predict(qgam(logarea.sim~s(logarea.t0,k=4),data=XH,qu=0.5))
  q.75<-predict(qgam(logarea.sim~s(logarea.t0,k=4),data=XH,qu=0.75)) 
  q.90<-predict(qgam(logarea.sim~s(logarea.t0,k=4),data=XH,qu=0.90))
  q.95<-predict(qgam(logarea.sim~s(logarea.t0,k=4),data=XH,qu=0.95))
  
  gau_mean[,i]<-Q.mean(q.25,q.50,q.75)
  gau_sd[,i]<-Q.sd(q.25,q.75)
  gau_skew[,i]<-Q.skewness(q.10,q.50,q.90)
  gau_kurt[,i]<-Q.kurtosis(q.05,q.25,q.75,q.95)
}

## and now the real data
q.05<-predict(qgam(logarea.t1~s(logarea.t0,k=4),data=XH,qu=0.05)) 
q.10<-predict(qgam(logarea.t1~s(logarea.t0,k=4),data=XH,qu=0.10)) 
q.25<-predict(qgam(logarea.t1~s(logarea.t0,k=4),data=XH,qu=0.25)) 
q.50<-predict(qgam(logarea.t1~s(logarea.t0,k=4),data=XH,qu=0.5))
q.75<-predict(qgam(logarea.t1~s(logarea.t0,k=4),data=XH,qu=0.75))
q.90<-predict(qgam(logarea.t1~s(logarea.t0,k=4),data=XH,qu=0.90))
q.95<-predict(qgam(logarea.t1~s(logarea.t0,k=4),data=XH,qu=0.95))

par(mfrow=c(2,2),mar=c(4,4,1,1))
plot(XH$logarea.t0,Q.mean(q.25,q.50,q.75),type="n",
     xlab="size t",ylab="mean size t1",ylim=c(min(sim_mean),max(sim_mean)))
for(i in 1:n_sim){
  points(XH$logarea.t0,gau_mean[,i],col=alpha("black",0.25),pch=".")
}
points(XH$logarea.t0,Q.mean(q.25,q.50,q.75),col="red",pch=".",cex=2)
legend("topleft",legend=c("Real data","Simulated from \nbest Gaussian"),
       lty=1,col=c("red","black"),lwd=c(2,1),cex=0.8,bty="n")

plot(XH$logarea.t0,Q.sd(q.25,q.75),type="n",
     xlab="size t",ylab="sd size t1",ylim=c(min(sim_sd),max(sim_sd)))
for(i in 1:n_sim){
  points(XH$logarea.t0,gau_sd[,i],col=alpha("black",0.25),pch=".")
}
points(XH$logarea.t0,Q.sd(q.25,q.75),col="red",pch=".",cex=2)

plot(XH$logarea.t0,Q.skewness(q.10,q.50,q.90),type="n",
     xlab="size t",ylab="skewness size t1",ylim=c(min(sim_skew),max(sim_skew)))
for(i in 1:n_sim){
  points(XH$logarea.t0,gau_skew[,i],col=alpha("black",0.25),pch=".")
}
points(XH$logarea.t0,Q.skewness(q.10,q.50,q.90),col="red",pch=".",cex=2)

plot(XH$logarea.t0,Q.kurtosis(q.05,q.25,q.75,q.95),type="n",
     xlab="size t",ylab="kurtosis size t1",ylim=c(min(sim_kurt),max(sim_kurt)))
for(i in 1:n_sim){
  points(XH$logarea.t0,gau_kurt[,i],col=alpha("black",0.25),pch=".")
}
points(XH$logarea.t0,Q.kurtosis(q.05,q.25,q.75,q.95),col="red",pch=".",cex=2)

## improved model: gam SHASH
fitSHASH <- gam(list(logarea.t1 ~ s(logarea.t0,k=4), # <- location 
                           ~ s(logarea.t0,k=4),   # <- log-scale
                           ~ s(logarea.t0,k=4),   # <- skewness
                           ~ s(logarea.t0,k=4)), # <- log-kurtosis
                      data = XH, 
                      family = shash,  
                      optimizer = "efs")
SHASH_pred<-predict(fitSHASH,type="response")

## simulate from fitted model and compare to Gaussian model
n_sim<-100
shash_mean<-shash_sd<-shash_skew<-shash_kurt<-matrix(NA,nrow=nrow(XH),ncol=n_sim)
for(i in 1:n_sim){
  ## add this iteration of sim data to real df
  XH$logarea.sim <- rSHASHo2(n=nrow(XH),
                                      mu=SHASH_pred[,1],
                                      sigma=exp(SHASH_pred[,2]),
                                      nu=SHASH_pred[,3],
                                      tau=exp(SHASH_pred[,4]))
  ## Qreg on sim data
  q.05<-predict(qgam(logarea.sim~s(logarea.t0,k=4),data=XH,qu=0.05)) 
  q.10<-predict(qgam(logarea.sim~s(logarea.t0,k=4),data=XH,qu=0.10)) 
  q.25<-predict(qgam(logarea.sim~s(logarea.t0,k=4),data=XH,qu=0.25))
  q.50<-predict(qgam(logarea.sim~s(logarea.t0,k=4),data=XH,qu=0.5))
  q.75<-predict(qgam(logarea.sim~s(logarea.t0,k=4),data=XH,qu=0.75)) 
  q.90<-predict(qgam(logarea.sim~s(logarea.t0,k=4),data=XH,qu=0.90))
  q.95<-predict(qgam(logarea.sim~s(logarea.t0,k=4),data=XH,qu=0.95))
  
  shash_mean[,i]<-Q.mean(q.25,q.50,q.75)
  shash_sd[,i]<-Q.sd(q.25,q.75)
  shash_skew[,i]<-Q.skewness(q.10,q.50,q.90)
  shash_kurt[,i]<-Q.kurtosis(q.05,q.25,q.75,q.95)
}

## and now the real data
q.05<-predict(qgam(logarea.t1~s(logarea.t0,k=4),data=XH,qu=0.05)) 
q.10<-predict(qgam(logarea.t1~s(logarea.t0,k=4),data=XH,qu=0.10)) 
q.25<-predict(qgam(logarea.t1~s(logarea.t0,k=4),data=XH,qu=0.25)) 
q.50<-predict(qgam(logarea.t1~s(logarea.t0,k=4),data=XH,qu=0.5))
q.75<-predict(qgam(logarea.t1~s(logarea.t0,k=4),data=XH,qu=0.75))
q.90<-predict(qgam(logarea.t1~s(logarea.t0,k=4),data=XH,qu=0.90))
q.95<-predict(qgam(logarea.t1~s(logarea.t0,k=4),data=XH,qu=0.95))

## combo GAU and SHASH figure
alpha_scale<-0.15
pdf("../manuscript/figures/coral_SHASH_fit.pdf",height = 6, width = 6,useDingbats = F)
par(mfrow=c(2,2),mar=c(4,4,1,1))
plot(XH$logarea.t0,Q.mean(q.25,q.50,q.75),type="n",
     xlab="size t",ylab="mean size t1",ylim=c(min(sim_mean),max(sim_mean)))
title(main="A",adj=0,font=3)
for(i in 1:n_sim){
  points(XH$logarea.t0,gau_mean[,i],col=alpha("tomato",alpha_scale),pch=".")
  points(XH$logarea.t0,shash_mean[,i],col=alpha("cornflowerblue",alpha_scale),pch=".")
}
points(XH$logarea.t0,Q.mean(q.25,q.50,q.75),col="black",pch=".",cex=2)
legend("topleft",legend=c("Real data",
                          "Simulated (Gaussian)",
                          "Simulated (SHASH)"),
       lty=1,col=c("black","tomato","cornflowerblue"),lwd=1,cex=0.8,bty="n")

plot(XH$logarea.t0,Q.sd(q.25,q.75),type="n",
     xlab="size t",ylab="sd size t1",ylim=c(min(sim_sd),max(sim_sd)))
title(main="B",adj=0,font=3)
for(i in 1:n_sim){
  points(XH$logarea.t0,gau_sd[,i],col=alpha("tomato",alpha_scale),pch=".")
  points(XH$logarea.t0,shash_sd[,i],col=alpha("cornflowerblue",alpha_scale),pch=".")
}
points(XH$logarea.t0,Q.sd(q.25,q.75),col="black",pch=".",cex=2)

plot(XH$logarea.t0,Q.skewness(q.10,q.50,q.90),type="n",
     xlab="size t",ylab="skewness size t1",ylim=c(min(sim_skew),max(sim_skew)))
title(main="C",adj=0,font=3)
for(i in 1:n_sim){
  points(XH$logarea.t0,gau_skew[,i],col=alpha("tomato",alpha_scale),pch=".")
  points(XH$logarea.t0,shash_skew[,i],col=alpha("cornflowerblue",alpha_scale),pch=".")
}
points(XH$logarea.t0,Q.skewness(q.10,q.50,q.90),col="black",pch=".",cex=2)

plot(XH$logarea.t0,Q.kurtosis(q.05,q.25,q.75,q.95),type="n",
     xlab="size t",ylab="kurtosis size t1",ylim=c(min(sim_kurt),max(sim_kurt)))
title(main="D",adj=0,font=3)
for(i in 1:n_sim){
  points(XH$logarea.t0,gau_kurt[,i],col=alpha("tomato",alpha_scale),pch=".")
  points(XH$logarea.t0,shash_kurt[,i],col=alpha("cornflowerblue",alpha_scale),pch=".")
}
points(XH$logarea.t0,Q.kurtosis(q.05,q.25,q.75,q.95),col="black",pch=".",cex=2)
dev.off()

## if one cared to know the AIC difference:
AIC(fitGAU);AIC(fitSHASH)##shash is clear winner
