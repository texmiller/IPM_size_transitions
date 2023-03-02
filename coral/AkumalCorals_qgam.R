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
require(tidyverse); require(maxLik); 

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
q.10<-predict(S.10);q.50<-predict(S.50);q.90<-predict(S.90)
NPS_hat = (q.10 + q.90 - 2*q.50)/(q.90 - q.10)

## NP kurtosis (relative to Gaussian)
q.05<-predict(S.05);q.25<-predict(S.25);q.75<-predict(S.75);q.95<-predict(S.95)
qN = qnorm(c(0.05,0.25,0.75,0.95))
KG = (qN[4]-qN[1])/(qN[3]-qN[2])
NPK_hat = ((q.95-q.05)/(q.75-q.25))/KG - 1

## view diagnostics of scaled residuals
par(mfrow=c(1,2),mar = c(5, 4, 2, 3), oma=c(0,0,0,4)) 
plot(XH$logarea.t0,XH$logarea.t1,pch=1,col=alpha("black",0.25),
     xlab="size t",ylab="size t1")
points(XH$logarea.t0,predict(fitGAU,type="response")[,1],col=alpha("red",0.25),pch=16,cex=.5)
par(new = TRUE)                           
plot(XH$logarea.t0,fitted_sd,col=alpha("blue",0.25),pch=16,cex=.5,
     axes = FALSE, xlab = "", ylab = "")
axis(side = 4, at = pretty(range(fitted_sd)),col="blue")
mtext("sigma", side = 4, line = 2,col="blue")

plot(XH$logarea.t0,XH$scaledResids,col=alpha("black",0.25),
     xlab="Size at time t",ylab="Scaled residuals of size at t+1")
points(XH$logarea.t0,q.05,col="black",pch=".")
points(XH$logarea.t0,q.10,col="black",pch=".")
points(XH$logarea.t0,q.25,col="black",pch=".")
points(XH$logarea.t0,q.50,col="black",pch=".")
points(XH$logarea.t0,q.75,col="black",pch=".")
points(XH$logarea.t0,q.90,col="black",pch=".")
points(XH$logarea.t0,q.95,col="black",pch=".")
par(new = TRUE)                           
plot(c(XH$logarea.t0,XH$logarea.t0),c(NPS_hat,NPK_hat),
     col=c(rep(alpha("blue",0.25),nrow(XH)),rep(alpha("red",0.25),nrow(XH))),
     pch=16,cex=.5, axes = FALSE, xlab = "", ylab = "")
abline(h=0,col="lightgray",lty=3)
axis(side = 4, at = pretty(range(c(NPS_hat,NPK_hat))),cex.axis=0.8)
mtext("Skewness", side = 4, line = 2,col="blue")
mtext("Excess Kurtosis", side = 4, line =3,col="red")

## compare this to the "old" way
z = rollMomentsNP(XH$logarea.t0,XH$scaledResids,windows=8,smooth=TRUE,scaled=TRUE,xlab="Initial log area") 
