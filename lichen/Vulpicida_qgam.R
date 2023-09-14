################################################################
# Lichen case study, data from Peterson (de Marche) et al. MEE 
################################################################

rm(list=ls(all=TRUE))

### move to the right local directory 
tom = "C:/Users/tm9/Dropbox/github/IPM_size_transitions"
steve = "c:/repos/IPM_size_transitions" 
home = ifelse(Sys.info()["user"] == "Ellner", steve, tom)
setwd(home); setwd("lichen"); 

require(car); require(zoo); require(moments); require(mgcv); 
require(gamlss); require(AICcmodavg); 
require(tidyverse); require(maxLik); require(qgam)
library(gamlss.dist)

source("../code/Diagnostics.R"); 
source("../code/fitChosenDists.R"); 
source("../code/variance_diagnostics.R"); 

PLOTTING = TRUE; 

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


XH = read.csv("Vulpicida raw data.csv"); 
XH = XH[XH$survival==1,]; 
e = order(XH$t0); XH = XH[e,]; 

fitGAU <- gam(list(t1~s(t0),~s(t0)), data=XH, gamma=1.4,family=gaulss())
summary(fitGAU); plot(fitGAU); 

## the mean looks almost linear; is there evidence against this? 
fitGAU0 <- gam(list(t1~t0,~s(t0)), data=XH, gamma=1.4, family=gaulss())
AIC(fitGAU); AIC(fitGAU0); # Somewhat, Delta AIC of about 4 in favor of the spline 

## the log(sigma) fit looks almost linear; is there evidence against this? 
fitGAU00 <- gam(list(t1~s(t0),~t0), data=XH, gamma=1.4,family=gaulss())
AIC(fitGAU); AIC(fitGAU00);  # Delta AIC about 37, so very strong evidence

## what about quadratic? 
fitGAU01 <- gam(list(t1~s(t0,k=12), ~t0 + I(t0^2)), data=XH, gamma=1.4,family=gaulss())
AIC(fitGAU); AIC(fitGAU01);  # Delta AIC about 1

## what about quadratic for the mean, as well? 
fitGAU22 <- gam(list(t1~t0 + I(t0^2), ~t0 + I(t0^2)), data=XH, gamma=1.4,family=gaulss())
AIC(fitGAU); AIC(fitGAU22);  # the quadratic-quadratic wins by a hair, Delta AIC \approx 1 .  

## proceeding with fitGAU22 as "best" Gaussian model
XH$fitted_sd <- 1/predict(fitGAU22,type="response")[,2]
XH$fitted = predict(fitGAU22,type="response")[,1]
XH$scaledResids=residuals(fitGAU22,type="pearson")

mean(XH$scaledResids); sd(XH$scaledResids); ## all good 

### Parametric SD function is nearly the same as the spline 
plot(XH$t0, log(XH$fitted_sd)); 
points(XH$t0, log( 1/predict(fitGAU,type="response")[,2]), col="red" )

####################### Save the best-fitting Gaussian 
fitGAU = fitGAU22; rm(fitGAU22); rm(fitGAU0); rm(fitGAU00);  

################# diagnostics on fitted std. dev. function 
stopCluster(c1); 
c1<- makeCluster(8); 
registerDoParallel(c1);
out = multiple_levene_test(XH$fitted, XH$scaledResids, 3, 8, 2000);
out$p_value; ## > 0.89

out = multiple_bs_test(XH$fitted, XH$scaledResids, 4, 8, 2000) 
out$p_value; ## > 0.78; 
stopCluster(c1); 

## quantile regressions on stand resids
S.05<-qgam(scaledResids~s(t0,k=6), data=XH,qu=0.05)
S.10<-qgam(scaledResids~s(t0,k=6), data=XH,qu=0.1)
S.25<-qgam(scaledResids~s(t0,k=6), data=XH,qu=0.25)
S.50<-qgam(scaledResids~s(t0,k=6), data=XH,qu=0.5) 
S.75<-qgam(scaledResids~s(t0,k=6), data=XH,qu=0.75)
S.90<-qgam(scaledResids~s(t0,k=6), data=XH,qu=0.9) 
S.95<-qgam(scaledResids~s(t0,k=6), data=XH,qu=0.95)

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
pdf("lichen_qgam_diagnostics.pdf",height = 5, width = 11,useDingbats = F)
par(mfrow=c(1,2),mar = c(5, 5, 2, 3), oma=c(0,0,0,2)) 
plot(XH$t0,XH$t1,pch=1,col=alpha("black",0.25),cex.axis=0.8,
     xlab="log area, time t",ylab="log area, time t+1")
points(XH$t0,predict(fitGAU,type="response")[,1],col=alpha("red",0.25),pch=16,cex=.5)
par(new = TRUE)                           
plot(XH$t0,XH$fitted_sd,col=alpha("blue",0.25),pch=16,cex=.5,
     axes = FALSE, xlab = "", ylab = "",cex.axis=0.8)
axis(side = 4, at = pretty(range(XH$fitted_sd)),cex.axis=0.8)
mtext("std dev", side = 4, line = 2)
legend("topleft",legend=c("Fitted mean","Fitted sd"),bg="white",pch=1,col=c("red","blue"),cex=0.8)
title("A",font=3,adj=0)

plot(XH$t0,XH$scaledResids,col=alpha("black",0.25),cex.axis=0.8,
     xlab="log area, time t",ylab="Scaled residuals of size at time t+1")
points(XH$t0,predict(S.05),col="black",pch=".")
points(XH$t0,predict(S.10),col="black",pch=".")
points(XH$t0,predict(S.25),col="black",pch=".")
points(XH$t0,predict(S.50),col="black",pch=".")
points(XH$t0,predict(S.75),col="black",pch=".")
points(XH$t0,predict(S.90),col="black",pch=".")
points(XH$t0,predict(S.95),col="black",pch=".")
par(new = TRUE)                           
matplot(cbind(XH$t0,XH$t0),cbind(NPS_hat,NPK_hat), type="l",
     col=c("blue","red"), lty=1, axes = FALSE, xlab = "", ylab = "",cex.axis=0.8)
abline(h=0,col="lightgray",lty=3)
axis(side = 4, at = pretty(range(c(NPS_hat,NPK_hat))),cex.axis=0.8)
mtext("Skewness and kurtosis", side = 4, line = 2)
legend("topleft",legend=c("Quantiles","NP skewness","NP excess kurtosis"),bg="white",pch=c(20,1,1),col=c("black","blue","red"),cex=0.8)
title("B",font=3,adj=0)
dev.off()

## improved model: gam SHASH
fitSHASH <- gam(list(t1 ~ s(t0), # <- location 
                        ~ s(t0),  # <- log-scale
                        ~ s(t0),  # <- skewness
                        ~ 1), # <- log-kurtosis
                      data = XH, gamma=1.4, family = shash,  optimizer = "efs")
SHASH_pred<-predict(fitSHASH,type="response")
plot(fitSHASH,scale=FALSE); 

AIC(fitGAU,fitSHASH); 

fitSHASH2 <- gam(list(t1 ~ t0 + I(t0^2), # <- location 
                        ~ s(t0),  # <- log-scale
                        ~ s(t0),  # <- skewness
                        ~ 1), # <- log-kurtosis
                      data = XH, gamma=1.4, family = shash,  optimizer = "efs")
SHASH_pred<-predict(fitSHASH,type="response")
plot(fitSHASH2,scale=FALSE); 

