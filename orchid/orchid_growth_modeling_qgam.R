## load libraries
library(tidyverse)
library(lme4)
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

## read in demographic data provided by Hans Jacquemyn (basis for Miller et al. 2012 PRSB)
orchid<-read_csv("orchid/Orchis_IPM_data.csv") %>% 
  ## there are two sites/populations, one in light and one in shade.
  ## for the purposes of this analysis, just take the light population
  filter(light=="L") %>% 
  mutate(log_area_t=log(total.leaf.area),
         log_area_t1=log(end.total.leaf.area)) 

orchid %>% 
  dplyr::select(log_area_t,log_area_t1,flowering,begin.year) %>% 
  drop_na() -> orchid_grow

## pilot Gaussian fit: different means for flowering and non-flowering, year rfx
orchid_m1<-lmer(log_area_t1~log_area_t*as.logical(flowering)+(1|begin.year),data=orchid_grow)
orchid_grow$fitted<-fitted(orchid_m1)
orchid_grow$resid<-residuals(orchid_m1)

plot(orchid_grow$fitted,orchid_grow$resid,col=alpha(orchid_grow$flowering+1,0.5))

## fit sd as a function of fitted value-- I know from previous work they have different variances
sdloglik = function(pars,resid,fitted) {
  dnorm(resid, mean=0, sd=pars[1]*exp(pars[2]*orchid_grow$fitted),log=TRUE)
}	
sdfit_flow0=maxLik(logLik=sdloglik,
                   resid=orchid_grow$resid[orchid_grow$flowering==0],
                   fitted=orchid_grow$fitted[orchid_grow$flowering==0],
                   start=c(sd(orchid_grow$resid[orchid_grow$flowering==0]),0)) 
sdfit_flow1=maxLik(logLik=sdloglik,
                   resid=orchid_grow$resid[orchid_grow$flowering==1],
                   fitted=orchid_grow$fitted[orchid_grow$flowering==1],
                   start=c(sd(orchid_grow$resid[orchid_grow$flowering==1]),0)) 

plot(orchid_grow$fitted[orchid_grow$flowering==0],
     sdfit_flow0$estimate[1]*exp(sdfit_flow0$estimate[2]*orchid_grow$fitted[orchid_grow$flowering==0]),
     xlim=range(orchid_grow$fitted),ylim=c(0,0.6),col=1)
points(orchid_grow$fitted[orchid_grow$flowering==1],
       sdfit_flow1$estimate[1]*exp(sdfit_flow1$estimate[2]*orchid_grow$fitted[orchid_grow$flowering==1]),col=2)

orchid_grow$stdresid <- NA
orchid_grow$stdresid[orchid_grow$flowering==0]<-orchid_grow$resid[orchid_grow$flowering==0]/(sdfit_flow0$estimate[1]*exp(sdfit_flow0$estimate[2]*orchid_grow$fitted[orchid_grow$flowering==0]))
orchid_grow$stdresid[orchid_grow$flowering==1]<-orchid_grow$resid[orchid_grow$flowering==1]/(sdfit_flow1$estimate[1]*exp(sdfit_flow1$estimate[2]*orchid_grow$fitted[orchid_grow$flowering==1]))

## fit qgams by flowering status
q.05<-predict(qgam(stdresid~s(fitted,k=4,by=flowering),data=orchid_grow,qu=0.05,argGam=list(gamma=2))) 
q.10<-predict(qgam(stdresid~s(fitted,k=4,by=flowering),data=orchid_grow,qu=0.10,argGam=list(gamma=2))) 
q.25<-predict(qgam(stdresid~s(fitted,k=4,by=flowering),data=orchid_grow,qu=0.25,argGam=list(gamma=2))) 
q.50<-predict(qgam(stdresid~s(fitted,k=4,by=flowering),data=orchid_grow,qu=0.50,argGam=list(gamma=2))) 
q.75<-predict(qgam(stdresid~s(fitted,k=4,by=flowering),data=orchid_grow,qu=0.75,argGam=list(gamma=2))) 
q.90<-predict(qgam(stdresid~s(fitted,k=4,by=flowering),data=orchid_grow,qu=0.90,argGam=list(gamma=2))) 
q.95<-predict(qgam(stdresid~s(fitted,k=4,by=flowering),data=orchid_grow,qu=0.95,argGam=list(gamma=2))) 

## replot std resids w quantiles
par(mfrow=c(2,1),mar = c(5, 4, 2, 3), oma=c(0,0,0,4))
plot(orchid_grow$fitted[orchid_grow$flowering==0],orchid_grow$stdresid[orchid_grow$flowering==0],
     col=alpha(1,0.5),xlim=range(orchid_grow$fitted),xlab="Fitted value",ylab="Standardized residual")
points(orchid_grow$fitted[orchid_grow$flowering==0],q.05[orchid_grow$flowering==0],col="black",pch=".")
points(orchid_grow$fitted[orchid_grow$flowering==0],q.10[orchid_grow$flowering==0],col="black",pch=".")
points(orchid_grow$fitted[orchid_grow$flowering==0],q.25[orchid_grow$flowering==0],col="black",pch=".")
points(orchid_grow$fitted[orchid_grow$flowering==0],q.50[orchid_grow$flowering==0],col="black",pch=".")
points(orchid_grow$fitted[orchid_grow$flowering==0],q.75[orchid_grow$flowering==0],col="black",pch=".")
points(orchid_grow$fitted[orchid_grow$flowering==0],q.90[orchid_grow$flowering==0],col="black",pch=".")
points(orchid_grow$fitted[orchid_grow$flowering==0],q.95[orchid_grow$flowering==0],col="black",pch=".")
par(new = TRUE)                           
plot(c(orchid_grow$fitted[orchid_grow$flowering==0],orchid_grow$fitted[orchid_grow$flowering==0]),
     c(Q.skewness(q.10[orchid_grow$flowering==0],q.50[orchid_grow$flowering==0],q.90[orchid_grow$flowering==0]),
       Q.kurtosis(q.05[orchid_grow$flowering==0],q.25[orchid_grow$flowering==0],q.75[orchid_grow$flowering==0],q.95[orchid_grow$flowering==0])),
     col=c(rep(alpha("blue",0.25),length(orchid_grow$fitted[orchid_grow$flowering==0])),
           rep(alpha("red",0.25),length(orchid_grow$fitted[orchid_grow$flowering==0]))),
     pch=16,cex=.5,axes = FALSE, xlab = "", ylab = "")
abline(h=0,col="lightgray",lty=3)
axis(side = 4,cex.axis=0.8,at = pretty(range(c(Q.skewness(q.10,q.50,q.90),Q.kurtosis(q.05,q.25,q.75,q.95)))))
mtext("Skewness", side = 4, line = 2,col="blue")
mtext("Excess Kurtosis", side = 4, line =3,col="red")

plot(orchid_grow$fitted[orchid_grow$flowering==1],orchid_grow$stdresid[orchid_grow$flowering==1],
       xlim=range(orchid_grow$fitted),xlab="Fitted value",ylab="Standardized residual",col=alpha(2,0.5))
points(orchid_grow$fitted[orchid_grow$flowering==1],q.05[orchid_grow$flowering==1],col=2,pch=".")
points(orchid_grow$fitted[orchid_grow$flowering==1],q.10[orchid_grow$flowering==1],col=2,pch=".")
points(orchid_grow$fitted[orchid_grow$flowering==1],q.25[orchid_grow$flowering==1],col=2,pch=".")
points(orchid_grow$fitted[orchid_grow$flowering==1],q.50[orchid_grow$flowering==1],col=2,pch=".")
points(orchid_grow$fitted[orchid_grow$flowering==1],q.75[orchid_grow$flowering==1],col=2,pch=".")
points(orchid_grow$fitted[orchid_grow$flowering==1],q.90[orchid_grow$flowering==1],col=2,pch=".")
points(orchid_grow$fitted[orchid_grow$flowering==1],q.95[orchid_grow$flowering==1],col=2,pch=".")
par(new = TRUE)                           
plot(c(orchid_grow$fitted[orchid_grow$flowering==1],orchid_grow$fitted[orchid_grow$flowering==1]),
     c(Q.skewness(q.10[orchid_grow$flowering==1],q.50[orchid_grow$flowering==1],q.90[orchid_grow$flowering==1]),
       Q.kurtosis(q.05[orchid_grow$flowering==1],q.25[orchid_grow$flowering==1],q.75[orchid_grow$flowering==1],q.95[orchid_grow$flowering==1])),
     col=c(rep(alpha("blue",0.25),length(orchid_grow$fitted[orchid_grow$flowering==1])),
           rep(alpha("red",0.25),length(orchid_grow$fitted[orchid_grow$flowering==1]))),
     pch=16,cex=.5,axes = FALSE, xlab = "", ylab = "")
abline(h=0,col="lightgray",lty=3)
axis(side = 4,cex.axis=0.8,at = pretty(range(c(Q.skewness(q.10,q.50,q.90),Q.kurtosis(q.05,q.25,q.75,q.95)))))
mtext("Skewness", side = 4, line = 2,col="blue")
mtext("Excess Kurtosis", side = 4, line =3,col="red")
