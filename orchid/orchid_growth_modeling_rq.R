## load libraries
library(tidyverse)
library(lme4)
library(scales)
library(quantreg)
library(gamlss.dist)
library(popbio)
library(moments)
library(maxLik)
library(bbmle)

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
orchid<- bind_rows(read_csv("orchid/orchis 2003-2013.csv"),
                   read_csv("orchid/orchis 2013-2015.csv")) %>% 
  ## there are two sites/populations, one in light and one in shade.
  ## for the purposes of this analysis, just take the light population
  filter(light=="L") %>% 
  mutate(log_area_t=log(total.leaf.area),
         log_area_t1=log(end.total.leaf.area)) 

orchid %>% 
  dplyr::select(log_area_t,log_area_t1,flowering,begin.year) %>% 
  drop_na() -> orchid_grow

## pilot gaussian model selection
orchid_GAU<-list()
orchid_GAU[[1]]<-lmer(log_area_t1~log_area_t +(1|begin.year),data=orchid_grow,REML=F)
orchid_GAU[[2]]<-lmer(log_area_t1~log_area_t + as.logical(flowering)+(1|begin.year),data=orchid_grow,REML=F)
orchid_GAU[[3]]<-lmer(log_area_t1~log_area_t*as.logical(flowering)+(1|begin.year),data=orchid_grow,REML=F)
AICtab(orchid_GAU)
orchid_GAU_best<-orchid_GAU[[which.min(AICctab(orchid_GAU,sort=F)$dAICc)]]
summary(orchid_GAU_best)

## fit sd as a function of initial size
orchid_grow$GAU_resids <- residuals(orchid_GAU_best)
sdloglik = function(pars) {
  dnorm(orchid_grow$GAU_resids, mean=0, sd=exp(pars[1]+pars[2]*as.logical(orchid_grow$flowering)
                               +(pars[3]+pars[4]*as.logical(orchid_grow$flowering))*orchid_grow$log_area_t),log=TRUE)
}	
GAU_sd_coef<-maxLik(logLik=sdloglik,start=c(exp(sd(orchid_grow$GAU_resids)),0,0,0))
orchid_grow$GAU_sd<-exp(GAU_sd_coef$estimate[1]+GAU_sd_coef$estimate[2]*as.logical(orchid_grow$flowering)
                        +(GAU_sd_coef$estimate[3]+GAU_sd_coef$estimate[4]*as.logical(orchid_grow$flowering))*orchid_grow$log_area_t)
orchid_grow$GAU_std_resids<-orchid_grow$GAU_resids/orchid_grow$GAU_sd
## should be mean zero unit variance
mean(orchid_grow$GAU_std_resids);sd(orchid_grow$GAU_std_resids) ##checks out

## fit quantile regression with quantreg
q.fit<-predict(rq(GAU_std_resids~log_area_t*as.logical(flowering), data=orchid_grow,tau=c(0.05,0.1,0.25,0.5,0.75,0.9,0.95)))
q.05<-q.fit[,1]
q.10<-q.fit[,2] 
q.25<-q.fit[,3] 
q.50<-q.fit[,4]
q.75<-q.fit[,5]
q.90<-q.fit[,6] 
q.95<-q.fit[,7]

par(mar = c(5, 4, 2, 3), oma=c(0,0,0,2),mfrow=c(2,2)) 
plot(orchid_grow$log_area_t,orchid_grow$log_area_t1,type="n",
     xlab="Size at t",ylab="Size at t+1")
points(orchid_grow$log_area_t[orchid_grow$flowering==0],
     orchid_grow$log_area_t1[orchid_grow$flowering==0],col=alpha("black",0.25))
points(orchid_grow$log_area_t[orchid_grow$flowering==0],
       fixef(orchid_GAU_best)[1]+fixef(orchid_GAU_best)[2]*orchid_grow$log_area_t[orchid_grow$flowering==0],
       col=alpha("red",0.25),pch=".",cex=2)
par(new = TRUE) 
plot(orchid_grow$log_area_t,orchid_grow$GAU_sd,
     type="n",cex=2,axes = FALSE, xlab = "", ylab = "")
points(orchid_grow$log_area_t[orchid_grow$flowering==0],
       orchid_grow$GAU_sd[orchid_grow$flowering==0],
       col=alpha("blue",0.25),pch=".",cex=2)
axis(side = 4,cex.axis=0.8,at = pretty(range(orchid_grow$GAU_sd)))
legend("topleft",legend=c("Mean","SD"),bg="white",pch=16,col=c("red","blue"),cex=0.8)
title("A",font=3,adj=0)

plot(orchid_grow$log_area_t,orchid_grow$log_area_t1,type="n",
     xlab="Size at t",ylab="Size at t+1")
points(orchid_grow$log_area_t[orchid_grow$flowering==1],
       orchid_grow$log_area_t1[orchid_grow$flowering==1],col=alpha("black",0.25))
points(orchid_grow$log_area_t[orchid_grow$flowering==1],
       fixef(orchid_GAU_best)[1]+fixef(orchid_GAU_best)[3]+(fixef(orchid_GAU_best)[2]+fixef(orchid_GAU_best)[4])*orchid_grow$log_area_t[orchid_grow$flowering==1],
       col=alpha("red",0.25),pch=".",cex=2)
par(new = TRUE) 
plot(orchid_grow$log_area_t,orchid_grow$GAU_sd,type="n",
     col=alpha("blue",0.25),pch=".",cex=2,axes = FALSE, xlab = "", ylab = "")
points(orchid_grow$log_area_t[orchid_grow$flowering==1],
       orchid_grow$GAU_sd[orchid_grow$flowering==1],
       col=alpha("blue",0.25),pch=".",cex=2)
axis(side = 4,cex.axis=0.8,at = pretty(range(orchid_grow$GAU_sd)))
mtext("Standard deviation", side = 4, line = 2)
title("B",font=3,adj=0)

plot(orchid_grow$log_area_t,orchid_grow$GAU_std_resids,type="n",
     xlab="Expected size at t+1",ylab="Scaled residuals of size at t+1")
points(orchid_grow$log_area_t[orchid_grow$flowering==0],
     orchid_grow$GAU_std_resids[orchid_grow$flowering==0],col=alpha("black",0.25),
     xlab="Expected size at t+1",ylab="Scaled residuals of size at t+1")
points(orchid_grow$log_area_t[orchid_grow$flowering==0],q.05[orchid_grow$flowering==0],col="black",pch=".")
points(orchid_grow$log_area_t[orchid_grow$flowering==0],q.10[orchid_grow$flowering==0],col="black",pch=".")
points(orchid_grow$log_area_t[orchid_grow$flowering==0],q.25[orchid_grow$flowering==0],col="black",pch=".")
points(orchid_grow$log_area_t[orchid_grow$flowering==0],q.50[orchid_grow$flowering==0],col="black",pch=".")
points(orchid_grow$log_area_t[orchid_grow$flowering==0],q.75[orchid_grow$flowering==0],col="black",pch=".")
points(orchid_grow$log_area_t[orchid_grow$flowering==0],q.90[orchid_grow$flowering==0],col="black",pch=".")
points(orchid_grow$log_area_t[orchid_grow$flowering==0],q.95[orchid_grow$flowering==0],col="black",pch=".")
par(new = TRUE)      
plot(c(orchid_grow$log_area_t,
       orchid_grow$log_area_t),
     c(Q.skewness(q.10,q.50,q.90),
       Q.kurtosis(q.05,q.25,q.75,q.95)),type="n",
     axes = FALSE, xlab = "", ylab = "")
points(c(orchid_grow$log_area_t[orchid_grow$flowering==0],
       orchid_grow$log_area_t[orchid_grow$flowering==0]),
     c(Q.skewness(q.10[orchid_grow$flowering==0],
                  q.50[orchid_grow$flowering==0],
                  q.90[orchid_grow$flowering==0]),
       Q.kurtosis(q.05[orchid_grow$flowering==0],
                  q.25[orchid_grow$flowering==0],
                  q.75[orchid_grow$flowering==0],
                  q.95[orchid_grow$flowering==0])),
     col=c(rep(alpha("blue",0.25),nrow(orchid_grow[orchid_grow$flowering==0,])),
           rep(alpha("red",0.25),nrow(orchid_grow[orchid_grow$flowering==0,]))),
     pch=16,cex=.5)
abline(h=0,col="lightgray",lty=3)
axis(side = 4,cex.axis=0.8,at = pretty(range(c(Q.skewness(q.10,q.50,q.90),Q.kurtosis(q.05,q.25,q.75,q.95)))))
#mtext("Skewness", side = 4, line = 2,col="blue")
#mtext("Excess Kurtosis", side = 4, line =3,col="red")
title("C",font=3,adj=0)

plot(orchid_grow$log_area_t,orchid_grow$GAU_std_resids,type="n",
     xlab="Expected size at t+1",ylab="Scaled residuals of size at t+1")
points(orchid_grow$log_area_t[orchid_grow$flowering==1],
     orchid_grow$GAU_std_resids[orchid_grow$flowering==1],col=alpha("black",0.25),
     xlab="Expected size at t+1",ylab="Scaled residuals of size at t+1")
points(orchid_grow$log_area_t[orchid_grow$flowering==1],q.05[orchid_grow$flowering==1],col="black",pch=".")
points(orchid_grow$log_area_t[orchid_grow$flowering==1],q.10[orchid_grow$flowering==1],col="black",pch=".")
points(orchid_grow$log_area_t[orchid_grow$flowering==1],q.25[orchid_grow$flowering==1],col="black",pch=".")
points(orchid_grow$log_area_t[orchid_grow$flowering==1],q.50[orchid_grow$flowering==1],col="black",pch=".")
points(orchid_grow$log_area_t[orchid_grow$flowering==1],q.75[orchid_grow$flowering==1],col="black",pch=".")
points(orchid_grow$log_area_t[orchid_grow$flowering==1],q.90[orchid_grow$flowering==1],col="black",pch=".")
points(orchid_grow$log_area_t[orchid_grow$flowering==1],q.95[orchid_grow$flowering==1],col="black",pch=".")
par(new = TRUE)      
plot(c(orchid_grow$log_area_t,
         orchid_grow$log_area_t),
       c(Q.skewness(q.10,q.50,q.90),
         Q.kurtosis(q.05,q.25,q.75,q.95)),type="n",
     axes = FALSE, xlab = "", ylab = "")
points(c(orchid_grow$log_area_t[orchid_grow$flowering==1],
       orchid_grow$log_area_t[orchid_grow$flowering==1]),
     c(Q.skewness(q.10[orchid_grow$flowering==1],
                  q.50[orchid_grow$flowering==1],
                  q.90[orchid_grow$flowering==1]),
       Q.kurtosis(q.05[orchid_grow$flowering==1],
                  q.25[orchid_grow$flowering==1],
                  q.75[orchid_grow$flowering==1],
                  q.95[orchid_grow$flowering==1])),
     col=c(rep(alpha("blue",0.25),nrow(orchid_grow[orchid_grow$flowering==1,])),
           rep(alpha("red",0.25),nrow(orchid_grow[orchid_grow$flowering==1,]))),
     pch=16,cex=.5)
abline(h=0,col="lightgray",lty=3)
axis(side = 4,cex.axis=0.8,at = pretty(range(c(Q.skewness(q.10,q.50,q.90),Q.kurtosis(q.05,q.25,q.75,q.95)))))
mtext("Skewness", side = 4, line = 2,col="blue")
mtext("Excess Kurtosis", side = 4, line =3,col="red")
title("D",font=3,adj=0)



dev.off()

