## load libraries
library(tidyverse)
library(lme4)
library(nlme)
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

## create a data subset for growth modeling
orchid %>% 
  dplyr::select(log_area_t,log_area_t1,flowering,begin.year) %>% 
  drop_na() -> orchid_grow

## pilot gaussian model selection
orchid_GAU<-list()
orchid_GAU[[1]]<-lme(log_area_t1~log_area_t,weights=varExp(form=~log_area_t),
                     random=~1|begin.year,data=orchid_grow,method="ML")
orchid_GAU[[2]]<-lme(log_area_t1~log_area_t + as.logical(flowering),weights=varExp(form=~log_area_t),
                     random=~1|begin.year,data=orchid_grow,method="ML")
orchid_GAU[[3]]<-lme(log_area_t1~log_area_t*as.logical(flowering),weights=varExp(form=~log_area_t),
                     random=~1|begin.year,data=orchid_grow,method="ML")
AICtab(orchid_GAU)
orchid_GAU_best<-orchid_GAU[[which.min(AICctab(orchid_GAU,sort=F)$dAICc)]]
orchid_grow$GAU_fitted <- fitted(orchid_GAU_best)
summary(orchid_GAU_best)

## fit sd as a function of initial size
orchid_grow$GAU_resids <- residuals(orchid_GAU_best)
## derive sd from fitted params (see basement)
orchid_grow$GAU_sd <- orchid_GAU_best$sigma*exp(orchid_GAU_best$modelStruct$varStruct[1]*orchid_grow$log_area_t)
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
     xlab="Size at t+1",ylab="Scaled residuals of size at t+1")
points(orchid_grow$log_area_t[orchid_grow$flowering==0],
     orchid_grow$GAU_std_resids[orchid_grow$flowering==0],col=alpha("black",0.25))
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
     xlab="Size at t+1",ylab="Scaled residuals of size at t+1")
points(orchid_grow$log_area_t[orchid_grow$flowering==1],
     orchid_grow$GAU_std_resids[orchid_grow$flowering==1],col=alpha("black",0.25))
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

## now try JSU and skewed t, allowing skew and kurtosis to differ by flowering status
JSULogLik=function(pars){
  dJSU(orchid_grow$log_area_t1, 
       mu=orchid_grow$GAU_fitted,
       sigma=orchid_grow$GAU_sd,
       nu = pars[1] + pars[2]*orchid_grow$log_area_t + pars[3]*as.logical(orchid_grow$flowering) + pars[4]*orchid_grow$log_area_t*as.logical(orchid_grow$flowering),
       tau = exp(pars[5] + pars[6]*orchid_grow$log_area_t + pars[7]*as.logical(orchid_grow$flowering) + pars[8]*orchid_grow$log_area_t*as.logical(orchid_grow$flowering)), 
       log=TRUE)
}
## SST is the gamlss parameterization for which mu is mean and sigma sd
## see documentation for explanation of link fns
SSTLogLik=function(pars){
  dSST(orchid_grow$log_area_t1, 
       mu=orchid_grow$GAU_fitted,
       sigma=orchid_grow$GAU_sd,
       nu = exp(pars[1] + pars[2]*orchid_grow$log_area_t + pars[3]*as.logical(orchid_grow$flowering) + pars[4]*orchid_grow$log_area_t*as.logical(orchid_grow$flowering)),
       tau = exp(pars[5] + pars[6]*orchid_grow$log_area_t + pars[7]*as.logical(orchid_grow$flowering) + pars[8]*orchid_grow$log_area_t*as.logical(orchid_grow$flowering))+2, 
       log=TRUE)
}

## starting parameters
p0<-c(0,0,0,0,0,0,0,0)
## pass this through several ML algorithms
JSUout=maxLik(logLik=JSULogLik,start=p0*exp(0.2*rnorm(length(p0))),method="BHHH",control=list(iterlim=5000,printLevel=2),finalHessian=FALSE) 
JSUout=maxLik(logLik=JSULogLik,start=JSUout$estimate,method="NM",control=list(iterlim=5000,printLevel=1),finalHessian=FALSE)
JSUout=maxLik(logLik=JSULogLik,start=JSUout$estimate,method="BHHH",control=list(iterlim=5000,printLevel=2),finalHessian=FALSE)

SSTout=maxLik(logLik=SSTLogLik,start=p0*exp(0.2*rnorm(length(p0))),method="BHHH",control=list(iterlim=5000,printLevel=2),finalHessian=FALSE) 
SSTout=maxLik(logLik=SSTLogLik,start=SSTout$estimate,method="NM",control=list(iterlim=5000,printLevel=1),finalHessian=FALSE)
SSTout=maxLik(logLik=SSTLogLik,start=SSTout$estimate,method="BHHH",control=list(iterlim=5000,printLevel=2),finalHessian=FALSE)


## simulate from fitted model
n_sim<-100
JSUsim_mean<-JSUsim_sd<-JSUsim_skew<-JSUsim_kurt<-matrix(NA,nrow=nrow(orchid_grow),ncol=n_sim)
SSTsim_mean<-SSTsim_sd<-SSTsim_skew<-SSTsim_kurt<-matrix(NA,nrow=nrow(orchid_grow),ncol=n_sim)
GAUsim_mean<-GAUsim_sd<-GAUsim_skew<-GAUsim_kurt<-matrix(NA,nrow=nrow(orchid_grow),ncol=n_sim)
for(i in 1:n_sim){
  orchid_grow$log_area_t1.simJSU <- rJSU(n=nrow(orchid_grow),
                                         mu=orchid_grow$GAU_fitted,
                                         sigma=orchid_grow$GAU_sd,
                                         nu=JSUout$estimate[1]+JSUout$estimate[2]*orchid_grow$log_area_t+
                                           JSUout$estimate[3]*as.logical(orchid_grow$flowering)+JSUout$estimate[4]*as.logical(orchid_grow$flowering)*orchid_grow$log_area_t,
                                         tau=exp(JSUout$estimate[5]+JSUout$estimate[6]*orchid_grow$log_area_t+
                                                   JSUout$estimate[7]*as.logical(orchid_grow$flowering)+JSUout$estimate[8]*as.logical(orchid_grow$flowering)*orchid_grow$log_area_t))
  q.fit<-predict(rq(log_area_t1.simJSU~log_area_t*as.logical(flowering),data=orchid_grow,tau=c(0.05,0.1,0.25,0.5,0.75,0.9,0.95)))
  JSUsim_mean[,i]<-Q.mean(q.fit[,3],q.fit[,4],q.fit[,5])
  JSUsim_sd[,i]<-Q.sd(q.fit[,3],q.fit[,5])
  JSUsim_skew[,i]<-Q.skewness(q.fit[,2],q.fit[,4],q.fit[,6])
  JSUsim_kurt[,i]<-Q.kurtosis(q.fit[,1],q.fit[,3],q.fit[,5],q.fit[,7])
  
  orchid_grow$log_area_t1.simSST <- rSST(n=nrow(orchid_grow),
                                         mu=orchid_grow$GAU_fitted,
                                         sigma=orchid_grow$GAU_sd,
                                         nu=exp(SSTout$estimate[1]+SSTout$estimate[2]*orchid_grow$log_area_t+
                                           SSTout$estimate[3]*as.logical(orchid_grow$flowering)+SSTout$estimate[4]*as.logical(orchid_grow$flowering)*orchid_grow$log_area_t),
                                         tau=exp(SSTout$estimate[5]+SSTout$estimate[6]*orchid_grow$log_area_t+
                                                   SSTout$estimate[7]*as.logical(orchid_grow$flowering)+SSTout$estimate[8]*as.logical(orchid_grow$flowering)*orchid_grow$log_area_t)+2)
  q.fit<-predict(rq(log_area_t1.simSST~log_area_t*as.logical(flowering),data=orchid_grow,tau=c(0.05,0.1,0.25,0.5,0.75,0.9,0.95)))
  SSTsim_mean[,i]<-Q.mean(q.fit[,3],q.fit[,4],q.fit[,5])
  SSTsim_sd[,i]<-Q.sd(q.fit[,3],q.fit[,5])
  SSTsim_skew[,i]<-Q.skewness(q.fit[,2],q.fit[,4],q.fit[,6])
  SSTsim_kurt[,i]<-Q.kurtosis(q.fit[,1],q.fit[,3],q.fit[,5],q.fit[,7])
  
  orchid_grow$log_area_t1.simGAU <- rnorm(n=nrow(orchid_grow),
                                          mean=orchid_grow$GAU_fitted,
                                          sd=orchid_grow$GAU_sd)
  q.fit<-predict(rq(log_area_t1.simGAU~log_area_t*as.logical(flowering),data=orchid_grow,tau=c(0.05,0.1,0.25,0.5,0.75,0.9,0.95)))
  GAUsim_mean[,i]<-Q.mean(q.fit[,3],q.fit[,4],q.fit[,5])
  GAUsim_sd[,i]<-Q.sd(q.fit[,3],q.fit[,5])
  GAUsim_skew[,i]<-Q.skewness(q.fit[,2],q.fit[,4],q.fit[,6])
  GAUsim_kurt[,i]<-Q.kurtosis(q.fit[,1],q.fit[,3],q.fit[,5],q.fit[,7])
}

## and now the real data
q.fit<-predict(rq(log_area_t1~log_area_t*as.logical(flowering),data=orchid_grow,tau=c(0.05,0.1,0.25,0.5,0.75,0.9,0.95)))
q.05<-q.fit[,1]
q.10<-q.fit[,2] 
q.25<-q.fit[,3] 
q.50<-q.fit[,4]
q.75<-q.fit[,5]
q.90<-q.fit[,6] 
q.95<-q.fit[,7]

pdf("./manuscript/figures/orchid_SST_fit.pdf",height = 6, width = 6,useDingbats = F)
par(mfrow=c(2,4),mar=c(4,4,1,1))
plot(orchid_grow$log_area_t,Q.mean(q.25,q.50,q.75),type="n",
     xlab="size t",ylab="mean size t1",
     ylim=c(min(c(GAUsim_mean,JSUsim_mean,SSTsim_mean)),max(c(GAUsim_mean,JSUsim_mean,SSTsim_mean))))
for(i in 1:n_sim){
  points(orchid_grow$log_area_t[orchid_grow$flowering==0],GAUsim_mean[orchid_grow$flowering==0,i],col=alpha("tomato",0.25),pch=".")
  #points(orchid_grow$log_area_t[orchid_grow$flowering==0],JSUsim_mean[orchid_grow$flowering==0,i],col=alpha("forestgreen",0.25),pch=".")
  points(orchid_grow$log_area_t[orchid_grow$flowering==0],SSTsim_mean[orchid_grow$flowering==0,i],col=alpha("cornflowerblue",0.25),pch=".")
}
points(orchid_grow$log_area_t[orchid_grow$flowering==0],Q.mean(q.25[orchid_grow$flowering==0],
                                     q.50[orchid_grow$flowering==0],
                                     q.75[orchid_grow$flowering==0]),col="black",pch=".",cex=2)
legend("topleft",legend=c("Real data","Gaussian simulation","SST simulation"),
       lty=1,col=c("black","tomato","cornflowerblue"),cex=0.8,bty="n")

plot(orchid_grow$log_area_t,Q.sd(q.25,q.75),type="n",
     xlab="size t",ylab="sd size t1",
     ylim=c(min(c(GAUsim_sd,JSUsim_sd,SSTsim_sd)),max(c(GAUsim_sd,JSUsim_sd,SSTsim_sd))))
for(i in 1:n_sim){
  points(orchid_grow$log_area_t[orchid_grow$flowering==0],GAUsim_sd[orchid_grow$flowering==0,i],col=alpha("tomato",0.25),pch=".")
  #points(orchid_grow$log_area_t[orchid_grow$flowering==0],JSUsim_sd[orchid_grow$flowering==0,i],col=alpha("forestgreen",0.25),pch=".")
  points(orchid_grow$log_area_t[orchid_grow$flowering==0],SSTsim_sd[orchid_grow$flowering==0,i],col=alpha("cornflowerblue",0.25),pch=".")
}
points(orchid_grow$log_area_t[orchid_grow$flowering==0],
       Q.sd(q.25[orchid_grow$flowering==0],
            q.75[orchid_grow$flowering==0]),col="black",pch=".",cex=2)

plot(orchid_grow$log_area_t,Q.skewness(q.10,q.50,q.90),type="n",
     xlab="size t",ylab="skewness size t1",
     ylim=c(min(c(GAUsim_skew,JSUsim_skew,SSTsim_skew)),max(c(GAUsim_skew,JSUsim_skew,SSTsim_skew))))
for(i in 1:n_sim){
  points(orchid_grow$log_area_t[orchid_grow$flowering==0],GAUsim_skew[orchid_grow$flowering==0,i],col=alpha("tomato",0.25),pch=".")
  #points(orchid_grow$log_area_t[orchid_grow$flowering==0],JSUsim_skew[orchid_grow$flowering==0,i],col=alpha("forestgreen",0.25),pch=".")
  points(orchid_grow$log_area_t[orchid_grow$flowering==0],SSTsim_skew[orchid_grow$flowering==0,i],col=alpha("cornflowerblue",0.25),pch=".")
}
points(orchid_grow$log_area_t[orchid_grow$flowering==0],
       Q.skewness(q.10[orchid_grow$flowering==0],
                  q.50[orchid_grow$flowering==0],
                  q.90[orchid_grow$flowering==0]),col="black",pch=".",cex=2)

plot(orchid_grow$log_area_t,Q.kurtosis(q.05,q.25,q.75,q.95),type="n",
     xlab="size t",ylab="kurtosis size t1",
     ylim=c(min(GAUsim_kurt,JSUsim_kurt,SSTsim_kurt),max(GAUsim_kurt,JSUsim_kurt,SSTsim_kurt)))
for(i in 1:n_sim){
  points(orchid_grow$log_area_t[orchid_grow$flowering==0],GAUsim_kurt[orchid_grow$flowering==0,i],col=alpha("tomato",0.25),pch=".")
  #points(orchid_grow$log_area_t[orchid_grow$flowering==0],JSUsim_kurt[orchid_grow$flowering==0,i],col=alpha("forestgreen",0.25),pch=".")
  points(orchid_grow$log_area_t[orchid_grow$flowering==0],SSTsim_kurt[orchid_grow$flowering==0,i],col=alpha("cornflowerblue",0.25),pch=".")
}
points(orchid_grow$log_area_t[orchid_grow$flowering==0],
       Q.kurtosis(q.05[orchid_grow$flowering==0],
                  q.25[orchid_grow$flowering==0],
                  q.75[orchid_grow$flowering==0],
                  q.95[orchid_grow$flowering==0]),col="black",pch=".",cex=2)

plot(orchid_grow$log_area_t,Q.mean(q.25,q.50,q.75),type="n",
     xlab="size t",ylab="mean size t1",
     ylim=c(min(c(GAUsim_mean,JSUsim_mean,SSTsim_mean)),max(c(GAUsim_mean,JSUsim_mean,SSTsim_mean))))
for(i in 1:n_sim){
  points(orchid_grow$log_area_t[orchid_grow$flowering==1],GAUsim_mean[orchid_grow$flowering==1,i],col=alpha("tomato",0.25),pch=".")
  #points(orchid_grow$log_area_t[orchid_grow$flowering==1],JSUsim_mean[orchid_grow$flowering==1,i],col=alpha("forestgreen",0.25),pch=".")
  points(orchid_grow$log_area_t[orchid_grow$flowering==1],SSTsim_mean[orchid_grow$flowering==1,i],col=alpha("cornflowerblue",0.25),pch=".")
}
points(orchid_grow$log_area_t[orchid_grow$flowering==1],Q.mean(q.25[orchid_grow$flowering==1],
                                                               q.50[orchid_grow$flowering==1],
                                                               q.75[orchid_grow$flowering==1]),col="black",pch=".",cex=2)
plot(orchid_grow$log_area_t,Q.sd(q.25,q.75),type="n",
     xlab="size t",ylab="sd size t1",
     ylim=c(min(c(GAUsim_sd,JSUsim_sd,SSTsim_sd)),max(c(GAUsim_sd,JSUsim_sd,SSTsim_sd))))
for(i in 1:n_sim){
  points(orchid_grow$log_area_t[orchid_grow$flowering==1],GAUsim_sd[orchid_grow$flowering==1,i],col=alpha("tomato",0.25),pch=".")
  #points(orchid_grow$log_area_t[orchid_grow$flowering==1],JSUsim_sd[orchid_grow$flowering==1,i],col=alpha("forestgreen",0.25),pch=".")
  points(orchid_grow$log_area_t[orchid_grow$flowering==1],SSTsim_sd[orchid_grow$flowering==1,i],col=alpha("cornflowerblue",0.25),pch=".")
}
points(orchid_grow$log_area_t[orchid_grow$flowering==1],
       Q.sd(q.25[orchid_grow$flowering==1],
            q.75[orchid_grow$flowering==1]),col="black",pch=".",cex=2)

plot(orchid_grow$log_area_t,Q.skewness(q.10,q.50,q.90),type="n",
     xlab="size t",ylab="skewness size t1",
     ylim=c(min(c(GAUsim_skew,JSUsim_skew,SSTsim_skew)),max(c(GAUsim_skew,JSUsim_skew,SSTsim_skew))))
for(i in 1:n_sim){
  points(orchid_grow$log_area_t[orchid_grow$flowering==1],GAUsim_skew[orchid_grow$flowering==1,i],col=alpha("tomato",0.25),pch=".")
  #points(orchid_grow$log_area_t[orchid_grow$flowering==1],JSUsim_skew[orchid_grow$flowering==1,i],col=alpha("forestgreen",0.25),pch=".")
  points(orchid_grow$log_area_t[orchid_grow$flowering==1],SSTsim_skew[orchid_grow$flowering==1,i],col=alpha("cornflowerblue",0.25),pch=".")
}
points(orchid_grow$log_area_t[orchid_grow$flowering==1],
       Q.skewness(q.10[orchid_grow$flowering==1],
                  q.50[orchid_grow$flowering==1],
                  q.90[orchid_grow$flowering==1]),col="black",pch=".",cex=2)

plot(orchid_grow$log_area_t,Q.kurtosis(q.05,q.25,q.75,q.95),type="n",
     xlab="size t",ylab="kurtosis size t1",
     ylim=c(min(GAUsim_kurt,JSUsim_kurt,SSTsim_kurt),max(GAUsim_kurt,JSUsim_kurt,SSTsim_kurt)))
for(i in 1:n_sim){
  points(orchid_grow$log_area_t[orchid_grow$flowering==1],GAUsim_kurt[orchid_grow$flowering==1,i],col=alpha("tomato",0.25),pch=".")
  #points(orchid_grow$log_area_t[orchid_grow$flowering==1],JSUsim_kurt[orchid_grow$flowering==1,i],col=alpha("forestgreen",0.25),pch=".")
  points(orchid_grow$log_area_t[orchid_grow$flowering==1],SSTsim_kurt[orchid_grow$flowering==1,i],col=alpha("cornflowerblue",0.25),pch=".")
}
points(orchid_grow$log_area_t[orchid_grow$flowering==1],
       Q.kurtosis(q.05[orchid_grow$flowering==1],
                  q.25[orchid_grow$flowering==1],
                  q.75[orchid_grow$flowering==1],
                  q.95[orchid_grow$flowering==1]),col="black",pch=".",cex=2)
dev.off()
#### The Basement ###################################
## make sure I know how to use lme's variance function
x<-runif(1000,0,10)
group<-rep(LETTERS[1:10],each=100)
group_means<-rep(rnorm(10,50,20),each=100)
y<-rnorm(1000,mean=group_means-3.5*x,sd=exp(1.5+0.14*x))    
plot(x,y)
test<-lme(y~x,
          weights=varExp(form=~x),
          random=~1|group,method="ML")
summary(test)
plot(x,exp(1.5+0.14*x))
points(x,test$sigma*exp(test$modelStruct$varStruct[1]*x),col="red")

