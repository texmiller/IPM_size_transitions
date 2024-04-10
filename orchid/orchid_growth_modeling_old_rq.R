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
library(qpdf)
# library(Rage)

### move to the right local directory 
tom = "C:/Users/tm9/Dropbox/github/IPM_size_transitions"
steve = "c:/repos/IPM_size_transitions" 
home = ifelse(Sys.info()["user"] == "Ellner", steve, tom)
setwd(home); 

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
orchid<- read_csv("orchid/Orchis_IPM_data.csv") %>% 
  ## there are two sites/populations, one in light and one in shade.
  ## for the purposes of this analysis, just take the light population
  filter(light=="L") %>% 
  mutate(log_area_t=log(total.leaf.area),
         log_area_t1=log(end.total.leaf.area)) 
seedlings<-read.csv("orchid/Orchis_seedlings.csv",T) %>% 
  filter(light=="L") %>% 
  mutate(log_area_t1=log(end.total.leaf.area)) 

## create a data subset for growth modeling
orchid %>% 
  dplyr::select(individual,log_area_t,log_area_t1,flowering,begin.year) %>% 
  drop_na() -> orchid_grow

## pilot gaussian model selection
orchid_GAU<-list()
orchid_GAU[[1]]<-lme(log_area_t1~log_area_t,weights=varExp(form=~log_area_t),
                     random=~1|begin.year,data=orchid_grow,method="ML")
orchid_GAU[[2]]<-lme(log_area_t1~log_area_t + as.logical(flowering),weights=varExp(form=~log_area_t),
                     random=~1|begin.year,data=orchid_grow,method="ML")
orchid_GAU[[3]]<-lme(log_area_t1~log_area_t*as.logical(flowering),weights=varExp(form=~log_area_t),
                     random=~1|begin.year,data=orchid_grow,method="ML")
## same models but sd differs between veg and flow
orchid_GAU[[4]]<-lme(log_area_t1~log_area_t,weights=varExp(form=~log_area_t*as.logical(flowering)),
                     random=~1|begin.year,data=orchid_grow,method="ML")
orchid_GAU[[5]]<-lme(log_area_t1~log_area_t + as.logical(flowering),weights=varExp(form=~log_area_t*as.logical(flowering)),
                     random=~1|begin.year,data=orchid_grow,method="ML")
orchid_GAU[[6]]<-lme(log_area_t1~log_area_t*as.logical(flowering),weights=varExp(form=~log_area_t*as.logical(flowering)),
                     random=~1|begin.year,data=orchid_grow,method="ML")
AICtab(orchid_GAU)
orchid_GAU_best<-orchid_GAU[[which.min(AICctab(orchid_GAU,sort=F)$dAICc)]]
orchid_grow$GAU_fitted <- fitted(orchid_GAU_best)
summary(orchid_GAU_best)

##refit best model with REML
orchid_GAU_best<-update(orchid_GAU_best,method="REML")
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

pdf("./manuscript/figures/orchid_diagnostics.pdf",height = 5, width = 6,useDingbats = F)
par(mar = c(4, 4, 1, 3), oma=c(0,0,0,0),mfrow=c(2,2)) 
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
legend("topleft",legend=c("Skewness","Kurtosis"),bg="white",pch=16,col=c("blue","red"),cex=0.8)
title("B",font=3,adj=0)

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
mtext("Standard deviation", side = 4, line = 2,cex=0.9)
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
mtext("Skewness or Kurtosis", side = 4, line = 2,cex=0.9)
title("D",font=3,adj=0)
dev.off()

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
  cat("Simulation ",i, "\n"); 
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

pdf("./manuscript/figures/orchid_SST_fit.pdf",height = 4, width = 8, useDingbats = F)
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
title("A",font=3,adj=0)
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
title("B",font=3,adj=0)

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
title("C",font=3,adj=0)

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
title("D",font=3,adj=0)

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
title("E",font=3,adj=0)
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
title("F",font=3,adj=0)

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
title("G",font=3,adj=0)

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
title("H",font=3,adj=0)
dev.off()


# IPM and life history analysis -------------------------------------------
## here are the additional components for the IPM
### survival (not affected by flowering)
surv<-glmer(survival~log_area_t+(1|begin.year),family="binomial",data=orchid)
### flowering probability
flower<-glmer(flowering~log_area_t+(1|begin.year),data=orchid,family="binomial")
### number of flowers produced by flowering plants
nflowers<-glmer(number.flowers~log_area_t+(1|begin.year),data=subset(orchid,flowering==1),family="poisson")
### proportion of fruits that set seed; Hans et al. used a size-dependent function
### but I get non-significant slopes in both environments, so I will use means
propfruit<-glmer(number.fruits/number.flowers~(1|begin.year),weights=number.flowers,data=subset(orchid,flowering==1),family="binomial")
### seedling size distributions
kidsize<-mean(seedlings$log_area_t1)
kidsd<-sd(seedlings$log_area_t1)
### seeds per fruit;see text
seeds<-6000	
### seed to protocorm transition; see text
eta<-mean(c(0.007,0.023))
### protocorm survival probability (from Eelke)
sigmap<-0.01485249
### tuber survival probability (from Eelke)
sigmat<-0.05940997
### size-dependent probability of entering dormancy
dormancy<-glmer(dormant~log_area_t+(1|begin.year),family="binomial",data=orchid)
### size distribution (mean and sd) of plants emerging from dormancy (assumed to be common across environments)
Dsize<-mean(orchid$log_area_t1[orchid$number.leaves==0],na.rm=T)
Dsizesd<-sd(orchid$log_area_t1[orchid$number.leaves==0],na.rm=T)
### observed size limits
minsize<-min(na.omit(c(orchid$log_area_t,orchid$log_area_t1)))
maxsize<-max(na.omit(c(orchid$log_area_t,orchid$log_area_t1)))
matsize<-100	##size of approximating matrix

##################################################################################
### now store parameters in vectors
###################################################################################

params<-c()
## survival params
params$surv.int <- fixef(surv)[1]
params$surv.size <- fixef(surv)[2]
## gaussian growth params
params$grow.int <- fixef(orchid_GAU_best)[1]
params$grow.size <- fixef(orchid_GAU_best)[2]
params$grow.flow <- fixef(orchid_GAU_best)[3]
params$grow.size.flow <- fixef(orchid_GAU_best)[4]
params$growsd.int <- orchid_GAU_best$sigma
params$growsd.size <- orchid_GAU_best$modelStruct$varStruct[1]
## SST growth params
params$growSST.nu.int<-SSTout$estimate[1]
params$growSST.nu.size<-SSTout$estimate[2]
params$growSST.nu.flow<-SSTout$estimate[3]
params$growSST.nu.size.flow<-SSTout$estimate[4]
params$growSST.tau.int<-SSTout$estimate[5]
params$growSST.tau.size<-SSTout$estimate[6]
params$growSST.tau.flow<-SSTout$estimate[7]
params$growSST.tau.size.flow<-SSTout$estimate[8]
## flowering params
params$flow.int <- fixef(flower)[1]
params$flow.size <- fixef(flower)[2]
## fruit production
params$flowers.int <- fixef(nflowers)[1]
params$flowers.size <- fixef(nflowers)[2]
## fruit set
params$fruits<-fixef(propfruit)[1]
## offspring
params$seeds<-seeds
params$kidsize<-kidsize
params$kidsize.sd<-kidsd
params$germ<-eta 
params$sigmap<-sigmap 
## dormancy
params$dorm.int<-fixef(dormancy)[1]
params$dorm.size<-fixef(dormancy)[2]
params$sigmat<-sigmat 
params$Dsize<-Dsize						
params$Dsizesd<-Dsizesd						
## matrix
params$minsize<-minsize
params$maxsize<-maxsize

#IPM source functions are here:
source("orchid/orchis.IPM.source.R")


params$matsize<-matsize


flowint<-seq(-30,0,0.1)
R0out.beta0.GAU<-R0out.beta0.SST<-vector("numeric",length=length(flowint))
for(i in 1:length(flowint)){
  params$flow.int<-flowint[i]
  R0out.beta0.GAU[i]<-returnR0(params,dist="GAU")$R0
  R0out.beta0.SST[i]<-returnR0(params,dist="SST")$R0
  if(i==length(flowint)){params$flow.int<-fixef(flower)[1]}
}

## compare mean and variance of life expectancies
## use Chrissy's variance function instead of RAGE. 
source("lifespan_moments_CMH.R"); 


mean.life.GAU<-mean.life.SST<-vector("numeric",length=params$matsize+3)
matU.GAU<-returnR0(params=params,dist="GAU")$T
var.life.GAU<-var.life.SST<-vector("numeric",length=params$matsize+3)
matU.SST<-returnR0(params=params,dist="SST")$T
meshpts<-returnR0(params=params,dist="GAU")$meshpts
for(i in 1:length(mean.life.GAU)){
  mixdist = rep(0,ncol(matU.GAU));  mixdist[i]=1; 
  
  mean.life.GAU[i]<-mean_lifespan(matU.GAU, mixdist)
  mean.life.SST[i]<-mean_lifespan(matU.SST, mixdist)

  #mean.life.GAU[i]<-life_expect_mean(matU = matU.GAU, start = i)
  #mean.life.SST[i]<-life_expect_mean(matU = matU.SST, start = i)


   var.life.GAU[i]<-var_lifespan(matU.GAU, mixdist)
   var.life.SST[i]<-var_lifespan(matU.SST, mixdist)
#  var.life.GAU[i]<-life_expect_var(matU = matU.GAU, start = i)
#  var.life.SST[i]<-life_expect_var(matU = matU.SST, start = i)

}

pdf("manuscript/figures/orchis_life_history.pdf",height=3,width=9)
par(mar = c(5, 4, 2, 1),mfrow=c(1,3)) 
plot(-flowint/params$flow.size,R0out.beta0.GAU,type="l",lwd=2,col="tomato",
     xlab="Flowering size (log leaf area)",ylab="R0",cex.lab=1.4)
lines(-flowint/params$flow.size,R0out.beta0.SST,lwd=2,col="cornflowerblue")
abline(v=-params$flow.int/params$flow.size,lty=2)
legend("topleft",legend=c("Gaussian","Skewed t"),title="Growth model:",
       lwd=2,col=c("tomato","cornflowerblue"),bty="n")
title("A",adj=0,font=3)

plot(c(-2.5,-2,-1.5,meshpts),mean.life.GAU,axes=F,
     xlab=" ",ylab=" ",type="n",lwd=2,ylim=c(range(c(mean.life.GAU,mean.life.SST))))
lines(meshpts,mean.life.GAU[-(1:3)],col="tomato",lwd=2)
lines(meshpts,mean.life.SST[-(1:3)],col="cornflowerblue",lwd=2)
points(c(-2.5,-2,-1.5),mean.life.GAU[1:3],pch=16,col="tomato",cex=1.5)
points(c(-2.5,-2,-1.5),mean.life.SST[1:3],pch=16,col="cornflowerblue",cex=1.5)
abline(v=-1)
box()
axis(side = 1,cex.axis=0.8,at = seq(-0.5,7.5,1))
axis(side = 1,labels=c("protocorm","tuber","dormant"),at=c(-2.5,-2,-1.5),las=2)
mtext("log(leaf area)", side = 1, line = 3,cex=.9)
mtext("Mean remaining lifespan (yrs)", side = 2, line = 3,cex=.9)
axis(side = 2,cex.axis=0.8,at = pretty(range(c(mean.life.GAU,mean.life.SST))))
title("B",adj=0,font=3)

plot(c(-2.5,-2,-1.5,meshpts),log(var.life.GAU),axes=F,
     xlab=" ",ylab=" ",type="n",lwd=2,ylim=c(range(c(log(var.life.GAU),log(var.life.SST)))))
lines(meshpts,log(var.life.GAU[-(1:3)]),col="tomato",lwd=2)
lines(meshpts,log(var.life.SST[-(1:3)]),col="cornflowerblue",lwd=2)
points(c(-2.5,-2,-1.5),log(var.life.GAU[1:3]),pch=16,col="tomato",cex=1.5)
points(c(-2.5,-2,-1.5),log(var.life.SST[1:3]),pch=16,col="cornflowerblue",cex=1.5)
abline(v=-1)
box()
axis(side = 1,cex.axis=0.8,at = seq(-0.5,7.5,1))
axis(side = 1,labels=c("protocorm","tuber","dormant"),at=c(-2.5,-2,-1.5),las=2)
mtext("log(leaf area)", side = 1, line = 3,cex=.9)
mtext("log(variance of remaining lifespan)", side = 2, line = 3,cex=.9)
axis(side = 2,cex.axis=0.8,at = pretty(range(c(log(var.life.GAU),log(var.life.SST)))))
title("C",adj=0,font=3)
dev.off()

## by how much do the variances differ?
maxdiff<-which.max(var.life.SST-var.life.GAU)
(var.life.SST-var.life.GAU)[maxdiff]/var.life.SST[maxdiff]

## check sensitivity to matrix dimensions -- DON'T USE THE VARIANCE RESULT
## 
dims<-c(50,100,150,200,250,300,350,400,500,600,800)
lambda.GAU<-lambda.SST<-c()
mean.life.GAU<-mean.life.SST<-c()
var.life.GAU<-var.life.SST<-c()
mature.age.GAU<-mature.age.SST<-c()
for(i in 1:length(dims)){
  params$matsize<-dims[i]
  
  R0.GAU<-returnR0(params=params,dist="GAU")
  lambda.GAU[i]<-lambda(R0.GAU$matrix)
  mean.life.GAU[i]<-life_expect_mean(matU = R0.GAU$Tmatrix, start = 10)
  var.life.GAU[i]<-life_expect_var(matU = R0.GAU$Tmatrix, start = 10)
  mature.age.GAU[i]<-mature_age(matU = R0.GAU$Tmatrix,matR = R0.GAU$Fmatrix,1)
  
  R0.SST<-returnR0(params=params,dist="SST")
  lambda.SST[i]<-lambda(R0.SST$matrix)
  mean.life.SST[i]<-life_expect_mean(matU = R0.SST$Tmatrix, start = 10)
  var.life.SST[i]<-life_expect_var(matU = R0.SST$Tmatrix, start = 10)
  mature.age.SST[i]<-mature_age(matU = R0.SST$Tmatrix,matR = R0.SST$Fmatrix,1)
}

win.graph()
par(mfrow=c(1,3))
plot(dims,lambda.GAU,col="tomato",ylim=range(c(lambda.GAU,lambda.SST)),
     xlab="Matrix dimension",ylab=expression(paste(lambda)),pch=16)
points(dims,lambda.SST,col="cornflowerblue",pch=16)
plot(dims,mean.life.GAU,col="tomato",ylim=range(c(mean.life.GAU,mean.life.SST)),
     xlab="Matrix dimension",ylab="Mean life expectancy",pch=16)
points(dims,mean.life.SST,col="cornflowerblue",pch=16)
## this one is real bad
plot(dims,var.life.GAU,col="tomato",ylim=range(c(var.life.GAU,var.life.SST)),
     xlab="Matrix dimension",ylab="Variance in life expectancy",pch=16)
points(dims,var.life.SST,col="cornflowerblue",pch=16)
legend("topright",legend=c("Gaussian","Skewed t"),title="Growth model:",
       lwd=2,col=c("tomato","cornflowerblue"),bty="n")

plot(dims,mature.age.GAU,col="tomato",ylim=range(c(mature.age.GAU,mature.age.SST)))
points(dims,mature.age.SST,col="cornflowerblue")


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

