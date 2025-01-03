## load libraries
library(tidyverse)
library(lme4)
library(scales)
library(qgam)
library(gamlss.dist)
library(popbio)
library(moments)
library(maxLik)
library(bbmle)
library(qpdf)
library(oce)
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
orchid_GAU[[1]]<-lmer(log_area_t1~log_area_t+(1|begin.year), data=orchid_grow,REML=F)
orchid_GAU[[2]]<-lmer(log_area_t1~log_area_t + as.logical(flowering) + (1|begin.year),data=orchid_grow,REML=F)
orchid_GAU[[3]]<-lmer(log_area_t1~log_area_t * as.logical(flowering) + (1|begin.year),data=orchid_grow,REML=F)
AICtab(orchid_GAU)
orchid_GAU_best<-orchid_GAU[[which.min(AICctab(orchid_GAU,sort=F)$dAICc)]]
orchid_grow$GAU_fitted <- fitted(orchid_GAU_best)

## now use iterative re-weighting to fit sd as function of expected value
sdloglik = function(pars,resids,fitted) {
  dnorm(resids, mean=0, sd=exp(pars[1]+pars[2]*fitted+pars[3]*fitted^2),log=TRUE)
}	
pars<-list()
for(mod in 1:length(orchid_GAU)) {
  err = 1; rep=0; 
  while(err > 0.000001) {
    rep=rep+1
    model = orchid_GAU[[mod]]
    fitted_vals = fitted(model)
    resids = residuals(model) 
    out=maxLik(logLik=sdloglik,start=c(sd(resids),0,0),resids=resids,fitted=fitted_vals)
    pars[[mod]]=out$estimate 
    new_sigma = exp(pars[[mod]][1] + pars[[mod]][2]*fitted_vals + pars[[mod]][3]*fitted_vals^2)
    new_weights = 1/((new_sigma)^2)
    new_weights = 0.5*(weights(model) + new_weights) # cautious update 
    new_model <- update(model,weights=new_weights) 
    err = weights(model)-weights(new_model)
    err=sqrt(mean(err^2))
    cat(mod,rep,err,"\n") # check on convergence of estimated weights 
    orchid_GAU[[mod]]<-new_model 
  }
}

##re-do model selection using the best weights for all models
best_weights<-weights(orchid_GAU[[which.min(AICctab(orchid_GAU,sort=F)$dAICc)]])
for(mod in 1:length(orchid_GAU)) {
  orchid_GAU[[mod]]<-update(orchid_GAU[[mod]],weights=best_weights)
}
AICtab(orchid_GAU)
orchid_GAU_best<-orchid_GAU[[which.min(AICctab(orchid_GAU,sort=F)$dAICc)]]
##refit best model with REML
orchid_GAU_best<-update(orchid_GAU_best,REML=T)
best_weights<-weights(orchid_GAU_best)
orchid_grow$GAU_fitted <- fitted(orchid_GAU_best)
orchid_grow$GAU_sd <- 1/sqrt(best_weights)
orchid_grow$GAU_resids <- residuals(orchid_GAU_best)
orchid_grow$GAU_scaled_resids <- orchid_grow$GAU_resids*sqrt(best_weights) ##sqrt(weights)=1/sd
## should be mean zero unit variance
mean(orchid_grow$GAU_scaled_resids);sd(orchid_grow$GAU_scaled_resids)
## store stdev coefs
stdev_coef <- out$estimate
##this will help with plotting lines
orchid_grow<-arrange(orchid_grow,GAU_fitted)


#############################################################
# Test for variance trend in scaled residuals 
#############################################################
source("code/variance_diagnostics.R"); 

stopCluster(c1); 
c1<- makeCluster(6); 
registerDoParallel(c1);
out_levene = multiple_levene_test(orchid_grow$GAU_fitted,orchid_grow$GAU_scaled_resids,3,8,2500); 
out_levene$p_value; # p=0.09
out_levene = multiple_levene_test(orchid_grow$log_area_t,orchid_grow$GAU_scaled_resids,3,8,2500); 
out_levene$p_value; # p=0.14

out_ss = ss_test(orchid_grow$GAU_fitted,orchid_grow$GAU_scaled_resids,2500); 
out_ss$p_value; # p=0.22
out_ss = ss_test(orchid_grow$log_area_t,orchid_grow$GAU_scaled_resids,2500); 
out_ss$p_value; # p=0.22
stopCluster(c1);

### No trend in mean  
mfit = rsq.smooth.spline(orchid_grow$GAU_fitted,orchid_grow$GAU_scaled_resids) 
mfit$rsq; mfit$adj.rsq
mfit = rsq.smooth.spline(orchig_grow$log_area_t,orchid_grow$GAU_scaled_resids) 
mfit$rsq; mfit$adj.rsq


### No trend in variance 
vfit = rsq.smooth.spline(orchid_grow$GAU_fitted,abs(orchid_grow$GAU_scaled_resids)) 
vfit$rsq; vfit$adj.rsq
vfit = rsq.smooth.spline(orchid_grow$log_area_t,abs(orchid_grow$GAU_scaled_resids)) 
vfit$rsq; vfit$adj.rsq



 

## quantile regressions on stand resids
k_param=4
q.05<-predict(qgam(GAU_scaled_resids~s(GAU_fitted,k=k_param), data=orchid_grow,qu=0.05))
q.10<-predict(qgam(GAU_scaled_resids~s(GAU_fitted,k=k_param), data=orchid_grow,qu=0.1))
q.25<-predict(qgam(GAU_scaled_resids~s(GAU_fitted,k=k_param), data=orchid_grow,qu=0.25))
q.50<-predict(qgam(GAU_scaled_resids~s(GAU_fitted,k=k_param), data=orchid_grow,qu=0.5)) 
q.75<-predict(qgam(GAU_scaled_resids~s(GAU_fitted,k=k_param), data=orchid_grow,qu=0.75))
q.90<-predict(qgam(GAU_scaled_resids~s(GAU_fitted,k=k_param), data=orchid_grow,qu=0.9)) 
q.95<-predict(qgam(GAU_scaled_resids~s(GAU_fitted,k=k_param), data=orchid_grow,qu=0.95))
## NP skewness
NPS_hat = Q.skewness(q.10,q.50,q.90)
## NP kurtosis (relative to Gaussian)
NPK_hat = Q.kurtosis(q.05,q.25,q.75,q.95)

## plot raw growth data, sd(fitted as inset), and skew and kurt against std resids
veg_size<-pretty(orchid_grow$log_area_t[orchid_grow$flowering==0],n=50)
flow_size<-pretty(orchid_grow$log_area_t[orchid_grow$flowering==1],n=50)

pdf("./manuscript/figures/orchid_diagnostics.pdf",height = 5, width = 10,useDingbats = F)
par(mar = c(4, 4, 1, 3), oma=c(0,0,0,0),mfrow=c(1,2)) 
plot(orchid_grow$log_area_t,orchid_grow$log_area_t1,type="n",
     xlab="Size at t",ylab="Size at t+1")
points(orchid_grow$log_area_t[orchid_grow$flowering==0],
       orchid_grow$log_area_t1[orchid_grow$flowering==0],col=alpha("red",0.1))
lines(veg_size,fixef(orchid_GAU_best)[1]+fixef(orchid_GAU_best)[2]*veg_size,
       col=alpha("red",0.75),lwd=2)
points(orchid_grow$log_area_t[orchid_grow$flowering==1],
       orchid_grow$log_area_t1[orchid_grow$flowering==1],col=alpha("blue",0.1))
lines(flow_size,fixef(orchid_GAU_best)[1]+fixef(orchid_GAU_best)[3]+(fixef(orchid_GAU_best)[2]+fixef(orchid_GAU_best)[4])*flow_size,
       col=alpha("blue",0.75),lwd=2)
plotInset(2.5, 0.25, 7, 2,
          expr = plot(orchid_grow$GAU_fitted,orchid_grow$GAU_sd,
                      type="l", xlab = "Fitted value",
                      ylab = "Std Dev",cex.lab=0.8,
                      cex.axis = 0.5, mgp = c(3/2, 1/2, 0)),
          mar = c(0, 3, 0, 0))
title("A",font=3,adj=0)

plot(orchid_grow$GAU_fitted,orchid_grow$GAU_scaled_resids,type="n",
     xlab="Fitted value",ylab="Scaled residuals of size at t+1")
points(orchid_grow$GAU_fitted,q.05,col="black",pch=".")
points(orchid_grow$GAU_fitted,q.10,col="black",pch=".")
points(orchid_grow$GAU_fitted,q.25,col="black",pch=".")
points(orchid_grow$GAU_fitted,q.50,col="black",pch=".")
points(orchid_grow$GAU_fitted,q.75,col="black",pch=".")
points(orchid_grow$GAU_fitted,q.90,col="black",pch=".")
points(orchid_grow$GAU_fitted,q.95,col="black",pch=".")
par(new = TRUE)      
plot(c(orchid_grow$GAU_fitted,
         orchid_grow$GAU_fitted),
       c(NPS_hat,NPK_hat),
       col=c(rep(alpha("blue",0.25),nrow(orchid_grow)),
             rep(alpha("red",0.25),nrow(orchid_grow))),
       pch=16,cex=.5,axes = FALSE, xlab = "", ylab = "")
abline(h=0,col="lightgray",lty=3)
axis(side = 4,cex.axis=0.8,at = pretty(range(c(NPS_hat,NPK_hat))))
legend("topleft",legend=c("Skewness","Kurtosis"),bg="white",pch=16,col=c("blue","red"),cex=0.8)
title("B",font=3,adj=0)
dev.off()

## save plot elements so these can be re-drawn with other species
orchid_out<-list(
orchid_grow = orchid_grow[,c("log_area_t","log_area_t1","flowering","GAU_fitted","GAU_sd","GAU_scaled_resids")],
veg_size = veg_size,
flow_size = flow_size,
orchid_GAU_best = orchid_GAU_best,
q.05 = q.05,
q.10 = q.10,
q.25 = q.25,
q.50 = q.50,
q.75 = q.75,
q.90 = q.90,
q.95 = q.95,
NPS_hat = NPS_hat,
NPK_hat = NPK_hat)
############################################################
## Improvement: skewness and excess kurtosis
## try JSU and skewed t

JSULogLik=function(pars){
  dJSU(orchid_grow$log_area_t1, 
       mu=orchid_grow$GAU_fitted,
       sigma=orchid_grow$GAU_sd,
       nu = pars[1] + pars[2]*orchid_grow$GAU_fitted,
       tau = exp(pars[3] + pars[4]*orchid_grow$GAU_fitted), 
       log=TRUE)
}
## SST is the gamlss parameterization for which mu is mean and sigma sd
## see documentation for explanation of link fns
SSTLogLik=function(pars){
  dSST(orchid_grow$log_area_t1, 
       mu=orchid_grow$GAU_fitted,
       sigma=orchid_grow$GAU_sd,
       nu = exp(pars[1] + pars[2]*orchid_grow$GAU_fitted),
       tau = exp(pars[3] + pars[4]*orchid_grow$GAU_fitted)+2, 
       log=TRUE)
}

## starting parameters
p0<-c(0,0,0,0)
## pass this through several ML algorithms
JSUout=maxLik(logLik=JSULogLik,start=p0*exp(0.2*rnorm(length(p0))),method="BHHH",control=list(iterlim=5000,printLevel=2),finalHessian=FALSE) 
JSUout=maxLik(logLik=JSULogLik,start=JSUout$estimate,method="NM",control=list(iterlim=5000,printLevel=1),finalHessian=FALSE)
JSUout=maxLik(logLik=JSULogLik,start=JSUout$estimate,method="BHHH",control=list(iterlim=5000,printLevel=2),finalHessian=FALSE)

SSTout=maxLik(logLik=SSTLogLik,start=p0*exp(0.2*rnorm(length(p0))),method="BHHH",control=list(iterlim=5000,printLevel=2),finalHessian=FALSE) 
SSTout=maxLik(logLik=SSTLogLik,start=SSTout$estimate,method="NM",control=list(iterlim=5000,printLevel=1),finalHessian=FALSE)
SSTout=maxLik(logLik=SSTLogLik,start=SSTout$estimate,method="BHHH",control=list(iterlim=5000,printLevel=2),finalHessian=FALSE)

AIC(JSUout);AIC(SSTout) ## SST favored by AIC


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
                                         nu=JSUout$estimate[1]+JSUout$estimate[2]*orchid_grow$GAU_fitted,
                                         tau=exp(JSUout$estimate[3]+JSUout$estimate[4]*orchid_grow$GAU_fitted))
  ## Qreg on sim data
  q.05<-predict(qgam(log_area_t1.simJSU~s(log_area_t,k=k_param), data=orchid_grow,qu=0.05)) 
  q.10<-predict(qgam(log_area_t1.simJSU~s(log_area_t,k=k_param), data=orchid_grow,qu=0.10)) 
  q.25<-predict(qgam(log_area_t1.simJSU~s(log_area_t,k=k_param), data=orchid_grow,qu=0.25))
  q.50<-predict(qgam(log_area_t1.simJSU~s(log_area_t,k=k_param), data=orchid_grow,qu=0.5))
  q.75<-predict(qgam(log_area_t1.simJSU~s(log_area_t,k=k_param), data=orchid_grow,qu=0.75)) 
  q.90<-predict(qgam(log_area_t1.simJSU~s(log_area_t,k=k_param), data=orchid_grow,qu=0.90))
  q.95<-predict(qgam(log_area_t1.simJSU~s(log_area_t,k=k_param), data=orchid_grow,qu=0.95))
  JSUsim_mean[,i]<-Q.mean(q.25,q.50,q.75)
  JSUsim_sd[,i]<-Q.sd(q.25,q.75)
  JSUsim_skew[,i]<-Q.skewness(q.10,q.50,q.90)
  JSUsim_kurt[,i]<-Q.kurtosis(q.05,q.25,q.75,q.95)
  
  orchid_grow$log_area_t1.simSST <- rSST(n=nrow(orchid_grow),
                                         mu=orchid_grow$GAU_fitted,
                                         sigma=orchid_grow$GAU_sd,
                                         nu=exp(SSTout$estimate[1]+SSTout$estimate[2]*orchid_grow$GAU_fitted),
                                         tau=exp(SSTout$estimate[3]+SSTout$estimate[4]*orchid_grow$GAU_fitted)+2)
  q.05<-predict(qgam(log_area_t1.simSST~s(log_area_t,k=k_param),data=orchid_grow,qu=0.05)) 
  q.10<-predict(qgam(log_area_t1.simSST~s(log_area_t,k=k_param), data=orchid_grow,qu=0.10)) 
  q.25<-predict(qgam(log_area_t1.simSST~s(log_area_t,k=k_param), data=orchid_grow,qu=0.25))
  q.50<-predict(qgam(log_area_t1.simSST~s(log_area_t,k=k_param), data=orchid_grow,qu=0.5))
  q.75<-predict(qgam(log_area_t1.simSST~s(log_area_t,k=k_param), data=orchid_grow,qu=0.75)) 
  q.90<-predict(qgam(log_area_t1.simSST~s(log_area_t,k=k_param), data=orchid_grow,qu=0.90))
  q.95<-predict(qgam(log_area_t1.simSST~s(log_area_t,k=k_param), data=orchid_grow,qu=0.95))
  SSTsim_mean[,i]<-Q.mean(q.25,q.50,q.75)
  SSTsim_sd[,i]<-Q.sd(q.25,q.75)
  SSTsim_skew[,i]<-Q.skewness(q.10,q.50,q.90)
  SSTsim_kurt[,i]<-Q.kurtosis(q.05,q.25,q.75,q.95)
  
  orchid_grow$log_area_t1.simGAU <- rnorm(n=nrow(orchid_grow),
                                          mean=orchid_grow$GAU_fitted,
                                          sd=orchid_grow$GAU_sd)
  q.05<-predict(qgam(log_area_t1.simGAU~s(log_area_t,k=k_param), data=orchid_grow,qu=0.05)) 
  q.10<-predict(qgam(log_area_t1.simGAU~s(log_area_t,k=k_param), data=orchid_grow,qu=0.10)) 
  q.25<-predict(qgam(log_area_t1.simGAU~s(log_area_t,k=k_param), data=orchid_grow,qu=0.25))
  q.50<-predict(qgam(log_area_t1.simGAU~s(log_area_t,k=k_param), data=orchid_grow,qu=0.5))
  q.75<-predict(qgam(log_area_t1.simGAU~s(log_area_t,k=k_param), data=orchid_grow,qu=0.75)) 
  q.90<-predict(qgam(log_area_t1.simGAU~s(log_area_t,k=k_param), data=orchid_grow,qu=0.90))
  q.95<-predict(qgam(log_area_t1.simGAU~s(log_area_t,k=k_param), data=orchid_grow,qu=0.95))
  GAUsim_mean[,i]<-Q.mean(q.25,q.50,q.75)
  GAUsim_sd[,i]<-Q.sd(q.25,q.75)
  GAUsim_skew[,i]<-Q.skewness(q.10,q.50,q.90)
  GAUsim_kurt[,i]<-Q.kurtosis(q.05,q.25,q.75,q.95)
}

## and now the real data
q.05<-predict(qgam(log_area_t1~s(log_area_t,k=k_param), data=orchid_grow,qu=0.05)) 
q.10<-predict(qgam(log_area_t1~s(log_area_t,k=k_param), data=orchid_grow,qu=0.10)) 
q.25<-predict(qgam(log_area_t1~s(log_area_t,k=k_param), data=orchid_grow,qu=0.25))
q.50<-predict(qgam(log_area_t1~s(log_area_t,k=k_param), data=orchid_grow,qu=0.5))
q.75<-predict(qgam(log_area_t1~s(log_area_t,k=k_param), data=orchid_grow,qu=0.75)) 
q.90<-predict(qgam(log_area_t1~s(log_area_t,k=k_param), data=orchid_grow,qu=0.90))
q.95<-predict(qgam(log_area_t1~s(log_area_t,k=k_param), data=orchid_grow,qu=0.95))



pdf("./manuscript/figures/orchid_SST_fit.pdf",height = 8, width = 8, useDingbats = F)
par(mfrow=c(2,2),mar=c(4,4,1,1))
plot(orchid_grow$log_area_t,Q.mean(q.25,q.50,q.75),type="n",
     xlab="size t",ylab="mean size t1",
     ylim=c(min(c(GAUsim_mean,JSUsim_mean,SSTsim_mean)),1+max(c(GAUsim_mean,JSUsim_mean,SSTsim_mean))))
matpoints(orchid_grow$log_area_t,GAUsim_mean,col=alpha("tomato",0.25),type="l",lty=1)
matpoints(orchid_grow$log_area_t,1+SSTsim_mean,col=alpha("cornflowerblue",0.25),type="l",lty=1)
lines(orchid_grow$log_area_t,Q.mean(q.25,q.50,q.75),col="black",lwd=2)
lines(orchid_grow$log_area_t,1+Q.mean(q.25,q.50,q.75),col="black",lwd=2)
title("A",font=3,adj=0)
legend("topleft",legend=c("Real data","Gaussian simulation","SST simulation + offset"),
       lty=1,col=c("black","tomato","cornflowerblue"),cex=0.8,bty="n")

plot(orchid_grow$log_area_t,Q.sd(q.25,q.75),type="n",
     xlab="size t",ylab="sd size t1",
     ylim=c(min(c(GAUsim_sd,JSUsim_sd,SSTsim_sd)),1+max(c(GAUsim_sd,JSUsim_sd,SSTsim_sd))))
matpoints(orchid_grow$log_area_t,GAUsim_sd,col=alpha("tomato",0.25),type="l",lty=1)
matpoints(orchid_grow$log_area_t,1+SSTsim_sd,col=alpha("cornflowerblue",0.25),type="l",lty=1)
points(orchid_grow$log_area_t,Q.sd(q.25,q.75),col="black",type="l",lwd=2)
points(orchid_grow$log_area_t,1+Q.sd(q.25,q.75),col="black",type="l",lwd=2)
title("B",font=3,adj=0)

plot(orchid_grow$log_area_t,Q.skewness(q.10,q.50,q.90),type="n",
     xlab="size t",ylab="skewness size t1",
     ylim=c(min(c(GAUsim_skew,JSUsim_skew,SSTsim_skew)),1+max(c(GAUsim_skew,JSUsim_skew,SSTsim_skew))))
matpoints(orchid_grow$log_area_t,GAUsim_skew,col=alpha("tomato",0.25),type="l",lty=1)
matpoints(orchid_grow$log_area_t,1+SSTsim_skew,col=alpha("cornflowerblue",0.25),type="l",lty=1)
points(orchid_grow$log_area_t,Q.skewness(q.10,q.50,q.90),col="black",type="l",lwd=2)
points(orchid_grow$log_area_t,1+Q.skewness(q.10,q.50,q.90),col="black",type="l",lwd=2)
title("C",font=3,adj=0)

plot(orchid_grow$log_area_t,Q.kurtosis(q.05,q.25,q.75,q.95),type="n",
     xlab="size t",ylab="kurtosis size t1",
     ylim=c(min(GAUsim_kurt,JSUsim_kurt,SSTsim_kurt),1+max(GAUsim_kurt,JSUsim_kurt,SSTsim_kurt)))
matpoints(orchid_grow$log_area_t,GAUsim_kurt,col=alpha("tomato",0.25),type="l",lty=1)
matpoints(orchid_grow$log_area_t,1+SSTsim_kurt,col=alpha("cornflowerblue",0.25),type="l",lty=1)
points(orchid_grow$log_area_t,Q.kurtosis(q.05,q.25,q.75,q.95),col="black",type="l",lwd=2)
points(orchid_grow$log_area_t,1+Q.kurtosis(q.05,q.25,q.75,q.95),col="black",type="l",lwd=2)
title("D",font=3,adj=0)

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

##################################################################################
### store parameters in vector
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
params$growsd.int <- stdev_coef[1]
params$growsd.fit <- stdev_coef[2]
params$growsd.fit2 <- stdev_coef[3]
## SST growth params
params$growSST.nu.int<-SSTout$estimate[1]
params$growSST.nu.fit<-SSTout$estimate[2]
params$growSST.tau.int<-SSTout$estimate[3]
params$growSST.tau.fit<-SSTout$estimate[4]
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
params$minsize<-minsize*0.95 #5% lower than observed
params$maxsize<-maxsize*1.05 #5% higher than observed
params$matsize<-100 ##explored below

#IPM source functions are here:
source("orchid/orchis.IPM.source.R")

##check eviction by setting survival to 100%
plot(colSums(returnR0(params,dist="GAU",lower.extend=5,upper.extend=2)$T))
plot(colSums(returnR0(params,dist="SST",lower.extend=5,upper.extend=2)$T))
## these are the extensions I need -- change survival back!

## test matrix dimensions
dims <- c(25,50,100,200,300,500)
R0_dims <- c()
for(i in 1:length(dims)){
  params$matsize<-dims[i]
  R0_dims[i]<-returnR0(params,dist="GAU",lower.extend=5,upper.extend=2)$R0
}
plot(dims,R0_dims)

flowint<-seq(-30,0,0.1)
R0out.beta0.GAU<-R0out.beta0.SST<-vector("numeric",length=length(flowint))
for(i in 1:length(flowint)){
  params$flow.int<-flowint[i]
  R0out.beta0.GAU[i]<-returnR0(params,dist="GAU",lower.extend=5,upper.extend=2)$R0
  R0out.beta0.SST[i]<-returnR0(params,dist="SST",lower.extend=5,upper.extend=2)$R0
  if(i==length(flowint)){params$flow.int<-fixef(flower)[1]}
}

pdf("./manuscript/figures/orchid_R0.pdf",height = 5, width = 5,useDingbats = F)
par(mar=c(5,5,1,1))
plot(-flowint/params$flow.size,R0out.beta0.GAU,type="l",lwd=3,col="#ef8a62",
     xlab="Flowering size (log leaf area)",ylab=expression(R[0]),cex.lab=1.4)
lines(-flowint/params$flow.size,R0out.beta0.SST,lwd=3,col="#67a9cf")
abline(v=-params$flow.int/params$flow.size,lty=2)
legend("topleft",legend=c("Gaussian","Skewed t"),title="Growth model:",
       lwd=3,col=c("#ef8a62","#67a9cf"),bty="n")
dev.off()

## write out matrices
orchid_out$mat_GAU<-returnR0(params,dist="GAU",lower.extend=5,upper.extend=2)
orchid_out$mat_SST<-returnR0(params,dist="SST",lower.extend=5,upper.extend=2)
write_rds(orchid_out,file="orchid/orchid_out.rds")

#### The Basement ###################################
## I thought I might do better fitting the sd as a gam--turns out not really

### I do not like the variance trend in the scaled residuals (again!)
## would we do better with a gam?
orchid_grow$year_rfx <- as.factor(orchid_grow$begin.year)
orchid_GAU[[4]]<-gam(list(log_area_t1~log_area_t * as.logical(flowering) + s(year_rfx,bs="re"),~s(log_area_t)),
                     data=orchid_grow,family="gaulss",method="ML",gamma=1.4)
AICtab(orchid_GAU) ## gam much (much?) better

## iteratively re-weight for sd as function of fitted
fit_gaulss = orchid_GAU[[4]]
fitted_all = predict(fit_gaulss,type="response",data=orchid_grow)                  
new_fitted_vals = fitted_all[,1]
orchid_grow$fitted_vals = new_fitted_vals 
weights = fitted_all[,2]
fit_gaulss = gam(list(log_area_t1~log_area_t * as.logical(flowering) + s(year_rfx,bs="re"),~s(fitted_vals,k=4)),
                 data=orchid_grow,family="gaulss",method="ML",gamma=1.4) 
err=100; k=0; 
while(err>10^(-6)) {
  orchid_grow$fitted_vals = new_fitted_vals 
  fit_gaulss<-update(fit_gaulss)
  fitted_all = predict(fit_gaulss,type="response",data=orchid_grow)   
  new_fitted_vals = fitted_all[,1]; new_weights = fitted_all[,2]
  err = weights - new_weights; err=sqrt(mean(err^2)); 
  weights = new_weights; 
  k=k+1; cat(k,err,"\n"); 
}   
orchid_GAU[[4]] = fit_gaulss 
AICtab(orchid_GAU)

##update lmer models with best weights from gam
orchid_GAU[[1]]<-update(orchid_GAU[[1]],weights=out[,2]^2)
orchid_GAU[[2]]<-update(orchid_GAU[[2]],weights=out[,2]^2)
orchid_GAU[[3]]<-update(orchid_GAU[[3]],weights=out[,2]^2)
AICtab(orchid_GAU)

##finally, re-fit best model with REML
orchid_GAU_best<-update(orchid_GAU[[4]],method="REML")
out = predict(orchid_GAU_best,type="response") 
orchid_grow$GAU_fitted <- out[,1]
orchid_grow$GAU_sd <- 1/out[,2]
## note the residuals are already "scaled" under the hood
orchid_grow$GAU_scaled_resids <- residuals(orchid_GAU_best)
## should be mean zero unit variance
mean(orchid_grow$GAU_scaled_resids);sd(orchid_grow$GAU_scaled_resids)
## plot fitted sd
out = predict(orchid_GAU[[4]],type="response") 
plot(orchid_grow$GAU_fitted,orchid_grow$GAU_sd)

## let's look again now at the variance trend
## quantile regressions on stand resids
k_param=4
S.05<-qgam(GAU_scaled_resids~s(GAU_fitted,k=k_param), data=orchid_grow,qu=0.05)
S.10<-qgam(GAU_scaled_resids~s(GAU_fitted,k=k_param), data=orchid_grow,qu=0.1)
S.25<-qgam(GAU_scaled_resids~s(GAU_fitted,k=k_param), data=orchid_grow,qu=0.25)
S.50<-qgam(GAU_scaled_resids~s(GAU_fitted,k=k_param), data=orchid_grow,qu=0.5) 
S.75<-qgam(GAU_scaled_resids~s(GAU_fitted,k=k_param), data=orchid_grow,qu=0.75)
S.90<-qgam(GAU_scaled_resids~s(GAU_fitted,k=k_param), data=orchid_grow,qu=0.9) 
S.95<-qgam(GAU_scaled_resids~s(GAU_fitted,k=k_param), data=orchid_grow,qu=0.95)
## NP skewness
NPS_hat = Q.skewness(q.10=predict(S.10),q.50=predict(S.50),q.90=predict(S.90))
## NP kurtosis (relative to Gaussian)
NPK_hat = Q.kurtosis(q.05=predict(S.05),q.25=predict(S.25),q.75=predict(S.75),q.95=predict(S.95))

plot(orchid_grow$GAU_fitted,orchid_grow$GAU_scaled_resids,type="n",
     xlab=expression(paste("E[", Size[t + 1], "]")),ylab="Scaled residuals of size at t+1")
points(orchid_grow$GAU_fitted,predict(S.05),col="black",pch=".")
points(orchid_grow$GAU_fitted,predict(S.10),col="black",pch=".")
points(orchid_grow$GAU_fitted,predict(S.25),col="black",pch=".")
points(orchid_grow$GAU_fitted,predict(S.50),col="black",pch=".")
points(orchid_grow$GAU_fitted,predict(S.75),col="black",pch=".")
points(orchid_grow$GAU_fitted,predict(S.90),col="black",pch=".")
points(orchid_grow$GAU_fitted,predict(S.95),col="black",pch=".")
par(new = TRUE)      
plot(c(orchid_grow$GAU_fitted,
       orchid_grow$GAU_fitted),
     c(Q.skewness(predict(S.10),
                  predict(S.50),
                  predict(S.90)),
       Q.kurtosis(predict(S.05),
                  predict(S.25),
                  predict(S.75),
                  predict(S.95))),
     col=c(rep(alpha("blue",0.25),nrow(orchid_grow)),
           rep(alpha("red",0.25),nrow(orchid_grow))),
     pch=16,cex=.5,axes = FALSE, xlab = "", ylab = "")
abline(h=0,col="lightgray",lty=3)
axis(side = 4,cex.axis=0.8,at = pretty(range(c(Q.skewness(predict(S.10),
                                                          predict(S.50),
                                                          predict(S.90)),
                                               Q.kurtosis(predict(S.05),
                                                          predict(S.25),
                                                          predict(S.75),
                                                          predict(S.95))))))
legend("topleft",legend=c("Skewness","Kurtosis"),bg="white",pch=16,col=c("blue","red"),cex=0.8)

