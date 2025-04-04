## Returning to cactus growth, now estimating skew and kurtosis by quantile regression

### move to the right local directory 
tom = "C:/Users/tm9/Dropbox/github/IPM_size_transitions"
steve = "c:/repos/IPM_size_transitions" 
home = ifelse(Sys.info()["user"] == "Ellner", steve, tom)
setwd(home); setwd("cactus"); 

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
library(wesanderson)

# function for converting cactus size measurements to volume
volume <- function(h, w, p){
  (1/3)*pi*h*(((w + p)/2)/2)^2
}
invlogit<-function(x){exp(x)/(1+exp(x))}

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

## read in cactus demography data
## these data are published on EDI: https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-sev.323.1
CYIM_full<-read_csv("cholla_demography_20042018_EDI.csv")%>% 
  ## filter out transplants and drop seed addition plots (which start with letter H)
  ## also drop Year_t==2018 because I don't have 2019 size data (not entered yet). 2018 data still included in 2017-2018 transition.
  filter(Transplant == 0,
         str_sub(Plot,1,1)!="H",
         Year_t!=2018) %>% 
  ## convert height, max width, perp width to volume of cone, take natural log
  mutate(vol_t = volume(Height_t,Width_t,Perp_t),
         vol_t1 = volume(Height_t1,Width_t1,Perp_t1),
         plot = as.factor(Plot),
         year_t = as.factor(Year_t),
         ID = interaction(TagID,plot)) %>%
  #select(ID,year_t,plot,vol_t,vol_t1,Survival_t1,Goodbuds_t1) %>% 
  ## sort by initial size
  arrange(vol_t) 

## In prelim analysis I inspected several unrealistic size transitions
## this file identifies plants to drop
CYIM_outliers<-read_csv("CYIM_outliers.csv")%>% 
  filter(FLAG==1) %>% dplyr::select(ID) %>% unique()
## drop outliers and create log size variables
CYIM_full %<>% 
  filter(!ID%in%CYIM_outliers$ID) %>% 
  mutate(logvol_t=log(vol_t),
         logvol_t1=log(vol_t1)) 

## pull out and na.omit size transitions for growth modeling
## NAs come from new recruits (missing year_t size) and mortality (missing year_t1 size)
CYIM_full %>% 
  dplyr::select(ID,year_t,plot,logvol_t,logvol_t1) %>% 
  ## drop rows with NAs
  drop_na() -> CYIM_grow

## use gam to fit Gaussian growth model with non-constant variance
CYIM_grow_m1 <- gam(list(logvol_t1 ~ s(logvol_t,k=4) + s(plot,bs="re") + s(year_t,bs="re"), ~s(logvol_t,k=4)), 
                    data=CYIM_grow, family=gaulss())
CYIM_gam_pred <- predict(CYIM_grow_m1,type="response",exclude=c("s(plot)","s(year_t)"))
CYIM_grow$fitted_norfx<-CYIM_gam_pred[,1]
## inspect residuals, scaled by sd; re-run predict now w/RFX
CYIM_grow$fitted_sd<-1/predict(CYIM_grow_m1,type="response")[,2]
CYIM_grow$scaledResids=residuals(CYIM_grow_m1,type="response")/CYIM_grow$fitted_sd

##are the standardized residuals gaussian? -- no
jarque.test(CYIM_grow$scaledResids) # normality test: FAILS, P < 0.001 
anscombe.test(CYIM_grow$scaledResids) # kurtosis: FAILS, P < 0.001 
agostino.test(CYIM_grow$scaledResids) # skewness: FAILS, P<0.001 

## fit qgam -- we will need several quantiles for skewness and kurtosis
S.05<-qgam(scaledResids~s(logvol_t,k=4), data=CYIM_grow,qu=0.05)#,argGam=list(gamma=gamma_param)) 
S.10<-qgam(scaledResids~s(logvol_t,k=4), data=CYIM_grow,qu=0.1)#,argGam=list(gamma=gamma_param)) 
S.25<-qgam(scaledResids~s(logvol_t,k=4), data=CYIM_grow,qu=0.25)#,argGam=list(gamma=gamma_param)) 
S.50<-qgam(scaledResids~s(logvol_t,k=4), data=CYIM_grow,qu=0.5)#,argGam=list(gamma=gamma_param)) 
S.75<-qgam(scaledResids~s(logvol_t,k=4), data=CYIM_grow,qu=0.75)#,argGam=list(gamma=gamma_param)) 
S.90<-qgam(scaledResids~s(logvol_t,k=4), data=CYIM_grow,qu=0.9)#,argGam=list(gamma=gamma_param)) 
S.95<-qgam(scaledResids~s(logvol_t,k=4), data=CYIM_grow,qu=0.95)#,argGam=list(gamma=gamma_param)) 

## NP skewness
q.10<-predict(S.10);q.50<-predict(S.50);q.90<-predict(S.90)
NPS_hat = (q.10 + q.90 - 2*q.50)/(q.90 - q.10)

## NP kurtosis (relative to Gaussian)
q.05<-predict(S.05);q.25<-predict(S.25);q.75<-predict(S.75);q.95<-predict(S.95)
qN = qnorm(c(0.05,0.25,0.75,0.95))
KG = (qN[4]-qN[1])/(qN[3]-qN[2])
NPK_hat = ((q.95-q.05)/(q.75-q.25))/KG - 1

## view diagnostics of scaled residuals
pdf("../manuscript/figures/cactus_qgam_diagnostics.pdf",height = 5, width = 11,useDingbats = F)
par(mfrow=c(1,2),mar = c(5, 5, 2, 3), oma=c(0,0,0,2)) 
plot(CYIM_grow$logvol_t,CYIM_grow$logvol_t1,pch=1,col=alpha("black",0.25),
     xlab="log volume, time t",ylab="log volume, time t+1")
points(CYIM_grow$logvol_t,CYIM_grow$fitted_norfx,col=alpha("red",0.25),pch=16,cex=.5)
par(new = TRUE)                           
plot(CYIM_grow$logvol_t,CYIM_grow$fitted_sd,col=alpha("blue",0.25),pch=16,cex=.5,
     axes = FALSE, xlab = "", ylab = "")
axis(side = 4, at = pretty(range(CYIM_grow$fitted_sd)))
mtext("std dev", side = 4, line = 2)
legend("topleft",legend=c("Fitted mean","Fitted sd"),bg="white",pch=1,col=c("red","blue"),cex=0.8)
title("A",font=3,adj=0)

plot(CYIM_grow$logvol_t,CYIM_grow$scaledResids,col=alpha("black",0.25),
     xlab="log volume, time t",ylab="Scaled residuals of size at t+1")
points(CYIM_grow$logvol_t,q.05,col="black",pch=".")
points(CYIM_grow$logvol_t,q.10,col="black",pch=".")
points(CYIM_grow$logvol_t,q.25,col="black",pch=".")
points(CYIM_grow$logvol_t,q.50,col="black",pch=".")
points(CYIM_grow$logvol_t,q.75,col="black",pch=".")
points(CYIM_grow$logvol_t,q.90,col="black",pch=".")
points(CYIM_grow$logvol_t,q.95,col="black",pch=".")
par(new = TRUE)                           
plot(c(CYIM_grow$logvol_t,CYIM_grow$logvol_t),c(NPS_hat,NPK_hat),
     col=c(rep(alpha("blue",0.25),nrow(CYIM_grow)),rep(alpha("red",0.25),nrow(CYIM_grow))),
           pch=16,cex=.5, axes = FALSE, xlab = "", ylab = "")
abline(h=0,col="lightgray",lty=3)
axis(side = 4, at = pretty(range(c(NPS_hat,NPK_hat))),cex.axis=0.8)
mtext("Skewness", side = 4, line = 2,col="blue")
mtext("Excess Kurtosis", side = 4, line =3,col="red")
title("B",font=3,adj=0)
dev.off()

plot(CYIM_grow$logvol_t,CYIM_grow$fitted)

## collect data objects to re-draw plot
cactus_out <- list(
  cactus_grow = CYIM_grow[,c("logvol_t","logvol_t1","fitted_norfx","fitted_sd","scaledResids")],
  q.05 = q.05,
  q.10 = q.10,
  q.25 = q.25,
  q.50 = q.50,
  q.75 = q.75,
  q.90 = q.90,
  q.95 = q.95,
  NPS_hat = NPS_hat,
  NPK_hat = NPK_hat
)

## now I need to fit a distribution with negative skew and positive excess kurtosis
## both skewness and kurtosis should be non-monotonic wrt size
## turns out mgcv can do this!
CYIM_gam_shash <- gam(list(logvol_t1 ~ s(logvol_t,k=4) + s(plot,bs="re") + s(year_t,bs="re"), # <- model for location 
                           ~ s(logvol_t,k=4),   # <- model for log-scale
                           ~ s(logvol_t,k=4),   # <- model for skewness
                           ~ s(logvol_t,k=4)), # <- model for log-kurtosis
                      data = CYIM_grow, 
                      family = shash,  
                      optimizer = "efs")
CYIM_shash_pred <- predict(CYIM_gam_shash,type="response",exclude=c("s(plot)","s(year_t)"))

## view parameter estimates
par(mfrow=c(2,2))
plot(CYIM_grow$logvol_t,CYIM_grow$logvol_t1,pch=1,col=alpha("black",0.25))
points(CYIM_grow$logvol_t,CYIM_shash_pred[,1],col=alpha("red",0.25),pch=16,cex=.5)
plot(CYIM_grow$logvol_t,exp(CYIM_shash_pred[,2]),col=alpha("red",0.25),pch=16,cex=.5)
plot(CYIM_grow$logvol_t,CYIM_shash_pred[,3],col=alpha("red",0.25),pch=16,cex=.5)
plot(CYIM_grow$logvol_t,exp(CYIM_shash_pred[,4]),col=alpha("red",0.25),pch=16,cex=.5)


## simulate data from fitted models
n_sim<-100
## four sets of simulated data!
NOsim_mean<-NOsim_sd<-NOsim_skew<-NOsim_kurt<-matrix(NA,nrow=nrow(CYIM_grow),ncol=n_sim)
SHASHsim_mean<-SHASHsim_sd<-SHASHsim_skew<-SHASHsim_kurt<-matrix(NA,nrow=nrow(CYIM_grow),ncol=n_sim)

for(i in 1:n_sim){
  cat("################# Simulation number ", i, "\n") 
  ## add this iteration of sim data to real df
  CYIM_grow$logvol_t1.sim.NO <- rnorm(n=nrow(CYIM_grow),
                                             mean=CYIM_gam_pred[,1],
                                             sd=1/CYIM_gam_pred[,2])
  CYIM_grow$logvol_t1.sim.SHASH <- rSHASHo2(n=nrow(CYIM_grow),
                          mu=CYIM_shash_pred[,1],
                          sigma=exp(CYIM_shash_pred[,2]),
                          nu=CYIM_shash_pred[,3],
                          tau=exp(CYIM_shash_pred[,4]))
  ## Qreg on sim data
  q.05.NO<-predict(qgam(logvol_t1.sim.NO~s(logvol_t,k=4),data=CYIM_grow,qu=0.05)) 
  q.10.NO<-predict(qgam(logvol_t1.sim.NO~s(logvol_t,k=4), data=CYIM_grow,qu=0.10)) 
  q.25.NO<-predict(qgam(logvol_t1.sim.NO~s(logvol_t,k=4), data=CYIM_grow,qu=0.25))
  q.50.NO<-predict(qgam(logvol_t1.sim.NO~s(logvol_t,k=4), data=CYIM_grow,qu=0.5))
  q.75.NO<-predict(qgam(logvol_t1.sim.NO~s(logvol_t,k=4), data=CYIM_grow,qu=0.75)) 
  q.90.NO<-predict(qgam(logvol_t1.sim.NO~s(logvol_t,k=4), data=CYIM_grow,qu=0.90))
  q.95.NO<-predict(qgam(logvol_t1.sim.NO~s(logvol_t,k=4), data=CYIM_grow,qu=0.95))
  NOsim_mean[,i]<-Q.mean(q.25.NO,q.50.NO,q.75.NO)
  NOsim_sd[,i]<-Q.sd(q.25.NO,q.75.NO)
  NOsim_skew[,i]<-Q.skewness(q.10.NO,q.50.NO,q.90.NO)
  NOsim_kurt[,i]<-Q.kurtosis(q.05.NO,q.25.NO,q.75.NO,q.95.NO)

  q.05.SHASH<-predict(qgam(logvol_t1.sim.SHASH~s(logvol_t,k=4),data=CYIM_grow,qu=0.05)) 
  q.10.SHASH<-predict(qgam(logvol_t1.sim.SHASH~s(logvol_t,k=4), data=CYIM_grow,qu=0.10)) 
  q.25.SHASH<-predict(qgam(logvol_t1.sim.SHASH~s(logvol_t,k=4), data=CYIM_grow,qu=0.25))
  q.50.SHASH<-predict(qgam(logvol_t1.sim.SHASH~s(logvol_t,k=4), data=CYIM_grow,qu=0.5))
  q.75.SHASH<-predict(qgam(logvol_t1.sim.SHASH~s(logvol_t,k=4), data=CYIM_grow,qu=0.75)) 
  q.90.SHASH<-predict(qgam(logvol_t1.sim.SHASH~s(logvol_t,k=4), data=CYIM_grow,qu=0.90))
  q.95.SHASH<-predict(qgam(logvol_t1.sim.SHASH~s(logvol_t,k=4), data=CYIM_grow,qu=0.95))
  SHASHsim_mean[,i]<-Q.mean(q.25.SHASH,q.50.SHASH,q.75.SHASH)
  SHASHsim_sd[,i]<-Q.sd(q.25.SHASH,q.75.SHASH)
  SHASHsim_skew[,i]<-Q.skewness(q.10.SHASH,q.50.SHASH,q.90.SHASH)
  SHASHsim_kurt[,i]<-Q.kurtosis(q.05.SHASH,q.25.SHASH,q.75.SHASH,q.95.SHASH)
}

## and now the real data
q.05<-predict(qgam(logvol_t1~s(logvol_t,k=4),data=CYIM_grow,qu=0.05)) 
q.10<-predict(qgam(logvol_t1~s(logvol_t,k=4), data=CYIM_grow,qu=0.10)) 
q.25<-predict(qgam(logvol_t1~s(logvol_t,k=4), data=CYIM_grow,qu=0.25)) 
q.50<-predict(qgam(logvol_t1~s(logvol_t,k=4), data=CYIM_grow,qu=0.5))
q.75<-predict(qgam(logvol_t1~s(logvol_t,k=4), data=CYIM_grow,qu=0.75))
q.90<-predict(qgam(logvol_t1~s(logvol_t,k=4), data=CYIM_grow,qu=0.90))
q.95<-predict(qgam(logvol_t1~s(logvol_t,k=4), data=CYIM_grow,qu=0.95))


## NO vs SHASH gam fit
pdf("../manuscript/figures/cactus_SHASH_fit.pdf",height = 8, width = 8, useDingbats = F)
par(mfrow=c(2,2),mar=c(4,4,1,1),cex.lab=1.25,mgp=c(2.1,1,0))
plot(CYIM_grow$logvol_t,Q.mean(q.25,q.50,q.75),type="n",
     xlab="Size(t)",ylab="Mean size(t+1)",ylim=c(min(NOsim_mean),max(SHASHsim_mean)+1))
title(main="A",adj=0,font=3)
  matpoints(CYIM_grow$logvol_t,NOsim_mean,col=alpha("tomato",0.25),type="l",lty=1)
  matpoints(CYIM_grow$logvol_t,1+SHASHsim_mean,col=alpha("cornflowerblue",0.25),type="l",lty=1)
lines(CYIM_grow$logvol_t,Q.mean(q.25,q.50,q.75),col="black",lwd=2)
lines(CYIM_grow$logvol_t,1+Q.mean(q.25,q.50,q.75),col="black",lwd=2)
legend("topleft",legend=c("Real data","Simulated from Gaussian",
                          "Simulated from SHASH + offset"),
       lty=1,col=c("black","tomato","cornflowerblue"),cex=0.9,bty="n")

plot(CYIM_grow$logvol_t,Q.sd(q.25,q.75),type="n",
     xlab="Size(t)",ylab="SD size(t+1)",ylim=c(min(NOsim_sd),1+max(SHASHsim_sd)))
title(main="B",adj=0,font=3)
  matpoints(CYIM_grow$logvol_t,NOsim_sd,col=alpha("tomato",0.25),type="l",lty=1)
  matpoints(CYIM_grow$logvol_t,1+SHASHsim_sd,col=alpha("cornflowerblue",0.25),type="l",lty=1)
lines(CYIM_grow$logvol_t,Q.sd(q.25,q.75),col="black",lwd=2)
lines(CYIM_grow$logvol_t,1+Q.sd(q.25,q.75),col="black",lwd=2)

plot(CYIM_grow$logvol_t,Q.skewness(q.10,q.50,q.90), type="n", 
     xlab="Size(t)",ylab="Skewness size(t+1)",ylim=c(min(NOsim_skew),1+max(SHASHsim_skew)))
title(main="C",adj=0,font=3)
matpoints(CYIM_grow$logvol_t,NOsim_skew,col=alpha("tomato",0.25),type="l",lty=1)
matpoints(CYIM_grow$logvol_t,1+SHASHsim_skew,col=alpha("cornflowerblue",0.25),type="l",lty=1)
lines(CYIM_grow$logvol_t,Q.skewness(q.10,q.50,q.90),col="black",lwd=2)
lines(CYIM_grow$logvol_t,1+Q.skewness(q.10,q.50,q.90),col="black",lwd=2)

plot(CYIM_grow$logvol_t,Q.kurtosis(q.05,q.25,q.75,q.95),type="n",
    xlab="Size(t)",ylab="Kurtosis size(t+1)",ylim=c(min(NOsim_kurt),1+max(SHASHsim_kurt)))
title(main="D",adj=0,font=3)
  matpoints(CYIM_grow$logvol_t,NOsim_kurt,col=alpha("tomato",0.25),type="l",lty=1)
  matpoints(CYIM_grow$logvol_t,1+SHASHsim_kurt,col=alpha("cornflowerblue",0.25),type="l",lty=1)
  lines(CYIM_grow$logvol_t,Q.kurtosis(q.05,q.25,q.75,q.95),col="black",lwd=2)
  lines(CYIM_grow$logvol_t,1+Q.kurtosis(q.05,q.25,q.75,q.95),col="black",lwd=2)
dev.off()


# compare IPM results between Gaussian and SHASH growth kernel ------------
# Here are the size-dependent functions for survival, flowering, and flowerbud production, 
# fit in mgcv with year and plot random effects.
surv_mod <- gam(Survival_t1 ~ s(logvol_t,k=4) + s(plot,bs="re") + s(year_t,bs="re"), family="binomial", data=CYIM_full)
flow_mod <- gam(Goodbuds_t>0 ~ s(logvol_t,k=4) + s(plot,bs="re") + s(year_t,bs="re"), family="binomial", data=CYIM_full)
fert_mod <- gam(Goodbuds_t ~ s(logvol_t,k=4) + s(plot,bs="re") + s(year_t,bs="re"), family="nb", data=subset(CYIM_full,Goodbuds_t>0))

## misc parameters for seeds and seedlings
seeds_per_fruit<-read.csv("JO_fruit_data_final_dropplant0.csv",T)  %>% drop_na() %>% 
  summarise(seeds_per_fruit = mean(seed_count))
seed_survival <- read.csv("FruitSurvival.csv",T) %>% drop_na() %>% mutate(fruit_frac = Fr.on.grnd.not.chewed/Fr.on.plant) %>% 
  summarise(seed_survival = mean(fruit_frac))
germination <- read.csv("Germination.csv") %>% drop_na() %>% 
  mutate(germ1_frac = Seedlings04/Input,
         germ2_frac = Seedlings05/(Input-Seedlings04)) %>% 
  summarise(germ1 = mean(germ1_frac), germ2 = mean(germ2_frac))
precensus_survival <- read.csv("PrecensusSurvival.csv") %>% dplyr::select(survive0405) %>% drop_na() %>% 
  summarise(precensus_survival = mean(survive0405))
seedling_size <- read_csv("cholla_demography_20042018_EDI.csv") %>% 
  ## subset seed germination plots
  filter(str_sub(Plot,1,1)=="H") %>% 
  mutate(vol_t = volume(Height_t,Width_t,Perp_t)) %>% 
  summarise(mean_size = mean(log(vol_t),na.rm=T),
            sd_size = sd(log(vol_t),na.rm=T))
## size bounds
min.size <- log(min(CYIM_full$vol_t,na.rm=T)) 
max.size <- log(max(CYIM_full$vol_t1,na.rm=T))

## IPM source functions
## SURVIVAL
sx<-function(x,exclude,year=2004,plot=1){
  xb=pmin(pmax(x,min.size),max.size)
  pred=predict(surv_mod,newdata = data.frame(logvol_t=xb,plot=plot,year_t=year),
                exclude=paste0(exclude))
  return(invlogit(pred))
}
## GROWTH - SHASH
gxy_SHASH<-function(x,y,exclude,year=2004,plot=1){
  xb=pmin(pmax(x,min.size),max.size) #Transforms all values below/above limits in min/max size
  pred=predict(CYIM_gam_shash,
               newdata = data.frame(logvol_t=xb,plot=plot,year_t=year),
               exclude=paste0(exclude))
  return(dSHASHo2(x=y, 
                  mu=pred[,1],
                  sigma = exp(pred[,2]), 
                  nu = pred[,3], 
                  tau = exp(pred[,4])))
}
## GROWTH - Gaussian
gxy_GAU<-function(x,y,exclude,year=2004,plot=1){
  xb=pmin(pmax(x,min.size),max.size) #Transforms all values below/above limits in min/max size
  pred = predict(CYIM_grow_m1,
                 newdata = data.frame(logvol_t=xb,plot=plot,year_t=year),
                 exclude=paste0(exclude))
  return(dnorm(y,mean=pred[,1],sd=exp(pred[,2])))
}
## COMBINED GROWTH_SURVIVAL
pxy <- function(x,y,dist,exclude,year=2004,plot=1){
  result <- sx(x,exclude,year,plot)*do.call(paste0("gxy_",dist),list(x,y,exclude,year,plot))
  return(result)
}
#PR FLOWERING
flow.x <- function(x,exclude,year=2004,plot=1){
  xb=pmin(pmax(x,min.size),max.size)
  pred=predict(flow_mod,
               newdata = data.frame(logvol_t=xb,plot=plot,year_t=year),
               exclude=paste0(exclude))
  return(invlogit(pred))
}
##FLOWERBUD PRODUCTION BY FLOWERING PLANTS
fert.x <- function(x,exclude,year=2004,plot=1){
  xb=pmin(pmax(x,min.size),max.size)
  pred=predict(fert_mod,
               newdata = data.frame(logvol_t=xb,plot=plot,year_t=year),
               exclude=paste0(exclude))
  return(exp(pred))
}
##SEED BANK CONTRIBUTION OF X-SIZED PLANTS
fx<-function(x,exclude,year=2004,plot=1){
  return(flow.x(x,exclude,year,plot)*fert.x(x,exclude,year,plot)*seeds_per_fruit$seeds_per_fruit*seed_survival$seed_survival)  
}
#SIZE DISTRIBUTION OF RECRUITS
recruit.size<-function(y){
  dnorm(x=y,mean=seedling_size$mean_size,sd=seedling_size$sd_size)
}

#PUT IT ALL TOGETHER
bigmatrix<-function(lower.extension = 0,upper.extension = 0,
                    mat.size,dist,exclude,year=2004,plot=1){
  n<-mat.size
  L<-min.size + lower.extension
  U<-max.size + upper.extension
  #these are the upper and lower integration limits
  h<-(U-L)/n                   #Bin size
  b<-L+c(0:n)*h;               #Lower boundaries of bins 
  y<-0.5*(b[1:n]+b[2:(n+1)]);  #Bins' midpoints
  
  # Fertility matrix
  Fmat<-matrix(0,(n+2),(n+2))
  
  # 1-yo banked seeds go in top row
  Fmat[1,3:(n+2)]<-fx(y,exclude,year,plot)
  
  # Growth/survival transition matrix
  Tmat<-matrix(0,(n+2),(n+2))
  
  # Graduation to 2-yo seed bank = pr(not germinating as 1-yo)
  Tmat[2,1]<-(1-germination$germ1)
  
  # Graduation from 1-yo bank to cts size = germination * size distn * pre-census survival
  Tmat[3:(n+2),1]<- germination$germ1 * precensus_survival$precensus_survival * recruit.size(y) * h   
  
  # Graduation from 2-yo bank to cts size = germination * size distn * pre-census survival
  Tmat[3:(n+2),2]<- germination$germ2 * precensus_survival$precensus_survival * recruit.size(y) * h  
  
  # Growth/survival transitions among cts sizes
  Tmat[3:(n+2),3:(n+2)]<-t(outer(y,y,pxy,dist=dist,exclude=exclude,year=year,plot=plot)) * h 
  
  # Put it all together
  IPMmat<-Fmat+Tmat     #Full Kernel is simply a summation ot fertility
  #and transition matrix
  return(list(IPMmat=IPMmat,Fmat=Fmat,Tmat=Tmat,meshpts=y))
}
# lambdaS function##########################################################
lambdaSim<-function(mat_list, ## a list of transition matrices, each corresponding to a study year
                    max_yrs=1000, ## how many years the simulation runs (arbitrarily large)
                    seed=NULL){
  ## grab the dimension of the projection matrix
  matdim<-dim(mat_list[[1]])[1]
  ## grab the number of study years / matrices we have available
  n_years <- length(mat_list)
  ## vector that will hold year-by-year growth rates
  rtracker <- rep(0,max_yrs)
  ## initial vector of population structure -- note that this sums to one, which will be convenient
  n0 <- rep(1/matdim,matdim)
  for(t in 1:max_yrs){ #Start loop
    ## for each year, randomly sample one of the matrices
    if(!is.null(seed)){set.seed(seed[t])}
    A_t <- mat_list[[sample.int(n=n_years,size=1)]]
    ## project the population one step forward
    n0 <- A_t %*% n0
    ## total population size after one year of growth
    N  <- sum(n0)
    ## calculate r as log(N_t+1 / N_t), note that here N_t=1
    rtracker[t]<-log(N)
    ## rescale population vector to sum to one, so the same trick works again next time step
    n0 <-n0/N
  }
  #discard first 10% of time series
  burnin    <- round(max_yrs*0.1)
  #Finish and return
  log_lambdaS <- mean(rtracker[-c(1:burnin)])
  lambdaS<-exp(log_lambdaS)
  return(list(log_lambdaS=log_lambdaS,lambdaS=lambdaS))
}

## create an array of projection kernels for each year
mat.size = 200
lower.extension = -1
upper.extension = 1.5

## store kernel outputs for average plot and year
exclude_plot_year<-c("s(plot)","s(year_t)")
cactus_out$mat_GAU<-bigmatrix(lower.extension = lower.extension, 
                              upper.extension = upper.extension,
                              mat.size = mat.size,exclude=exclude_plot_year,
                              dist="GAU")
cactus_out$mat_SHASH<-bigmatrix(lower.extension = lower.extension, 
                              upper.extension = upper.extension,
                              mat.size = mat.size,exclude=exclude_plot_year,
                              dist="SHASH")
write_rds(cactus_out,file="cactus_out.rds")


## store year-specific K's and lambda's
studyyears<-sort(unique(CYIM_full$Year_t))
K_t_SHASH<-K_t_GAU<-vector("list",length=length(studyyears))
lambda_t_SHASH<-lambda_t_GAU<-vector("numeric",length=length(studyyears))
## loop over years for an average plot
exclude_plot<-exclude_plot_year[1]
for(t in 1:length(studyyears)){
  K_t_SHASH[[t]]<-bigmatrix(lower.extension = lower.extension, 
                      upper.extension = upper.extension,
                      mat.size = mat.size,exclude=exclude_plot,
                      dist="SHASH",year=studyyears[t])$IPMmat
  lambda_t_SHASH[t]<-lambda(K_t_SHASH[[t]])
  K_t_GAU[[t]]<-bigmatrix(lower.extension = lower.extension, 
                            upper.extension = upper.extension,
                            mat.size = mat.size,exclude=exclude_plot,
                            dist="GAU",year=studyyears[t])$IPMmat
  lambda_t_GAU[t]<-lambda(K_t_GAU[[t]])
}
sd(lambda_t_GAU);sd(lambda_t_SHASH)

## now the same for plots
## store plot-specific K's and lambda's
studyplots<-unique(CYIM_full$plot)
K_p_SHASH<-K_p_GAU<-vector("list",length=length(studyplots))
lambda_p_SHASH<-lambda_p_GAU<-vector("numeric",length=length(studyplots))
## loop over plots for an average year
exclude_year<-exclude_plot_year[2]
for(p in 1:length(studyplots)){
  K_p_SHASH[[p]]<-bigmatrix(lower.extension = lower.extension, 
                            upper.extension = upper.extension,
                            mat.size = mat.size,exclude=exclude_year,
                            dist="SHASH",plot=studyplots[p])$IPMmat
  lambda_p_SHASH[p]<-lambda(K_p_SHASH[[p]])
  K_p_GAU[[p]]<-bigmatrix(lower.extension = lower.extension, 
                          upper.extension = upper.extension,
                          mat.size = mat.size,exclude=exclude_year,
                          dist="GAU",plot=studyplots[p])$IPMmat
  lambda_p_GAU[p]<-lambda(K_p_GAU[[p]])
}
sd(lambda_p_GAU);sd(lambda_p_SHASH)

## compare plot-specific asypmtotic growth rates
## compare year-specific asypmtotic growth rates
cols<-wes_palette("Zissou1")
pdf("../manuscript/figures/cactus_lambda_years_plots.pdf",height = 4, width = 8,useDingbats = F)
par(mfrow=c(1,2),mar=c(4,6,1,1))
plot(studyyears,lambda_t_SHASH,type="b",col=alpha(cols[1],0.85),pch=16,
     xlab="Year",ylab=expression(paste(lambda[t])),cex=1.5,cex.lab=1.8)
lines(studyyears,lambda_t_GAU,type="b",col=alpha(cols[5],0.85),pch=16,cex=1.5)
title("A",adj=0,font=3)
legend("bottomright",legend=c("SHASH","Gaussian"),title="Growth function:",
       pch=16,col=cols[c(1,5)])

plot(as.numeric(studyplots),lambda_p_SHASH[as.numeric(studyplots)],type="n",cex.lab=1.8,
    ylim=range(lambda_t_SHASH),xlab="Plot",ylab=expression(paste(lambda[p])))
axis(1,at=1:11,labels=1:11)
points(as.numeric(studyplots),lambda_p_SHASH[as.numeric(studyplots)],col=alpha(cols[1],0.85),pch=16,type="p",cex=1.5)
points(as.numeric(studyplots),lambda_p_GAU[as.numeric(studyplots)],col=alpha(cols[5],0.85),pch=16,type="p",cex=1.5)
title("B",adj=0,font=3)
box()
dev.off()


## compare stochastic growth rates for average plot
## for fair comparison, run them through the same environments
max_yrs<-50000
seed<-sample.int(max_yrs)
lambdaSim(mat_list = K_t_SHASH,max_yrs=max_yrs,seed=seed)$lambdaS
lambdaSim(mat_list = K_t_GAU,max_yrs=max_yrs,seed=seed)$lambdaS
## effectively identical, probably because the growth difference is swamped by 
## size structure transients

##lastly, life table response experiment for both plot and year variation
## make a data frame of lambda_t and time-varying vital rate parameters
## first, find the rfx in the parameter vectors
surv_years<-which(names(surv_mod$coefficients)=="s(year_t).1"):which(names(surv_mod$coefficients)=="s(year_t).13")
growGAU_years<-which(names(CYIM_grow_m1$coefficients)=="s(year_t).1"):which(names(CYIM_grow_m1$coefficients)=="s(year_t).13")
growSHASH_years<-which(names(CYIM_gam_shash$coefficients)=="s(year_t).1"):which(names(CYIM_gam_shash$coefficients)=="s(year_t).13")
flow_years<-which(names(flow_mod$coefficients)=="s(year_t).1"):which(names(flow_mod$coefficients)=="s(year_t).13")
fert_years<-which(names(fert_mod$coefficients)=="s(year_t).1"):which(names(fert_mod$coefficients)=="s(year_t).13")

year_ltre_df<-data.frame(
  lambda_t_GAU=lambda_t_GAU,
  lambda_t_SHASH=lambda_t_SHASH,
  surv_t=surv_mod$coefficients[surv_years],
  grow_t_GAU=CYIM_grow_m1$coefficients[growGAU_years],
  grow_t_SHASH=CYIM_gam_shash$coefficients[growSHASH_years],
  flow_t=flow_mod$coefficients[flow_years],
  fert_t=fert_mod$coefficients[fert_years])

#fit linear models for years
fitGAU_year<-lm(lambda_t_GAU~surv_t+grow_t_GAU+flow_t+fert_t,data=year_ltre_df)
fitSHASH_year<-lm(lambda_t_SHASH~surv_t+grow_t_SHASH+flow_t+fert_t,data=year_ltre_df)
#use th fitted slopes to approximate a sensitivity matrix
sensGAU_year<-outer(coef(fitGAU_year)[2:5],coef(fitGAU_year)[2:5],FUN="*")
sensSHASH_year<-outer(coef(fitSHASH_year)[2:5],coef(fitSHASH_year)[2:5],FUN="*")
#VcoV of vital rates
paramcovGAU_year<-cov(year_ltre_df[,c("surv_t","grow_t_GAU","flow_t","fert_t")])
paramcovSHASH_year<-cov(year_ltre_df[,c("surv_t","grow_t_SHASH","flow_t","fert_t")])
#Calculate 1st order approx
ltre_approxGAU_year <-sensGAU_year*paramcovGAU_year
ltre_approxSHASH_year <-sensSHASH_year*paramcovSHASH_year
#compare var(almbda) to the sum of these terms
var(year_ltre_df$lambda_t_GAU);sum(ltre_approxGAU_year)#pretty good!
var(year_ltre_df$lambda_t_SHASH);sum(ltre_approxSHASH_year)
#calculate contributions (thank you Mark and Steve)
(ltre_contGAU_year <- sort(apply(ltre_approxGAU_year,1,sum)/sum(ltre_approxGAU_year)))
(ltre_contSHASH_year <- sort(apply(ltre_approxSHASH_year,1,sum)/sum(ltre_approxSHASH_year)))

barplot(rbind(ltre_contGAU_year,ltre_contSHASH_year),beside=T)

##same now for plots
surv_plots<-which(names(surv_mod$coefficients)=="s(plot).1"):which(names(surv_mod$coefficients)=="s(plot).11")
growGAU_plots<-which(names(CYIM_grow_m1$coefficients)=="s(plot).1"):which(names(CYIM_grow_m1$coefficients)=="s(plot).11")
growSHASH_plots<-which(names(CYIM_gam_shash$coefficients)=="s(plot).1"):which(names(CYIM_gam_shash$coefficients)=="s(plot).11")
flow_plots<-which(names(flow_mod$coefficients)=="s(plot).1"):which(names(flow_mod$coefficients)=="s(plot).11")
fert_plots<-which(names(fert_mod$coefficients)=="s(plot).1"):which(names(fert_mod$coefficients)=="s(plot).11")

plot_ltre_df<-data.frame(
  lambda_p_GAU=lambda_p_GAU,
  lambda_p_SHASH=lambda_p_SHASH,
  surv_p=surv_mod$coefficients[surv_plots],
  grow_p_GAU=CYIM_grow_m1$coefficients[growGAU_plots],
  grow_p_SHASH=CYIM_gam_shash$coefficients[growSHASH_plots],
  flow_p=flow_mod$coefficients[flow_plots],
  fert_p=fert_mod$coefficients[fert_plots])

#fit linear models for plots
fitGAU_plot<-lm(lambda_p_GAU~surv_p+grow_p_GAU+flow_p+fert_p,data=plot_ltre_df)
fitSHASH_plot<-lm(lambda_p_SHASH~surv_p+grow_p_SHASH+flow_p+fert_p,data=plot_ltre_df)
#use th fitted slopes to approximate a sensitivity matrix
sensGAU_plot<-outer(coef(fitGAU_plot)[2:5],coef(fitGAU_plot)[2:5],FUN="*")
sensSHASH_plot<-outer(coef(fitSHASH_plot)[2:5],coef(fitSHASH_plot)[2:5],FUN="*")
#VcoV of vital rates
paramcovGAU_plot<-cov(plot_ltre_df[,c("surv_p","grow_p_GAU","flow_p","fert_p")])
paramcovSHASH_plot<-cov(plot_ltre_df[,c("surv_p","grow_p_SHASH","flow_p","fert_p")])
#Calculate 1st order approx
ltre_approxGAU_plot <-sensGAU_plot*paramcovGAU_plot
ltre_approxSHASH_plot <-sensSHASH_plot*paramcovSHASH_plot
#compare var(almbda) to the sum of these terms
var(plot_ltre_df$lambda_p_GAU);sum(ltre_approxGAU_plot)#pretty good!
var(plot_ltre_df$lambda_p_SHASH);sum(ltre_approxSHASH_plot)
#calculate contributions (thank you Mark and Steve)
(ltre_contGAU_plot <- sort(apply(ltre_approxGAU_plot,1,sum)/sum(ltre_approxGAU_plot)))
(ltre_contSHASH_plot <- sort(apply(ltre_approxSHASH_plot,1,sum)/sum(ltre_approxSHASH_plot)))

barplot(rbind(ltre_contGAU_plot,ltre_contSHASH_plot),beside=T)
## not sure I believe this -- the approximation does not work well
## and there is not much variance in the first place

# the basement ------------------------------------------------------------

## test that I understand how to get sigma from the gaulss fit
x <- runif(500,-1,1)
y <- rnorm(500,mean=8.6+3*x,sd=exp(0.5+0.8*x))
plot(x,y)
testgam<-gam(list(y~s(x),~s(x)),family=gaulss())
pred.response <- predict(testgam,type="response")
plot(exp(0.5+0.8*x),1/pred.response[,2]);abline(0,1)

pred.lpmat <- predict(testgam,type="lpmatrix")
grow_sd_index <- which(as.factor(names(coef(testgam)))=="(Intercept).1") ## this is where the sd coefficients start
gam_coef_length <- length(coef(testgam))
plot(exp(0.5+0.8*x),
exp(pred.lpmat[, grow_sd_index:length(coef(testgam))] %*% coef(testgam)[grow_sd_index:length(coef(testgam))]));abline(0,1)


## now curious to look at the outliers
CYIM_grow %>% 
  filter(scaledResids < quantile(scaledResids,probs=0.025) |
           scaledResids >  quantile(scaledResids,probs=0.975) ) %>% 
  select(ID,year_t) %>% 
  left_join(.,CYIM_full,
            by=c("ID","year_t"))-> outliers
write_csv(outliers,"cactus/CYIM_outliers.csv")

separate(as.character(CYIM_grow$ID))

df <- data.frame(x = c(NA, "x.y", "x.z", "y.z"))
df %>% separate(x, c("A", "B"))

CYIM_test<-read_csv("cactus/cholla_demography_20042018_EDI.csv")
CYIM_test %>% 
  filter(Plot=="3",TagID=="45") %>% View()

# And finally for IPM results. The mean model (averaging across years and plots) predicts very different growth rates:
tibble(Growth_Dist = c("SHASH","Gaussian"),
       lambda = c(round(lambda(kernel_SHASH$IPMmat),4),round(lambda(kernel_GAU$IPMmat),4)))

#Also, the SSD's are very different and the observed distribution is maybe in between 
# the two of them (in the data plot, colored bars are years). The Gaussian ssd is bi-modal but the 
# lower mode disappears with the skewed $t$ growth kernel, probably because the big reproductive plants 
# are rare in the ssd so you lose the signal of recruitment. 
ssd_SHASH <- stable.stage(kernel_SHASH$IPMmat)[3:(mat.size+2)] / sum(stable.stage(kernel_SHASH$IPMmat)[3:(mat.size+2)])
ssd_norm <- stable.stage(kernel_GAU$IPMmat)[3:(mat.size+2)] / sum(stable.stage(kernel_GAU$IPMmat)[3:(mat.size+2)])
empirical_sd <- density(na.omit(log(CYIM_full$vol_t1)),n=mat.size)

pdf("manuscript/figures/cactus_ssd.pdf",height = 6,width = 6,useDingbats = F)
plot(kernel_GAU$meshpts,ssd_norm,type="l",lwd=3,
     xlab="log volume",ylab="Density",col=alpha("blue",0.5))
lines(kernel_SHASH$meshpts,ssd_SHASH,type="l",lwd=3,col=alpha("red",0.5))
lines(empirical_sd$x,empirical_sd$y/sum(empirical_sd$y),lwd=2,lty=3)
text(0,0.015,expression(paste(lambda[GAU],"=0.98788")),col=alpha("blue",0.5))
text(0,0.01,expression(paste(lambda[SHASH],"=0.98818")),col=alpha("red",0.5))
legend("topleft",c("Gaussian SSD","SHASH SSD", "Empirical SD"),lty=c(1,1,2),
       col=c(alpha("blue",0.5),alpha("red",0.5),"black"),lwd=3,bty="n")
dev.off()

