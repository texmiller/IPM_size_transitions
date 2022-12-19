## Returning to cactus growth, now estimating skew and kurtosis by quantile regression

## setwd
setwd("C:/Users/tm9/Dropbox/github/IPM_size_transitions")

## load libraries
library(tidyverse)
library(mgcv)
library(scales)
library(qgam)
library(gamlss.dist)

# function for converting cactus size measurements to volume
volume <- function(h, w, p){
  (1/3)*pi*h*(((w + p)/2)/2)^2
}

## quartile-based estimates of mean and sd
## see https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/1471-2288-14-135
Q.mean<-function(q.25,q.50,q.75){(q.25+q.50+q.75)/3}
Q.sd<-function(q.25,q.75){(q.75-q.25)/1.35}
Q.skewness<-function(q.10,q.50,q.90){(q.10 + q.90 - 2*q.50)/(q.90 - q.10)}
Q.kurtosis<-function(q.05,q.25,q.75,q.95){
  qN = qnorm(c(0.05,0.25,0.75,0.95))
  KG = (qN[4]-qN[1])/(qN[3]-qN[2])
  return(((q.95-q.05)/(q.75-q.25))/KG - 1)
}

## read in cactus demography data
## these data are published on EDI: https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-sev.323.1
CYIM_full<-read_csv("cactus/cholla_demography_20042018_EDI.csv")%>% 
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
CYIM_outliers<-read_csv("cactus/CYIM_outliers.csv") %>% 
  filter(FLAG==1) %>% select(ID) %>% unique()

## pull out and na.omit size transitions for growth modeling
## NAs come from new recruits (missing year_t size) and mortality (missing year_t1 size)
CYIM_full %>% 
  filter(!ID%in%CYIM_outliers$ID) %>% 
  mutate(logvol_t=log(vol_t),
         logvol_t1=log(vol_t1)) %>% 
  select(ID,year_t,plot,logvol_t,logvol_t1) %>% 
  ## drop rows with NAs
  drop_na() -> CYIM_grow

## use gam to fit Gaussian growth model with non-constant variance
CYIM_grow_m1 <- gam(list(logvol_t1 ~ s(logvol_t,k=4) + s(plot,bs="re") + s(year_t,bs="re"), ~s(logvol_t,k=4)), 
                    data=CYIM_grow, family=gaulss())
CYIM_gam_pred <- predict(CYIM_grow_m1,type="response",exclude=c("s(plot)","s(year_t)"))

## inspect residuals, scaled by sd; re-run predict now w/RFX
fitted_sd<-1/predict(CYIM_grow_m1,type="response")[,2]
CYIM_grow$scaledResids=residuals(CYIM_grow_m1,type="response")/fitted_sd

## fit qgam -- we will need several quantiles for skewness and kurtosis
gamma_param<-2
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
plot(CYIM_grow$logvol_t,CYIM_grow$logvol_t1,pch=1,col=alpha("black",0.25),
     xlab="size t",ylab="size t1")
points(CYIM_grow$logvol_t,CYIM_gam_pred[,1],col=alpha("red",0.25),pch=16,cex=.5)
par(new = TRUE)                           
plot(CYIM_grow$logvol_t,1/CYIM_gam_pred[,2],col=alpha("blue",0.25),pch=16,cex=.5,
     axes = FALSE, xlab = "", ylab = "")
axis(side = 4, at = pretty(range(1/CYIM_gam_pred[,2])),col="blue")
mtext("sigma", side = 4, line = 2,col="blue")

pdf("./manuscript/figures/cactus_qgam_diagnostics.pdf",height = 6, width = 8,useDingbats = F)
par(mar = c(5, 4, 2, 3), oma=c(0,0,0,4)) 
plot(CYIM_grow$logvol_t,CYIM_grow$scaledResids,col=alpha("black",0.25),
     xlab="Size at time t",ylab="Scaled residuals of size at t+1")
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
dev.off()

plot(CYIM_grow$logvol_t,NPK_hat,col=alpha("blue",0.25),pch=16,cex=.5)

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

## simulate data from fitted model and compare to real data
n_sim<-100
sim_mean<-sim_sd<-sim_skew<-sim_kurt<-matrix(NA,nrow=nrow(CYIM_grow),ncol=n_sim)
for(i in 1:n_sim){
  ## add this iteration of sim data to real df
  CYIM_grow$logvol_t1.sim <- rSHASHo2(n=nrow(CYIM_grow),
                          mu=CYIM_shash_pred[,1],
                          sigma=exp(CYIM_shash_pred[,2]),
                          nu=CYIM_shash_pred[,3],
                          tau=exp(CYIM_shash_pred[,4]))
  ## Qreg on sim data
  q.05<-predict(qgam(logvol_t1.sim~s(logvol_t,k=4),data=CYIM_grow,qu=0.05)) 
  q.10<-predict(qgam(logvol_t1.sim~s(logvol_t,k=4), data=CYIM_grow,qu=0.10)) 
  q.25<-predict(qgam(logvol_t1.sim~s(logvol_t,k=4), data=CYIM_grow,qu=0.25))
  q.50<-predict(qgam(logvol_t1.sim~s(logvol_t,k=4), data=CYIM_grow,qu=0.5))
  q.75<-predict(qgam(logvol_t1.sim~s(logvol_t,k=4), data=CYIM_grow,qu=0.75)) 
  q.90<-predict(qgam(logvol_t1.sim~s(logvol_t,k=4), data=CYIM_grow,qu=0.90))
  q.95<-predict(qgam(logvol_t1.sim~s(logvol_t,k=4), data=CYIM_grow,qu=0.95))
  
  sim_mean[,i]<-Q.mean(q.25,q.50,q.75)
  sim_sd[,i]<-Q.sd(q.25,q.75)
  sim_skew[,i]<-Q.skewness(q.10,q.50,q.90)
  sim_kurt[,i]<-Q.kurtosis(q.05,q.25,q.75,q.95)
}

## and now the real data
q.05<-predict(qgam(logvol_t1~s(logvol_t,k=4),data=CYIM_grow,qu=0.05)) 
q.10<-predict(qgam(logvol_t1~s(logvol_t,k=4), data=CYIM_grow,qu=0.10)) 
q.25<-predict(qgam(logvol_t1~s(logvol_t,k=4), data=CYIM_grow,qu=0.25)) 
q.50<-predict(qgam(logvol_t1~s(logvol_t,k=4), data=CYIM_grow,qu=0.5))
q.75<-predict(qgam(logvol_t1~s(logvol_t,k=4), data=CYIM_grow,qu=0.75))
q.90<-predict(qgam(logvol_t1~s(logvol_t,k=4), data=CYIM_grow,qu=0.90))
q.95<-predict(qgam(logvol_t1~s(logvol_t,k=4), data=CYIM_grow,qu=0.95))

pdf("./manuscript/figures/cactus_SHASH_fit.pdf",height = 6, width = 6,useDingbats = F)
par(mfrow=c(2,2),mar=c(4,4,1,1))
plot(CYIM_grow$logvol_t,Q.mean(q.25,q.50,q.75),type="n",
     xlab="size t",ylab="mean size t1",ylim=c(min(sim_mean),max(sim_mean)))
for(i in 1:n_sim){
  points(CYIM_grow$logvol_t,sim_mean[,i],col=alpha("black",0.25))
}
points(CYIM_grow$logvol_t,Q.mean(q.25,q.50,q.75),col="red",lwd=2)
legend("topleft",legend=c("Real data","Simulated from \nfitted SHASH gam"),
       lty=1,col=c("red","black"),lwd=c(2,1),cex=0.8,bty="n")

plot(CYIM_grow$logvol_t,Q.sd(q.25,q.75),type="n",
     xlab="size t",ylab="sd size t1",ylim=c(min(sim_sd),max(sim_sd)))
for(i in 1:n_sim){
  points(CYIM_grow$logvol_t,sim_sd[,i],col=alpha("black",0.25))
}
points(CYIM_grow$logvol_t,Q.sd(q.25,q.75),col="red",lwd=2)

plot(CYIM_grow$logvol_t,Q.skewness(q.10,q.50,q.90),type="n",
     xlab="size t",ylab="skewness size t1",ylim=c(min(sim_skew),max(sim_skew)))
for(i in 1:n_sim){
  points(CYIM_grow$logvol_t,sim_skew[,i],col=alpha("black",0.25))
}
points(CYIM_grow$logvol_t,Q.skewness(q.10,q.50,q.90),col="red",lwd=2)

plot(CYIM_grow$logvol_t,Q.kurtosis(q.05,q.25,q.75,q.95),type="n",
     xlab="size t",ylab="kurtosis size t1",ylim=c(min(sim_kurt),max(sim_kurt)))
for(i in 1:n_sim){
  points(CYIM_grow$logvol_t,sim_kurt[,i],col=alpha("black",0.25))
}
points(CYIM_grow$logvol_t,Q.kurtosis(q.05,q.25,q.75,q.95),col="red",lwd=2)

dev.off()


# compare IPM results between Gaussian and SHASH growth kernel ------------

# Here are the size-dependent functions for survival, flowering, and flowerbud production, 
# fit in lme4 with year and plot random effects.
surv_mod <- glmer(Survival_t1 ~ log(vol_t) + (1|year_t) + (1|plot), family="binomial", data=CYIM_full)
flow_mod <- glmer(Goodbuds_t1>0 ~ log(vol_t1) + (1|year_t) + (1|plot), family="binomial", data=CYIM_full)
fert_mod <- glmer(Goodbuds_t1 ~ log(vol_t1) + (1|year_t) + (1|plot), family="poisson", data=subset(CYIM_full,Goodbuds_t1>0))

par(mfrow=c(2,2),mar=c(4,4,2,1))
plot(log(CYIM_full$vol_t),CYIM_full$Survival_t1,xlab="log Size_t",ylab="Survival",col=alpha("gray",0.5))
lines(size_dum,invlogit(fixef(surv_mod)[1]+fixef(surv_mod)[2]*size_dum),lwd=3)
plot(log(CYIM_full$vol_t),log(CYIM_full$vol_t1),xlab="log Size_t",ylab="log Size_t+1",col=alpha("gray",0.5))
lines(size_dum,mean(coefs[c(years,plots)]) + coefs["log(vol_t)"] * size_dum + coefs["I(log(vol_t)^2)"] * size_dum^2,lwd=3)
lines(size_dum,fixef(CYIM_lmer_best)[1]+fixef(CYIM_lmer_best)[2]*size_dum+fixef(CYIM_lmer_best)[3]*size_dum^2,col="red",lwd=3)
legend("topleft",legend=c("SHASH location","lmer mean"),lwd=c(2,1),col=c("black","red"))
plot(log(CYIM_full$vol_t),CYIM_full$Goodbuds_t1>0,xlab="log Size_t",ylab="Flowering",col=alpha("gray",0.5))
lines(size_dum,invlogit(fixef(flow_mod)[1]+fixef(flow_mod)[2]*size_dum),lwd=3)
plot(log(CYIM_full$vol_t[CYIM_full$Goodbuds_t1>0]),CYIM_full$Goodbuds_t1[CYIM_full$Goodbuds_t1>0],xlab="log Size_t",ylab="Fertility",col=alpha("gray",0.5))
lines(size_dum,exp(fixef(fert_mod)[1]+fixef(fert_mod)[2]*size_dum),lwd=3)

source("cactus_IPM_source_fns.R")

seeds_per_fruit<-read.csv("JO_fruit_data_final_dropplant0.csv",T)  %>% drop_na() %>% 
  summarise(seeds_per_fruit = mean(seed_count))
seed_survival <- read.csv("FruitSurvival.csv",T) %>% drop_na() %>% mutate(fruit_frac = Fr.on.grnd.not.chewed/Fr.on.plant) %>% 
  summarise(seed_survival = mean(fruit_frac))
germination <- read.csv("Germination.csv") %>% drop_na() %>% 
  mutate(germ1_frac = Seedlings04/Input,
         germ2_frac = Seedlings05/(Input-Seedlings04)) %>% 
  summarise(germ1 = mean(germ1_frac), germ2 = mean(germ2_frac))
precensus_survival <- read.csv("PrecensusSurvival.csv") %>% select(survive0405) %>% drop_na() %>% 
  summarise(precensus_survival = mean(survive0405))
seedling_size <- read_csv("cholla_demography_20042018_EDI.csv") %>% 
  ## subset seed germination plots
  filter(str_sub(Plot,1,1)=="H") %>% 
  mutate(vol_t = volume(Height_t,Width_t,Perp_t)) %>% 
  summarise(mean_size = mean(log(vol_t),na.rm=T),
            sd_size = sd(log(vol_t),na.rm=T))

cactus_params<-c()
## skewed t growth params
## intecept is the mean across years and plots, corresponds to lme4 coefs for other vital rates
cactus_params$grow.mu <- mean(coefs[c(years,plots)]) 
cactus_params$grow.bsize <- coefs["log(vol_t)"]
cactus_params$grow.bsize2 <- coefs["I(log(vol_t)^2)"]
cactus_params$sigma_b0 <- coefs["sigma_b0"]
cactus_params$sigma_b1 <- coefs["sigma_b1"]
cactus_params$sigma_b2 <- coefs["sigma_b2"]
cactus_params$nu_b0 <- coefs["nu_b0"]
cactus_params$nu_b1 <- coefs["nu_b1"]
cactus_params$nu_b2 <- coefs["nu_b2"]
cactus_params$tau_b0 <- coefs["tau_b0"]
cactus_params$tau_b1 <- coefs["tau_b1"]
cactus_params$tau_b2 <- coefs["tau_b2"]
## Gaussian growth from best lme4 model
cactus_params$grow.mu.norm <- fixef(CYIM_lmer_best)[1]
cactus_params$grow.bsize.norm <- fixef(CYIM_lmer_best)[2]
cactus_params$grow.bsize2.norm <- fixef(CYIM_lmer_best)[3]
cactus_params$grow.sd.b0 <- pars[[best_model]][1]
cactus_params$grow.sd.b1 <- pars[[best_model]][2]
cactus_params$grow.sd.b2 <- pars[[best_model]][3]
## survival
cactus_params$surv.mu <- fixef(surv_mod)[1]
cactus_params$surv.bsize <- fixef(surv_mod)[2]
## flowering
cactus_params$flow.mu <- fixef(flow_mod)[1]
cactus_params$flow.bsize <- fixef(flow_mod)[2]
## fruit production
cactus_params$fert.mu <- fixef(fert_mod)[1]
cactus_params$fert.bsize <- fixef(fert_mod)[2]
## seeds per fruit
cactus_params$mu_spf <- seeds_per_fruit$seeds_per_fruit
## seed survival
cactus_params$seedsurv <- seed_survival$seed_survival
## germination rates
cactus_params$germ1 <- germination$germ1
cactus_params$germ2 <- germination$germ2
## precensus survival
cactus_params$precenus_surv <- precensus_survival$precensus_survival
## seedling size
cactus_params$mu_sdlgsize <- seedling_size$mean_size
cactus_params$sigma_sdlgsize <- seedling_size$sd_size
## size bounds
cactus_params$min.size <- log(min(CYIM$vol_t)) 
cactus_params$max.size <- log(max(CYIM$vol_t1))

mat.size = 200
lower.extension = -1
upper.extension = 1.5
kernel_SHASH <- bigmatrix(params = cactus_params,
                          lower.extension = lower.extension, 
                          upper.extension = upper.extension,
                          mat.size = mat.size,
                          dist="SHASH")
kernel_norm <- bigmatrix(params = cactus_params,
                         lower.extension = lower.extension, 
                         upper.extension = upper.extension,
                         mat.size = mat.size,
                         dist="norm")


# And finally for IPM results. The mean model (averaging across years and plots) predicts very different growth rates:
tibble(Growth_Dist = c("SHASH","Gaussian"),
       lambda = c(round(lambda(kernel_SHASH$IPMmat),4),round(lambda(kernel_norm$IPMmat),4)))

#Also, the SSD's are very different and the observed distribution is maybe in between 
# the two of them (in the data plot, colored bars are years). The Gaussian ssd is bi-modal but the 
# lower mode disappears with the skewed $t$ growth kernel, probably because the big reproductive plants 
# are rare in the ssd so you lose the signal of recruitment. 
ssd_SHASH <- stable.stage(kernel_SHASH$IPMmat)[3:(mat.size+2)] / sum(stable.stage(kernel_SHASH$IPMmat)[3:(mat.size+2)])
ssd_norm <- stable.stage(kernel_norm$IPMmat)[3:(mat.size+2)] / sum(stable.stage(kernel_norm$IPMmat)[3:(mat.size+2)])
empirical_sd <- density(log(CYIM$vol_t1),n=mat.size)

pdf("../manuscript/figures/cactus_ssd.pdf",height = 6,width = 6,useDingbats = F)
plot(kernel_norm$meshpts,ssd_norm,type="l",lty=2,lwd=3,xlab="log volume",ylab="Density")
lines(kernel_SHASH$meshpts,ssd_SHASH,type="l",lwd=3)
lines(empirical_sd$x,empirical_sd$y/sum(empirical_sd$y),col="red",lwd=3)
legend("topleft",c("Gaussian SSD","SHASH SSD", "Empirical SD"),lty=c(2,1,1),col=c("black","black","red"),lwd=3,bty="n")
dev.off()

# The difference in ssd's is pretty striking, so I just wanted to have a closer look at the two 
# growth kernels. The skewed $t$ is "peak-ier" than the Gaussian and this causes the Gaussian 
# to allow for large increases in size that don't happen with the skewed $t$. I think that is 
# where the big difference in SSD comes from. 
x = c(-2,5,12)
pdf("../manuscript/figures/cactus_growth_compare.pdf",height = 5,width = 5,useDingbats = F)
plot(size_dum,gxy_SHASH(x=5,y=size_dum,params=cactus_params),
     xlab="Future size",ylab="Density",type="n",ylim=c(0,1.6))
lines(size_dum,gxy_SHASH(x=x[1],y=size_dum,params=cactus_params),lwd=2,col="tomato")
lines(size_dum,gxy_norm(x=x[1],y=size_dum,params=cactus_params),lty=2,lwd=2,col="tomato")
abline(v=x[1],lty=3,col="tomato")
lines(size_dum,gxy_SHASH(x=x[2],y=size_dum,params=cactus_params),lwd=2,col="darkgrey")
lines(size_dum,gxy_norm(x=x[2],y=size_dum,params=cactus_params),lty=2,lwd=2,col="darkgrey")
abline(v=x[2],lty=3,col="darkgrey")
lines(size_dum,gxy_SHASH(x=x[3],y=size_dum,params=cactus_params),lwd=2,col="dodgerblue")
lines(size_dum,gxy_norm(x=x[3],y=size_dum,params=cactus_params),lty=2,lwd=2,col="dodgerblue")
abline(v=x[3],lty=3,col="dodgerblue")
legend("topleft",legend=c("SHASH","Gaussian"),lty=c(1,2),bty="n")
dev.off()




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
