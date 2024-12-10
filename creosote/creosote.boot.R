### move to the right local directory 
tom = "C:/Users/tm9/Dropbox/github/IPM_size_transitions"
steve = "c:/repos/IPM_size_transitions" 
home = ifelse(Sys.info()["user"] == "Ellner", steve, tom)
setwd(home); setwd("creosote"); 

library(lme4)
library(mgcv)
library(qgam)
library(tidyverse)
library(quantreg)
library(gamlss.dist)
library(maxLik)
library(bbmle)
library(popbio)
library(moments)
library(minpack.lm)
library(sqldf)
library(SuppDists)
library(Rage)
invlogit<-function(x){exp(x)/(1+exp(x))}

## grab the creosote demography data from github
source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/04_CDataPrep.R")
## source IPM functions
source("creosote_IPM_source_fns.R")
# functions for life history metrics
source("../code/metaluck_fns_CMH.R")
# functions for bias-corrected bootstrap
source("../code/bca.R")

## define k parameter (basis number) and gamma for global use
k_param<-6
gamma = 1.4

## bootstrap settings and storage
nboot = 501
traits_G = traits_J = matrix(NA,nboot,5) ## hold life history outputs 

## start bootstrap here
for(i in 1:nboot){
  ## bootstrap the main data sources
  CData_boot<-CData[sample(1:nrow(CData),nrow(CData),replace=T),]
  CData.Transplants_boot<-CData.Transplants[sample(1:nrow(CData.Transplants),nrow(CData.Transplants),replace=T),]
  ## run the real data through iteration 1
  if(i==1){CData_boot<-CData; CData.Transplants_boot<-CData.Transplants}
# data prep ---------------------------------------------------------------
# all the pieces are derivatives of CData
# Growth data
LATR_full <- CData_boot %>% 
  mutate(unique.transect = interaction(transect, site)) %>%
  mutate(log_volume_t = log(volume_t),
         log_volume_t1 = log(volume_t1),
         dens_scaled = weighted.dens/100)
# Prepare a data subset for growth that drops rows missing either t or t1 size data
## update to the growth data -- dropping a few unbelievable outliers following additional QA/QC
outliers<-c("MOD.2.50.3.2016","MOD.3.200.1.2015","MOD.3.200.1.2014","PDC.2.0.5.2014")
LATR_grow <- LATR_full %>% 
  ## this ID will help us drop outliers below
  mutate(ID=interaction(site,transect,actual.window,plant,year_t)) %>% 
  drop_na(volume_t, volume_t1) %>% 
  ##need to scale weighted density because 1st and 2nd order variables were hugely different in range
  filter(!ID%in%outliers)
e = order(LATR_grow$log_volume_t)
LATR_grow = LATR_grow[e,]

# Reproduction data
# Populate year t of 2017-2018 transition year
# There are no 2018 data but this way we get all four years in the reproduction models
# Do this by creating the 2017-18 data as a stand-alone df then bind rows
LATR_dat_201718 <- LATR_full[LATR_full$year_t == 2016 & LATR_full$survival_t1 == 1, ]
# These are the 2017 survivors; make their year t demography last year's data
LATR_dat_201718$year_t <- 2017
LATR_dat_201718$year_t1 <- 2018
LATR_dat_201718$log_volume_t <- LATR_dat_201718$log_volume_t1
LATR_dat_201718$flowers_t <- LATR_dat_201718$flowers_t1
LATR_dat_201718$fruits_t <- LATR_dat_201718$fruits_t1
LATR_dat_201718$reproductive_fraction_t <- LATR_dat_201718$reproductive_fraction_t1
LATR_dat_201718$total.reproduction_t <- LATR_dat_201718$total.reproduction_t1
# Now set all the t1 data to NA
LATR_dat_201718$log_volume_t1 <- NA
LATR_dat_201718$flowers_t1 <- NA
LATR_dat_201718$fruits_t1 <- NA
LATR_dat_201718$reproductive_fraction_t1 <- NA
LATR_dat_201718$total.reproduction_t1 <- NA
# Bind rows and create log_vol as new variables (easier for GAMs)
LATR_flow_dat <- bind_rows(LATR_full,LATR_dat_201718) %>% 
  dplyr::select(unique.transect,log_volume_t,total.reproduction_t,dens_scaled) %>% drop_na()
# Create new df with plants that have produced at least one reproductive structure
LATR_fruits_dat <- subset(LATR_flow_dat, total.reproduction_t > 0)

# Survival data
# Combine transplants with large shrubs; keep only location info, survival, volume, and density
CData.Transplants_boot %>% 
  dplyr::select("site", "transect", "actual.window", 
                "spring_survival_t1", "volume_t", "weighted.dens", "transplant") %>% 
  rename("survival_t1" = "spring_survival_t1") %>% 
  mutate(unique.transect = interaction(transect, site)) %>% 
  rbind(dplyr::select(LATR_full, "site", "transect", "actual.window", 
                      "survival_t1", "volume_t", "weighted.dens", "transplant","unique.transect")) %>% 
  mutate(log_volume_t = log(volume_t),
         dens_scaled = weighted.dens/100) %>% 
  drop_na() -> LATR_surv_dat

# Filter out seedlings and get their sizes
LATR_recruit_size <- LATR_full %>% 
  filter(seedling_t1 == 1) %>% 
  mutate(log_volume = log(volume_t1)) %>% 
  arrange(dens_scaled)

# Create maximum and minimum size bounds for the IPM
LATR_size_bounds <- data.frame(min_size = log(min(LATR_full$volume_t, LATR_full$volume_t1[LATR_full$transplant == FALSE], na.rm = TRUE)),
                               max_size = log(max(LATR_full$volume_t, LATR_full$volume_t1[LATR_full$transplant == FALSE], na.rm = TRUE)))

# fit models --------------------------------------------------------------
# "best" models follow model selection from original analysis

# fit gaussian growth model with mgcv
LATR_GAU <- gam(list(log_volume_t1~s(log_volume_t,k=k_param) + s(dens_scaled,k=k_param) + s(unique.transect,bs="re"),~1), 
                  family="gaulss", data=LATR_grow, method="ML",gamma=gamma)

fitted_GAU = predict(LATR_GAU,type="response",data=LATR_grow)                  
new_fitted_vals = fitted_GAU[,1]
LATR_grow$fitted_vals = new_fitted_vals
weights = fitted_GAU[,2] #1/sigma_i
#re-fit with sd~f(fitted)
if (i==1) {
  LATR_GAU <- gam(list(log_volume_t1~s(log_volume_t,k=k_param) + s(dens_scaled,k=k_param) + s(unique.transect,bs="re"),~s(fitted_vals,k=k_param)), 
                  family="gaulss", data=LATR_grow, method="ML",gamma=gamma)
} else {
  LATR_GAU <- gam(list(log_volume_t1~s(log_volume_t,k=k_param) + s(dens_scaled,k=k_param) + s(unique.transect,bs="re"),~s(fitted_vals,k=k_param)), 
                  family="gaulss", data=LATR_grow, method="ML",gamma=gamma,sp=LATR_GAU_best_sp)
}
err=100; k=0; 
while(err>10^(-6)) {
  LATR_grow$fitted_vals = new_fitted_vals; 
  LATR_GAU<-update(LATR_GAU)
  fitted_all = predict(LATR_GAU,type="response",data=LATR_grow)   
  new_fitted_vals = fitted_all[,1]; new_weights = fitted_all[,2]
  err = weights - new_weights; err=sqrt(mean(err^2)) 
  weights = new_weights; 
  k=k+1
}
LATR_GAU_best<-LATR_GAU
LATR_grow$GAU_mean<- fitted(LATR_GAU_best)[,1]
LATR_grow$GAU_sd<- 1/fitted(LATR_GAU_best)[,2]
if(i==1){LATR_GAU_best_sp<-LATR_GAU_best$sp}

## fit JSU using mean and sd of LATR_GAU_best
JSULogLik=function(pars){
  dJSU(LATR_grow$log_volume_t1, 
       mu=LATR_grow$GAU_mean,
       sigma=LATR_grow$GAU_sd,
       nu = pars[1]+pars[2]*LATR_grow$GAU_mean,
       tau = exp(pars[3]+pars[4]*LATR_grow$GAU_mean), log=TRUE)
}
## starting parameters
p0<-c(0,0,0,0)
## pass this through several ML algorithms
JSUout=maxLik(logLik=JSULogLik,start=p0*exp(0.2*rnorm(length(p0))),method="BHHH",control=list(iterlim=5000,printLevel=0),finalHessian=FALSE) 
JSUout=maxLik(logLik=JSULogLik,start=JSUout$estimate,method="NM",control=list(iterlim=5000,printLevel=0),finalHessian=FALSE)
JSUout=maxLik(logLik=JSULogLik,start=JSUout$estimate,method="BHHH",control=list(iterlim=5000,printLevel=0),finalHessian=FALSE)

# probability of flowering, fruit production, survival
if (i==1) {
  LATR_flower_best <- gam(total.reproduction_t > 0 ~ s(log_volume_t) + s(dens_scaled) + ti(log_volume_t,dens_scaled) + s(unique.transect, bs = "re"),
                          data = LATR_flow_dat, gamma = gamma, family = "binomial")
  LATR_fruits_best <- gam(total.reproduction_t ~ s(log_volume_t) + s(dens_scaled) + s(unique.transect, bs = "re"),
                          data = LATR_fruits_dat, gamma = gamma, family = "nb")
  LATR_surv_best <- gam(survival_t1 ~ s(log_volume_t,by=as.factor(transplant)) + s(dens_scaled,by=as.factor(transplant))  + s(unique.transect, bs = "re"),
                        data = LATR_surv_dat, gamma = gamma, family = "binomial")
  # store smoothing params
  LATR_flower_best_sp<-LATR_flower_best$sp
  LATR_fruits_best_sp<-LATR_fruits_best$sp
  LATR_surv_best_sp<-LATR_surv_best$sp
} else {
  LATR_flower_best <- gam(total.reproduction_t > 0 ~ s(log_volume_t) + s(dens_scaled) + ti(log_volume_t,dens_scaled) + s(unique.transect, bs = "re"),
                          data = LATR_flow_dat, gamma = gamma, family = "binomial",sp=LATR_flower_best_sp)
  LATR_fruits_best <- gam(total.reproduction_t ~ s(log_volume_t) + s(dens_scaled) + s(unique.transect, bs = "re"),
                          data = LATR_fruits_dat, gamma = gamma, family = "nb",sp=LATR_fruits_best_sp)
  LATR_surv_best <- gam(survival_t1 ~ s(log_volume_t,by=as.factor(transplant)) + s(dens_scaled,by=as.factor(transplant))  + s(unique.transect, bs = "re"),
                        data = LATR_surv_dat, gamma = gamma, family = "binomial",sp=LATR_surv_best_sp)
}

## don't bootstrap any of the recruit stuff bc it's too small a data set
if(i==1){
# Recruitment data
## number of seeds per fruit
seeds_per_fruit <- 5
# Create subset df of recruits
LATR_recruits <- LATR_full %>% 
  mutate(unique.transect = interaction(transect, site)) %>% 
  group_by(year_t1, unique.transect, actual.window) %>% 
  filter(seedling_t1 == 1)
suppressMessages(summarise(LATR_recruits, recruits = n())) %>% 
  rename(window = actual.window) -> LATR_recruits
# Estimate total seeds produced in each window
# This is computed using the known plant sizes and the fitted flowering and fruiting models
LATR_transects <- Cdata.Transects.Windows %>% 
  mutate(unique.transect = interaction(transect, site),
         log_volume_t = log(volume),
         dens_scaled = weighted.dens/100)
LATR_transects$seeds = ceiling(invlogit(predict.gam(LATR_flower_best, newdata = LATR_transects))* 
                                 seeds_per_fruit*exp(predict.gam(LATR_fruits_best, newdata = LATR_transects)))
LATR_transects <- group_by(LATR_transects, unique.transect, window)
suppressMessages(summarise(LATR_transects, total_seeds = sum(seeds),
                           dens_scaled = unique(dens_scaled))) -> LATR_transects
# Take three copies of this df, assigning each one to a different year and assigning recruits to zero (for now)
LATR_recruitment <- bind_rows(LATR_transects %>% filter(unique.transect == "1.FPS" | unique.transect == "2.FPS" | unique.transect == "3.FPS") %>% 
                                mutate(year_t1 = 2014, recruits = 0), ## only FPS for 2013-2014
                              LATR_transects %>% mutate(year_t1 = 2015, recruits = 0),
                              LATR_transects %>% mutate(year_t1 = 2016, recruits = 0),
                              LATR_transects %>% mutate(year_t1 = 2017, recruits = 0)) %>% 
  left_join(., LATR_recruits, by = c("year_t1", "unique.transect", "window")) %>% 
  mutate(recruits.y = replace_na(recruits.y, 0),
         recruits = pmax(recruits.x, recruits.y, na.rm = T)) %>% 
  drop_na()
## pr recruitment 
LATR_recruit_best <- gam(cbind(recruits,total_seeds - recruits) ~ s(unique.transect, bs = "re"),
                           data = LATR_recruitment, gamma = gamma, family = "binomial")
## recruit size distribution
LATR_recruitsize_best <- gam(list(log_volume ~ s(dens_scaled) + s(unique.transect,bs = "re"),~1), 
                               data = LATR_recruit_size, gamma = gamma, family = gaulss())
recruit_size_sd_index <- which(as.factor(names(coef(LATR_recruitsize_best)))=="(Intercept).1") ## this is where the sd coefficients start
recruit_size_coef_length <- length(coef(LATR_recruitsize_best))
}

# Finally - construct transition matrix for minimum weighted density (zero)
mat_GAU <- ApproxMatrix(dens = 0, dist="GAU",mat.size=200,ext.lower=-8,ext.upper=2)
mat_JSU <- ApproxMatrix(dens = 0, dist="JSU",mat.size=200,ext.lower=-8,ext.upper=2)

## calculate and store trait values
## define mixing distribution based on recruit size distribution -- need to add the two seed banks
c0 = recruitsize(y=mat_GAU$meshpts,d=0)
c0 = c0/sum(c0)
traits_G[i,] <-c(
  Re(eigen(mat_GAU$IPMmat)$values[1]), 
  mean_lifespan(mat_GAU$Pmat, mixdist=c0),
  mean_LRO(mat_GAU$Pmat,mat_GAU$Fmat,mixdist=c0),
  mean_age_repro(mat_GAU$Pmat,mat_GAU$Fmat,mixdist=c0),
  gen_time_mu1_v(mat_GAU$Pmat,mat_GAU$Fmat))
traits_J[i,] <-c(
  Re(eigen(mat_JSU$IPMmat)$values[1]), 
  mean_lifespan(mat_JSU$Pmat, mixdist=c0),
  mean_LRO(mat_JSU$Pmat,mat_JSU$Fmat,mixdist=c0),
  mean_age_repro(mat_JSU$Pmat,mat_JSU$Fmat,mixdist=c0),
  gen_time_mu1_v(mat_JSU$Pmat,mat_JSU$Fmat))
print(i)
}#end bootstrap loop

# start here
load("creosote_boot.Rdata")

################################################### 
## Output results: GAUSSIAN 
###################################################
traits_G_true = traits_G[1,]
traits_G_boot = traits_G[-1,]
xbar = apply(traits_G_boot,2,mean)
xsd = apply(traits_G_boot,2,var)^0.5

### Compute BCA intervals 
CI_G = matrix(NA,2,5)
for(j in 1:5) {
  CI_G[1:2,j]=bca(theta = traits_G_boot[,j], theta_hat = traits_G_true[j], a = 0, conf.level = 0.95) 
}

cat("GAUSSIAN", "\n")
cat("point    ", signif(traits_G_true,3),"\n")
cat("boot mean", signif(xbar,3),"\n")
cat("boot sd  ", signif(xsd,3), "\n")
cat("BCA 95% confidence intervals", "\n") 
print(signif(CI_G,4))

par(mfrow=c(2,3))
hist(traits_G_boot[,1])
abline(v=CI_G[,1],col="red",lwd=2,lty=2)
abline(v=quantile(traits_G_boot[,1],c(0.025,0.975)),col="blue",lwd=2,lty=2)
hist(traits_G_boot[,2])
abline(v=CI_G[,2],col="red",lwd=2,lty=2)
abline(v=quantile(traits_G_boot[,2],c(0.025,0.975)),col="blue",lwd=2,lty=2)
hist(traits_G_boot[,3])
abline(v=CI_G[,3],col="red",lwd=2,lty=2)
abline(v=quantile(traits_G_boot[,3],c(0.025,0.975)),col="blue",lwd=2,lty=2)
hist(traits_G_boot[,4])
abline(v=CI_G[,4],col="red",lwd=2,lty=2)
abline(v=quantile(traits_G_boot[,4],c(0.025,0.975)),col="blue",lwd=2,lty=2)
hist(traits_G_boot[,5])
abline(v=CI_G[,5],col="red",lwd=2,lty=2)
abline(v=quantile(traits_G_boot[,5],c(0.025,0.975)),col="blue",lwd=2,lty=2)

###################################################
##### Output results: JSU
###################################################
traits_J_true = traits_J[1,]
traits_J_boot = traits_J[-1,]
xbar = apply(traits_J_boot,2,mean) 
xsd = apply(traits_J_boot,2,var)^0.5 

### Compute BC intervals 
CI_J = matrix(NA,2,5) 
for(j in 1:5) {
  CI_J[1:2,j]=bca(theta  = traits_J_boot[,j], theta_hat = traits_J_true[j], a = 0, conf.level = 0.95) 
}

cat("JSU", "\n")
cat("point    ", signif(traits_J_true,3),"\n")
cat("boot mean", signif(xbar,3),"\n")
cat("boot sd  ", signif(xsd,3), "\n")
cat("BC 95% confidence intervals", "\n") 
print(signif(CI_J,3))

#save.image(file="creosote_boot.Rdata")

