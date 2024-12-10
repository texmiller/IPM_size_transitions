### move to the right local directory 
tom = "C:/Users/tm9/Dropbox/github/IPM_size_transitions"
steve = "c:/repos/IPM_size_transitions" 
home = ifelse(Sys.info()["user"] == "Ellner", steve, tom)
setwd(home);setwd("pike")

## load libraries
library(tidyverse)
library(lme4)
library(mgcv)
library(gamlss.dist)

## source IPM functions
source("pike_IPM_source_fns.R")
# functions for life history metrics
source("../code/metaluck_fns_CMH.R")
# functions for bias-corrected bootstrap
source("../code/bca.R")

## read in data
pike_dat <- read.csv("data/PikeGrowthData1944_1995.csv")
## this will need some data manipulation, starting with wide format
pike_wide <- pike_dat %>% 
  arrange(Year,Ind) %>% dplyr::select(-RowID) %>% 
  pivot_wider(names_from = Year,values_from = Length)
##now stack transition years
pike_trans_year <- pike_wide[,(1:3)]
pike_trans_year$year <- as.numeric(names(pike_wide)[2])
names(pike_trans_year)[2:3]<-c("t0","t1")
for(i in 3:(ncol(pike_wide)-1)){
  hold <- pike_wide[,c(1,i,(i+1))]
  hold$year <- as.numeric(names(pike_wide)[i])
  names(hold)[2:3]<-c("t0","t1")
  pike_trans_year <- bind_rows(pike_trans_year,hold)
}

## survival data and gam fit
pike_surv <- read_csv("data/PikeSurvivalData1953_1990.csv") %>% 
  mutate(log_t0 = log(Length)) %>% 
  arrange(log_t0)

## fertility data and gam fit
pike_fert <- read_csv("data/FecundityData1963_2002.csv") %>% 
  mutate(log_t0 = log(Length),
         log_eggs = log(Eggs)) %>% 
  arrange(log_t0)

## bootstrap settings and life history outputs 
nboot = 501
traits_G = traits_S = matrix(NA,nboot,5) 

## start bootstrap
for(i in 1:nboot){
  ## bootstrap the main data sources
  pike_boot<-pike_trans_year[sample(1:nrow(pike_trans_year),nrow(pike_trans_year),replace=T),]
  pike_surv_boot<-pike_surv[sample(1:nrow(pike_surv),nrow(pike_surv),replace=T),]
  pike_fert_boot<-pike_fert[sample(1:nrow(pike_fert),nrow(pike_fert),replace=T),]
  if(i==1){pike_boot<-pike_trans_year;pike_surv_boot<-pike_surv;pike_fert_boot<-pike_fert}
  
  pike_boot %>% filter(!is.na(t0) & !is.na(t1)) %>% 
    mutate(log_t0 = log(t0),
         log_t1 = log(t1)) %>% 
    arrange(log_t0,log_t1)-> pike_final
  pike_boot %>% filter(is.na(t0) & !is.na(t1)) %>% 
    mutate(log_t1 = log(t1)) %>% 
    arrange(log_t1)-> pike_recruits
  
## growth
  if (i==1) {
    pike_gau<-gam(list(log_t1 ~ s(log_t0,k=5), ~s(log_t0,k=5)), data=pike_final, family=gaulss()) 
    pike_shash <- gam(list(log_t1 ~ s(log_t0,k=5), # <- model for location 
                               ~ s(log_t0,k=5),   # <- model for log-scale
                               ~ s(log_t0,k=5),   # <- model for skewness
                               ~ s(log_t0,k=5)), # <- model for log-kurtosis
                          data = pike_final,family = shash,optimizer = "efs")
    pike_surv_gam <- gam(Survival ~ s(log_t0,k=5),data = pike_surv_boot,family = binomial)
    pike_fert_gam <- gam(Eggs ~ s(log_t0,k=5),data = pike_fert_boot,family = nb)
    pike_recruitsize_shash <- gam(list(log_t1 ~ 1, # <- model for location 
                                       ~ 1,   # <- model for log-scale
                                       ~ 1,   # <- model for skewness
                                       ~ 1), # <- model for log-kurtosis
                                  data = pike_recruits,family = shash,optimizer = "efs")
    ## store smoothing parameters for boot iterations (no smooths for recruit size)
    pike_gau_sp<-pike_gau$sp
    pike_shash_sp<-pike_shash$sp
    pike_surv_sp<-pike_surv_gam$sp
    pike_fert_sp<-pike_fert_gam$sp
  } else {
    pike_gau<-gam(list(log_t1 ~ s(log_t0,k=5), ~s(log_t0,k=5)), data=pike_final, family=gaulss(), sp=pike_gau_sp)
    pike_shash <- gam(list(log_t1 ~ s(log_t0,k=5), # <- model for location 
                           ~ s(log_t0,k=5),   # <- model for log-scale
                           ~ s(log_t0,k=5),   # <- model for skewness
                           ~ s(log_t0,k=5)), # <- model for log-kurtosis
                      data = pike_final,family = shash,optimizer = "efs",sp=pike_shash_sp)
    pike_surv_gam <- gam(Survival ~ s(log_t0,k=5),data = pike_surv_boot,family = binomial,sp=pike_surv_sp)
    pike_fert_gam <- gam(Eggs ~ s(log_t0,k=5),data = pike_fert_boot,family = nb, sp=pike_fert_sp)
    pike_recruitsize_shash <- gam(list(log_t1 ~ 1, # <- model for location 
                                       ~ 1,   # <- model for log-scale
                                       ~ 1,   # <- model for skewness
                                       ~ 1), # <- model for log-kurtosis
                                  data = pike_recruits,family = shash,optimizer = "efs")
  }

## approximating matrices
mat_GAU<-ApproxMatrix(dist="GAU",mat.size=1000,lower=min(pike_final$log_t0),upper=max(pike_final$log_t0))
mat_SHASH<-ApproxMatrix(dist="SHASH",mat.size=1000,lower=min(pike_final$log_t0),upper=max(pike_final$log_t0))

## calculate and store trait values
## define mixing distribution based on recruit size distribution -- need to add the two seed banks
c0 = recruit.size(y=mat_GAU$meshpts)
c0 = c0/sum(c0)
traits_G[i,] <-c(
  Re(eigen(mat_GAU$IPMmat)$values[1]), 
  mean_lifespan(mat_GAU$Pmat, mixdist=c0),
  mean_LRO(mat_GAU$Pmat,mat_GAU$Fmat,mixdist=c0),
  mean_age_repro(mat_GAU$Pmat,mat_GAU$Fmat,mixdist=c0),
  gen_time_mu1_v(mat_GAU$Pmat,mat_GAU$Fmat))
traits_S[i,] <-c(
  Re(eigen(mat_SHASH$IPMmat)$values[1]), 
  mean_lifespan(mat_SHASH$Pmat, mixdist=c0),
  mean_LRO(mat_SHASH$Pmat,mat_SHASH$Fmat,mixdist=c0),
  mean_age_repro(mat_SHASH$Pmat,mat_SHASH$Fmat,mixdist=c0),
  gen_time_mu1_v(mat_SHASH$Pmat,mat_SHASH$Fmat))
print(i)
}

# or start here
load(file="pike_boot.Rdata")

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
print(signif(CI_G,3))


################################################### 
## Output results: SHASH
###################################################
traits_S_true = traits_S[1,]
traits_S_boot = traits_S[-1,]
xbar = apply(traits_S_boot,2,mean)
xsd = apply(traits_S_boot,2,var)^0.5

### Compute BCA intervals 
CI_S = matrix(NA,2,5)
for(j in 1:5) {
  CI_S[1:2,j]=bca(theta = traits_S_boot[,j], theta_hat = traits_G_true[j], a = 0, conf.level = 0.95) 
}

cat("SHASH", "\n")
cat("point    ", signif(traits_S_true,3),"\n")
cat("boot mean", signif(xbar,3),"\n")
cat("boot sd  ", signif(xsd,3), "\n")
cat("BCA 95% confidence intervals", "\n") 
print(signif(CI_S,3))

#save.image(file="pike_boot.Rdata")