## purpose: bootstrap resampling to estimate confidence limits on the cactus model

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

# functions for life history metrics
source("../code/metaluck_fns_CMH.R")
source("../code/bca.R")
# functions for cactus IPM
source("cactus_IPM_source_fns.R")

## read in cactus demography data
## these data are published on EDI: https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-sev.323.1
CYIM_full<-read_csv("cholla_demography_20042018_EDI.csv")%>% 
  ## filter out transplants and drop seed addition plots (which start with letter H)
  ## also drop Year_t==2018 because I don't have 2019 size data. 2018 data still included in 2017-2018 transition.
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

## misc parameters for seeds and seedlings
## there is not a lot of data here so I am not bootstrapping this
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

## fit growth gams to get smoothing parameters
grow_GAU <- gam(list(logvol_t1 ~ s(logvol_t,k=4) + s(plot,bs="re") + s(year_t,bs="re"), ~s(logvol_t,k=6)), 
                data=CYIM_full, family=gaulss())
grow_GAU_sp<-grow_GAU$sp
grow_SHASH <- gam(list(logvol_t1 ~ s(logvol_t,k=4) + s(plot,bs="re") + s(year_t,bs="re"), # <- model for location 
                       ~ s(logvol_t,k=4),   # <- model for log-scale
                       ~ s(logvol_t,k=4),   # <- model for skewness
                       ~ s(logvol_t,k=4)), # <- model for log-kurtosis
                  data = CYIM_full, 
                  family = shash,  
                  optimizer = "efs")
grow_SHASH_sp<-grow_SHASH$sp
## settings for approximating  matrix
mat.size = 200
lower.extension = -2.5 #-1 (I think I was evicting tiny seedling, though it did not matter much)
upper.extension = 1.5
exclude_plot_year<-c("s(plot)","s(year_t)")

## bootstrap settings and storage
nboot = 501
traits_G = traits_S = matrix(NA,nboot,5) ## hold life history outputs 

## start bootstrap
for(i in 1:nboot){
  cactus_boot<-CYIM_full[sample(1:nrow(CYIM_full),nrow(CYIM_full),replace=T),]
  ## first boot iteration, use the real data
  if(i==1){cactus_boot<-CYIM_full}
  
  ## gaussian and shash growth
  grow_GAU <- gam(list(logvol_t1 ~ s(logvol_t,k=4) + s(plot,bs="re") + s(year_t,bs="re"), ~s(logvol_t,k=6)), 
                      data=cactus_boot, family=gaulss(), sp=grow_GAU_sp)
  grow_SHASH <- gam(list(logvol_t1 ~ s(logvol_t,k=4) + s(plot,bs="re") + s(year_t,bs="re"), # <- model for location 
                             ~ s(logvol_t,k=4),   # <- model for log-scale
                             ~ s(logvol_t,k=4),   # <- model for skewness
                             ~ s(logvol_t,k=4)), # <- model for log-kurtosis
                        data = cactus_boot,family = shash,optimizer = "efs",sp=grow_SHASH_sp)
  surv_mod <- gam(Survival_t1 ~ s(logvol_t,k=4) + s(plot,bs="re") + s(year_t,bs="re"), family="binomial", data=cactus_boot)
  flow_mod <- gam(Goodbuds_t>0 ~ s(logvol_t,k=4) + s(plot,bs="re") + s(year_t,bs="re"), family="binomial", data=cactus_boot)
  fert_mod <- gam(Goodbuds_t ~ s(logvol_t,k=4) + s(plot,bs="re") + s(year_t,bs="re"), family="nb", data=subset(cactus_boot,Goodbuds_t>0))
  ## size bounds
  min.size <- log(min(cactus_boot$vol_t,na.rm=T)) 
  max.size <- log(max(cactus_boot$vol_t1,na.rm=T))
  
  mat_GAU<-bigmatrix(lower.extension = lower.extension, 
                                upper.extension = upper.extension,
                                mat.size = mat.size,exclude=exclude_plot_year,
                                dist="GAU")
  mat_SHASH<-bigmatrix(lower.extension = lower.extension, 
                                  upper.extension = upper.extension,
                                  mat.size = mat.size,exclude=exclude_plot_year,
                                  dist="SHASH")
  ## define mixing distribution based on recruit size distribution -- need to add the two seed banks
  cactus_c0 = c(0,0,recruit.size(mat_GAU$meshpts))
  cactus_c0 = cactus_c0/sum(cactus_c0)
  traits_G[i,]<-c(
    Re(eigen(mat_GAU$IPMmat)$values[1]),
    mean_lifespan(mat_GAU$Tmat, mixdist=cactus_c0),
    mean_LRO(mat_GAU$Tmat,mat_GAU$Fmat,mixdist=cactus_c0),
    mean_age_repro(mat_GAU$Tmat,mat_GAU$Fmat,mixdist=cactus_c0),
    gen_time_mu1_v(mat_GAU$Tmat,mat_GAU$Fmat)
  )
  traits_S[i,]<-c(
    Re(eigen(mat_SHASH$IPMmat)$values[1]),
    mean_lifespan(mat_SHASH$Tmat, mixdist=cactus_c0),
    mean_LRO(mat_SHASH$Tmat,mat_SHASH$Fmat,mixdist=cactus_c0),
    mean_age_repro(mat_SHASH$Tmat,mat_SHASH$Fmat,mixdist=cactus_c0),
    gen_time_mu1_v(mat_SHASH$Tmat,mat_SHASH$Fmat)
  )
  print(i)
}

save.image(file="cactus_boot.Rdata")

traits_G_true = traits_G[1,]
traits_S_true = traits_S[1,]


################################################### 
## Output results: GAUSSIAN 
###################################################
traits_G_boot = traits_G[-1,]
xbar = apply(traits_G_boot,2,mean)
xsd = apply(traits_G_boot,2,var)^0.5

### Compute BCA intervals 
CI_G = matrix(NA,2,5)
for(j in 1:5) {
  CI_G[1:2,j]=bca(traits_G_boot[,j], conf.level = 0.95) 
}

cat("GAUSSIAN", "\n")
cat("point    ", signif(traits_G_true,3),"\n")
cat("boot mean", signif(xbar,3),"\n")
cat("boot sd  ", signif(xsd,3), "\n")
cat("BCA 95% confidence intervals", "\n") 
print(signif(CI_G,3))

## these intervals seem bizarrely wide; see how they look on the histograms
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
## Output results: SHASH
###################################################
traits_S_boot = traits_S[-1,]
xbar = apply(traits_S_boot,2,mean)
xsd = apply(traits_S_boot,2,var)^0.5

### Compute BCA intervals 
CI_S = matrix(NA,2,5)
for(j in 1:5) {
  CI_S[1:2,j]=bca(traits_S_boot[,j], conf.level = 0.95) 
}

cat("SHASH", "\n")
cat("point    ", signif(traits_S_true,3),"\n")
cat("boot mean", signif(xbar,3),"\n")
cat("boot sd  ", signif(xsd,3), "\n")
cat("BCA 95% confidence intervals", "\n") 
print(signif(CI_S,3))

par(mfrow=c(2,3))
hist(traits_S_boot[,1])
abline(v=CI_S[,1],col="red",lwd=2,lty=2)
abline(v=quantile(traits_S_boot[,1],c(0.025,0.975)),col="blue",lwd=2,lty=2)
hist(traits_S_boot[,2])
abline(v=CI_S[,2],col="red",lwd=2,lty=2)
abline(v=quantile(traits_S_boot[,2],c(0.025,0.975)),col="blue",lwd=2,lty=2)
hist(traits_S_boot[,3])
abline(v=CI_S[,3],col="red",lwd=2,lty=2)
abline(v=quantile(traits_S_boot[,3],c(0.025,0.975)),col="blue",lwd=2,lty=2)
hist(traits_S_boot[,4])
abline(v=CI_S[,4],col="red",lwd=2,lty=2)
abline(v=quantile(traits_S_boot[,4],c(0.025,0.975)),col="blue",lwd=2,lty=2)
hist(traits_S_boot[,5])
abline(v=CI_S[,5],col="red",lwd=2,lty=2)
abline(v=quantile(traits_S_boot[,5],c(0.025,0.975)),col="blue",lwd=2,lty=2)

