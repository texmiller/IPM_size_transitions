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
library(doParallel) 

# functions for life history metrics
source("../code/metaluck_fns_CMH.R")
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

## setting for approximating  matrix
mat.size = 200
lower.extension = -1
upper.extension = 1.5
exclude_plot_year<-c("s(plot)","s(year_t)")

## bootstrap settings and storage
nboot = 500
traits_G = traits_S = matrix(NA,nboot,5) ## hold life history outputs 

## start bootstrap
for(i in 1:nboot){
  cactus_boot<-CYIM_full[sample(1:nrow(CYIM_full),nrow(CYIM_full),replace=T),]
  ## gaussian and shash growth
  grow_GAU <- gam(list(logvol_t1 ~ s(logvol_t,k=4) + s(plot,bs="re") + s(year_t,bs="re"), ~s(logvol_t,k=6)), 
                      data=cactus_boot, family=gaulss())
  grow_SHASH <- gam(list(logvol_t1 ~ s(logvol_t,k=4) + s(plot,bs="re") + s(year_t,bs="re"), # <- model for location 
                             ~ s(logvol_t,k=4),   # <- model for log-scale
                             ~ s(logvol_t,k=4),   # <- model for skewness
                             ~ s(logvol_t,k=4)), # <- model for log-kurtosis
                        data = cactus_boot, 
                        family = shash,  
                        optimizer = "efs")
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
  cactus_c0 = rep(0,nrow(mat_GAU$Tmat)); cactus_c0[1]=1
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

## add "true" values as final row
grow_GAU <- gam(list(logvol_t1 ~ s(logvol_t,k=4) + s(plot,bs="re") + s(year_t,bs="re"), ~s(logvol_t,k=6)),data=CYIM_full, family=gaulss())
grow_SHASH <- gam(list(logvol_t1 ~ s(logvol_t,k=4) + s(plot,bs="re") + s(year_t,bs="re"), # <- model for location 
                       ~ s(logvol_t,k=4),   # <- model for log-scale
                       ~ s(logvol_t,k=4),   # <- model for skewness
                       ~ s(logvol_t,k=4)), # <- model for log-kurtosis
                  data = CYIM_full,family = shash,optimizer = "efs")
surv_mod <- gam(Survival_t1 ~ s(logvol_t,k=4) + s(plot,bs="re") + s(year_t,bs="re"), family="binomial", data=CYIM_full)
flow_mod <- gam(Goodbuds_t>0 ~ s(logvol_t,k=4) + s(plot,bs="re") + s(year_t,bs="re"), family="binomial", data=CYIM_full)
fert_mod <- gam(Goodbuds_t ~ s(logvol_t,k=4) + s(plot,bs="re") + s(year_t,bs="re"), family="nb", data=subset(CYIM_full,Goodbuds_t>0))
## size bounds
min.size <- log(min(CYIM_full$vol_t,na.rm=T)) 
max.size <- log(max(CYIM_full$vol_t1,na.rm=T))

mat_GAU<-bigmatrix(lower.extension = lower.extension, 
                   upper.extension = upper.extension,
                   mat.size = mat.size,exclude=exclude_plot_year,
                   dist="GAU")
mat_SHASH<-bigmatrix(lower.extension = lower.extension, 
                     upper.extension = upper.extension,
                     mat.size = mat.size,exclude=exclude_plot_year,
                     dist="SHASH")

traits_GAU<-rbind(traits_G,c(
  Re(eigen(mat_GAU$IPMmat)$values[1]),
  mean_lifespan(mat_GAU$Tmat, mixdist=cactus_c0),
  mean_LRO(mat_GAU$Tmat,mat_GAU$Fmat,mixdist=cactus_c0),
  mean_age_repro(mat_GAU$Tmat,mat_GAU$Fmat,mixdist=cactus_c0),
  gen_time_mu1_v(mat_GAU$Tmat,mat_GAU$Fmat)
))
traits_SHASH<-rbind(traits_S,c(
  Re(eigen(mat_SHASH$IPMmat)$values[1]),
  mean_lifespan(mat_SHASH$Tmat, mixdist=cactus_c0),
  mean_LRO(mat_SHASH$Tmat,mat_SHASH$Fmat,mixdist=cactus_c0),
  mean_age_repro(mat_SHASH$Tmat,mat_SHASH$Fmat,mixdist=cactus_c0),
  gen_time_mu1_v(mat_SHASH$Tmat,mat_SHASH$Fmat)
))
colnames(traits_GAU)<-c("lambda","meanlife","meanLRO","meanagerepro","gentime")
colnames(traits_SHASH)<-c("lambda","meanlife","meanLRO","meanagerepro","gentime")
write.csv(traits_GAU,"cactus_traits_GAU_boot.csv",row.names=F)
write.csv(traits_SHASH,"cactus_traits_SHASH_boot.csv",row.names=F)
