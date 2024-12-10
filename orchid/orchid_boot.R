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
setwd(home);setwd("orchid") 

#IPM source functions are here:
source("orchid/orchis.IPM.source.R")
# functions for life history metrics
source("../code/metaluck_fns_CMH.R")
# functions for bias-corrected bootstrap
source("../code/bca.R")
invlogit<-function(x){exp(x)/(1+exp(x))}

## read in demographic data provided by Hans Jacquemyn (basis for Miller et al. 2012 PRSB)
orchid<- read_csv("Orchis_IPM_data.csv") %>% 
  ## there are two sites/populations, one in light and one in shade.
  ## for the purposes of this analysis, just take the light population
  filter(light=="L") %>% 
  mutate(log_area_t=log(total.leaf.area),
         log_area_t1=log(end.total.leaf.area)) 
seedlings<-read.csv("Orchis_seedlings.csv",T) %>% 
  filter(light=="L") %>% 
  mutate(log_area_t1=log(end.total.leaf.area))

## bootstrap settings and storage
nboot = 501
traits_G = traits_S = matrix(NA,nboot,5) ## hold life history outputs 

## start boot here
for(i in 1:nboot){
  
orchid_boot<-orchid[sample(1:nrow(orchid),nrow(orchid),replace=T),]
seedlings_boot<-seedlings[sample(1:nrow(seedlings),nrow(seedlings),replace=T),]
## run the real data through iteration 1
if(i==1){orchid_boot<-orchid; seedlings_boot<-seedlings}

orchid_boot %>% 
  dplyr::select(individual,log_area_t,log_area_t1,flowering,begin.year) %>% 
  drop_na() -> orchid_grow
orchid_GAU_best<-lmer(log_area_t1~log_area_t * as.logical(flowering) + (1|begin.year),data=orchid_grow,REML=F)
orchid_grow$GAU_fitted <- fitted(orchid_GAU_best)
## now use iterative re-weighting to fit sd as function of expected value
sdloglik = function(pars,resids,fitted) {
  dnorm(resids, mean=0, sd=exp(pars[1]+pars[2]*fitted+pars[3]*fitted^2),log=TRUE)
}	
err = 1; rep=0; 
while(err > 0.000001) {
    rep=rep+1
    model = orchid_GAU_best
    fitted_vals = fitted(model)
    resids = residuals(model) 
    out=maxLik(logLik=sdloglik,start=c(sd(resids),0,0),resids=resids,fitted=fitted_vals)
    pars=out$estimate 
    new_sigma = exp(pars[1] + pars[2]*fitted_vals + pars[3]*fitted_vals^2)
    new_weights = 1/((new_sigma)^2)
    new_weights = 0.5*(weights(model) + new_weights) # cautious update 
    new_model <- update(model,weights=new_weights) 
    err = weights(model)-weights(new_model)
    err=sqrt(mean(err^2))
    orchid_GAU_best<-new_model 
  }
orchid_grow$GAU_fitted <- fitted(orchid_GAU_best)
orchid_grow$GAU_sd <- 1/sqrt(weights(orchid_GAU_best))
stdev_coef <- out$estimate

SSTLogLik=function(pars){
  dSST(orchid_grow$log_area_t1, 
       mu=orchid_grow$GAU_fitted,
       sigma=orchid_grow$GAU_sd,
       nu = exp(pars[1] + pars[2]*orchid_grow$GAU_fitted),
       tau = exp(pars[3] + pars[4]*orchid_grow$GAU_fitted)+2, 
       log=TRUE)
}
p0<-c(0,0,0,0)
SSTout=maxLik(logLik=SSTLogLik,start=p0*exp(0.2*rnorm(length(p0))),method="BHHH",control=list(iterlim=5000,printLevel=0),finalHessian=FALSE) 
SSTout=maxLik(logLik=SSTLogLik,start=SSTout$estimate,method="NM",control=list(iterlim=5000,printLevel=0),finalHessian=FALSE)
SSTout=maxLik(logLik=SSTLogLik,start=SSTout$estimate,method="BHHH",control=list(iterlim=5000,printLevel=0),finalHessian=FALSE)
## other vital rates
surv<-glmer(survival~log_area_t+(1|begin.year),family="binomial",data=orchid_boot)
flower<-glmer(flowering~log_area_t+(1|begin.year),data=orchid_boot,family="binomial")
nflowers<-glmer(number.flowers~log_area_t+(1|begin.year),data=subset(orchid_boot,flowering==1),family="poisson")
propfruit<-glmer(number.fruits/number.flowers~(1|begin.year),weights=number.flowers,data=subset(orchid_boot,flowering==1),family="binomial")
kidsize<-mean(seedlings_boot$log_area_t1)
kidsd<-sd(seedlings_boot$log_area_t1)
seeds<-6000	
eta<-mean(c(0.007,0.023))
sigmap<-0.01485249
sigmat<-0.05940997
dormancy<-glm(dormant~log_area_t,family="binomial",data=orchid_boot)
Dsize<-mean(orchid_boot$log_area_t1[orchid_boot$number.leaves==0],na.rm=T)
Dsizesd<-sd(orchid_boot$log_area_t1[orchid_boot$number.leaves==0],na.rm=T)
minsize<-min(na.omit(c(orchid_boot$log_area_t,orchid_boot$log_area_t1)))
maxsize<-max(na.omit(c(orchid_boot$log_area_t,orchid_boot$log_area_t1)))

### store parameters in vector
params<-c()
params$surv.int <- fixef(surv)[1]
params$surv.size <- fixef(surv)[2]
params$grow.int <- fixef(orchid_GAU_best)[1]
params$grow.size <- fixef(orchid_GAU_best)[2]
params$grow.flow <- fixef(orchid_GAU_best)[3]
params$grow.size.flow <- fixef(orchid_GAU_best)[4]
params$growsd.int <- stdev_coef[1]
params$growsd.fit <- stdev_coef[2]
params$growsd.fit2 <- stdev_coef[3]
params$growSST.nu.int<-SSTout$estimate[1]
params$growSST.nu.fit<-SSTout$estimate[2]
params$growSST.tau.int<-SSTout$estimate[3]
params$growSST.tau.fit<-SSTout$estimate[4]
params$flow.int <- fixef(flower)[1]
params$flow.size <- fixef(flower)[2]
params$flowers.int <- fixef(nflowers)[1]
params$flowers.size <- fixef(nflowers)[2]
params$fruits<-fixef(propfruit)[1]
params$seeds<-seeds
params$kidsize<-kidsize
params$kidsize.sd<-kidsd
params$germ<-eta 
params$sigmap<-sigmap 
params$dorm.int<-coef(dormancy)[1]
params$dorm.size<-coef(dormancy)[2]
params$sigmat<-sigmat 
params$Dsize<-Dsize						
params$Dsizesd<-Dsizesd						
params$minsize<-minsize*0.95 #5% lower than observed
params$maxsize<-maxsize*1.05 #5% higher than observed
params$matsize<-100 ##explored below

mat_GAU<-returnR0(params,dist="GAU",lower.extend=5,upper.extend=2)
mat_SST<-returnR0(params,dist="SST",lower.extend=5,upper.extend=2)

## calculate and store trait values
## mixing distribution is dominated by a cohort of newborns (protocorms)
c0 = rep(0,nrow(mat_GAU$matrix)); c0[1]=1
traits_G[i,] <-c(
  Re(eigen(mat_GAU$matrix)$values[1]), 
  mean_lifespan(mat_GAU$Tmatrix, mixdist=c0),
  mean_LRO(mat_GAU$Tmatrix,mat_GAU$Fmatrix,mixdist=c0),
  mean_age_repro(mat_GAU$Tmatrix,mat_GAU$Fmatrix,mixdist=c0),
  gen_time_mu1_v(mat_GAU$Tmatrix,mat_GAU$Fmatrix))
traits_S[i,] <-c(
  Re(eigen(mat_SST$matrix)$values[1]), 
  mean_lifespan(mat_SST$Tmatrix, mixdist=c0),
  mean_LRO(mat_SST$Tmatrix,mat_SST$Fmatrix,mixdist=c0),
  mean_age_repro(mat_SST$Tmatrix,mat_SST$Fmatrix,mixdist=c0),
  gen_time_mu1_v(mat_SST$Tmatrix,mat_SST$Fmatrix))
print(i)
}

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
cat("point    ", signif(traits_G_true,4),"\n")
cat("boot mean", signif(xbar,3),"\n")
cat("boot sd  ", signif(xsd,3), "\n")
cat("BCA 95% confidence intervals", "\n") 
print(signif(CI_G,3))

################################################### 
## Output results: SST 
###################################################
traits_S_true = traits_S[1,]
traits_S_boot = traits_S[-1,]
xbar = apply(traits_S_boot,2,mean)
xsd = apply(traits_S_boot,2,var)^0.5

### Compute BCA intervals 
CI_S = matrix(NA,2,5)
for(j in 1:5) {
  CI_S[1:2,j]=bca(theta = traits_S_boot[,j], theta_hat = traits_S_true[j], a = 0, conf.level = 0.95) 
}

cat("SST", "\n")
cat("point    ", signif(traits_S_true,4),"\n")
cat("boot mean", signif(xbar,3),"\n")
cat("boot sd  ", signif(xsd,3), "\n")
cat("BCA 95% confidence intervals", "\n") 
print(signif(CI_S,3))

## save bootstrap data
#save.image(file="orchid_boot.Rdata")