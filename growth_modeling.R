## Model size transitions using data list prepared in dryad_data_prep.R
library(tidyverse)
library(moments)
library(zoo)
library(mgcv)
library(scales)
library(gamlss)
library(gamlss.tr)
library(bbmle)
library(lme4)

home = "c:/repos" ## edit as needed 
setwd(home); setwd("IPM_size_transitions"); 

## load data
growth_dat <- readRDS("growth_dat.rds")

## Steve's diagnostics code
### Scatterplot smoothing functions
spline.scatter.smooth=function(x,y,...) {
  fit=gam(y~s(x),gamma=10,method="REML")
  plot(x,y,type="p",...);
  out=predict(fit,type="response"); 
  points(x,out,type="l",lwd=2)
  #fit2 = lm(y~x+I(x^2));
  #points(x,fit2$fitted,type="l",col="red",lwd=2,lty=2);
}

par(mfrow=c(4,4),bty="l",mar=c(4,4,1,1),mgp=c(2,1,0),cex.axis=1.3,cex.lab=1.3);
for(i in 1:4){
z0 = growth_dat[[i]]$t0; z1=growth_dat[[i]]$t1; 
e = order(z0); z0=z0[e]; z1 = z1[e]; 

rollmean<-rollapply(z0,width=20,mean,na.rm=T)
rollmean1<-rollapply(z1,width=20,mean,na.rm=T)
rollsd <-rollapply(z1,width=20,sd,na.rm=T); 
rollskew<-rollapply(z1,width=20,skewness,na.rm=T)
rollkurt<-rollapply(z1,width=20,kurtosis,na.rm=T)/3 - 1
spline.scatter.smooth(rollmean,rollmean1,col="grey50",xlab="z0",ylab="z1");title(main=names(growth_dat[i])) 
spline.scatter.smooth(rollmean,rollsd,col="grey50",xlab="z0",ylab="SD"); 
spline.scatter.smooth(rollmean,rollskew,col="grey50",xlab="z0",ylab="Skew"); 
spline.scatter.smooth(rollmean,rollkurt,col="grey50",xlab="z0",ylab="Excess Kurtosis"); 

#scatter.smooth(rollmean,rollmean1,col="grey50",xlab="z0",ylab="z1",degree=2);title(main=names(growth_dat[i])) 
#scatter.smooth(rollmean,rollsd,col="grey50",xlab="z0",ylab="SD",degree=2); 
#scatter.smooth(rollmean,rollskew,col="grey50",xlab="z0",ylab="Skew",degree=2); 
#scatter.smooth(rollmean,rollkurt,col="grey50",xlab="z0",ylab="Excess Kurtosis",degree=2); 
}


# coral data --------------------------------------------------------------
coral <- growth_dat[["coral"]]
## there are sites, plots, and years here - but mostly sites are visited in separate years
table(coral$Plot,coral$Year,coral$Site)

## a simple linear regression
coral_lm <- lm(t1 ~ t0, data = coral)
plot(t1 ~ t0, data = coral)
abline(coef(coral_lm))

## simluate data with fitted model
layout.matrix <- matrix(c(1, 2, 3, 1, 4, 5), nrow = 3, ncol = 2)
layout(mat = layout.matrix, heights = c(2, 1, 1))
layout.show(5)
plot(density(coral$t1),type="n")
n_sim <- 500
coral_sim <- matrix(NA,dim(coral)[1],n_sim); coral_mean <- coral_sd <- coral_skew <- coral_kurt <- vector("numeric",n_sim)
for(i in 1:n_sim){
  coral_sim[,i] <- rnorm(n = nrow(coral_sim), mean = coef(coral_lm)[1] + coef(coral_lm)[2] * coral$t0, sd = sigma(coral_lm))
  lines(density(coral_sim[,i]),col=alpha("gray",0.5))
  coral_mean[i] <- mean(coral_sim[,i])
  coral_sd[i] <- sd(coral_sim[,i])
  coral_skew[i] <- skewness(coral_sim[,i])
  coral_kurt[i] <- kurtosis(coral_sim[,i])
}
lines(density(coral$t1),lwd=4)
plot(density(coral_mean),lwd=2,col="gray",xlim=c(min(c(coral_mean,mean(coral$t1))),max(c(coral_mean,mean(coral$t1)))))
abline(v=mean(coral$t1),lwd=4)
plot(density(coral_sd),lwd=2,col="gray",xlim=c(min(c(coral_sd,sd(coral$t1))),max(c(coral_sd,sd(coral$t1)))))
abline(v=sd(coral$t1),lwd=4)
plot(density(coral_skew),lwd=2,col="gray",xlim=c(min(c(coral_skew,skewness(coral$t1))),max(c(coral_skew,skewness(coral$t1)))))
abline(v=skewness(coral$t1),lwd=4)
plot(density(coral_kurt),lwd=2,col="gray",xlim=c(min(c(coral_kurt,kurtosis(coral$t1))),max(c(coral_kurt,kurtosis(coral$t1)))))
abline(v=kurtosis(coral$t1),lwd=4)

## compare to more complex models with gamlss
## c.crit affects when the algorithm thinks it's converged, and n.cyc is the number of iterations
## all of these include size-dependence in sigma
con.opts <- gamlss.control(c.crit=0.001, n.cyc = 80)
## normal
coral_NO <- gamlss(t1 ~ t0, sigma.formula = ~t0, data = coral, family="NO",control=con.opts)
## skewed normal
coral_SN <- gamlss(t1 ~ t0, sigma.formula = ~t0, data = coral, family="SN1",control=con.opts)
refit(coral_SN)
## SHASH
coral_SHASH <- gamlss(t1 ~ t0, sigma.formula = ~t0, data = coral, family="SHASH",control=con.opts)
## skewed t
coral_ST <- gamlss(t1 ~ t0, sigma.formula = ~t0, data = coral, family="ST1",control=con.opts)
## AIC ranks
AICctab(coral_NO,coral_SN,coral_SHASH,coral_ST) ## something weird with the SHASH

## Can the SHASH do a better job if nu and tau are functions of size0?
coral_SHASH_all_sizedep <- gamlss(t1 ~ t0, sigma.formula = ~t0, nu.formula = ~t0 + I(t0^2), tau.formula =  ~t0 + I(t0^2),
                                  data = coral, family="SHASH",control=con.opts)
AICtab(coral_SHASH,coral_SHASH_all_sizedep)

## simulate data from best model
coral_sim_SHASH<-matrix(NA,nrow=nrow(coral),ncol=n_sim)
for(i in 1:n_sim){
  coral_sim_SHASH[,i] <- rSHASH(n = nrow(coral_NO_sim), mu = coral_SHASH_all_sizedep$mu.coefficients[1] + coral_SHASH_all_sizedep$mu.coefficients[2] * coral$t0, 
                                            sigma = exp(coral_SHASH_all_sizedep$sigma.coefficients[1] + coral_SHASH_all_sizedep$sigma.coefficients[2] * coral$t0),
                                            nu = exp(coral_SHASH_all_sizedep$nu.coefficients[1] + coral_SHASH_all_sizedep$nu.coefficients[2]*coral$t0 +
                                                       coral_SHASH_all_sizedep$nu.coefficients[3]*(coral$t0^2)),
                                            tau = exp(coral_SHASH_all_sizedep$tau.coefficients[1] + coral_SHASH_all_sizedep$tau.coefficients[2]*coral$t0 +
                                                        coral_SHASH_all_sizedep$tau.coefficients[3]*(coral$t0^2)))
}

## moments of the data by size bin
n_bins = 10
alpha_scale = 0.7
coral_moments <- coral %>% 
  arrange(t0) %>% 
  mutate(size_bin = cut_number(t0,n=n_bins)) %>% 
  group_by(size_bin) %>% 
  summarise(mean_t1 = mean(t1),
            sd_t1 = sd(t1),
            skew_t1 = skewness(t1),
            kurt_t1 = kurtosis(t1),
            bin_mean = mean(t0),
            bin_n = n()) 

## plot moments of real data and of many simulations from fitted model
par(mfrow=c(2,2))
plot(coral_moments$bin_mean,coral_moments$mean_t1,type="n",ylim=c(0,50),xlab="Mean size t0",ylab="mean(Size t1)")
for(i in 1:n_sim){
  sim_moments <- bind_cols(coral,data.frame(sim=coral_sim_SHASH[,i])) %>% 
    arrange(t0) %>% 
    mutate(size_bin = cut_number(t0,n=n_bins)) %>% 
    group_by(size_bin) %>% 
    summarise(mean_t1 = mean(sim),
              bin_mean = mean(t0))
  points(sim_moments$bin_mean,sim_moments$mean_t1,col=alpha("gray",0.5))
}
points(coral_moments$bin_mean,coral_moments$mean_t1,cex=3,pch=16,col=alpha("red",alpha_scale))

plot(coral_moments$bin_mean,coral_moments$sd_t1,type="n",ylim=c(0,20),xlab="Mean size t0",ylab="sd(Size t1)")
for(i in 1:n_sim){
  sim_moments <- bind_cols(coral,data.frame(sim=coral_sim_SHASH[,i])) %>% 
    arrange(t0) %>% 
    mutate(size_bin = cut_number(t0,n=n_bins)) %>% 
    group_by(size_bin) %>% 
    summarise(sd_t1 = sd(sim),
              bin_mean = mean(t0))
  points(sim_moments$bin_mean,sim_moments$sd_t1,col=alpha("gray",0.5))
}
points(coral_moments$bin_mean,coral_moments$sd_t1,cex=3,pch=16,col=alpha("red",alpha_scale))

plot(coral_moments$bin_mean,coral_moments$skew_t1,type="n",ylim=c(-5,2),xlab="Mean size t0",ylab="skewness(Size t1)")
for(i in 1:n_sim){
  sim_moments <- bind_cols(coral,data.frame(sim=coral_sim_SHASH[,i])) %>% 
    arrange(t0) %>% 
    mutate(size_bin = cut_number(t0,n=n_bins)) %>% 
    group_by(size_bin) %>% 
    summarise(skew_t1 = skewness(sim),
              bin_mean = mean(t0))
  points(sim_moments$bin_mean,sim_moments$skew_t1,col=alpha("gray",0.5))
}
points(coral_moments$bin_mean,coral_moments$skew_t1,cex=3,pch=16,col=alpha("red",alpha_scale))

plot(coral_moments$bin_mean,coral_moments$kurt_t1,type="n",ylim=c(0,30),xlab="Mean size t0",ylab="kurtosis(Size t1)")
for(i in 1:n_sim){
  sim_moments <- bind_cols(coral,data.frame(sim=coral_sim_SHASH[,i])) %>% 
    arrange(t0) %>% 
    mutate(size_bin = cut_number(t0,n=n_bins)) %>% 
    group_by(size_bin) %>% 
    summarise(kurt_t1 = kurtosis(sim),
              bin_mean = mean(t0))
  points(sim_moments$bin_mean,sim_moments$kurt_t1,col=alpha("gray",0.5))
}
points(coral_moments$bin_mean,coral_moments$kurt_t1,cex=3,pch=16,col=alpha("red",alpha_scale))



# Cactus data -------------------------------------------------------------

## compare random-effect fits of gamlss vs lme4
cactus <- growth_dat[["cactus"]] %>% 
  select(t0,t1,Year_t,Plot)

cactus_rfx_lme4 <- lmer(t1 ~ t0 + (1|Year_t), data=cactus, REML = F)
cactus_rfx_gamlss <- gamlss(t1 ~ t0 + random(as.factor(Year_t)), data = cactus, family="NO",control=con.opts)

fixef(cactus_rfx_lme4)
coef(cactus_rfx_gamlss)

plot(ranef(cactus_rfx_lme4)$Year_t[,1], getSmo(cactus_rfx_gamlss)$coef,type="n",
     xlab="lme4 random year effects",ylab="gamlss random year effects")
text(ranef(cactus_rfx_lme4)$Year_t[,1], getSmo(cactus_rfx_gamlss)$coef,
     labels = rownames(ranef(cactus_rfx_lme4)$Year_t))
abline(0,1)
