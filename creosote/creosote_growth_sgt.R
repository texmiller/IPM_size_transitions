rm(list=ls(all=TRUE));

# setwd("c:/repos/IPM_size_transitions/creosote"); #Steve
setwd("C:/Users/tm9/Desktop/git local/IPM_size_transitions/creosote"); #Tom

require(car); require(lme4); require(zoo); require(moments); require(mgcv); 
require(gamlss); require(gamlss.tr); require(AICcmodavg); 
require(lmerTest); require(tidyverse); require(maxLik); 
require(actuar); require(lattice); require(grid); require(scales);
require(sgt); require(formatR); require(popbio); require(bbmle)

# misc functions
invlogit <- function(x){exp(x)/(1+exp(x))}

# log kurtosis function for diagnostics
Lkurtosis=function(x) log(kurtosis(x)); 

# Steve's diagnostics functions
source("../Diagnostics.R")
# function for choosing distribution family
source("../fitChosenDists.R")

## read in data for Larrea tridentata (LATR) -- derived data frame generated here: https://github.com/TrevorHD/LTEncroachment/blob/master/04_CDataPrep.R (line 261)
LATR_full <- read.csv("CData.csv") %>% 
  #calculate volume
  mutate(#standardize weighted density to mean zero -- dropping this bc different data subsets woudl have different values
    #d.stand = (weighted.dens - mean(weighted.dens, na.rm = TRUE)) / sd(weighted.dens, na.rm = TRUE),
    #create unique transect as interaction of transect and site
    unique.transect = interaction(transect, site))
## prep a data subset for growth that drops rows missing either t or t1 size data
LATR_grow <- LATR_full  %>% drop_na(volume_t,volume_t1) %>% 
  ## also create log_volume as a new variable because gam doesn't like functions of variables as variables
  mutate(log_volume_t = log(volume_t),
         log_volume_t1 = log(volume_t1))

# first look at size transitions
plot(log(LATR_grow$volume_t),log(LATR_grow$volume_t1))

############################################################################
# Gaussian fits and model selection using mgcv 
############################################################################
## three candidate models for the mean: size only, size + density, or size, density, and size:density
## three candidates for variance: size only, size+density, fitted value (all the covariates plus rfx)
LATR_gam_models=list()
## Pilot fits, where sigma depends on initial size only
LATR_gam_models[[1]] <- gam(list(log_volume_t1 ~s(log_volume_t) + s(unique.transect,bs="re"), ~s(log_volume_t)), 
                            data=LATR_grow, gamma=1.4, family=gaulss())
LATR_gam_models[[2]] <- gam(list(log_volume_t1 ~s(log_volume_t) + s(weighted.dens) + s(unique.transect,bs="re"), ~s(log_volume_t)), 
                            data=LATR_grow, gamma=1.4, family=gaulss())                
LATR_gam_models[[3]] <- gam(list(log_volume_t1 ~s(log_volume_t) + s(weighted.dens) + weighted.dens:log_volume_t + s(unique.transect,bs="re"), ~s(log_volume_t)), 
                            data=LATR_grow, gamma=1.4, family=gaulss())  
## Fits where sigma depends on both initial size and density
LATR_gam_models[[4]] <- gam(list(log_volume_t1 ~s(log_volume_t) + s(unique.transect,bs="re"), ~s(log_volume_t) + s(weighted.dens)), 
                            data=LATR_grow, gamma=1.4, family=gaulss())
LATR_gam_models[[5]] <- gam(list(log_volume_t1 ~s(log_volume_t) + s(weighted.dens) + s(unique.transect,bs="re"), ~s(log_volume_t) + s(weighted.dens)), 
                            data=LATR_grow, gamma=1.4, family=gaulss())                
LATR_gam_models[[6]] <- gam(list(log_volume_t1 ~s(log_volume_t) + s(weighted.dens) + weighted.dens:log_volume_t + s(unique.transect,bs="re"), ~s(log_volume_t) + s(weighted.dens)), 
                            data=LATR_grow, gamma=1.4, family=gaulss()) 
## these models will be iterated to fit sigma as f(fitted value)
LATR_grow$fitted_vals = LATR_grow$log_volume_t 
LATR_gam_models[[7]] <- gam(list(log_volume_t1 ~s(log_volume_t) + s(unique.transect,bs="re"), ~s(fitted_vals)), 
                            data=LATR_grow, gamma=1.4, family=gaulss())
LATR_gam_models[[8]] <- gam(list(log_volume_t1 ~s(log_volume_t) + s(weighted.dens) + s(unique.transect,bs="re"), ~s(fitted_vals)), 
                            data=LATR_grow, gamma=1.4, family=gaulss())                
LATR_gam_models[[9]] <- gam(list(log_volume_t1 ~s(log_volume_t) + s(weighted.dens) + weighted.dens:log_volume_t + s(unique.transect,bs="re"), ~s(fitted_vals)), 
                            data=LATR_grow, gamma=1.4, family=gaulss())  
for(mod in 7:9) {
  fitGAU = LATR_gam_models[[mod]]
  fitted_all = predict(fitGAU,type="response",data=LATR);                  
  fitted_vals = new_fitted_vals = fitted_all[,1]; 
  weights = fitted_all[,2]; # what I call "weights" here are 1/sigma values; see ?gaulss for details.
  
  err=100; k=0; 
  while(err>10^(-6)) {
    LATR_grow$fitted_vals = new_fitted_vals; 
    fitGAU <- update(fitGAU); 
    fitted_all = predict(fitGAU,type="response",data=LATR_grow);   
    new_fitted_vals = fitted_all[,1]; new_weights = fitted_all[,2];
    err = weights - new_weights; err=sqrt(mean(err^2)); 
    weights = new_weights; 
    k=k+1; cat(k,err,"\n"); 
  }   
  LATR_gam_models[[mod]] =  fitGAU;
}

grow_aic <- AICtab(LATR_gam_models,base=T,sort=F) 
## Model 5 is the winner: mean ~ s(size) + s(density), sd ~ s(size) + s(density)

## define models 5 as our best Gaussian model
LATR_gam_model <- LATR_gam_models[[which.min(grow_aic$AIC)]]
saveRDS(LATR_gam_model,"LATR_grow_best.rds")
LATR_grow$fitted_vals = new_fitted_vals
## extract the linear predictor for the mean -- we'll use this later
LATR_Xp <- predict.gam(LATR_gam_model,type="lpmatrix")
## fitted coefficients
LATR_beta <- coef(LATR_gam_model)

##################################################################  
# Extract values of the fitted splines to explore their properties 
##################################################################
fitted_terms = predict(LATR_gam_model,type="terms")  

##### effect of initial size on mean of final size 
plot(LATR_grow$log_volume_t, fitted_terms[,"s(log_volume_t)"]) 
##### effect of weighted.dens on mean of final size 
plot(LATR_grow$weighted.dens, fitted_terms[,"s(weighted.dens)"]); ## complicated 

##### sigma versus size - presumably log scale
plot(LATR_grow$log_volume_t, fitted_terms[,"s.1(log_volume_t)"]) 
##### sigma versus density 
plot(LATR_grow$weighted.dens, fitted_terms[,"s.1(weighted.dens)"]) 

##################################################################  
# Inspect scaled residuals to evaluate the pilot model: FAILS 
##################################################################
fitted_all = predict(LATR_gam_model,type="response")  
sigma.hat = 1/fitted_all[,2]
scaledResids = residuals(LATR_gam_model,type="response")/sigma.hat;  # note the 'type' argument is needed
par(mfrow=c(1,2))
plot(fitted_all[,1], scaledResids,xlab="Fitted values", ylab="Scaled residuals") 

qqPlot(scaledResids) # really bad in both tails
jarque.test(scaledResids) # normality test: FAILS, P < 0.001 
anscombe.test(scaledResids) # kurtosis: FAILS, P < 0.001 
agostino.test(scaledResids) # skewness: FAILS, P<0.001 

######## rolling NP moments diagnostic: skew is small but variable, tails are fat. 
px = LATR_grow$fitted_vals; py=scaledResids; 
z = rollMomentsNP(px,py,windows=8,smooth=TRUE,scaled=TRUE) 

###########################################################
## fit sgt using design matrices from normal gam
###########################################################

sgtLogLik=function(pars,response){
  val = dsgt(x = response, 
             mu=LATR_Xp[,1:31]%*%pars[1:31],
             sigma=exp(LATR_Xp[,32:50]%*%pars[32:50]),
             lambda=-invlogit(pars[51]+pars[52]*LATR_grow$log_volume_t),
             p=exp(pars[53]),
             q=exp(pars[54]),
             mean.cent=T,
             var.adj=T,
             log=T) 
  return(val); 
}
## initial parameter values
p0=c(LATR_beta,-10,0,2,2) 
out=maxLik(logLik=sgtLogLik,start=p0*exp(0.2*rnorm(length(p0))), response=LATR_grow$log_volume_t1,
           method="BHHH",control=list(iterlim=5000,printLevel=2),finalHessian=FALSE); 

out=maxLik(logLik=sgtLogLik,start=out$estimate,response=LATR_grow$log_volume_t1,
           method="NM",control=list(iterlim=5000,printLevel=1),finalHessian=FALSE); 

out=maxLik(logLik=sgtLogLik,start=out$estimate,response=LATR_grow$log_volume_t1,
           method="BHHH",control=list(iterlim=5000,printLevel=2),finalHessian=FALSE); 

out=maxLik(logLik=sgtLogLik,start=out$estimate,response=LATR_grow$log_volume_t1,
           method="BHHH",control=list(iterlim=5000,printLevel=2),finalHessian=TRUE) 

## compare to original (gaussian) gam parameter estimates
plot(LATR_beta,out$estimate[1:50])
abline(0,1)
## compare the expected value of the two models
plot(LATR_Xp[,1:31]%*%LATR_beta[1:31],LATR_Xp[,1:31]%*%out$estimate[1:31])
abline(0,1)
## I'm satisfied that the sgt can recover the same expected value as the gaussian gam()

########################################################
# compare simulated and real data -----------------------------------------
###########################################################
# Simulate data from normal and sgt models
n_sim <- 500
LATR_sim_NO<-LATR_sim_SGT<-matrix(NA,nrow=nrow(LATR_grow),ncol=n_sim)
for(i in 1:n_sim){
  print(i)
  LATR_sim_SGT[,i] <- rsgt(n = nrow(LATR_grow), 
                           mu = LATR_Xp[,1:31]%*%out$estimate[1:31], 
                           sigma = exp(LATR_Xp[,32:50]%*%out$estimate[32:50]),
                           lambda=-invlogit(out$estimate[51]+out$estimate[52]*LATR_grow$log_volume_t),
                           p=exp(out$estimate[53]),
                           q=exp(out$estimate[54]),
                           mean.cent=T,
                           var.adj=T)
  LATR_sim_NO[,i] <- rnorm(n = nrow(LATR_grow),
                           mean = fitted_all[,1],
                           sd = 1/fitted_all[,2])
}

n_bins = 12
alpha_scale = 0.7
LATR_moments <- LATR_grow %>% 
  arrange(log_volume_t) %>% 
  mutate(size_bin = cut_interval(log_volume_t,n=n_bins)) %>% 
  group_by(size_bin) %>% 
  summarise(mean_t1 = mean(log_volume_t1),
            sd_t1 = sd(log_volume_t1),
            skew_t1 = NPskewness(log_volume_t1),
            kurt_t1 = NPkurtosis(log_volume_t1),
            bin_mean = mean(log_volume_t),
            bin_n = n()) 

par(mfrow=c(2,2),mar=c(4,4,2,1),cex.axis=1.3,cex.lab=1.3,mgp=c(2,1,0),bty="l"); 
sim_bin_means=sim_moment_means=sim_moment_means_norm = matrix(NA,n_bins,n_sim); 
for(i in 1:n_sim){
  sim_moments <- bind_cols(LATR_grow,data.frame(sim=LATR_sim_SGT[,i],
                                                sim_norm=LATR_sim_NO[,i])) %>% 
    arrange(log_volume_t) %>% 
    mutate(size_bin = cut_interval(log_volume_t,n=n_bins)) %>% 
    group_by(size_bin) %>% 
    summarise(mean_t1 = mean(sim),
              mean_t1_norm = mean(sim_norm),
              bin_mean = mean(log_volume_t))
  sim_bin_means[,i]=sim_moments$bin_mean; 
  sim_moment_means[,i]=sim_moments$mean_t1; sim_moment_means_norm[,i]=sim_moments$mean_t1_norm;		  
}
matplot(LATR_moments$bin_mean, sim_moment_means,col=alpha("gray",0.5),pch=16,xlab="Mean size t0",ylab="Mean(Size t1)",cex=1,
        xlim=c(min(LATR_moments$bin_mean),max(LATR_moments$bin_mean)+0.4))
matplot(LATR_moments$bin_mean+0.4, sim_moment_means_norm,col=alpha("cornflowerblue",0.5),pch=16,add=T)
points(LATR_moments$bin_mean+0.2, LATR_moments$mean_t1,pch=16,lwd=2,col=alpha("red",alpha_scale),cex=1.6)
points(LATR_moments$bin_mean, apply(sim_moment_means,1,median),pch=1,lwd=2,col=alpha("black",alpha_scale),cex=1.6)
points(LATR_moments$bin_mean+0.4, apply(sim_moment_means_norm,1,median),pch=1,lwd=2,col=alpha("black",alpha_scale),cex=1.6)
legend("topleft",legend=c("SGT","Gaussian","Data"),
       col=c("gray","cornflowerblue","red"),pch=16,bty="n",cex=1,pt.lwd=2,pt.cex = 1.2) 
add_panel_label("a")

for(i in 1:n_sim){
  sim_moments <- bind_cols(LATR_grow,data.frame(sim=LATR_sim_SGT[,i],
                                                sim_norm=LATR_sim_NO[,i])) %>% 
    arrange(log_volume_t) %>% 
    mutate(size_bin = cut_interval(log_volume_t,n=n_bins)) %>% 
    group_by(size_bin) %>% 
    summarise(mean_t1 = sd(sim),
              mean_t1_norm = sd(sim_norm),
              bin_mean = mean(log_volume_t))
  sim_bin_means[,i]=sim_moments$bin_mean; 
  sim_moment_means[,i]=sim_moments$mean_t1; sim_moment_means_norm[,i]=sim_moments$mean_t1_norm;		  
}
matplot(LATR_moments$bin_mean, sim_moment_means,col=alpha("gray",0.5),pch=16,xlab="Mean size t0",ylab="SD(Size t1)",cex=1,
        xlim=c(min(LATR_moments$bin_mean),max(LATR_moments$bin_mean)+0.4)) 
matplot(LATR_moments$bin_mean+0.4, sim_moment_means_norm,col=alpha("cornflowerblue",0.5),pch=16,add=T)
points(LATR_moments$bin_mean+0.2, LATR_moments$sd_t1,pch=16,lwd=2,col=alpha("red",alpha_scale),cex=1.6)
points(LATR_moments$bin_mean, apply(sim_moment_means,1,median),pch=1,lwd=2,col=alpha("black",alpha_scale),cex=1.6)
points(LATR_moments$bin_mean+0.4, apply(sim_moment_means_norm,1,median),pch=1,lwd=2,col=alpha("black",alpha_scale),cex=1.6)
add_panel_label("b")

for(i in 1:n_sim){
  sim_moments <- bind_cols(LATR_grow,data.frame(sim=LATR_sim_SGT[,i],
                                                sim_norm=LATR_sim_NO[,i])) %>% 
    arrange(log_volume_t) %>% 
    mutate(size_bin = cut_interval(log_volume_t,n=n_bins)) %>% 
    group_by(size_bin) %>% 
    summarise(mean_t1 = NPskewness(sim),
              mean_t1_norm = NPskewness(sim_norm),
              bin_mean = mean(log_volume_t))
  sim_bin_means[,i]=sim_moments$bin_mean; 
  sim_moment_means[,i]=sim_moments$mean_t1; sim_moment_means_norm[,i]=sim_moments$mean_t1_norm;	  
}
matplot(LATR_moments$bin_mean, sim_moment_means,col=alpha("gray",0.5),pch=16,xlab="Mean size t0",ylab="Skew(Size t1)",cex=1,
        xlim=c(min(LATR_moments$bin_mean),max(LATR_moments$bin_mean)+0.4))
matplot(LATR_moments$bin_mean+0.4, sim_moment_means_norm,col=alpha("cornflowerblue",0.5),pch=16,add=T)
points(LATR_moments$bin_mean+0.2, LATR_moments$skew_t1,pch=16,lwd=2,col=alpha("red",alpha_scale),cex=1.6)
points(LATR_moments$bin_mean, apply(sim_moment_means,1,median),pch=1,lwd=2,col=alpha("black",alpha_scale),cex=1.6)
points(LATR_moments$bin_mean+0.4, apply(sim_moment_means_norm,1,median),pch=1,lwd=2,col=alpha("black",alpha_scale),cex=1.6)
add_panel_label("c")

for(i in 1:n_sim){
  sim_moments <- bind_cols(LATR_grow,data.frame(sim=LATR_sim_SGT[,i],
                                                sim_norm=LATR_sim_NO[,i])) %>% 
    arrange(log_volume_t) %>% 
    mutate(size_bin = cut_interval(log_volume_t,n=n_bins)) %>% 
    group_by(size_bin) %>% 
    summarise(mean_t1 = NPkurtosis(sim),
              mean_t1_norm = NPkurtosis(sim_norm),
              bin_mean = mean(log_volume_t))
  sim_bin_means[,i]=sim_moments$bin_mean; 
  sim_moment_means[,i]=sim_moments$mean_t1; sim_moment_means_norm[,i]=sim_moments$mean_t1_norm;	  
}
matplot(LATR_moments$bin_mean, sim_moment_means,col=alpha("gray",0.5),pch=16,xlab="Mean size t0",ylab="Kurtosis(Size t1)",cex=1,
        xlim=c(min(LATR_moments$bin_mean),max(LATR_moments$bin_mean)+0.4))
matplot(LATR_moments$bin_mean+0.4, sim_moment_means_norm,col=alpha("cornflowerblue",0.5),pch=16,add=T)
points(LATR_moments$bin_mean+0.2, LATR_moments$kurt_t1,pch=16,lwd=2,col=alpha("red",alpha_scale),cex=1.6)
points(LATR_moments$bin_mean, apply(sim_moment_means,1,median),pch=1,lwd=2,col=alpha("black",alpha_scale),cex=1.6)
points(LATR_moments$bin_mean+0.4, apply(sim_moment_means_norm,1,median),pch=1,lwd=2,col=alpha("black",alpha_scale),cex=1.6)
add_panel_label("d")

