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

################################################################################
## Finding a better distribution. Must at least allow positive excess kurtosis 
#################################################################################

cand_dist=c("NO","GT","LO","TF","ST1","ST2") ##NET, SEP1, and EGB2 are all highly unstable
n_bins <- 8
select_dist <- tibble(init_size = LATR_grow$log_volume_t,
                      scale_resid = scaledResids) %>% 
  mutate(bin = as.integer(cut_number(init_size,n_bins)),
         best_dist = NA,
         secondbest_dist = NA,
         aic_margin = NA) 
for(b in 1:n_bins){
  bin_fit <- gamlssMaxlik(y=select_dist$scale_resid[select_dist$bin==b],DIST=cand_dist)
  select_dist$best_dist[select_dist$bin==b] <- cand_dist[which(bin_fit$aics==sort(bin_fit$aics)[1])]#names(bin_fit$fits[1])
  select_dist$secondbest_dist[select_dist$bin==b] <- cand_dist[which(bin_fit$aics==sort(bin_fit$aics)[2])]#names(bin_fit$fits[2])
  select_dist$aic_margin[select_dist$bin==b] <- sort(bin_fit$aics)[2] - sort(bin_fit$aics)[1]
}
## keep in mind here that we are ignoring density variation, which can affect variance
select_dist %>% 
  group_by(bin) %>% 
  summarise(n_bin = n(),
            best_dist = unique(best_dist),
            secondbest_dist = unique(secondbest_dist),
            aic_margin = unique(aic_margin))
## Going ahead with ST1, which can hopefully get both the skew and kurtosis

## visualize ST1 parameters in relation to fitted value by bin
LATR_bin_fit <-LATR_grow %>% 
  mutate(init_size = log_volume_t,
         bin = as.integer(cut_number(init_size,n_bins)),
         dens_bin = as.integer(cut_number(weighted.dens,2))) %>% 
  mutate(mu_lowdens=NA,sigma_lowdens=NA,nu_lowdens=NA,tau_lowdens=NA,
         mu_highdens=NA,sigma_highdens=NA,nu_highdens=NA,tau_highdens=NA)
for(b in 1:n_bins){
  bin_fit_lowdens <- gamlssMaxlik(y=LATR_bin_fit$log_volume_t1[LATR_bin_fit$bin==b & LATR_bin_fit$dens_bin==1],DIST="ST1")
  bin_fit_highdens <- gamlssMaxlik(y=LATR_bin_fit$log_volume_t1[LATR_bin_fit$bin==b & LATR_bin_fit$dens_bin==2],DIST="ST1")
  ## ST1 links are identity-log-identity-log
  LATR_bin_fit$mu_lowdens[LATR_bin_fit$bin==b] <- bin_fit_lowdens$out[[1]]$estimate["eta.mu"] 
  LATR_bin_fit$sigma_lowdens[LATR_bin_fit$bin==b] <- exp(bin_fit_lowdens$out[[1]]$estimate["eta.sigma"])
  LATR_bin_fit$nu_lowdens[LATR_bin_fit$bin==b] <- (bin_fit_lowdens$out[[1]]$estimate["eta.nu"]) 
  LATR_bin_fit$tau_lowdens[LATR_bin_fit$bin==b] <- exp(bin_fit_lowdens$out[[1]]$estimate["eta.tau"]) 
  LATR_bin_fit$mu_highdens[LATR_bin_fit$bin==b] <- bin_fit_highdens$out[[1]]$estimate["eta.mu"] 
  LATR_bin_fit$sigma_highdens[LATR_bin_fit$bin==b] <- exp(bin_fit_highdens$out[[1]]$estimate["eta.sigma"]) 
  LATR_bin_fit$nu_highdens[LATR_bin_fit$bin==b] <- (bin_fit_highdens$out[[1]]$estimate["eta.nu"]) 
  LATR_bin_fit$tau_highdens[LATR_bin_fit$bin==b] <- exp(bin_fit_highdens$out[[1]]$estimate["eta.tau"]) 
  }
LATR_bin_fit %>% 
  group_by(bin) %>% 
  summarise(N = n(),
            mean_size = mean(init_size),
            mu_lowdens = unique(mu_lowdens),
            sigma_lowdens=unique(sigma_lowdens),
            nu_lowdens=unique(nu_lowdens),
            tau_lowdens=unique(tau_lowdens),
            mu_highdens = unique(mu_highdens),
            sigma_highdens=unique(sigma_highdens),
            nu_highdens=unique(nu_highdens),
            tau_highdens=unique(tau_highdens)) -> LATR_bin_fit

par(mfrow=c(2,2),bty="l",mar=c(4,4,2,1),mgp=c(2.2,1,0),cex.axis=1.4,cex.lab=1.4);
plot(LATR_bin_fit$mean_size,LATR_bin_fit$mu_lowdens,xlab="Initial size",ylab=expression(paste("ST parameter  ", mu )),type="b",pch=1)
points(LATR_bin_fit$mean_size,LATR_bin_fit$mu_highdens,type="b",pch=16)

plot(LATR_bin_fit$mean_size,LATR_bin_fit$sigma_lowdens,xlab="Initial size",
     ylab=expression(paste("ST  parameter  ", sigma)),type="b",pch=1)
points(LATR_bin_fit$mean_size,LATR_bin_fit$sigma_highdens,type="b",pch=16)

plot(LATR_bin_fit$mean_size,LATR_bin_fit$nu_lowdens,xlab="Initial size",
     ylab=expression(paste("ST  parameter  ", nu)),type="b",pch=1,ylim=c(-20,20))
points(LATR_bin_fit$mean_size,LATR_bin_fit$nu_highdens,type="b",pch=16)

plot(LATR_bin_fit$mean_size,LATR_bin_fit$tau_highdens,xlab="Initial size",
     ylab=expression(paste("ST  parameter  ", tau)),type="b",pch=16,ylim=c(-20,50))
points(LATR_bin_fit$mean_size,LATR_bin_fit$tau_lowdens,xlab="Initial size",
     ylab=expression(paste("ST  parameter  ", tau)),type="b",pch=1)
## nu and tau have some weirdly unstable values

# Fitting the final model -------------------------------------------------
## we will create basis functions for size and density, 
## for now going with a basis dimension of 10 for each (can revisit with AIC selection)
## sigma and nu also look like they should have size and density dependence
## tau is causing trouble in the binned fits, so I will fit a constant value

## ST1 likelihood -- with size and density dependence in sigma and nu
LogLik=function(pars,response,U){
  pars1 = pars[1:ncol(U)]; pars2=pars[-(1:ncol(U))];
  mu = U%*%pars1;  
  val = dST1(x = response, 
             mu=mu,
             sigma = exp(pars2[1] + pars2[2]*LATR_grow$log_volume_t + pars2[3]*LATR_grow$weighted.dens),
             nu = pars2[4] + pars2[5]*LATR_grow$log_volume_t + pars2[6]*LATR_grow$weighted.dens,
             tau = exp(pars2[7]),
             log=T) 
  return(val); 
}

## design matrix for random effects on the mean
R <- model.matrix(~LATR_grow$unique.transect-1)
## design matrix for size effect on mean (smooth term) -- can and will choose k by model selection. but not now. 
S_size <- smoothCon(s(log_volume_t,k=10),data=LATR_grow,absorb.cons=TRUE)[[1]]$X
## design matrix for density effect on mean (smooth term) -- can and will choose k by model selection. but not now. 
S_dens <- smoothCon(s(weighted.dens,k=10),data=LATR_grow,absorb.cons=TRUE)[[1]]$X
## bind together in a single design matrix for the mean
Xb = cbind(R,S_size,S_dens)

# Starting values for smooths and fac levels -- starting "blind" here
fixed_start = rep(0,times=ncol(Xb))
## starting values for sigma, nu, tau
fit_sigma = lm(log(sigma_highdens)~mean_size, data=LATR_bin_fit)
fit_nu = lm(nu_highdens~mean_size, data=LATR_bin_fit)
fit_tau = lm(log(tau_highdens)~1, data=LATR_bin_fit)
## bundle coefficients for mu and sigma
p0=c(fixed_start,c(coef(fit_sigma),0),c(coef(fit_nu),0),5) 

paranoid_iter <- 3
coefs = list(paranoid_iter); LL=numeric(paranoid_iter);  
for(j in 1:paranoid_iter) {
  out=maxLik(logLik=LogLik,start=p0*exp(0.2*rnorm(length(p0))), response=LATR_grow$log_volume_t1,U=Xb,
             method="BHHH",control=list(iterlim=5000,printLevel=2),finalHessian=FALSE); 
  
  out=maxLik(logLik=LogLik,start=out$estimate,response=LATR_grow$log_volume_t1,U=Xb,
             method="NM",control=list(iterlim=5000,printLevel=1),finalHessian=FALSE); 
  
  out=maxLik(logLik=LogLik,start=out$estimate,response=LATR_grow$log_volume_t1,U=Xb,
             method="BHHH",control=list(iterlim=5000,printLevel=2),finalHessian=FALSE); 
  
  coefs[[j]] = out$estimate; LL[j] = out$maximum;
  cat(j, "#--------------------------------------#",out$maximum,"\n"); 
}

j = min(which(LL==max(LL))) 
out=maxLik(logLik=LogLik,start=coefs[[j]],response=LATR_grow$log_volume_t1,U=Xb,
           method="BHHH",control=list(iterlim=5000,printLevel=2),finalHessian=TRUE) 

# AIC 
(AIC_ST1 <- 2*length(coefs) - 2*out$maximum)

######### save results of ML fit.  
names(out$estimate)<-c(colnames(Xb),"sigma_b0","sigma_b1","sigma_b2","nu_b0","nu_b1","nu_b2","tau_b0")
coefs=out$estimate # parameters
V = vcov(out); SEs = sqrt(diag(V)); # standard errors

############# RFX shrinkage -- NEW, following the theory and code in Steve's appendix
n_transects <- length(LATR_gam_model$coefficients[paste0("s(unique.transect).",1:12)])
transects=1:n_transects ## these are the indices of transect estimates in the parameter estimate vector
# Variance-covariance matrices for random intercepts
V1 = V[transects,transects]; 
# Extract year-specific intercepts, center them to zero
fixed.fx = coefs[transects]; fixed.fx = fixed.fx-mean(fixed.fx);
# Estimate sigma^2
var.hat = mean(fixed.fx^2) - mean(diag(V1)) + (sum(V1)-sum(diag(V1)))/(2*n_transects*(n_transects-1)); ## still comes out negative
# Shrink deviations from the mean
shrinkRanIntercept = fixed.fx*sqrt(var.hat/(var.hat + diag(V1)));
## compare to gam re's
plot(LATR_gam_model$coefficients[paste0("s(unique.transect).",1:12)],
     shrinkRanIntercept);abline(0,1)

# compare simulated and real data -----------------------------------------
# Simulate data from fitted ST1 model
MLmu = Xb%*%coefs[1:ncol(Xb)] 
n_sim <- 500
LATR_sim_NO<-LATR_sim_ST1<-matrix(NA,nrow=nrow(LATR_grow),ncol=n_sim)
for(i in 1:n_sim){
  print(i)
  LATR_sim_ST1[,i] <- rST1(n = nrow(LATR_grow), 
                           mu = MLmu, 
                           sigma = exp(coefs["sigma_b0"] + coefs["sigma_b1"]*LATR_grow$log_volume_t  + coefs["sigma_b2"]*LATR_grow$weighted.dens),
                           nu=coefs["nu_b0"] + coefs["nu_b1"]*LATR_grow$log_volume_t  + coefs["nu_b2"]*LATR_grow$weighted.dens,
                           tau=exp(coefs["tau_b0"]))
  LATR_sim_NO[,i] <- rnorm(n = nrow(LATR_grow),
                           mean = fitted_all[,1],
                           sd = 1/fitted_all[,2])
}

n_bins = 10
alpha_scale = 0.7
LATR_moments <- LATR_grow %>% 
  arrange(log_volume_t) %>% 
  mutate(size_bin = cut_number(log_volume_t,n=n_bins)) %>% 
  group_by(size_bin) %>% 
  summarise(mean_t1 = mean(log_volume_t1),
            sd_t1 = sd(log_volume_t1),
            skew_t1 = NPskewness(log_volume_t1),
            kurt_t1 = NPkurtosis(log_volume_t1),
            bin_mean = mean(log_volume_t),
            bin_n = n()) 

pdf("../manuscript/figures/creosote_sim_moments_ST1.pdf",height = 10,width = 10,useDingbats = F)
par(mfrow=c(2,2),mar=c(4,4,2,1),cex.axis=1.3,cex.lab=1.3,mgp=c(2,1,0),bty="l"); 
sim_bin_means=sim_moment_means=sim_moment_means_norm = matrix(NA,n_bins,n_sim); 
for(i in 1:n_sim){
  sim_moments <- bind_cols(LATR_grow,data.frame(sim=LATR_sim_ST1[,i],
                                           sim_norm=LATR_sim_NO[,i])) %>% 
    arrange(log_volume_t) %>% 
    mutate(size_bin = cut_number(log_volume_t,n=n_bins)) %>% 
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
legend("topleft",legend=c("Skewed t model","Gaussian model","Data"),
       col=c("gray","cornflowerblue","red"),pch=16,bty="n",cex=1.4,pt.lwd=2,pt.cex = 1.6) 
add_panel_label("a")

for(i in 1:n_sim){
  sim_moments <- bind_cols(LATR_grow,data.frame(sim=LATR_sim_ST1[,i],
                                           sim_norm=LATR_sim_NO[,i])) %>% 
    arrange(log_volume_t) %>% 
    mutate(size_bin = cut_number(log_volume_t,n=n_bins)) %>% 
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
  sim_moments <- bind_cols(LATR_grow,data.frame(sim=LATR_sim_ST1[,i],
                                           sim_norm=LATR_sim_NO[,i])) %>% 
    arrange(log_volume_t) %>% 
    mutate(size_bin = cut_number(log_volume_t,n=n_bins)) %>% 
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
  sim_moments <- bind_cols(LATR_grow,data.frame(sim=LATR_sim_ST1[,i],
                                           sim_norm=LATR_sim_NO[,i])) %>% 
    arrange(log_volume_t) %>% 
    mutate(size_bin = cut_number(log_volume_t,n=n_bins)) %>% 
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
dev.off()

# IPM ---------------------------------------------------------------------
## the improvement of the skewed t (or the logistic for that matter) over the Gaussian is not particularly
## striking in the binned simulated moment figures. But let's see how much this changes the IPM.

## First, we need functions for survival and reproduction, and I will use gam() here following the models above

# Flowering ---------------------------------------------------------------
# populate year t of 2017-2018 transition year (there are no 2018 data but this way we get all four years in the reproduction models)
# I'll do this by creating the 2017-18 data as a stand-alone df then bind rows
LATR_dat_201718 <- LATR_full[LATR_full$year_t==2016 & LATR_full$survival_t1==1,]
## these are the 2017 survivors. Make their year t demography last year's data
LATR_dat_201718$year_t<-2017;LATR_dat_201718$year_t1<-2018
LATR_dat_201718$max.ht_t<-LATR_dat_201718$max.ht_t1
LATR_dat_201718$max.w_t<-LATR_dat_201718$max.w_t1
LATR_dat_201718$volume_t<-LATR_dat_201718$volume_t1
LATR_dat_201718$perp.w_t<-LATR_dat_201718$perp.w_t1
LATR_dat_201718$flowers_t<-LATR_dat_201718$flowers_t1
LATR_dat_201718$fruits_t<-LATR_dat_201718$fruits_t1
LATR_dat_201718$reproductive_fraction_t<-LATR_dat_201718$reproductive_fraction_t1
LATR_dat_201718$total.reproduction_t<-LATR_dat_201718$total.reproduction_t1
## now set all the t1 data to NA
LATR_dat_201718$max.ht_t1<-NA
LATR_dat_201718$max.w_t1<-NA
LATR_dat_201718$volume_t1<-NA
LATR_dat_201718$perp.w_t1<-NA
LATR_dat_201718$flowers_t1<-NA
LATR_dat_201718$fruits_t1<-NA
LATR_dat_201718$reproductive_fraction_t1<-NA
LATR_dat_201718$total.reproduction_t1<-NA
## bind rows
LATR_flow_dat <- bind_rows(LATR_full,LATR_dat_201718) %>% 
  select(unique.transect,volume_t,total.reproduction_t,weighted.dens) %>% drop_na()
## create log_vol as new variables (easier for gams)
LATR_flow_dat$log_volume_t <- log(LATR_flow_dat$volume_t)

LATR_flower <- list()
LATR_flower[[1]] <-  gam(total.reproduction_t>0 ~ s(log_volume_t) + s(unique.transect,bs="re"),
                   data=LATR_flow_dat, gamma=1.4, family="binomial")
LATR_flower[[2]] <-  gam(total.reproduction_t>0 ~ s(log_volume_t) + s(weighted.dens) + s(unique.transect,bs="re"),
                         data=LATR_flow_dat, gamma=1.4, family="binomial")
LATR_flower[[3]] <-  gam(total.reproduction_t>0 ~ s(log_volume_t) + s(weighted.dens) + weighted.dens:log_volume_t  + s(unique.transect,bs="re"),
                         data=LATR_flow_dat, gamma=1.4, family="binomial")
flower_aic<-AICtab(LATR_flower,base=T,sort=F)
LATR_flower_best <- LATR_flower[[which.min(flower_aic$AIC)]]
saveRDS(LATR_flower_best,"LATR_flower_best.rds")
LATR_flower_fitted_terms = predict(LATR_flower_best,type="terms") 
LATR_flow_dat$pred = predict.gam(LATR_flower_best,newdata = LATR_flow_dat, exclude = "s(unique.transect)")

##### effect of size on pr(flower) -- fairly linear
plot(LATR_flow_dat$log_volume_t,LATR_flower_fitted_terms[,"s(log_volume_t)"]) 
#### effect of d.stand on pr(flower) -- linear 
plot(LATR_flow_dat$weighted.dens,LATR_flower_fitted_terms[,"s(weighted.dens)"]) 

## visualize data + model
n_cuts_dens <- 6
n_cuts_size <- 4
LATR_flow_dat %>% 
  mutate(size_bin = as.integer(cut_number(log_volume_t,n_cuts_size)),
         dens_bin = as.integer(cut_number(weighted.dens,n_cuts_dens))) %>% 
  group_by(size_bin,dens_bin) %>% 
  summarise(mean_size = mean(log_volume_t),
         mean_density = mean(weighted.dens),
         mean_flower = mean(total.reproduction_t > 0),
         pred_flower = mean(pred),
         bin_n = n()) -> LATR_flow_dat_plot

## generate predictions for plotting
size_means_flow <- LATR_flow_dat_plot %>% group_by(size_bin) %>% summarise(mean_size=mean(mean_size))
LATR_flow_pred <- data.frame(
  weighted.dens = rep(seq(min(LATR_flow_dat$weighted.dens),max(LATR_flow_dat$weighted.dens),length.out = 20),times=n_cuts_size),
  log_volume_t = rep(size_means_flow$mean_size,each=20),
  unique.transect=1,
  size_bin = rep(size_means_flow$size_bin,each=20)
)
LATR_flow_pred$pred <- predict.gam(LATR_flower_best,newdata = LATR_flow_pred, exclude = "s(unique.transect)")

plot(LATR_flow_dat_plot$mean_density,LATR_flow_dat_plot$mean_flower,type="n",ylim=c(0,1),
     xlab="Weighted density",ylab="Pr(Flowering)")
for(i in 1:n_cuts_size){
  points(LATR_flow_dat_plot$mean_density[LATR_flow_dat_plot$size_bin==i],
         LATR_flow_dat_plot$mean_flower[LATR_flow_dat_plot$size_bin==i],pch=16,col=i,
         cex=(LATR_flow_dat_plot$bin_n[LATR_flow_dat_plot$size_bin==i]/max(LATR_flow_dat_plot$bin_n))*3)
  lines(LATR_flow_pred$weighted.dens[LATR_flow_pred$size_bin==i],
      invlogit(LATR_flow_pred$pred[LATR_flow_pred$size_bin==i]),col=i)
}


# Fruit production --------------------------------------------------------
LATR_fruits_dat <- subset(LATR_flow_dat,total.reproduction_t>0)
LATR_fruits <- list()
LATR_fruits[[1]] <-  gam(total.reproduction_t ~ s(log_volume_t) + s(unique.transect,bs="re"),
                         data=LATR_fruits_dat, gamma=1.4, family="nb")
LATR_fruits[[2]] <-  gam(total.reproduction_t ~ s(log_volume_t) + s(weighted.dens) + s(unique.transect,bs="re"),
                         data=LATR_fruits_dat, gamma=1.4, family="nb")
LATR_fruits[[3]] <-  gam(total.reproduction_t ~ s(log_volume_t) + s(weighted.dens) + weighted.dens:log_volume_t  + s(unique.transect,bs="re"),
                         data=LATR_fruits_dat, gamma=1.4, family="nb")
fruits_aic<-AICtab(LATR_fruits,base=T,sort=F)
LATR_fruits_best <- LATR_fruits[[which.min(fruits_aic$AIC)]]
saveRDS(LATR_fruits_best,"LATR_fruits_best.rds")
LATR_fruits_fitted_terms = predict(LATR_fruits_best,type="terms") 
LATR_fruits_dat$pred = predict.gam(LATR_fruits_best,newdata = LATR_fruits_dat,exclude="s(unique.transect)")

##### effect of size on pr(flower) -- linear, positive
plot(LATR_fruits_dat$log_volume_t,LATR_fruits_fitted_terms[,"s(log_volume_t)"]) 
#### effect of d.stand on pr(flower) -- negative, a bit of a kink
plot(LATR_fruits_dat$weighted.dens,LATR_fruits_fitted_terms[,"s(weighted.dens)"]) 

LATR_fruits_dat %>% 
  mutate(size_bin = as.integer(cut_number(log_volume_t,n_cuts_size)),
         dens_bin = as.integer(cut_number(weighted.dens,n_cuts_dens))) %>% 
  group_by(size_bin,dens_bin) %>% 
  summarise(mean_size = mean(log_volume_t),
         mean_density = mean(weighted.dens),
         mean_fruits = mean(total.reproduction_t),
         pred_fruits = mean(pred),
         bin_n = n()) -> LATR_fruits_dat_plot

## new data set for gam prediction
size_means_fruit <- LATR_fruits_dat_plot %>% group_by(size_bin) %>% summarise(mean_size=mean(mean_size))
LATR_fruit_pred <- data.frame(
  weighted.dens = rep(seq(min(LATR_fruits_dat$weighted.dens),max(LATR_fruits_dat$weighted.dens),length.out = 20),times=n_cuts_size),
  log_volume_t = rep(size_means_fruit$mean_size,each=20),
  unique.transect=1,
  size_bin = rep(size_means_fruit$size_bin,each=20)
)
LATR_fruit_pred$pred <- predict.gam(LATR_fruits_best,newdata = LATR_fruit_pred, exclude = "s(unique.transect)")


plot(LATR_fruits_dat_plot$mean_density,LATR_fruits_dat_plot$mean_fruits,type="n",
     xlab="Weighted density",ylab="Flowers and Fruits")
for(i in 1:n_cuts_size){
  points(LATR_fruits_dat_plot$mean_density[LATR_fruits_dat_plot$size_bin==i],
         LATR_fruits_dat_plot$mean_fruits[LATR_fruits_dat_plot$size_bin==i],pch=16,col=i,
         cex=(LATR_fruits_dat_plot$bin_n[LATR_fruits_dat_plot$size_bin==i]/max(LATR_fruits_dat_plot$bin_n))*3)
  lines(LATR_fruit_pred$weighted.dens[LATR_fruit_pred$size_bin==i],
        exp(LATR_fruit_pred$pred[LATR_fruit_pred$size_bin==i]),col=i)
}


# Survival ----------------------------------------------------------------
# read in transplant experiment and merge with obs data
# Combine transplants with large shrubs for later survival analysis
# Keep only location info, survival, volume, and density
CData.Transplants<-read.csv("CData.Transplants.csv")%>% 
  select("site", "transect", "actual.window", 
       "spring_survival_t1", "volume_t", "weighted.dens", "transplant") %>% 
  rename("survival_t1" = "spring_survival_t1") %>% 
  mutate(unique.transect = interaction(transect, site)) %>% 
  rbind(select(LATR_full, "site", "transect", "actual.window", 
               "survival_t1", "volume_t", "weighted.dens", "transplant","unique.transect")) %>% 
  mutate(log_volume_t = log(volume_t)) %>% 
  drop_na()-> LATR_surv_dat

## how much size overlap do we have between transplant experiment and observational census?
hist(log(LATR_surv_dat$volume_t[LATR_surv_dat$transplant==F]))
hist(log(LATR_surv_dat$volume_t[LATR_surv_dat$transplant==T]),add=T,col=alpha("gray",0.5))

plot(log(LATR_surv_dat$volume_t[LATR_surv_dat$transplant==F]),
     LATR_surv_dat$survival_t1[LATR_surv_dat$transplant==F])
points(log(LATR_surv_dat$volume_t[LATR_surv_dat$transplant==T]),
       LATR_surv_dat$survival_t1[LATR_surv_dat$transplant==T]-0.025,pch=2)

LATR_surv <- list()
LATR_surv[[1]] <-  gam(survival_t1 ~ s(log_volume_t) + transplant + s(unique.transect,bs="re"),
                         data=LATR_surv_dat, gamma=1.4, family="binomial")
LATR_surv[[2]] <-  gam(survival_t1 ~ s(log_volume_t) + s(weighted.dens)  + transplant + s(unique.transect,bs="re"),
                       data=LATR_surv_dat, gamma=1.4, family="binomial")
LATR_surv[[3]] <-  gam(survival_t1 ~ s(log_volume_t) + s(weighted.dens) + transplant + weighted.dens:log_volume_t + s(unique.transect,bs="re"),
                       data=LATR_surv_dat, gamma=1.4, family="binomial")
surv_aic<-AICtab(LATR_surv,base=T,sort=F)
LATR_surv_best <- LATR_surv[[which.min(surv_aic$AIC)]]
saveRDS(LATR_surv_best,"LATR_surv_best.rds")
LATR_surv_fitted_terms = predict(LATR_surv_best,type="terms") 
LATR_surv_dat$pred = predict.gam(LATR_surv_best,newdata = LATR_surv_dat,exclude="s(unique.transect)")

##### effect of size on pr(survival) -- linear, positive
plot(LATR_surv_dat$log_volume_t,LATR_surv_fitted_terms[,"s(log_volume_t)"]) 
#### effect of d.stand on pr(survival) -- negativ, a bit of a kink
plot(LATR_surv_dat$weighted.dens,LATR_surv_fitted_terms[,"s(weighted.dens)"]) 

## visualize data + model -- this is for the natural census
n_cuts_dens <- 4
n_cuts_size <- 4
LATR_surv_dat %>% 
  filter(transplant==F) %>% 
  mutate(size_bin = as.integer(cut_interval(log_volume_t,n_cuts_size)),
         dens_bin = as.integer(cut_interval(weighted.dens,n_cuts_dens))) %>% 
  group_by(size_bin,dens_bin) %>% 
  summarise(mean_size = mean(log_volume_t),
            mean_density = mean(weighted.dens),
            mean_surv = mean(survival_t1),
            pred_surv = mean(pred),
            bin_n = n()) -> LATR_surv_nat_plot

## generate predictions for plotting
size_means_surv_nat <- LATR_surv_nat_plot %>% group_by(size_bin) %>% summarise(mean_size=mean(mean_size))
LATR_surv_nat_pred <- data.frame(
  weighted.dens = rep(seq(min(LATR_surv_nat_plot$mean_density),max(LATR_surv_nat_plot$mean_density),length.out = 20),times=n_cuts_size),
  log_volume_t = rep(size_means_surv_nat$mean_size,each=20),
  unique.transect=1,
  transplant=F,
  size_bin = rep(size_means_surv_nat$size_bin,each=20)
)
LATR_surv_nat_pred$pred <- predict.gam(LATR_surv_best,newdata = LATR_surv_nat_pred, exclude = "s(unique.transect)")

plot(LATR_surv_nat_plot$mean_density,LATR_surv_nat_plot$mean_surv,type="n",ylim=c(0,1),
     xlab="Weighted density",ylab="Pr(Survival)")
for(i in 1:n_cuts_size){
  points(LATR_surv_nat_plot$mean_density[LATR_surv_nat_plot$size_bin==i],
         LATR_surv_nat_plot$mean_surv[LATR_surv_nat_plot$size_bin==i],pch=16,col=i,
         cex=(LATR_surv_nat_plot$bin_n[LATR_surv_nat_plot$size_bin==i]/max(LATR_surv_nat_plot$bin_n))*3)
  lines(LATR_surv_nat_pred$weighted.dens[LATR_surv_nat_pred$size_bin==i],
        invlogit(LATR_surv_nat_pred$pred[LATR_surv_nat_pred$size_bin==i]),col=i)
}

## now transplants
LATR_surv_dat %>% 
  filter(transplant==T) %>% 
  mutate(dens_bin = as.integer(cut_interval(weighted.dens,n_cuts_dens)),
         mean_size = mean(log_volume_t)) %>% 
  group_by(dens_bin) %>% 
  summarise(mean_size = unique(mean_size),
            mean_density = mean(weighted.dens),
            mean_surv = mean(survival_t1),
            bin_n = n()) -> LATR_surv_exp_plot
LATR_surv_exp_pred <- data.frame(
  weighted.dens = seq(min(LATR_surv_exp_plot$mean_density),max(LATR_surv_exp_plot$mean_density),length.out = 20),
  log_volume_t = LATR_surv_exp_plot$mean_size[1],
  unique.transect=1,
  transplant=T
)
LATR_surv_exp_pred$pred <- predict.gam(LATR_surv_best,newdata = LATR_surv_exp_pred, exclude = "s(unique.transect)")

points(LATR_surv_exp_plot$mean_density,LATR_surv_exp_plot$mean_surv,ylim=c(0,1))
lines(LATR_surv_exp_pred$weighted.dens,invlogit(LATR_surv_exp_pred$pred))


# Recruitment -------------------------------------------------------------
## estimate per-seed recruitment probability by estimating total seeds per window and total recruits per window

LATR_recruits <- LATR_full %>% 
  mutate(unique.transect = interaction(transect, site)) %>% 
  group_by(year_t1,unique.transect,actual.window) %>% 
  filter(seedling_t1==1) %>% 
  summarise(recruits = n()) %>% 
  rename(window=actual.window)

## now estimate total seeds produced in each window using the known plant sizes and the fitted flowering and fruiting models
LATR_transects <- read.csv("CData.Transects.Windows.csv") %>% 
  mutate(unique.transect = interaction(transect, site),
         log_volume_t = log(volume))
LATR_transects$seeds = ceiling(invlogit(predict.gam(LATR_flower_best,newdata = LATR_transects)) * 
           exp(predict.gam(LATR_fruits_best,newdata = LATR_transects))) 
LATR_transects %>% 
  group_by(unique.transect,window) %>% 
  summarise(total_seeds=sum(seeds),
            weighted.dens = unique(weighted.dens)) -> LATR_transects

## now do something weird. take three copies of this df, assigning each one to a different year and assigning recruits to zero (for now)
LATR_recruitment <- bind_rows(LATR_transects %>% filter(unique.transect=="1.FPS"|unique.transect=="2.FPS"|unique.transect=="3.FPS") %>% 
                                mutate(year_t1=2014,recruits=0), ## only FPS for 2013-2014
          LATR_transects %>% mutate(year_t1=2015,recruits=0),
          LATR_transects %>% mutate(year_t1=2016,recruits=0),
          LATR_transects %>% mutate(year_t1=2017,recruits=0)) %>% 
  left_join(.,LATR_recruits,by=c("year_t1","unique.transect","window")) %>% 
  mutate(recruits.y=replace_na(recruits.y,0),
         recruits = pmax(recruits.x,recruits.y,na.rm=T)) %>% 
  drop_na()

LATR_recruit <- list()
LATR_recruit[[1]] <-  gam(cbind(recruits,total_seeds-recruits) ~ s(unique.transect,bs="re"),
                       data=LATR_recruitment, gamma=1.4, family="binomial")
LATR_recruit[[2]] <-  gam(cbind(recruits,total_seeds-recruits) ~ s(weighted.dens) + s(unique.transect,bs="re"),
                          data=LATR_recruitment, gamma=1.4, family="binomial")
recruit_aic<-AICtab(LATR_recruit,base=T,sort=F)
LATR_recruit_best <- LATR_recruit[[which.min(recruit_aic$AIC)]]
saveRDS(LATR_recruit_best,"LATR_recruit_best.rds")

plot(LATR_recruitment$weighted.dens,LATR_recruitment$recruits/LATR_recruitment$total_seeds)
LATR_recruitment$pred = predict.gam(LATR_recruit_best,newdata = LATR_recruitment,exclude="s(unique.transect)")
points(LATR_recruitment$weighted.dens,invlogit(LATR_recruitment$pred),col="red",pch=".")


# Seedlings size distribution ---------------------------------------------
LATR_recruit_size <- LATR_full %>% 
  filter(seedling_t1==1) %>% 
  mutate(log_volume = log(volume_t1))

hist(LATR_recruit_size$log_volume)
saveRDS(data.frame(recruit_mean = mean(LATR_recruit_size$log_volume),
                   recruit_sd = sd(LATR_recruit_size$log_volume)),"LATR_recruit_size.rds")


# Size bounds -------------------------------------------------------------
saveRDS(data.frame(min_size = log(min(LATR_full$volume_t,LATR_full$volume_t1[LATR_full$transplant==F],na.rm=T)),
                   max_size = log(max(LATR_full$volume_t,LATR_full$volume_t1[LATR_full$transplant==F],na.rm=T))),
        "LATR_size_bounds.rds")


# IPM results -------------------------------------------------------------
source("creosote_IPM_source_fns.R")

## lambda over density variation
density_dummy <- seq(min(LATR_full$weighted.dens,na.rm=T),max(LATR_full$weighted.dens,na.rm=T),length.out = 30)
lambda_density <- c()
for(d in 1:length(density_dummy)){
  print(d)
  lambda_density[d] <- lambda(bigmatrix(dens=density_dummy[d], mat.size=200)$IPMmat)
}

plot(density_dummy,lambda_density,type="l",lwd=3,xlab="Weighted density",ylab="lambda")
abline(h=1,lty=3)

test <- bigmatrix(dens=0,mat.size = 200)
## check the P matrix
plot(test$Pmat)
colSums(test$Pmat)
lambda(test$Pmat)
## checks out

## check the F matrix
plot(test$Fmat)


lambda(test$IPMmat)
