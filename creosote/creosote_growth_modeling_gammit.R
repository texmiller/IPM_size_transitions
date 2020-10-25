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
LATR <- read.csv("CData.csv") %>% 
  #calculate volume
  mutate(#standardize weighted density to mean zero
         d.stand = (weighted.dens - mean(weighted.dens, na.rm = TRUE)) / sd(weighted.dens, na.rm = TRUE),
         #create unique transect as interaction of transect and site
         unique.transect = interaction(transect, site)) %>% 
  drop_na(volume_t,volume_t1)

# first look at size transitions
plot(log(LATR$volume_t),log(LATR$volume_t1))

############################################################################
# Gaussian fits and model selection using mgcv 
############################################################################
## three candidate models for the mean: size only, size + density, or size, density, and size:density
## three candidates for variance: size only, size+density, fitted value (all the covariates plus rfx)
LATR_gam_models=list()
## Pilot fits, where sigma depends on initial size only
LATR_gam_models[[1]] <- gam(list(log(vol_t1) ~s(log(vol_t)) + s(unique.transect,bs="re"), ~s(log(vol_t))), 
                            data=LATR, gamma=1.4, family=gaulss())
LATR_gam_models[[2]] <- gam(list(log(vol_t1) ~s(log(vol_t)) + s(d.stand) + s(unique.transect,bs="re"), ~s(log(vol_t))), 
                            data=LATR, gamma=1.4, family=gaulss())                
LATR_gam_models[[3]] <- gam(list(log(vol_t1) ~s(log(vol_t)) + s(d.stand) + d.stand:log(vol_t) + s(unique.transect,bs="re"), ~s(log(vol_t))), 
                            data=LATR, gamma=1.4, family=gaulss())  
## Fits where sigma depends on both initial size and density
LATR_gam_models[[4]] <- gam(list(log(vol_t1) ~s(log(vol_t)) + s(unique.transect,bs="re"), ~s(log(vol_t)) + s(d.stand)), 
                            data=LATR, gamma=1.4, family=gaulss())
LATR_gam_models[[5]] <- gam(list(log(vol_t1) ~s(log(vol_t)) + s(d.stand) + s(unique.transect,bs="re"), ~s(log(vol_t)) + s(d.stand)), 
                            data=LATR, gamma=1.4, family=gaulss())                
LATR_gam_models[[6]] <- gam(list(log(vol_t1) ~s(log(vol_t)) + s(d.stand) + d.stand:log(vol_t) + s(unique.transect,bs="re"), ~s(log(vol_t)) + s(d.stand)), 
                            data=LATR, gamma=1.4, family=gaulss()) 
## these models will be iterated to fit sigma as f(fitted value)
LATR$fitted_vals = log(LATR$vol_t); 
LATR_gam_models[[7]] <- gam(list(log(vol_t1) ~s(log(vol_t)) + s(unique.transect,bs="re"), ~s(fitted_vals)), 
                            data=LATR, gamma=1.4, family=gaulss())
LATR_gam_models[[8]] <- gam(list(log(vol_t1) ~s(log(vol_t)) + s(d.stand) + s(unique.transect,bs="re"), ~s(fitted_vals)), 
                            data=LATR, gamma=1.4, family=gaulss())                
LATR_gam_models[[9]] <- gam(list(log(vol_t1) ~s(log(vol_t)) + s(d.stand) + d.stand:log(vol_t) + s(unique.transect,bs="re"), ~s(fitted_vals)), 
                            data=LATR, gamma=1.4, family=gaulss())  
for(mod in 7:9) {
  fitGAU = LATR_gam_models[[mod]]
  fitted_all = predict(fitGAU,type="response",data=LATR);                  
  fitted_vals = new_fitted_vals = fitted_all[,1]; 
  weights = fitted_all[,2]; # what I call "weights" here are 1/sigma values; see ?gaulss for details.

  err=100; k=0; 
  while(err>10^(-6)) {
    LATR$fitted_vals = new_fitted_vals; 
    fitGAU <- update(fitGAU); 
    fitted_all = predict(fitGAU,type="response",data=LATR);   
    new_fitted_vals = fitted_all[,1]; new_weights = fitted_all[,2];
    err = weights - new_weights; err=sqrt(mean(err^2)); 
    weights = new_weights; 
    k=k+1; cat(k,err,"\n"); 
  }   
  LATR_gam_models[[mod]] =  fitGAU;
}

AICtab(LATR_gam_models) 
## Model 2 is the winner: mean ~ s(size) + s(density), sd ~ s(size)
## It surprises me that initial size was the better covariate than initial value, because I expected the
## additional effect of density to affect sigma. But maybe size and density are correlated.
plot(LATR$d.stand,log(LATR$vol_t)) ## yes, they are - onward

## define models 2 as our best Gaussian model
LATR_gam_model <- LATR_gam_models[[5]]
LATR$fitted_vals = new_fitted_vals
## extract the linear predictor for the mean -- we'll use this later
LATR_Xp <- predict.gam(LATR_gam_model,type="lpmatrix")
##################################################################  
# Extract values of the fitted splines to explore their properties 
##################################################################
fitted_terms = predict(LATR_gam_model,type="terms")  

##### effect of initial size on mean of final size 
plot(log(LATR$vol_t), fitted_terms[,"s(log(vol_t))"]) 
##### effect of d.stand on mean of final size 
plot(LATR$d.stand, fitted_terms[,"s(d.stand)"]); ## complicated 

##### sigma versus size - presumably log scale
plot(log(LATR$vol_t), fitted_terms[,"s.1(log(vol_t))"]) 
##### sigma versus density 
plot(LATR$d.stand, fitted_terms[,"s.1(d.stand)"]) 

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
px = LATR$fitted_vals; py=scaledResids; 
z = rollMomentsNP(px,py,windows=8,smooth=TRUE,scaled=TRUE) 

################################################################################
## Finding a better distribution. Must at least allow positive excess kurtosis 
#################################################################################

cand_dist=c("NO","GT","LO","TF","ST1","ST2") ##NET, SEP1, and EGB2 are all highly unstable
n_bins <- 8
select_dist <- tibble(init_size = log(LATR$vol_t),
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
## TF and LO pop up a lot, but for the smallest bin it seems like a skewed distribution is best
## Going ahead with ST1, which can hopefully get both the skew and kurtosis
## Noteworthy that the normal is not supported (except last bin)

## visualize ST1 parameters in relation to fitted value by bin
LATR_bin_fit <-LATR %>% 
  mutate(init_size = log(LATR$vol_t),
         bin = as.integer(cut_number(init_size,n_bins)),
         dens_bin = as.integer(cut_number(d.stand,2))) %>% 
  mutate(mu_lowdens=NA,sigma_lowdens=NA,nu_lowdens=NA,tau_lowdens=NA,
         mu_highdens=NA,sigma_highdens=NA,nu_highdens=NA,tau_highdens=NA)
for(b in 1:n_bins){
  bin_fit_lowdens <- gamlssMaxlik(y=log(LATR_bin_fit$vol_t1[LATR_bin_fit$bin==b & LATR_bin_fit$dens_bin==1]),DIST="ST1")
  bin_fit_highdens <- gamlssMaxlik(y=log(LATR_bin_fit$vol_t1[LATR_bin_fit$bin==b & LATR_bin_fit$dens_bin==2]),DIST="ST1")
  LATR_bin_fit$mu_lowdens[LATR_bin_fit$bin==b] <- bin_fit_lowdens$out[[1]]$estimate["eta.mu"] ## identity link
  LATR_bin_fit$sigma_lowdens[LATR_bin_fit$bin==b] <- exp(bin_fit_lowdens$out[[1]]$estimate["eta.sigma"]) ## log link
  LATR_bin_fit$nu_lowdens[LATR_bin_fit$bin==b] <- (bin_fit_lowdens$out[[1]]$estimate["eta.nu"]) ## identity link
  LATR_bin_fit$tau_lowdens[LATR_bin_fit$bin==b] <- exp(bin_fit_lowdens$out[[1]]$estimate["eta.tau"]) ## log link
  LATR_bin_fit$mu_highdens[LATR_bin_fit$bin==b] <- bin_fit_highdens$out[[1]]$estimate["eta.mu"] ## identity link
  LATR_bin_fit$sigma_highdens[LATR_bin_fit$bin==b] <- exp(bin_fit_highdens$out[[1]]$estimate["eta.sigma"]) ## log link
  LATR_bin_fit$nu_highdens[LATR_bin_fit$bin==b] <- (bin_fit_highdens$out[[1]]$estimate["eta.nu"]) ## identity link
  LATR_bin_fit$tau_highdens[LATR_bin_fit$bin==b] <- exp(bin_fit_highdens$out[[1]]$estimate["eta.tau"]) ## log link
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

pdf("../manuscript/figures/creosote_binned_LO.pdf",height = 5,width = 10,useDingbats = F)
par(mfrow=c(2,2),bty="l",mar=c(4,4,2,1),mgp=c(2.2,1,0),cex.axis=1.4,cex.lab=1.4);
plot(LATR_bin_fit$mean_size,LATR_bin_fit$mu_lowdens,xlab="Initial size",ylab=expression(paste("ST parameter  ", mu )),type="b",pch=1)
points(LATR_bin_fit$mean_size,LATR_bin_fit$mu_highdens,type="b",pch=16)

plot(LATR_bin_fit$mean_size,LATR_bin_fit$sigma_highdens,xlab="Initial size",
     ylab=expression(paste("ST  parameter  ", sigma)),type="b",pch=16)
points(LATR_bin_fit$mean_size,LATR_bin_fit$sigma_lowdens,type="b",pch=1)

plot(LATR_bin_fit$mean_size,LATR_bin_fit$nu_highdens,xlab="Initial size",
     ylab=expression(paste("ST  parameter  ", nu)),type="b",pch=16)
points(LATR_bin_fit$mean_size,LATR_bin_fit$nu_lowdens,type="b",pch=1)

plot(LATR_bin_fit$mean_size,LATR_bin_fit$tau_highdens,xlab="Initial size",
     ylab=expression(paste("ST  parameter  ", tau)),type="b",pch=16)
points(LATR_bin_fit$mean_size,LATR_bin_fit$tau_lowdens,xlab="Initial size",
     ylab=expression(paste("ST  parameter  ", tau)),type="b",pch=1)

dev.off()


# Fitting the final model -------------------------------------------------
## we will use the design matrix from gam() to fit the linear predictor
## for the mean exactly as it was fit for the gaussian model (except now it's location)
## gam() also fit smooths for sigma but I think it is not the same sigma as ST1's...?
## sigma and nu also look like they should have size and density dependence
## tau is causing trouble in the binned fits, so I will fit a constant value

## I am going to separate the design matrix into matrices for the mean and sigma
## sigma terms start with "(Intercept).1"
U <- LATR_Xp[,1:(which(names(coef(LATR_gam_model))=="(Intercept).1")-1)]
U_sigma <- LATR_Xp[,(which(names(coef(LATR_gam_model))=="(Intercept).1"):length(coef(LATR_gam_model)))]

## ST1 likelihood -- with design matrix from gam()
LogLik=function(pars,response,U){
  pars1 = pars[1:ncol(U)]; pars2=pars[-(1:ncol(U))];
  mu = U%*%pars1;  
  val = dST1(x = response, 
            mu=mu,
            sigma = exp(pars2[1] + pars2[2]*log(LATR$vol_t) + pars2[3]*LATR$d.stand),
            nu = pars2[4] + pars2[5]*log(LATR$vol_t) + pars2[6]*LATR$d.stand,
            tau = exp(pars2[7]),
            log=T) 
  return(val); 
}

paranoid_iter <- 3
coefs = list(paranoid_iter); LL=numeric(paranoid_iter);  

# Starting values from the gam fit  
fixed_start = coef(LATR_gam_model)[1:ncol(LATR_Xp_mu)]

## make sure the dimensions line up
length(fixed_start);ncol(U);colnames(U) 
## starting values for sigma, nu, tau
fit_sigma = lm(log(sigma_highdens)~mean_size, data=LATR_bin_fit)
fit_nu = lm(nu_highdens~mean_size, data=LATR_bin_fit)
fit_tau = lm(log(tau_highdens)~1, data=LATR_bin_fit)

## bundle coefficients for mu and sigma
p0=c(fixed_start,c(coef(fit_sigma),0),c(coef(fit_nu),0),coef(fit_tau)) 

for(j in 1:paranoid_iter) {
  out=maxLik(logLik=LogLik,start=p0*exp(0.2*rnorm(length(p0))), response=log(LATR$vol_t1),U=U,
             method="BHHH",control=list(iterlim=5000,printLevel=2),finalHessian=FALSE); 
  
  out=maxLik(logLik=LogLik,start=out$estimate,response=log(LATR$vol_t1),U=U,
             method="NM",control=list(iterlim=5000,printLevel=1),finalHessian=FALSE); 
  
  out=maxLik(logLik=LogLik,start=out$estimate,response=log(LATR$vol_t1),U=U,
             method="BHHH",control=list(iterlim=5000,printLevel=2),finalHessian=FALSE); 
  
  coefs[[j]] = out$estimate; LL[j] = out$maximum;
  cat(j, "#--------------------------------------#",out$maximum,"\n"); 
}

j = min(which(LL==max(LL))) ## they actually all land on the same likelihood-that's good!
out=maxLik(logLik=LogLik,start=coefs[[j]],response=log(LATR$vol_t1),U=U,
           method="BHHH",control=list(iterlim=5000,printLevel=2),finalHessian=TRUE) 

# AIC improvement over gaussian
(AIC_ST <- 2*length(coefs) - 2*out$maximum)
(AIC_norm <- AIC(LATR_gam_model) ) 

## How do the param estimates compare to original gam()?
## keep in mind I have constant variance, skew, kurtosis for now
plot(fixed_start,out$estimate[1:ncol(LATR_Xp_mu)]);abline(0,1)

######### save results of ML fit.  
names(out$estimate)<-c(colnames(U),"sigma_b0","sigma_b1","sigma_b2","nu_b0","nu_b1","nu_b2","tau_b0")
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

# Simulate data from fitted LO model
MLmu = U%*%coefs[1:ncol(U)] 
n_sim <- 500
LATR_sim_LO<-LATR_sim_NO<-LATR_sim_ST1<-matrix(NA,nrow=nrow(LATR),ncol=n_sim)
for(i in 1:n_sim){
  #LATR_sim_LO[,i] <- rLO(n = nrow(LATR), 
  #                       mu = MLmu, 
  #                       sigma = exp(coefs["sigma_b0"] + coefs["sigma_b1"]*MLmu  + coefs["sigma_b2"]*MLmu^2))
  ## is this a fair comparison? The linear predictor for the mean is different, because the LO approximates the gam
  LATR_sim_ST1[,i] <- rST1(n = nrow(LATR), 
                           mu = MLmu, 
                           sigma = exp(coefs["sigma_b0"] + coefs["sigma_b1"]*MLmu  + coefs["sigma_b2"]*MLmu^2),
                           nu=coefs["nu_b0"] + coefs["nu_b1"]*MLmu  + coefs["nu_b2"]*MLmu^2,
                           tau=exp(coefs["tau_b0"]))
  LATR_sim_NO[,i] <- rnorm(n = nrow(LATR),
                           mean = fitted_all[,1],
                           sd = 1/fitted_all[,2])
}


n_bins = 10
alpha_scale = 0.7
LATR_moments <- LATR %>% 
  arrange(log(vol_t)) %>% 
  mutate(size_bin = cut_number(log(vol_t),n=n_bins)) %>% 
  group_by(size_bin) %>% 
  summarise(mean_t1 = mean(log(vol_t1)),
            sd_t1 = sd(log(vol_t1)),
            skew_t1 = NPskewness(log(vol_t1)),
            kurt_t1 = NPkurtosis(log(vol_t1)),
            bin_mean = mean(log(vol_t)),
            bin_n = n()) 

pdf("../manuscript/figures/creosote_sim_moments_ST1.pdf",height = 10,width = 10,useDingbats = F)
par(mfrow=c(2,2),mar=c(4,4,2,1),cex.axis=1.3,cex.lab=1.3,mgp=c(2,1,0),bty="l"); 
sim_bin_means=sim_moment_means=sim_moment_means_norm = matrix(NA,n_bins,n_sim); 
for(i in 1:n_sim){
  sim_moments <- bind_cols(LATR,data.frame(sim=LATR_sim_ST1[,i],
                                           sim_norm=LATR_sim_NO[,i])) %>% 
    arrange(log(vol_t)) %>% 
    mutate(size_bin = cut_number(log(vol_t),n=n_bins)) %>% 
    group_by(size_bin) %>% 
    summarise(mean_t1 = mean(sim),
              mean_t1_norm = mean(sim_norm),
              bin_mean = mean(log(vol_t)))
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
  sim_moments <- bind_cols(LATR,data.frame(sim=LATR_sim_ST1[,i],
                                           sim_norm=LATR_sim_NO[,i])) %>% 
    arrange(log(vol_t)) %>% 
    mutate(size_bin = cut_number(log(vol_t),n=n_bins)) %>% 
    group_by(size_bin) %>% 
    summarise(mean_t1 = sd(sim),
              mean_t1_norm = sd(sim_norm),
              bin_mean = mean(log(vol_t)))
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
  sim_moments <- bind_cols(LATR,data.frame(sim=LATR_sim_ST1[,i],
                                           sim_norm=LATR_sim_NO[,i])) %>% 
    arrange(log(vol_t)) %>% 
    mutate(size_bin = cut_number(log(vol_t),n=n_bins)) %>% 
    group_by(size_bin) %>% 
    summarise(mean_t1 = NPskewness(sim),
              mean_t1_norm = NPskewness(sim_norm),
              bin_mean = mean(log(vol_t)))
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
  sim_moments <- bind_cols(LATR,data.frame(sim=LATR_sim_ST1[,i],
                                           sim_norm=LATR_sim_NO[,i])) %>% 
    arrange(log(vol_t)) %>% 
    mutate(size_bin = cut_number(log(vol_t),n=n_bins)) %>% 
    group_by(size_bin) %>% 
    summarise(mean_t1 = NPkurtosis(sim),
              mean_t1_norm = NPkurtosis(sim_norm),
              bin_mean = mean(log(vol_t)))
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
LATR_flow_dat <- LATR %>% select(vol_t,total.reproduction_t,d.stand,unique.transect) %>% drop_na() %>% 
  mutate(log_vol_t = log(vol_t))
LATR_flower <- list()
LATR_flower[[1]] <-  gam(total.reproduction_t>0 ~ s(log_vol_t) + s(unique.transect,bs="re"),
                   data=LATR_flow_dat, gamma=1.4, family="binomial")
LATR_flower[[2]] <-  gam(total.reproduction_t>0 ~ s(log_vol_t) + s(d.stand) + s(unique.transect,bs="re"),
                         data=LATR_flow_dat, gamma=1.4, family="binomial")
LATR_flower[[3]] <-  gam(total.reproduction_t>0 ~ s(log_vol_t) + s(d.stand) + d.stand:log(vol_t)  + s(unique.transect,bs="re"),
                         data=LATR_flow_dat, gamma=1.4, family="binomial")
AICtab(LATR_flower)
LATR_flower_best <- LATR_flower[[2]]
LATR_flower_fitted_terms = predict(LATR_flower_best,type="terms") 
LATR_flow_dat$pred = predict.gam(LATR_flower_best,newdata = LATR_flow_dat, exclude = "s(unique.transect)")

##### effect of size on pr(flower) -- simple linear
plot(LATR_flow_dat$log_vol_t,LATR_flower_fitted_terms[,"s(log_vol_t)"]) 
gam_size_smooth <- lm(LATR_flower_fitted_terms[,"s(log_vol_t)"]~LATR_flow_dat$log_vol_t)
abline(coef(gam_size_smooth)[1],coef(gam_size_smooth)[2])
#### effect of d.stand on pr(flower) -- linear POSITIVE effect of denity
plot(LATR_flow_dat$d.stand,LATR_flower_fitted_terms[,"s(d.stand)"]) 
gam_dens_smooth <- lm(LATR_flower_fitted_terms[,"s(d.stand)"]~LATR_flow_dat$d.stand)
abline(coef(gam_dens_smooth)[1],coef(gam_dens_smooth)[2])

## here is a frankenstein model that I will derive from the gam for prediction
## easy to do because the smooths are linear
gam_predict_flower <- function(size,dens){
  invlogit(coef(LATR_flower_best)[1] + 
             coef(gam_size_smooth)[1] + coef(gam_size_smooth)[2]*size + 
             coef(gam_dens_smooth)[1] + coef(gam_dens_smooth)[2]*dens)
}
## see how it compares to gam predict --  this works
plot(LATR_flow_dat$pred,logit(gam_predict_flower(LATR_flow_dat$log_vol_t,LATR_flow_dat$d.stand)))

n_cuts_dens <- 8
n_cuts_size <- 4
LATR_flow_dat %>% 
  mutate(size_bin = as.integer(cut_number(log_vol_t,n_cuts_size)),
         dens_bin = as.integer(cut_number(d.stand,n_cuts_dens))) %>% 
  group_by(size_bin,dens_bin) %>% 
  mutate(mean_size = mean(log_vol_t),
         mean_density = mean(d.stand),
         mean_flower = mean(total.reproduction_t > 0),
         pred_flower = mean(pred),
         bin_n = n()) -> LATR_flow_dat_plot

d_dummy <- seq(min(LATR_flow_dat_plot$d.stand),max(LATR_flow_dat_plot$d.stand),0.1)
plot(LATR_flow_dat_plot$mean_density,LATR_flow_dat_plot$total.reproduction_t>0,type="n",
     xlab="Weighted density",ylab="Pr(Flowering)")
for(i in 1:n_cuts_size){
  points(LATR_flow_dat_plot$mean_density[LATR_flow_dat_plot$size_bin==i],
         LATR_flow_dat_plot$mean_flower[LATR_flow_dat_plot$size_bin==i],pch=16,col=i,
         cex=(LATR_flow_dat_plot$bin_n[LATR_flow_dat_plot$size_bin==i]/max(LATR_flow_dat_plot$bin_n))*3)
  lines(d_dummy,
      gam_predict_flower(mean(LATR_flow_dat_plot$log_vol_t[LATR_flow_dat_plot$size_bin==i]),d_dummy),col=i)
}

# Fruit production --------------------------------------------------------
LATR_fruits_dat <- subset(LATR_flow_dat,total.reproduction_t>0)
LATR_fruits <- list()
LATR_fruits[[1]] <-  gam(total.reproduction_t ~ s(log_vol_t) + s(unique.transect,bs="re"),
                         data=LATR_fruits_dat, gamma=1.4, family="nb")
LATR_fruits[[2]] <-  gam(total.reproduction_t ~ s(log_vol_t) + s(d.stand) + s(unique.transect,bs="re"),
                         data=LATR_fruits_dat, gamma=1.4, family="nb")
LATR_fruits[[3]] <-  gam(total.reproduction_t ~ s(log_vol_t) + s(d.stand) + d.stand:log(vol_t)  + s(unique.transect,bs="re"),
                         data=LATR_fruits_dat, gamma=1.4, family="nb")
AICtab(LATR_fruits)
LATR_fruits_best <- LATR_fruits[[2]]
LATR_fruits_fitted_terms = predict(LATR_fruits_best,type="terms") 
LATR_fruits_dat$pred = predict.gam(LATR_fruits_best,newdata = LATR_fruits_dat,exclude="s(unique.transect)")

##### effect of size on pr(flower) -- a bit of a kink
plot(LATR_fruits_dat$log_vol_t,LATR_fruits_fitted_terms[,"s(log_vol_t)"]) 
#### effect of d.stand on pr(flower) -- funky
plot(LATR_fruits_dat$d.stand,LATR_fruits_fitted_terms[,"s(d.stand)"]) 

LATR_fruits_dat %>% 
  mutate(size_bin = as.integer(cut_number(log_vol_t,n_cuts_size)),
         dens_bin = as.integer(cut_number(d.stand,n_cuts_dens))) %>% 
  group_by(size_bin,dens_bin) %>% 
  mutate(mean_size = mean(log_vol_t),
         mean_density = mean(d.stand),
         mean_fruits = mean(total.reproduction_t),
         pred_fruits = mean(pred),
         bin_n = n()) -> LATR_fruits_dat_plot

## new data set for gam prediction
fruits_newd <- data.frame(
  size_bin = rep(1:n_cuts_size,each=length(d_dummy)),
  log_vol_t = rep(aggregate(log_vol_t ~ size_bin, data=LATR_fruits_dat_plot, mean)$log_vol_t,
                  each=length(d_dummy)),
  d.stand = rep(d_dummy,times=n_cuts_size),
  unique.transect = "transect"
)
fruits_newd$pred <- predict.gam(LATR_fruits_best,newdata = fruits_newd, exclude = "s(unique.transect)")

plot(LATR_fruits_dat_plot$mean_density,LATR_fruits_dat_plot$mean_fruits,type="n",
     xlab="Weighted density",ylab="Flowers and Fruits")
for(i in 1:n_cuts_size){
  points(LATR_fruits_dat_plot$mean_density[LATR_fruits_dat_plot$size_bin==i],
         LATR_fruits_dat_plot$mean_fruits[LATR_fruits_dat_plot$size_bin==i],pch=16,col=i,
         cex=(LATR_fruits_dat_plot$bin_n[LATR_fruits_dat_plot$size_bin==i]/max(LATR_fruits_dat_plot$bin_n))*3)
  lines(fruits_newd$d.stand[fruits_newd$size_bin==i],
        exp(fruits_newd$pred[fruits_newd$size_bin==i]),col=i)
}


# Survival ----------------------------------------------------------------
# read in transplant experiment and merge with obs data
# Combine transplants with large shrubs for later survival analysis
# Keep only location info, survival, volume, and density
CData.Transplants<-read.csv("CData.Transplants.csv") %>% 
  select("site", "transect", "actual.window", 
       "spring_survival_t1", "volume_t", "weighted.dens", "transplant") %>% 
  rename("survival_t1" = "spring_survival_t1") %>% 
  rbind(select(LATR, "site", "transect", "actual.window", 
               "survival_t1", "volume_t", "weighted.dens", "transplant")) -> LATR_surv_dat

## how much size overlap do we have between transplant experiment and observational census?
hist(log(LATR_surv_dat$volume_t[LATR_surv_dat$transplant==F]))
hist(log(LATR_surv_dat$volume_t[LATR_surv_dat$transplant==T]),add=T,col=alpha("gray",0.5))

plot(log(LATR_surv_dat$volume_t[LATR_surv_dat$transplant==F]),
     LATR_surv_dat$survival_t1[LATR_surv_dat$transplant==F])
