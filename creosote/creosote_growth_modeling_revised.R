rm(list=ls(all=TRUE));

# setwd("c:/repos/IPM_size_transitions/creosote"); #Steve
setwd("C:/Users/tm9/Desktop/git local/IPM_size_transitions/creosote"); #Tom

require(car); require(lme4); require(zoo); require(moments); require(mgcv); 
require(gamlss); require(gamlss.tr); require(AICcmodavg); 
require(lmerTest); require(tidyverse); require(maxLik); 
require(actuar); require(lattice); require(grid); require(scales);
require(sgt); require(formatR); require(popbio); require(bbmle)

# misc functions
volume <- function(h, w, p){
  (1/3)*pi*h*(((w + p)/2)/2)^2
}

invlogit <- function(x){exp(x)/(1+exp(x))}
 
# log kurtosis function for diagnostics
Lkurtosis=function(x) log(kurtosis(x)); 

# Steve's diagnostics functions
source("../Diagnostics.R")
# function for choosing distribution family
source("../fitChosenDists.R")

## read in data for Larrea tridentata (LATR)
LATR <- read.csv("creosote_growth_density.csv") %>% 
  #calculate volume
  mutate(vol_t = volume(max.ht_t,max.w_t,perp.w_t),
         vol_t1 = volume(max.ht_t1,max.w_t1,perp.w_t1),
         #standardize weighted density to mean zero
         d.stand = (weighted.dens - mean(weighted.dens, na.rm = TRUE)) / sd(weighted.dens, na.rm = TRUE),
         #create unique transect as interaction of transect and site
         unique.transect = interaction(transect, site)) %>% 
  drop_na(vol_t,vol_t1)

# closer look at the outliers
LATR %>% mutate(change=(vol_t1)-(vol_t)) %>% 
  filter(change > quantile(change,0.99,na.rm=T) | change < quantile(change,0.01,na.rm=T)) %>% 
  select(X,site,transect,designated.window,plant,year_t,max.ht_t,max.w_t,perp.w_t,max.ht_t1,max.w_t1,perp.w_t1)
## FPS-3-500 and MOD-3-200 both look like errors. These are the suspect row numbers
outliers <- c(617,684,686,688,882)
LATR %>% filter(!(X %in% outliers)) -> LATR

# first look at size transitions
plot(log(LATR$vol_t),log(LATR$vol_t1))

############################################################################
# Gaussian fits and model selection using mgcv 
############################################################################
LATR_gam_models=list()
## Pilot fits, where sigma depends on initial size only
LATR$fitted_vals = log(LATR$vol_t); 
LATR_gam_models[[1]] <- gam(list(log(vol_t1) ~s(log(vol_t)) + s(unique.transect,bs="re"), ~s(fitted_vals)), 
                            data=LATR, gamma=1.4, family=gaulss())
LATR_gam_models[[2]] <- gam(list(log(vol_t1) ~s(log(vol_t)) + s(d.stand) + s(unique.transect,bs="re"), ~s(fitted_vals)), 
                            data=LATR, gamma=1.4, family=gaulss())                
LATR_gam_models[[3]] <- gam(list(log(vol_t1) ~s(log(vol_t)) + s(d.stand) + d.stand:log(vol_t) + s(unique.transect,bs="re"), ~s(fitted_vals)), 
                            data=LATR, gamma=1.4, family=gaulss())   
## these models will be iterated to fit sigma as f(fitted value)
LATR_gam_models[[4]] <- gam(list(log(vol_t1) ~s(log(vol_t)) + s(unique.transect,bs="re"), ~s(fitted_vals)), 
                            data=LATR, gamma=1.4, family=gaulss())
LATR_gam_models[[5]] <- gam(list(log(vol_t1) ~s(log(vol_t)) + s(d.stand) + s(unique.transect,bs="re"), ~s(fitted_vals)), 
                            data=LATR, gamma=1.4, family=gaulss())                
LATR_gam_models[[6]] <- gam(list(log(vol_t1) ~s(log(vol_t)) + s(d.stand) + d.stand:log(vol_t) + s(unique.transect,bs="re"), ~s(fitted_vals)), 
                            data=LATR, gamma=1.4, family=gaulss())  
for(mod in 4:6) {
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
## Models 2 and 3 are top and basically equivalent. I''ll proceed with 2, the simpler one
## It surprises me that initial size was the better covariated than initial value, because I expected the
## additional effect of density to affect sigma. But maybe size and density are correlated.
plot(LATR$d.stand,log(LATR$vol_t)) ## yes, they are - onward
LATR_gam_model <- LATR_gam_models[[2]]; 
LATR$fitted_vals = new_fitted_vals; 


##################################################################  
# Extract values of the fitted splines to explore their properties 
##################################################################
fitted_terms = predict(LATR_gam_model,type="terms");   

##### effect of initial size on mean of final size 
plot(log(LATR$vol_t), fitted_terms[,1]); 
mean_fit1 = lm(fitted_terms[,1]~log(LATR$vol_t)); ## linear initial size is a perfect fit
abline(coef(mean_fit1)[1],coef(mean_fit1)[2],col="red")

##### effect of d.stand on mean of final size 
plot(LATR$d.stand, fitted_terms[,2]); ## complicated 
d.stand_fit1 = lm(fitted_terms[,2]~LATR$d.stand); 
d.stand_fit2 = lm(fitted_terms[,2]~LATR$d.stand + I(LATR$d.stand^2)); # R^2 = 0.88 
d.stand_fit3 = lm(fitted_terms[,2]~LATR$d.stand + I(LATR$d.stand^2) + I(LATR$d.stand^3)); # R^2 = 0.99
points(LATR$d.stand,fitted(d.stand_fit3), col="red"); 

##### sigma versus fitted values 
fitted_all = predict(LATR_gam_model,type="response");   
sigma.hat = 1/fitted_all[,2]; 
plot(LATR$fitted_vals, log(sigma.hat)) ; 
sigma_fit1 = lm(log(sigma.hat)~LATR$fitted_vals); ## linear in fitted value is a perfect fit 
abline(coef(sigma_fit1)[1],coef(sigma_fit1)[2],col="red")

##################################################################  
# Inspect scaled residuals to evaluate the pilot model: FAILS 
##################################################################
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

cand_dist=c("GT","LO","TF","ST1","ST2") ##NET, SEP1, and EGB2 are all highly unstable
n_bins <- 10
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

select_dist %>% 
  group_by(bin) %>% 
  summarise(n_bin = n(),
            best_dist = unique(best_dist),
            secondbest_dist = unique(secondbest_dist),
            aic_margin = unique(aic_margin))
## TF and LO pop up a lot, but for the smallest bin it seems like a skewed distribution is best
## But the skewed t's are giving me problems so I am going ahead with the logistic and we'll see how it does

## visualize TF parameters in relation to fitted value by bin
LATR_bin_fit <-LATR %>% 
  mutate(init_size = log(LATR$vol_t),
         bin = as.integer(cut_number(init_size,n_bins))) %>% 
  mutate(mu=NA, sigma=NA,nu=NA)
for(b in 1:n_bins){
  bin_fit <- gamlssMaxlik(y=log(LATR_bin_fit$vol_t1[LATR_bin_fit$bin==b]),DIST="ST1")
  LATR_bin_fit$mu[LATR_bin_fit$bin==b] <- bin_fit$out[[1]]$estimate["eta.mu"] ## identity link
  LATR_bin_fit$sigma[LATR_bin_fit$bin==b] <- exp(bin_fit$out[[1]]$estimate["eta.sigma"]) ## log link
  LATR_bin_fit$nu[LATR_bin_fit$bin==b] <- (bin_fit$out[[1]]$estimate["eta.nu"]) ## identity link
  LATR_bin_fit$tau[LATR_bin_fit$bin==b] <- exp(bin_fit$out[[1]]$estimate["eta.tau"]) ## log link
}
LATR_bin_fit %>% 
  group_by(bin) %>% 
  summarise(N = n(),
            mean_size = mean(init_size),
            mu = unique(mu),
            sigma=unique(sigma),
            nu=unique(nu),
            tau=unique(tau)) -> LATR_bin_fit

pdf("../manuscript/figures/creosote_binned_LO.pdf",height = 5,width = 10,useDingbats = F)
par(mfrow=c(2,2),bty="l",mar=c(4,4,2,1),mgp=c(2.2,1,0),cex.axis=1.4,cex.lab=1.4);
plot(LATR_bin_fit$mean_size,LATR_bin_fit$mu,xlab="Initial size",ylab=expression(paste("ST parameter  ", mu )),type="b",pch=16)
plot(LATR_bin_fit$mu,LATR_bin_fit$sigma,xlab=expression(paste("Fitted value  ", mu )),
     ylab=expression(paste("ST  parameter  ", sigma)),type="b",pch=16)
plot(LATR_bin_fit$mu,LATR_bin_fit$nu,xlab=expression(paste("Fitted value  ", mu )),
     ylab=expression(paste("ST  parameter  ", nu)),type="b",pch=16)
plot(LATR_bin_fit$mu,LATR_bin_fit$tau,xlab=expression(paste("Fitted value  ", mu )),
     ylab=expression(paste("ST  parameter  ", tau)),type="b",pch=16)
dev.off()


# Fitting the final model -------------------------------------------------
## Now we will define the U matrix based on what we found in the pilot gaussian fit.
## That fit was based on gam() but we have approximations to the gam terms that we will use here. 
## This includes a cubic term for density

U=model.matrix(~  0 + unique.transect + log(vol_t) + d.stand + I(d.stand^2) + I(d.stand^3), data=LATR)

# Likelihood function with sigma as a quadratic function of mu

## This was the logistic likelihood function, but I am trying to improve upon it
#LogLik=function(pars,response,U){
#  pars1 = pars[1:ncol(U)]; pars2=pars[-(1:ncol(U))];
#  mu = U%*%pars1;  
#  val = dLO(x = response, 
#            mu=mu,
#            sigma = exp(pars2[1] + pars2[2]*mu + pars2[3]*mu^2),
#            log=T) 
#  return(val); 
#}

## ST1 likelihood
LogLik=function(pars,response,U){
  pars1 = pars[1:ncol(U)]; pars2=pars[-(1:ncol(U))];
  mu = U%*%pars1;  
  val = dST1(x = response, 
            mu=mu,
            sigma = exp(pars2[1] + pars2[2]*mu + pars2[3]*mu^2),
            nu = pars2[4] + pars2[5]*mu + pars2[6]*mu^2,
            tau = exp(pars2[7]),
            log=T) 
  return(val); 
}

paranoid_iter <- 3
coefs = list(paranoid_iter); LL=numeric(paranoid_iter);  

# Starting values from the pilot model are jittered to do multi-start optimization). 
# Using good starting values really speeds up convergence in the ML fits  
fixed_start = c(LATR_gam_model$coefficients["(Intercept)"] + LATR_gam_model$coefficients[paste0("s(unique.transect).",1:12)],##transect estimates
                coef(mean_fit1)[2],#size slope
                coef(d.stand_fit3)[2:4])
## make sure the dimensions line up
length(fixed_start);ncol(U);colnames(U) 
## starting coefficients for sigma(mu)
fit_sigma = lm(log(sigma)~mu + I(mu^2), data=LATR_bin_fit)
fit_nu = lm(nu~mu + I(mu^2), data=LATR_bin_fit)
fit_tau = lm(log(tau)~1, data=LATR_bin_fit)

## bundle coefficients for mu and sigma
p0=c(fixed_start,coef(fit_sigma),coef(fit_nu),coef(fit_tau)) 

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
(AIC_norm <- AIC(LATR_gam_model) ) ## the logistic is quite an improvement

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

pdf("../manuscript/figures/creosote_sim_moments_ST1.pdf",height = 8,width = 8,useDingbats = F)
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
