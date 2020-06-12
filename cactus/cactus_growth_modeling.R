##############################################################################
# Growth kernel modeling for Cylindropuntia imbricata - modeled off of Steve's PSSPgrowthModeling-ML.R script. 
#
# Data collection described in Miller et al. 2009, Ohm and Miller 2014, Czachura and Miller 2020, and elsewhere  
# This data set includes 2018 data, which has not been previously published. 
#
# Last update: 8 June, 2020
#
##############################################################################

rm(list=ls(all=TRUE));
setwd("./cactus"); 

require(car); require(lme4); require(zoo); require(moments); require(mgcv); 
require(gamlss); require(gamlss.tr); require(AICcmodavg); 
require(lmerTest); require(tidyverse); require(maxLik); 
require(actuar); require(lattice); require(grid); require(scales);
require(sgt)

# Steve's diagnostics functions
source("../Diagnostics.R") 
# function for converting cactus size measurements to volume
volume <- function(h, w, p){
  (1/3)*pi*h*(((w + p)/2)/2)^2
}


##############################################################
# 1. Read in the data and pare down data set to just size transitions
##############################################################

CYIM<-read_csv("cholla_demography_20042018_EDI.csv") %>% 
  ## filter out transplants and drop seed addition plots (which start with letter H)
  ## also drop Year_t==2018 because I don't have 2019 size data (not entered yet). 2018 data still included in 2017-2018 transition
    filter(Transplant == 0,
         str_sub(Plot,1,1)!="H",
         Year_t!=2018) %>% 
  ## convert height, max width, perp width to volume of cone, take natural log
  mutate(vol_t = volume(Height_t,Width_t,Perp_t),
         vol_t1 = volume(Height_t1,Width_t1,Perp_t1),
         plot = as.factor(Plot),
         year_t = as.factor(Year_t)) %>%
  select(year_t,plot,vol_t,vol_t1) %>% 
  filter() %>% 
  ## sort by initial size
  arrange(vol_t) %>% 
  ## drop rows with NAs
  drop_na()

## how much coverage do I have across years and plots?
table(CYIM$plot,CYIM$year_t) ## first four years has fewer plants from fewer plots

########################################################################## 
## 2. Gaussian fits: log transformation, constant variance, test quadratic term for initial size
## There are no other fixed effects screaming out to be included in the models, but we could test for climate drivers
## and ant defense effects, as we have done elsewhere
########################################################################## 
CYIM_lmer_models <- list() ## drop plot rfx for now bc they create a headache in the shrinkage method later
## random effects are intercept-only - convergence troubles otherwise
CYIM_lmer_models[[1]] <- lmer(log(vol_t1) ~ log(vol_t) + (1|year_t) + (1|plot), data=CYIM,REML=F,control=lmerControl(optimizer="bobyqa"))
CYIM_lmer_models[[2]] <- lmer(log(vol_t1) ~ log(vol_t) + I(log(vol_t)^2) + (1|year_t) + (1|plot), data=CYIM,REML=F,control=lmerControl(optimizer="bobyqa"))

########################################################################## 
## 3. Use iterative re-weighting to fit with nonconstant variance,  
## then AIC for model selection.  
##########################################################################

## NegLogLik function to fit variance model for residuals 
varPars = function(pars) {
  return(-sum(dnorm(resids, mean=0, sd=exp(pars[1] + pars[2]*fitted_vals),log=TRUE)))
}	

##update to save sigma pars
pars<-list()
for(mod in 1:length(CYIM_lmer_models)) {
  err = 1; rep=0; 
  while(err > 0.000001) {
    rep=rep+1; model = CYIM_lmer_models[[mod]];
    fitted_vals = fitted(model);resids = residuals(model); 
    out=optim(c(sd(resids),0),varPars,control=list(maxit=5000)); 
    pars[[mod]]=out$par; 
    new_sigma = exp(pars[[mod]][1] + pars[[mod]][2]*fitted_vals); new_weights = 1/((new_sigma)^2)
    new_weights = 0.5*(weights(model) + new_weights); # cautious update 
    new_model <- update(model,weights=new_weights); 
    err = weights(model)-weights(new_model); err=sqrt(mean(err^2)); 
    cat(mod,rep,err,"\n") # check on convergence of estimated weights 
    CYIM_lmer_models[[mod]]<-new_model; 
  }}
aictab(CYIM_lmer_models) ## no strong support for quadratic term

######### For a fair AIC comparison, fit all models with the same weights 
aics = unlist(lapply(CYIM_lmer_models,AIC)); best_model=which(aics==min(aics)); 
best_weights=weights(CYIM_lmer_models[[best_model]]); 
for(mod in 1:length(CYIM_lmer_models)) {CYIM_lmer_models[[mod]] <- update(CYIM_lmer_models[[mod]],weights=best_weights)}
aictab(CYIM_lmer_models); # stronger support for model 1

######### Here's the best Gaussian model ########################################
CYIM_lmer_best = CYIM_lmer_models[[1]]; 
summary(CYIM_lmer_best); 

### refit with REML, as recommended for estimating random effects 
best_weights = weights(CYIM_lmer_best)
CYIM_lmer_best <- lmer(log(vol_t1) ~ log(vol_t) + (1|year_t) + (1|plot), data=CYIM,weights=best_weights,REML=TRUE) 
## this is as good as we can do with a Gaussian model - here's what it looks like

## Visualize the best Gaussian model
##dummy variable for initial size
size_dim <- 100
size_dum <- seq(min(log(CYIM$vol_t)),max(log(CYIM$vol_t)),length.out = size_dim)
##make a polygon for the Gaussian kernel with non-constant variance
CYIM_lmer_best_kernel <- matrix(NA,size_dim,size_dim)
for(i in 1:size_dim){
  mu_size <- fixef(CYIM_lmer_best)[1] + fixef(CYIM_lmer_best)[2] * size_dum[i]
  CYIM_lmer_best_kernel[,i] <- dnorm(size_dum,
                                     mean = mu_size,
                                     sd = exp(pars[[best_model]][1] + pars[[best_model]][2]*mu_size))
}

graphics.off(); dev.new(width=6,height=5); 
levelplot(CYIM_lmer_best_kernel,row.values = size_dum, column.values = size_dum,cuts=30,
          col.regions=rainbow(30),xlab="log Size t",ylab="log Size t+1",
          panel = function(...) {
            panel.levelplot(...)
            grid.points(log(CYIM$vol_t), log(CYIM$vol_t1), pch = ".",gp = gpar(cex=3,col=alpha("black",0.5)))
          }) 
dev.copy2pdf(file="../manuscript/figures/cactus_Gaussian_kernel_data.pdf") 

######################################################################
# Interrogate the scaled residuals - Gaussian? NO. 
###################################################################### 
plot(fitted(CYIM_lmer_best), 1/sqrt(weights(CYIM_lmer_best)),xlab="Fitted",ylab="Estimated residual Std Dev"); 
lines(fitted(CYIM_lmer_best),exp(pars[[best_model]][1] + pars[[best_model]][2]*fitted(CYIM_lmer_best)),col="red")
scaledResids = residuals(CYIM_lmer_best)*sqrt(weights(CYIM_lmer_best))
plot(fitted(CYIM_lmer_best), scaledResids); 
qqPlot(scaledResids); # really bad in both tails

jarque.test(scaledResids) # normality test: FAILS, P < 0.001 
anscombe.test(scaledResids) # kurtosis: FAILS, P < 0.001 
agostino.test(scaledResids) # skewness: FAILS, P<0.001 

########################################################################
## Rollapply diagnostics on the scaled residuals 
########################################################################
px = fitted(CYIM_lmer_best); py=scaledResids; 

graphics.off(); dev.new(width=8,height=6); 
par(mfrow=c(2,2),bty="l",mar=c(4,4,2,1),mgp=c(2.2,1,0),cex.axis=1.4,cex.lab=1.4);   
z = rollMoments(px,py,windows=10,smooth=TRUE,scaled=TRUE) 
dev.copy2pdf(file="../manuscript/figures/RollingMomentsCYIM.pdf") 
## scaled residuals look about mean zero and unit variance (but maybe some weird trends with fitted values)
## clear negative skew and excess kurtosis

########################################################################
## what is the likely distribution of the residuals?
## because there are trends wrt to fitted value, I will do this in slices. 
########################################################################
## I should re-write this code as a loop because it runs fitDist 3 times more than necessary
n_bins <- 8
select_dist <- tibble(fit_best = fitted(CYIM_lmer_best),
                      scale_resid = residuals(CYIM_lmer_best)*sqrt(weights(CYIM_lmer_best))) %>% 
  mutate(bin = as.integer(cut_number(fit_best,n_bins))) %>% 
  group_by(bin) %>% 
  mutate(best_dist = names(fitDist(scale_resid,type="realline")$fits[1]),
         secondbest_dist = names(fitDist(scale_resid,type="realline")$fits[2]),
         aic_margin = fitDist(scale_resid,type="realline")$fits[2] - fitDist(scale_resid,type="realline")$fits[1]) %>% 
  summarise(n_bin = n(),
            best_dist = unique(best_dist),
            secondbest_dist = unique(secondbest_dist),
            aic_margin = unique(aic_margin))
# a little bit of everything but some version of the skewed t is commonly favored, though the SHASH is strongly favored for the smaller bins

########################################################################
## Fit parameters of the skewed t to binned size data -- to get a visual sense of relationships between mu and other params
########################################################################
CYIM_bin_fit <-CYIM %>% 
  mutate(fitted = fitted(CYIM_lmer_best),
         bin = as.integer(cut_number(fitted,n_bins))) %>% 
  mutate(mu=NA, sigma=NA,nu=NA,tau=NA)
for(b in 1:n_bins){
  ## I get the most stable tau estimates with ST3
  bin_fit <- gamlssML(log(CYIM_bin_fit$vol_t1[CYIM_bin_fit$bin==b]) ~ 1,family="ST3")
  CYIM_bin_fit$mu[CYIM_bin_fit$bin==b] <- bin_fit$mu
  CYIM_bin_fit$sigma[CYIM_bin_fit$bin==b] <- bin_fit$sigma
  CYIM_bin_fit$nu[CYIM_bin_fit$bin==b] <- bin_fit$nu
  CYIM_bin_fit$tau[CYIM_bin_fit$bin==b] <- bin_fit$tau
}
CYIM_bin_fit %>% 
  group_by(bin) %>% 
  summarise(N = n(),
            mean_fitted = mean(fitted),
            mu = unique(mu),
            sigma=unique(sigma),
            nu=unique(nu),
            tau=unique(tau)) -> CYIM_bin_fit

par(mfrow=c(2,2),bty="l",mar=c(4,4,2,1),mgp=c(2.2,1,0),cex.axis=1.4,cex.lab=1.4);
## Steve's spline.scatter.smooth() function not working for me and I did not both trying to figure out why
plot(CYIM_bin_fit$mean_fitted,CYIM_bin_fit$mu,xlab="Fitted value",ylab=expression(paste("Location parameter  ", mu )),type="b")
plot(CYIM_bin_fit$mu,CYIM_bin_fit$sigma,xlab=expression(paste("Location parameter  ", mu )),
                      ylab=expression(paste("Scale parameter  ", sigma)),type="b")
plot(CYIM_bin_fit$mu,CYIM_bin_fit$nu,xlab=expression(paste("Location parameter  ", mu )),
                      ylab=expression(paste("Skewness parameter  ", nu )),type="b")
plot(CYIM_bin_fit$mu,CYIM_bin_fit$tau,xlab=expression(paste("Location parameter  ", mu )),
                      ylab=expression(paste("Kurtosis parameter  ", tau)),type="b") 

## The skewed t fits are a little unstable for some groups and ST3 is the only one that converges.
## and I am getting garbage tau estimate for group 1, wondering if I would do better with the skewed generalized t
## NegLogLik function to fit variance model for residuals 
sgt_fit = function(pars,dat) {
  return(-sum(dsgt(dat, mu=pars[1], sigma=pars[2], lambda=pars[3], p=pars[4], q=pars[5], log=TRUE, mean.cent = F)))
}	

CYIM_bin_sgt_fit <-CYIM %>% 
  mutate(fitted = fitted(CYIM_lmer_best),
         resids = residuals(CYIM_lmer_best),
         bin = as.integer(cut_number(fitted,n_bins))) 
for(b in 1:n_bins){
  bin_fit <- optim(c(mean(CYIM_bin_sgt_fit$fitted[CYIM_bin_sgt_fit$bin==b]),
                     sd(CYIM_bin_sgt_fit$resids[CYIM_bin_sgt_fit$bin==b]),0,2,2),
                   sgt_fit,dat=log(CYIM_bin_sgt_fit$vol_t1[CYIM_bin_sgt_fit$bin==b]),control=list(maxit=5000))
  CYIM_bin_sgt_fit$mu[CYIM_bin_sgt_fit$bin==b] <- bin_fit$par[1]
  CYIM_bin_sgt_fit$sigma[CYIM_bin_sgt_fit$bin==b] <- bin_fit$par[2]
  CYIM_bin_sgt_fit$lambda[CYIM_bin_sgt_fit$bin==b] <- bin_fit$par[3]
  CYIM_bin_sgt_fit$p[CYIM_bin_sgt_fit$bin==b] <- bin_fit$par[4]
  CYIM_bin_sgt_fit$q[CYIM_bin_sgt_fit$bin==b] <- bin_fit$par[5]
  CYIM_bin_sgt_fit$converge[CYIM_bin_sgt_fit$bin==b] <- bin_fit$convergence
}
CYIM_bin_sgt_fit %>% 
  group_by(bin) %>% 
  summarise(N = n(),
            mean_fitted = mean(fitted),
            mu = unique(mu),
            sigma=unique(sigma),
            lambda=unique(lambda),
            p=unique(p),
            q=unique(q),
            converge=unique(converge)) -> CYIM_bin_sgt_fit

par(mfrow=c(2,3),bty="l",mar=c(4,4,2,1),mgp=c(2.2,1,0),cex.axis=1.4,cex.lab=1.4);
## Steve's spline.scatter.smooth() function not working for me and I did not both trying to figure out why
plot(CYIM_bin_sgt_fit$mean_fitted,CYIM_bin_sgt_fit$mu,xlab="Fitted value",type="b")
plot(CYIM_bin_sgt_fit$mu,CYIM_bin_sgt_fit$sigma,type="b")
plot(CYIM_bin_sgt_fit$mu,CYIM_bin_sgt_fit$lambda,type="b")
plot(CYIM_bin_sgt_fit$mu,CYIM_bin_sgt_fit$p,type="b") 
plot(CYIM_bin_sgt_fit$mu,CYIM_bin_sgt_fit$q,type="b") 
## weird q estimates for the smallest group (like I saw in ST nu param). As q->Inf, SGT becomes skewed normal, among other things
################################################################################################################
# Proceed to fitting gamlss ST3 to log size data (may also try SGT just for fun). 
## Use the fixed effect structure corresponding to the best lmer fit,
# and use moment diagnostics to guide functions for other parameeters
#
# Fitting is done by maximum likelihood using maxLik (easier than mle) 
# 
# The random effects terms in the lmer fit are fitted here as fixed effects, then adjusted by shrinkage. 
# RFX were (1|year) and (1|plot) so there will be year- and plot-specific fixed effects
#################################################################################################################

# as a reminder, here is best model: lmer(log(vol_t1) ~ log(vol_t) + (1|year_t) + (1|plot), data=CYIM,REML=F,control=lmerControl(optimizer="bobyqa"))
# Model matrix for the fixed and random effects, specified so that each year gets its own coefficient
# dropping plot for now because I am not sure how to suppress intercept for it, and this will create some issues later on
U=model.matrix(~  log(vol_t) + year_t + plot, data=CYIM)

LogLik=function(pars,response,U){
  pars1 = pars[1:ncol(U)]; pars2=pars[-(1:ncol(U))];
  mu = U%*%pars1;  
  val = dST3(x = response, 
             mu=mu,
             ## maybe sigma response to mu is non-monotonic?
             sigma=exp(pars2[1]+pars2[2]*mu+pars2[3]*mu^2),
             ## nu could be accelerating -- I added a log link here bc I got an error in the maxLik fit that nu must be positive
             ## weird, because gamlss uses identity linl
             nu = exp(pars2[4]+pars2[5]*mu + pars2[6]*mu^2),
            tau = exp(pars2[7]+pars2[8]*mu), log=TRUE)
  return(val); 
}

############ Fit the data -- five separate times (Steve's paranoia)
coefs = list(5); LL=numeric(5);  

# Starting values from the pilot model are jittered to do multi-start optimization). 
# Using good starting values really speeds up convergence in the ML fits  

# Linear predictor coefficients extracted from the lmer model 
fixed_start = c(fixef(CYIM_lmer_best)[1] + unlist(ranef(CYIM_lmer_best)$year_t)[1] + unlist(ranef(CYIM_lmer_best)$plot)[1], #intercept adjusted by the year-1/ plot-1 random effects
                fixef(CYIM_lmer_best)[2], # size slope
                unlist(ranef(CYIM_lmer_best)$year_t)[-1], # all the other years
                unlist(ranef(CYIM_lmer_best)$plot)[-1]) # all the other plots
## make sure the dimensions line up
length(fixed_start);ncol(U) #nail

# Shape and scale coefficients from the rollaply diagnostic plots 
fit_sigma = lm(log(sigma)~mu + I(mu^2), data=CYIM_bin_fit)
fit_nu = lm(nu~mu + I(mu^2), data=CYIM_bin_fit)
fit_tau = lm(log(tau)~mu, data=CYIM_bin_fit)
p0=c(fixed_start, coef(fit_sigma), coef(fit_nu),coef(fit_tau))

for(j in 1:5) {
  out=maxLik(logLik=LogLik,start=p0*exp(0.2*rnorm(length(p0))), response=log(CYIM$vol_t1),U=U,
             method="BHHH",control=list(iterlim=5000,printLevel=2),finalHessian=FALSE); 
  
  out=maxLik(logLik=LogLik,start=out$estimate,response=log(CYIM$vol_t1),U=U,
             method="NM",control=list(iterlim=5000,printLevel=1),finalHessian=FALSE); 
  
  out=maxLik(logLik=LogLik,start=out$estimate,response=log(CYIM$vol_t1),U=U,
             method="BHHH",control=list(iterlim=5000,printLevel=2),finalHessian=FALSE); 
  
  coefs[[j]] = out$estimate; LL[j] = out$maximum;
  cat(j, "#--------------------------------------#",out$maximum,"\n"); 
}

j = min(which(LL==max(LL)));
out=maxLik(logLik=LogLik,start=coefs[[j]],response=log(CYIM$vol_t1),U=U,
           method="BHHH",control=list(iterlim=5000,printLevel=2),finalHessian=TRUE); 

######### save results of ML fit.  
names(out$estimate)<-colnames(U); coefs=out$estimate; SEs = sqrt(diag(vcov(out))); 

