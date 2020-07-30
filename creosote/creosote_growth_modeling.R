rm(list=ls(all=TRUE));

setwd("c:/repos/IPM_size_transitions/creosote"); #Steve
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


# Gaussian fits and fixed effect model selection --------------------------
LATR_lmer_models <- list() 
## simple model of size dependence
LATR_lmer_models[[1]] <- lmer(log(vol_t1) ~ log(vol_t) + (1|unique.transect), data=LATR,REML=F,control=lmerControl(optimizer="bobyqa"))
## quadratic term for size
LATR_lmer_models[[2]] <- lmer(log(vol_t1) ~ log(vol_t) + I(log(vol_t)^2) + (1|unique.transect), data=LATR,REML=F,control=lmerControl(optimizer="bobyqa"))
## size and density
LATR_lmer_models[[3]] <- lmer(log(vol_t1) ~ log(vol_t) + d.stand + (1|unique.transect), data=LATR,REML=F,control=lmerControl(optimizer="bobyqa"))
## size^2 and density
LATR_lmer_models[[4]] <- lmer(log(vol_t1) ~ log(vol_t) + I(log(vol_t)^2) + d.stand + (1|unique.transect), data=LATR,REML=F,control=lmerControl(optimizer="bobyqa"))
## size-density interaction
LATR_lmer_models[[5]] <- lmer(log(vol_t1) ~ log(vol_t) * d.stand + (1|unique.transect), data=LATR,REML=F,control=lmerControl(optimizer="bobyqa"))
## size and density interaction plus size2 (convergence issues with size2*density)
LATR_lmer_models[[6]] <- lmer(log(vol_t1) ~ log(vol_t) * d.stand + I(log(vol_t)^2) + (1|unique.transect), data=LATR,REML=F,control=lmerControl(optimizer="bobyqa"))
## quadratic term for density
LATR_lmer_models[[7]] <- lmer(log(vol_t1) ~ log(vol_t) + d.stand + I(d.stand^2) + (1|unique.transect), data=LATR,REML=F,control=lmerControl(optimizer="bobyqa"))
## quadratic terms for density and size plus size-density interaction
LATR_lmer_models[[8]] <- lmer(log(vol_t1) ~ log(vol_t)*d.stand + I(log(vol_t)^2) + I(d.stand^2) + (1|unique.transect), data=LATR,REML=F,control=lmerControl(optimizer="bobyqa"))
## quadratic terms for density and size
## this model will fit but the iterative re-weighting gets caught where the weights do not improve much, so skipping it
#LATR_lmer_models[[9]] <- lmer(log(vol_t1) ~ log(vol_t) + I(log(vol_t)^2) + d.stand + I(d.stand^2) + (1|unique.transect), data=LATR,REML=F,control=lmerControl(optimizer="bobyqa"))

## for now making variance a simple linear function of fitted value
varPars = function(pars) {
  return(-sum(dnorm(resids, mean=0, sd=exp(pars[1] + pars[2]*fitted_vals),log=TRUE)))
}	
## iterative re-weighting
pars<-list()
for(mod in 1:length(LATR_lmer_models)) {
  err = 1; rep=0; 
  while(err > 0.000001) {
    rep=rep+1; model = LATR_lmer_models[[mod]];
    fitted_vals = fitted(model);resids = residuals(model); 
    out=optim(c(sd(resids),0),varPars,control=list(maxit=5000)); 
    pars[[mod]]=out$par; 
    new_sigma = exp(pars[[mod]][1] + pars[[mod]][2]*fitted_vals); new_weights = 1/((new_sigma)^2)
    new_weights = 0.5*(weights(model) + new_weights); # cautious update 
    new_model <- update(model,weights=new_weights); 
    err = weights(model)-weights(new_model); err=sqrt(mean(err^2)); 
    cat(mod,rep,err,"\n") # check on convergence of estimated weights 
    LATR_lmer_models[[mod]]<-new_model; 
  }}

######### For a fair AIC comparison, fit all models with the same weights 
aics = unlist(lapply(LATR_lmer_models,AIC)); best_model=which(aics==min(aics)); 
best_weights=weights(LATR_lmer_models[[best_model]]); 
for(mod in 1:length(LATR_lmer_models)) {
  LATR_lmer_models[[mod]] <- update(LATR_lmer_models[[mod]],weights=best_weights)
}
AICtab(LATR_lmer_models)

# Model 7 is favored (density and density^2, no size*density). Last step is to re-fit with REML=T.
aics = unlist(lapply(LATR_lmer_models,AIC)); best_model=which(aics==min(aics));
LATR_lmer_best = LATR_lmer_models[[best_model]] 
best_weights = weights(LATR_lmer_best)
LATR_lmer_best <- lmer(log(vol_t1) ~ log(vol_t) + d.stand + I(d.stand^2) + (1|unique.transect), data=LATR,REML=F,control=lmerControl(optimizer="bobyqa"))
##refit the residuals as a function of mean
fitted_vals = fitted(LATR_lmer_best);resids = residuals(LATR_lmer_best)
best_pars <- optim(c(sd(resids),0),varPars,control=list(maxit=5000))

##### Inspect scaled residuals
scaledResids = residuals(LATR_lmer_best)*sqrt(weights(LATR_lmer_best))
par(mfrow=c(1,2))
plot(fitted(LATR_lmer_best), scaledResids) 
qqPlot(scaledResids) # really bad in both tails
jarque.test(scaledResids) # normality test: FAILS, P < 0.001 
anscombe.test(scaledResids) # kurtosis: FAILS, P < 0.001 
agostino.test(scaledResids) # skewness: FAILS, P<0.001 

px = fitted(LATR_lmer_best); py=scaledResids; 
par(mfrow=c(2,2),bty="l",mar=c(4,4,2,1),mgp=c(2.2,1,0),cex.axis=1.4,cex.lab=1.4);   
##### Alternatively, use nonparametric measures of skew and excess kurtosis. 
z = rollMomentsNP(px,py,windows=8,smooth=TRUE,scaled=TRUE) 
## there is still a size trend in the stdev of the residuals. hmm.
plot(LATR$d.stand,scaledResids) #-- there is clearly greater variance at low density, which is what I would expect
plot(log(LATR$vol_t),scaledResids) 

# visualize kernel
##dummy variable for initial size
size_dim <- 100
size_dum <- seq(min(log(LATR$vol_t)),max(log(LATR$vol_t)),length.out = size_dim)
##make a heatmap for the Gaussian kernel with non-constant variance -- updated with the quadratic terms for mean and variance
## Because there is density dependence, I will plot a low-density kernel
LATR_lmer_best_kernel_median_density <- matrix(NA,size_dim,size_dim)
for(i in 1:size_dim){
  mu_size <- fixef(LATR_lmer_best)[1] + fixef(LATR_lmer_best)[2] * size_dum[i] + fixef(LATR_lmer_best)[3]*median(LATR$d.stand) + fixef(LATR_lmer_best)[4]*median(LATR$d.stand)^2
  LATR_lmer_best_kernel_median_density[i,] <- dnorm(size_dum,
                                     mean = mu_size,
                                     sd = exp(best_pars$par[1] + best_pars$par[2]*mu_size))
}

levelplot(LATR_lmer_best_kernel_median_density,row.values = size_dum, column.values = size_dum,cuts=30,
          col.regions=rainbow(30),xlab="log Size t",ylab="log Size t+1",main="Gaussian, non-constant variance",
          panel = function(...) {
            panel.levelplot(...)
            grid.points(log(LATR$vol_t), log(LATR$vol_t1), pch = ".",gp = gpar(cex=3,col=alpha("black",0.5)))
          }) 
## hmmm, this does not look like a great fit

# Finding a better distribution -------------------------------------------
n_bins <- 6
select_dist <- tibble(fit_best = fitted(LATR_lmer_best),
                      scale_resid = residuals(LATR_lmer_best)*sqrt(weights(LATR_lmer_best))) %>% 
  mutate(bin = as.integer(cut_number(fit_best,n_bins)),
         best_dist = NA,
         secondbest_dist = NA,
         aic_margin = NA) 
for(b in 1:n_bins){
  bin_fit <- fitDist(select_dist$scale_resid[select_dist$bin==b],type="realline")
  select_dist$best_dist[select_dist$bin==b] <- names(bin_fit$fits[1])
  select_dist$secondbest_dist[select_dist$bin==b] <- names(bin_fit$fits[2])
  select_dist$aic_margin[select_dist$bin==b] <- bin_fit$fits[2] - bin_fit$fits[1]
}
select_dist %>% 
  group_by(bin) %>% 
  summarise(n_bin = n(),
            best_dist = unique(best_dist),
            secondbest_dist = unique(secondbest_dist),
            aic_margin = unique(aic_margin))
## TF, LO,  NET show up a lot, and this makes sense because the roll moments plot showed that skewness is not bad but kurtosis is a problem
## fit the TF by size bin (ignore density variation for now)
LATR_bin_fit <-LATR %>% 
  mutate(fitted = fitted(LATR_lmer_best),
         bin = as.integer(cut_number(fitted,n_bins))) %>% 
  mutate(mu=NA, sigma=NA,nu=NA,tau=NA)
for(b in 1:n_bins){
  bin_fit <- gamlssML(log(LATR_bin_fit$vol_t1[LATR_bin_fit$bin==b]) ~ 1,family="TF")
  LATR_bin_fit$mu[LATR_bin_fit$bin==b] <- bin_fit$mu
  LATR_bin_fit$sigma[LATR_bin_fit$bin==b] <- bin_fit$sigma
  #LATR_bin_fit$nu[LATR_bin_fit$bin==b] <- bin_fit$nu
  #LATR_bin_fit$tau[LATR_bin_fit$bin==b] <- bin_fit$tau
}
LATR_bin_fit %>% 
  group_by(bin) %>% 
  summarise(N = n(),
            mean_fitted = mean(fitted),
            mu = unique(mu),
            sigma=unique(sigma),
            nu=unique(nu),
            tau=unique(tau)) -> LATR_bin_fit

par(mfrow=c(2,2),bty="l",mar=c(4,4,2,1),mgp=c(2.2,1,0),cex.axis=1.4,cex.lab=1.4);
## Steve's spline.scatter.smooth() function not working for me and I did not both trying to figure out why
plot(LATR_bin_fit$mean_fitted,LATR_bin_fit$mu,xlab="Fitted value",ylab=expression(paste("Location parameter  ", mu )),type="b")
plot(LATR_bin_fit$mu,LATR_bin_fit$sigma,xlab=expression(paste("Location parameter  ", mu )),
     ylab=expression(paste("Scale parameter  ", sigma)),type="b")


# -------------------------------------------------------------------------
## as an alternative to fitDIst, try using Steve's improved fitDist
source("../fitChosenDists.R")
## these are gamlss' "realline" distributions
tryDists <- c("NO", "GU", "RG" ,"LO", "NET", "TF", "TF2", "PE","PE2", "SN1", "SN2", "exGAUS", "SHASH", "SHASHo","SHASHo2", "EGB2", "JSU", "JSUo", "SEP1", "SEP2", "SEP3", "SEP4", "ST1", "ST2", "ST3", "ST4", "ST5", "SST", "GT")
tryDensities <- paste("d",realline,sep="")
tryDensities <- list(dNO,dGU)

tryDists=c("EGB2","GT","JSU", "SHASHo","SEP1","SEP2","SEP3","SEP4"); 
tryDensities=list(dEGB2, dGT, dJSU, dSHASHo, dSEP1, dSEP2, dSEP3, dSEP4); 

bins = 1:n_bins
maxVals = matrix(NA,n_bins,length(tryDists))
for(j in 1:length(bins)){
  for(k in 1:length(tryDists)) {
    fitj = gamlssMaxlik(y=select_dist$scale_resid[select_dist$bin==j],
                        DIST=tryDists[k],
                        density=tryDensities[[k]]) 
    maxVals[j,k] = fitj$maximum
    cat("Finished ", tryDists[k]," ",j,k, fitj$maximum,"\n") 
  }
}

## best two for each bin 
for(j in 1:length(bins)){
  e = order(-maxVals[j,]); 
  cat(j, tryDists[e][1:2],"\n"); 
}	

# overall ranking 
e = order(-colSums(maxVals)); 
rbind(tryDists[e],round(colSums(maxVals)[e],digits=3)); 








# Fitting the final model -------------------------------------------------
# Now we can fit a custom model via maximum likelihood, matching the structure of the best lmer 
# model but fitting the random effects as fixed instead and estimating the corresponding variances 
# with Steve's shrinkage methods. 

#First we will defined the linear predictor for the location parameter mu 
#(which is not necessarily the expected value). Note that this with parameterzation, the intercept is year1/plot1.
U=model.matrix(~  0 + unique.transect + log(vol_t) + d.stand + I(d.stand^2), data=LATR)

# Next define a likelihood function using this linear predictor for the location. I am including 
# quadratic terms for sigma and nu based on the binned fits above.
LogLik=function(pars,response,U){
  pars1 = pars[1:ncol(U)]; pars2=pars[-(1:ncol(U))];
  mu = U%*%pars1;  
  val = dLO(x = response, 
               mu=mu,
               sigma = exp(pars2[1] + pars2[2]*mu + pars2[3]*mu^2), log=T) 
  return(val); 
}


paranoid_iter <- 1
coefs = list(paranoid_iter); LL=numeric(paranoid_iter);  

# Starting values from the pilot model are jittered to do multi-start optimization). 
# Using good starting values really speeds up convergence in the ML fits  
# Linear predictor coefficients extracted from the lmer model 
fixed_start = c(unlist(ranef(LATR_lmer_best)$unique.transect),
                fixef(LATR_lmer_best)[2],fixef(LATR_lmer_best)[3],fixef(LATR_lmer_best)[4])
## make sure the dimensions line up
length(fixed_start);ncol(U);colnames(U) 

fit_sigma = lm(log(sigma)~mu + I(mu^2), data=LATR_bin_fit)
p0=c(fixed_start, coef(fit_sigma))

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

######### save results of ML fit.  
names(out$estimate)<-c(colnames(U),"sigma_b0","sigma_b1","sigma_b2")
coefs=out$estimate
## these are the indices of plot and year effects, to make the next steps a little more intuitive
transects=1:12
SEs = sqrt(diag(vcov(out))) 
AIC_LO <- 2*length(coefs) - 2*out$maximum
AIC_norm <- AIC(LATR_lmer_best)

############# RFX shrinkage
transect_fixed.fx = coefs[transects] - mean(coefs[transects])
transect_fixed.se = SEs[transects]
transect_sigma2.hat = mean(transect_fixed.fx^2)-mean(transect_fixed.se^2)
transect_shrunk.fx = transect_fixed.fx*sqrt(transect_sigma2.hat/(transect_sigma2.hat + transect_fixed.se^2)) 
# lmer random effects for (1|year) 
transect_ran.fx = ranef(LATR_lmer_best)["unique.transect"]

plot(year_ran.fx$year_t$`(Intercept)`,year_shrunk.fx,xlab="lmer year random effects",ylab="Shrunk year fixed effects",type="n")
text(year_ran.fx$year_t$`(Intercept)`,year_shrunk.fx,labels=rownames(year_ran.fx$year_t))
abline(0,1,col="blue",lty=2);

tibble(sd_estimate = c(sd(year_fixed.fx),sd(year_shrunk.fx),sd(year_ran.fx$year_t$`(Intercept)`)),
       method = c("fixed","shrunk","lme4"))

