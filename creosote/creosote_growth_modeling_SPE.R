rm(list=ls(all=TRUE));

setwd("c:/repos/IPM_size_transitions/creosote"); #Steve
# setwd("C:/Users/tm9/Desktop/git local/IPM_size_transitions/creosote"); #Tom

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

############################################################################
# Gaussian fits and fixed effect model selection 
############################################################################
LATR_gam_models=list()

LATR_gam_models[[1]] <- gam(list(log(vol_t1) ~s(log(vol_t)) + s(unique.transect,bs="re"), ~s(log(vol_t))), 
                data=LATR, family=gaulss())
LATR_gam_models[[2]] <- gam(list(log(vol_t1) ~s(log(vol_t)) + s(d.stand) + s(unique.transect,bs="re"), ~s(log(vol_t))), 
                data=LATR, family=gaulss())                
LATR_gam_models[[3]] <- gam(list(log(vol_t1) ~s(log(vol_t)) + s(d.stand) + d.stand:log(vol_t) + s(unique.transect,bs="re"), ~s(log(vol_t))), 
                data=LATR, family=gaulss())                

for(mod in 1:3) {
  fitGAU = LATR_gam_models[[mod]]
  fitted_all = predict(fitGAU,type="response",data=LATR);                  
  fitted_vals = new_fitted_vals = fitted_all[,1]; weights = fitted_all[,2]; 

  err=100; k=0; 
  while(err>0.000001) {
    LATR$fitted_vals = new_fitted_vals; 
    if(mod==1) fitGAU <- gam(list(log(vol_t1) ~s(log(vol_t)) + s(unique.transect,bs="re"), ~s(fitted_vals)), 
                data=LATR, family=gaulss())
    if(mod==2)  fitGAU <- gam(list(log(vol_t1) ~s(log(vol_t)) + s(d.stand) + s(unique.transect,bs="re"), ~s(fitted_vals)), 
                data=LATR, family=gaulss())                   
    if(mod==3) fitGAU <- gam(list(log(vol_t1) ~s(log(vol_t)) + s(d.stand) + d.stand:log(vol_t) + s(unique.transect,bs="re"), ~s(fitted_vals)), 
                data=LATR, family=gaulss())  
    fitted_all = predict(fitGAU,type="response",data=LATR);   
    new_fitted_vals = fitted_all[,1]; new_weights = fitted_all[,2];
    matplot(new_fitted_vals,cbind(weights,new_weights),type="p"); 
    err = weights - new_weights; err=sqrt(mean(err^2)); 
    k=k+1; cat(k,err,"\n"); 
    weights = new_weights; 
  }   
  LATR_gam_models[[mod]] =  fitGAU;
}

AIC(LATR_gam_models[[1]]); AIC(LATR_gam_models[[2]]); AIC(LATR_gam_models[[3]]); # model 2. 
LATR_gam_model <- LATR_gam_models[[2]]; 
LATR$fitted_vals <- fitted(LATR_gam_model); 

plot(LATR_gam_model); 

####  Extract values of the fitted splines to explore their properties 
fitted_all = predict(LATR_gam_model,type="response");  

##### Mean vs. initial size 
logz = log(LATR$vol_t); 
mean_fit1 = lm(fitted_all[,1]~logz); 


#### log(sigma) is fitted well be a quadratic (spline has df just above 2) 
sigma_hat = 1/fitted_vals[,2]; 
sd_fit1 = lm(log(sigma_hat)~z_vals); # R^2 = 0.97 
sd_fit2 = lm(log(sigma_hat)~z_vals+I(z_vals^2)); # R^2 = 0.999 






##### Inspect scaled residuals
scaledResids = residuals(LATR_lmer_best)*sqrt(best_weights) ## here was a problem: weights of LATR_lmer_best were all 1. 
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
## there is still a size trend in the stdev of the residuals. hmm. ### Fixed! 

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
n_bins <- 8
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
## going with LO and we'll see if a 2-param distribution can do the job
LATR_bin_fit <-LATR %>% 
  mutate(fitted = fitted(LATR_lmer_best),
         bin = as.integer(cut_number(fitted,n_bins))) %>% 
  mutate(mu=NA, sigma=NA,nu=NA,tau=NA)
for(b in 1:n_bins){
  bin_fit <- gamlssML(log(LATR_bin_fit$vol_t1[LATR_bin_fit$bin==b]) ~ 1,family="LO")
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
## maybe a quadratic term for sigma


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
               sigma = exp(pars2[1] + pars2[2]*mu + pars2[3]*mu^2),log=T) 
  return(val); 
}

paranoid_iter <- 3
coefs = list(paranoid_iter); LL=numeric(paranoid_iter);  

# Starting values from the pilot model are jittered to do multi-start optimization). 
# Using good starting values really speeds up convergence in the ML fits  
# Linear predictor coefficients extracted from the lmer model 
fixed_start = c(unlist(ranef(LATR_lmer_best)$unique.transect),fixef(LATR_lmer_best)[2:4])
## make sure the dimensions line up
length(fixed_start);ncol(U);colnames(U) 

fit_sigma = lm(log(sigma)~mu + I(mu^2), data=LATR_bin_fit)
#fit_nu = lm(log(nu)~mu, data=LATR_bin_fit)
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

# AIC improvement over gaussian
(AIC_LO <- 2*length(coefs) - 2*out$maximum)
(AIC_norm <- AIC(LATR_lmer_best) ) ## the logistic is quite an improvement
######### save results of ML fit.  
names(out$estimate)<-c(colnames(U),"sigma_b0","sigma_b1","sigma_b2")
coefs=out$estimate # parameters
V = vcov(out); SEs = sqrt(diag(V)); # standard errors

############# RFX shrinkage -- NEW, following the theory and code in Steve's appendix
n_transects <- length(unlist(ranef(LATR_lmer_best)$unique.transect))
transects=1:n_transects ## these are the indices of transect estimates in the parameter estimate vector
# Variance-covariance matrices for random intercepts
V1 = V[transects,transects]; 
# Extract year-specific intercepts, center them to zero
fixed.fx = coefs[transects]; fixed.fx = fixed.fx-mean(fixed.fx);
# Estimate sigma^2
var.hat = mean(fixed.fx^2) - mean(diag(V1)) + (sum(V1)-sum(diag(V1)))/(2*n_transects*(n_transects-1)); ## still comes out negative
# Shrink deviations from the mean
shrinkRanIntercept = fixed.fx*sqrt(var.hat/(var.hat + diag(V1)));

############# RFX shrinkage -- OLD, the shrunk variance is different, ask Steve why
#transects=1:length(unlist(ranef(LATR_lmer_best)$unique.transect))
#SEs = sqrt(diag(vcov(out))) 
#transect_fixed.fx = coefs[transects] - mean(coefs[transects])
#transect_fixed.se = SEs[transects]
#transect_sigma2.hat = mean(transect_fixed.fx^2)-mean(transect_fixed.se^2)
#transect_shrunk.fx = transect_fixed.fx*sqrt(transect_sigma2.hat/(transect_sigma2.hat + transect_fixed.se^2)) 
# lmer random effects for (1|year) 
#transect_ran.fx = ranef(LATR_lmer_best)["unique.transect"]


# visualize logistic model and moment diagnostics --------------------------------------
LATR_LO_kernel_median_density <- matrix(NA,size_dim,size_dim)
for(i in 1:size_dim){
  mu_size <- mean(coefs[transects]) + coefs["log(vol_t)"]*size_dum[i] + coefs["d.stand"]*median(LATR$d.stand) + coefs["I(d.stand^2)"]*median(LATR$d.stand)^2
  LATR_LO_kernel_median_density[i,] <- dNET(x = size_dum, 
                                  mu=mu_size,
                                  sigma = exp(coefs["sigma_b0"] + coefs["sigma_b1"]*mu_size + coefs["sigma_b2"]*mu_size^2)) 
}


levelplot(LATR_LO_kernel_median_density,row.values = size_dum, column.values = size_dum,cuts=30,
          col.regions=rainbow(30),xlab="log Size t",ylab="log Size t+1",main="Logistic",
          panel = function(...) {
            panel.levelplot(...)
            grid.points(log(LATR$vol_t), log(LATR$vol_t1), pch = ".",gp = gpar(cex=3,col=alpha("black",0.5)))
          })

# Simulate data from fitted LO model
MLmu = U%*%coefs[1:ncol(U)] 
n_sim <- 500
LATR_sim_LO<-LATR_sim_NO<-matrix(NA,nrow=nrow(LATR),ncol=n_sim)
for(i in 1:n_sim){
  LATR_sim_LO[,i] <- rLO(n = nrow(LATR), 
                           mu = MLmu, 
                           sigma = exp(coefs["sigma_b0"] + coefs["sigma_b1"]*MLmu  + coefs["sigma_b2"]*MLmu^2))
  LATR_sim_NO[,i] <- rnorm(n = nrow(LATR),
                               mean = predict(LATR_lmer_best),
                               sd = exp(pars[[best_model]][1] + pars[[best_model]][2]*predict(LATR_lmer_best) + pars[[best_model]][3]*predict(LATR_lmer_best)^2))
}
n_bins = 8
alpha_scale = 0.9
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

## visualize how well final model describes real moments of the data
par(mfrow=c(2,2),mar=c(4,4,2,1),cex.axis=1.3,cex.lab=1.3,mgp=c(2,1,0),bty="l"); 
sim_bin_means=sim_moment_means_LO=sim_moment_means_NO = matrix(NA,n_bins,n_sim); 
for(i in 1:n_sim){
  sim_moments <- bind_cols(LATR,data.frame(sim_LO=LATR_sim_LO[,i],
                                           sim_NO=LATR_sim_NO[,i])) %>% 
    arrange(log(vol_t)) %>% 
    mutate(size_bin = cut_number(log(vol_t),n=n_bins)) %>% 
    group_by(size_bin) %>% 
    summarise(mean_LO = mean(sim_LO),
              mean_NO = mean(sim_NO),
              bin_mean = mean(log(vol_t)))
  sim_bin_means[,i]=sim_moments$bin_mean; 
  sim_moment_means_LO[,i]=sim_moments$mean_LO; sim_moment_means_NO[,i]=sim_moments$mean_NO;		  
}

matplot(LATR_moments$bin_mean, sim_moment_means_LO,col=alpha("gray",0.5),pch=16,xlab="Mean size t0",ylab="mean(Size t1)",cex=1.4)
points(LATR_moments$bin_mean, apply(sim_moment_means_LO,1,median),pch=1,lwd=2,col="grey40",cex=1.4)
matplot(LATR_moments$bin_mean+0.2, sim_moment_means_NO,col=alpha("cornflowerblue",0.5),pch=16,add=T)
points(LATR_moments$bin_mean+0.2, apply(sim_moment_means_NO,1,median),pch=1,lwd=2,col="darkblue",cex=1.4)
points(LATR_moments$bin_mean+0.1, LATR_moments$mean_t1,pch=16,lwd=2,col="red",cex=1.4)
legend("topleft",legend=c("Logistic","Gaussian","Data"),
       col=c(alpha("gray",alpha_scale),alpha("cornflowerblue",alpha_scale),alpha("red",alpha_scale)),pch=1,lwd=2,bty="n"); 
add_panel_label("a")

for(i in 1:n_sim){
  sim_moments <- bind_cols(LATR,data.frame(sim_LO=LATR_sim_LO[,i],
                                           sim_NO=LATR_sim_NO[,i])) %>% 
    arrange(log(vol_t)) %>% 
    mutate(size_bin = cut_number(log(vol_t),n=n_bins)) %>% 
    group_by(size_bin) %>% 
    summarise(mean_LO = sd(sim_LO),
              mean_NO = sd(sim_NO),
              bin_mean = mean(log(vol_t)))
  sim_bin_means[,i]=sim_moments$bin_mean; 
  sim_moment_means_LO[,i]=sim_moments$mean_LO; sim_moment_means_NO[,i]=sim_moments$mean_NO;		  
}

matplot(LATR_moments$bin_mean, sim_moment_means_LO,col=alpha("gray",0.5),pch=16,xlab="Mean size t0",ylab="SD(Size t1)",cex=1.4)
points(LATR_moments$bin_mean, apply(sim_moment_means_LO,1,median),pch=1,lwd=2,col="grey40",cex=1.4)
matplot(LATR_moments$bin_mean+0.2, sim_moment_means_NO,col=alpha("cornflowerblue",0.5),pch=16,add=T)
points(LATR_moments$bin_mean+0.2, apply(sim_moment_means_NO,1,median),pch=1,lwd=2,col="darkblue",cex=1.4)
points(LATR_moments$bin_mean+0.1, LATR_moments$sd_t1,pch=16,lwd=2,col=alpha("red",alpha_scale),cex=1.4)
add_panel_label("b")

for(i in 1:n_sim){
  sim_moments <- bind_cols(LATR,data.frame(sim_LO=LATR_sim_LO[,i],
                                           sim_NO=LATR_sim_NO[,i])) %>% 
    arrange(log(vol_t)) %>% 
    mutate(size_bin = cut_number(log(vol_t),n=n_bins)) %>% 
    group_by(size_bin) %>% 
    summarise(mean_LO = NPskewness(sim_LO),
              mean_NO = NPskewness(sim_NO),
              bin_mean = mean(log(vol_t)))
  sim_bin_means[,i]=sim_moments$bin_mean; 
  sim_moment_means_LO[,i]=sim_moments$mean_LO; sim_moment_means_NO[,i]=sim_moments$mean_NO;		  
}

matplot(LATR_moments$bin_mean, sim_moment_means_LO,col=alpha("gray",0.5),pch=16,xlab="Mean size t0",ylab="Skew(Size t1)",cex=1.4)
points(LATR_moments$bin_mean, apply(sim_moment_means_LO,1,median),pch=1,lwd=2,col="grey40",cex=1.4)
matplot(LATR_moments$bin_mean+0.2, sim_moment_means_NO,col=alpha("cornflowerblue",0.5),pch=16,add=T)
points(LATR_moments$bin_mean+0.2, apply(sim_moment_means_NO,1,median),pch=1,lwd=2,col="darkblue",cex=1.4)
points(LATR_moments$bin_mean+0.1, LATR_moments$skew_t1,pch=16,lwd=2,col=alpha("red",alpha_scale),cex=1.4)
add_panel_label("c")

for(i in 1:n_sim){
  sim_moments <- bind_cols(LATR,data.frame(sim_LO=LATR_sim_LO[,i],
                                           sim_NO=LATR_sim_NO[,i])) %>% 
    arrange(log(vol_t)) %>% 
    mutate(size_bin = cut_number(log(vol_t),n=n_bins)) %>% 
    group_by(size_bin) %>% 
    summarise(mean_LO = NPkurtosis(sim_LO),
              mean_NO = NPkurtosis(sim_NO),
              bin_mean = mean(log(vol_t)))
  sim_bin_means[,i]=sim_moments$bin_mean; 
  sim_moment_means_LO[,i]=sim_moments$mean_LO; sim_moment_means_NO[,i]=sim_moments$mean_NO;		  
}

matplot(LATR_moments$bin_mean, sim_moment_means_LO,col=alpha("gray",0.5),pch=16,xlab="Mean size t0",ylab="Kurtosis(Size t1)",cex=1.4)
points(LATR_moments$bin_mean, apply(sim_moment_means_LO,1,median),pch=1,lwd=2,col="grey40",cex=1.4)
matplot(LATR_moments$bin_mean+0.2, sim_moment_means_NO,col=alpha("cornflowerblue",0.5),pch=16,add=T)
points(LATR_moments$bin_mean+0.2, apply(sim_moment_means_NO,1,median),pch=1,lwd=2,col="darkblue",cex=1.4)
points(LATR_moments$bin_mean+0.1, LATR_moments$kurt_t1,pch=16,lwd=2,col=alpha("red",alpha_scale),cex=1.4)
add_panel_label("d")

