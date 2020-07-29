rm(list=ls(all=TRUE));

setwd("c:/repos/IPM_size_transitions/cactus"); #Steve
setwd("C:/Users/tm9/Desktop/git local/IPM_size_transitions/cactus"); #Tom

require(car); require(lme4); require(zoo); require(moments); require(mgcv); 
require(gamlss); require(gamlss.tr); require(AICcmodavg); 
require(lmerTest); require(tidyverse); require(maxLik); 
require(actuar); require(lattice); require(grid); require(scales);
require(sgt); require(formatR); require(popbio)

# function for converting cactus size measurements to volume
volume <- function(h, w, p){
  (1/3)*pi*h*(((w + p)/2)/2)^2
}

invlogit <- function(x){exp(x)/(1+exp(x))}

# log kurtosis function for diagnostics
Lkurtosis=function(x) log(kurtosis(x)); 

# Steve's diagnostics functions
source("../Diagnostics.R")

# read in data
CYIM_full<-read_csv("cholla_demography_20042018_EDI.csv") %>% 
  ## filter out transplants and drop seed addition plots (which start with letter H)
  ## also drop Year_t==2018 because I don't have 2019 size data (not entered yet). 2018 data still included in 2017-2018 transition.
  ## For now I am dropping the early years from plots T1-T3 because these cause problems with shrinkage parameter estimates.
    filter(Transplant == 0,
         str_sub(Plot,1,1)!="H",
         Year_t > 2007,
         Year_t!=2018) %>% 
  ## convert height, max width, perp width to volume of cone, take natural log
  mutate(vol_t = volume(Height_t,Width_t,Perp_t),
         vol_t1 = volume(Height_t1,Width_t1,Perp_t1),
         plot = as.factor(Plot),
         year_t = as.factor(Year_t),
         ID = interaction(TagID,plot)) %>%
  select(ID,year_t,plot,vol_t,vol_t1,Survival_t1,Goodbuds_t1) %>% 
  ## sort by initial size
  arrange(vol_t) 

## pull out and na.omit size transitions for growth modeling
CYIM_full %>% 
  select(ID,year_t,plot,vol_t,vol_t1) %>% 
  ## drop rows with NAs
  drop_na() -> CYIM

table(CYIM$plot,CYIM$year_t)


# Gaussian fits -----------------------------------------------------------
CYIM_lmer_models <- list() 
CYIM_lmer_models[[1]] <- lmer(log(vol_t1) ~ log(vol_t) + (1|year_t) + (1|plot), data=CYIM,REML=F,control=lmerControl(optimizer="bobyqa"))
CYIM_lmer_models[[2]] <- lmer(log(vol_t1) ~ log(vol_t) + I(log(vol_t)^2) + (1|year_t) + (1|plot), data=CYIM,REML=F,control=lmerControl(optimizer="bobyqa"))

## Now use the iterative re-weighting approach to re-fit these with non-constant variance:
## NegLogLik function to fit variance model for residuals 
## SPE: based on the residual diagnostic plot for a linear model of log(sigma), a quadratic term is added. 
varPars = function(pars) {
  return(-sum(dnorm(resids, mean=0, sd=exp(pars[1] + pars[2]*fitted_vals + pars[3]*fitted_vals^2),log=TRUE)))
}	

pars<-list()
for(mod in 1:length(CYIM_lmer_models)) {
  err = 1; rep=0; 
  while(err > 0.000001) {
    rep=rep+1; model = CYIM_lmer_models[[mod]];
    fitted_vals = fitted(model);resids = residuals(model); 
    out=optim(c(sd(resids),0,0),varPars,control=list(maxit=5000)); 
    pars[[mod]]=out$par; 
    new_sigma = exp(pars[[mod]][1] + pars[[mod]][2]*fitted_vals+ pars[[mod]][3]*fitted_vals^2); new_weights = 1/((new_sigma)^2)
    new_weights = 0.5*(weights(model) + new_weights); # cautious update 
    new_model <- update(model,weights=new_weights); 
    err = weights(model)-weights(new_model); err=sqrt(mean(err^2)); 
    cat(mod,rep,err,"\n") # check on convergence of estimated weights 
    CYIM_lmer_models[[mod]]<-new_model; 
  }}


######### For a fair AIC comparison, fit all models with the same weights 
aics = unlist(lapply(CYIM_lmer_models,AIC)); best_model=which(aics==min(aics)); 
best_weights=weights(CYIM_lmer_models[[best_model]]); 
for(mod in 1:length(CYIM_lmer_models)) {
  CYIM_lmer_models[[mod]] <- update(CYIM_lmer_models[[mod]],weights=best_weights)
  }
AIC(CYIM_lmer_models[[1]],CYIM_lmer_models[[2]])


# The quadratic model is favored. Last step is to re-fit with REML=T.
aics = unlist(lapply(CYIM_lmer_models,AIC)); best_model=which(aics==min(aics));
CYIM_lmer_best = CYIM_lmer_models[[best_model]] 
best_weights = weights(CYIM_lmer_best)
CYIM_lmer_best <- lmer(log(vol_t1) ~ log(vol_t) + I(log(vol_t)^2) + (1|year_t) + (1|plot), data=CYIM,weights=best_weights,REML=TRUE) 
##refit the residuals as a function of mean
fitted_vals = fitted(CYIM_lmer_best);resids = residuals(CYIM_lmer_best)
best_pars <- optim(c(sd(resids),0,0),varPars,control=list(maxit=5000))

# Now visualize this kernel. This model is the *very best* we could do in a Gaussian framework 
# (assuming there are no major fixed effects we are missing).
##dummy variable for initial size
size_dim <- 100
size_dum <- seq(min(log(CYIM$vol_t)),max(log(CYIM$vol_t)),length.out = size_dim)
##make a heatmap for the Gaussian kernel with non-constant variance -- updated with the quadratic terms for mean and variance
CYIM_lmer_best_kernel <- matrix(NA,size_dim,size_dim)
for(i in 1:size_dim){
  mu_size <- fixef(CYIM_lmer_best)[1] + fixef(CYIM_lmer_best)[2] * size_dum[i] + fixef(CYIM_lmer_best)[3]*size_dum[i]^2
  CYIM_lmer_best_kernel[,i] <- dnorm(size_dum,
                                     mean = mu_size,
                                     sd = exp(best_pars$par[1] + best_pars$par[2]*mu_size + best_pars$par[3]*mu_size^2))
}

levelplot(CYIM_lmer_best_kernel,row.values = size_dum, column.values = size_dum,cuts=30,
          col.regions=rainbow(30),xlab="log Size t",ylab="log Size t+1",main="Gaussian, non-constant variance",
          panel = function(...) {
            panel.levelplot(...)
            grid.points(log(CYIM$vol_t), log(CYIM$vol_t1), pch = ".",gp = gpar(cex=3,col=alpha("black",0.5)))
          }) 


scaledResids = residuals(CYIM_lmer_best)*sqrt(weights(CYIM_lmer_best))
par(mfrow=c(1,2))
plot(fitted(CYIM_lmer_best), scaledResids) 
qqPlot(scaledResids) # really bad in both tails
jarque.test(scaledResids) # normality test: FAILS, P < 0.001 
anscombe.test(scaledResids) # kurtosis: FAILS, P < 0.001 
agostino.test(scaledResids) # skewness: FAILS, P<0.001 

# One last look at the standardized residuals. They are roughly mean zero and unit variance -- 
# so that checks out. But there is negative skew and excess kurtosis, especially at large sizes. 
px = fitted(CYIM_lmer_best); py=scaledResids; 
par(mfrow=c(2,2),bty="l",mar=c(4,4,2,1),mgp=c(2.2,1,0),cex.axis=1.4,cex.lab=1.4);   
z = rollMoments(px,py,windows=8,smooth=TRUE,scaled=TRUE) 

##### Alternatively, use nonparametric measures of skew and excess kurtosis. 
z = rollMomentsNP(px,py,windows=8,smooth=TRUE,scaled=TRUE) 

# Finding a better distribution -------------------------------------------
# We now know the Gaussian provides a poor fit to the residual variance. We would like to know which 
# distribution provides a better - ideally _good_ - fit. Because there is size-dependence in skew and 
# kurtosis, we cannot marginalize over the entire distribution of residuals (this may point 
# us in the wrong direction). Instead, we can slice up the data into bins of expected 
# value and find the best distribution for each bin using gamlss' fitDist(). 
n_bins <- 8
select_dist <- tibble(fit_best = fitted(CYIM_lmer_best),
                      scale_resid = residuals(CYIM_lmer_best)*sqrt(weights(CYIM_lmer_best))) %>% 
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

# It's a little bit of everything and convergence is not great (I've suppressed those warnings 
# in this output). I am proceeding with the SHASH, which is support for 3/8 groups
CYIM_bin_fit <-CYIM %>% 
  mutate(fitted = fitted(CYIM_lmer_best),
         bin = as.integer(cut_number(fitted,n_bins))) %>% 
  mutate(mu=NA, sigma=NA,nu=NA,tau=NA)
for(b in 1:n_bins){
  bin_fit <- gamlssML(log(CYIM_bin_fit$vol_t1[CYIM_bin_fit$bin==b]) ~ 1,family="SHASH")
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


# Fitting the final model -------------------------------------------------
# Now we can fit a custom model via maximum likelihood, matching the structure of the best lmer 
# model but fitting the random effects as fixed instead and estimating the corresponding variances 
# with Steve's shrinkage methods. 

#First we will defined the linear predictor for the location parameter mu 
#(which is not necessarily the expected value). Note that this with parameterzation, the intercept is year1/plot1.
U=model.matrix(~  0 + year_t + plot + log(vol_t)+ I(log(vol_t)^2), data=CYIM)

# Next define a likelihood function using this linear predictor for the location. I am including 
# quadratic terms for sigma and nu based on the binned fits above.
LogLik=function(pars,response,U){
  pars1 = pars[1:ncol(U)]; pars2=pars[-(1:ncol(U))];
  mu = U%*%pars1;  
  val = dSHASH(x = response, 
             mu=mu,
             sigma = exp(pars2[1] + pars2[2]*mu + pars2[3]*mu^2), 
             nu = exp(pars2[4] + pars2[5]*mu + pars2[6]*mu^2), 
             tau = exp(pars2[7] + pars2[8]*mu), log=T) 
  return(val); 
}

# Now fit it, using the lmer() fit for the mean (location) and lm() fits for the higher moments as starting parameter values. 
# This comes straight from Steve's PSSP code, including his convergence paranoia, which is probably warranted. 
# In preliminary fits, all the iterations converged on the same likelihoods, and this takes forever, 
# so I am skimping on the paranoia here, but we might return to it to make this a good example to follow. 
# This could also be a model selection step; easy to pull AIC from the maxLik fit. 

paranoid_iter <- 1
coefs = list(paranoid_iter); LL=numeric(paranoid_iter);  

# Starting values from the pilot model are jittered to do multi-start optimization). 
# Using good starting values really speeds up convergence in the ML fits  
# Linear predictor coefficients extracted from the lmer model 
fixed_start = c(unlist(ranef(CYIM_lmer_best)$year_t) + unlist(ranef(CYIM_lmer_best)$plot)[1], #year estimates, conditioned on plot 1
                unlist(ranef(CYIM_lmer_best)$plot)[-1],
                fixef(CYIM_lmer_best)[2],fixef(CYIM_lmer_best)[3])
## make sure the dimensions line up
length(fixed_start);ncol(U);colnames(U) 

# Shape and scale coefficients from the rollaply diagnostic plots 
fit_sigma = lm(log(sigma)~mu + I(mu^2), data=CYIM_bin_fit)
fit_nu = lm(log(nu)~mu + I(mu^2), data=CYIM_bin_fit)
fit_tau = lm(log(tau)~mu, data=CYIM_bin_fit)
p0=c(fixed_start, coef(fit_sigma), coef(fit_nu),coef(fit_tau))

for(j in 1:paranoid_iter) {
  out=maxLik(logLik=LogLik,start=p0*exp(0.2*rnorm(length(p0))), response=log(CYIM$vol_t1),U=U,
             method="BHHH",control=list(iterlim=5000,printLevel=2),finalHessian=FALSE); 
  
  out=maxLik(logLik=LogLik,start=out$estimate,response=log(CYIM$vol_t1),U=U,
             method="NM",control=list(iterlim=5000,printLevel=1),finalHessian=FALSE); 
  
  out=maxLik(logLik=LogLik,start=out$estimate,response=log(CYIM$vol_t1),U=U,
             method="BHHH",control=list(iterlim=5000,printLevel=2),finalHessian=FALSE); 
  
  coefs[[j]] = out$estimate; LL[j] = out$maximum;
  cat(j, "#--------------------------------------#",out$maximum,"\n"); 
}

j = min(which(LL==max(LL))) ## they actually all land on the same likelihood-that's good!
out=maxLik(logLik=LogLik,start=coefs[[j]],response=log(CYIM$vol_t1),U=U,
           method="BHHH",control=list(iterlim=5000,printLevel=2),finalHessian=TRUE) 

######### save results of ML fit.  
names(out$estimate)<-c(colnames(U),"sigma_b0","sigma_b1","sigma_b2","nu_b0","nu_b1","nu_b2","tau_b0","tau_b1")
coefs=out$estimate
## these are the indices of plot and year effects, to make the next steps a little more intuitive
years=1:9
plots=10:16
SEs = sqrt(diag(vcov(out))) 
AIC_SHASH <- 2*length(coefs) - 2*out$maximum
AIC_norm <- AIC(CYIM_lmer_best)

# The AIC comparison is a blowout, though I am not sure if it is ok to 
# use the lmer AIC. Might be better to re-fit with plot and year as fixed. 
tibble(growth_function=c("Gaussian","SHASH"),AIC = c(AIC_norm,AIC_SHASH))

# Now use Steve's shrinkage code to get the random effect variances and compare these to the lmer estimates. 
# I will need Steve to explain to me what the shrinkage step is doing, and what is the theory for this (or I should do my own homework). 
# shrinkage random effects for (1|year) 
year_fixed.fx = coefs[years] - mean(coefs[years])
year_fixed.se = SEs[years]
year_sigma2.hat = mean(year_fixed.fx^2)-mean(year_fixed.se^2)
year_shrunk.fx = year_fixed.fx*sqrt(year_sigma2.hat/(year_sigma2.hat + year_fixed.se^2)) 
# lmer random effects for (1|year) 
year_ran.fx = ranef(CYIM_lmer_best)["year_t"]

plot(year_ran.fx$year_t$`(Intercept)`,year_shrunk.fx,xlab="lmer year random effects",ylab="Shrunk year fixed effects",type="n")
text(year_ran.fx$year_t$`(Intercept)`,year_shrunk.fx,labels=rownames(year_ran.fx$year_t))
abline(0,1,col="blue",lty=2);

tibble(sd_estimate = c(sd(year_fixed.fx),sd(year_shrunk.fx),sd(year_ran.fx$year_t$`(Intercept)`)),
       method = c("fixed","shrunk","lme4"))

# The code for plots is more complicated than years, because plot effects were parameterized as contrasts. 
# What is the SE of the plot 1 effect?? Is it the mean of the year SEs? -- probably not! 
# I am going to drop plot 1 from the shrunk effects, because I don't know how to calculate its SE. 
# This will bias the variance, hopefully not too much. 

# These estimates don't look as good as the year effects but, again, the shrinkage variance corresponds well to the lme4 esimate.
# shrinkage random effects for (1|plot) 
plot_coefs <- mean(coefs[years])+coefs[plots]
plot_fixed.fx = plot_coefs - mean(plot_coefs)
plot_fixed.se = SEs[plots]
plot_sigma2.hat = mean(plot_fixed.fx^2)-mean(plot_fixed.se^2)
plot_shrunk.fx = plot_fixed.fx*sqrt(plot_sigma2.hat/(plot_sigma2.hat + plot_fixed.se^2)) 
# lmer random effects for (1|plot) 
plot_ran.fx = ranef(CYIM_lmer_best)["plot"]

plot(plot_ran.fx$plot$`(Intercept)`,c(NA,plot_shrunk.fx),xlab="lmer year random effects",ylab="Shrunk year fixed effects",type="n")
text(plot_ran.fx$plot$`(Intercept)`,c(NA,plot_shrunk.fx),labels=rownames(plot_ran.fx$plot))
abline(0,1,col="blue",lty=2);

tibble(sd_estimate = c(sd(plot_fixed.fx),sd(plot_shrunk.fx),sd(plot_ran.fx$plot$`(Intercept)`)),
       method = c("fixed","shrunk","lme4"))

# Steve had an interesting idea for an alternative way to fit two random effects, with one used as an offset for the other. 
# Might swing back to that, especially since the reference-level problem is nagging at me. 

## Visual diagnostics of model fit via shrinkage
## Here is the top-level view of the skewed t kernel. This will be the mean kernel, averaged over years and plots.
CYIM_SHASH_kernel <- matrix(NA,size_dim,size_dim)
for(i in 1:size_dim){
  mu_size <- mean(coefs[c(years,plots)]) + coefs["log(vol_t)"] * size_dum[i] + coefs["I(log(vol_t)^2)"] * size_dum[i]^2
  CYIM_SHASH_kernel[,i] <- dSHASH(x = size_dum, 
             mu=mu_size,
             sigma = exp(coefs["sigma_b0"] + coefs["sigma_b1"]*mu_size + coefs["sigma_b2"]*mu_size^2), 
             nu = exp(coefs["nu_b0"] + coefs["nu_b1"]*mu_size + coefs["nu_b2"]*mu_size^2), 
             tau = exp(coefs["tau_b0"] + coefs["tau_b1"]*mu_size)) 
}


levelplot(CYIM_SHASH_kernel,row.values = size_dum, column.values = size_dum,cuts=30,
          col.regions=rainbow(30),xlab="log Size t",ylab="log Size t+1",main="SHASH",
          panel = function(...) {
            panel.levelplot(...)
            grid.points(log(CYIM$vol_t), log(CYIM$vol_t1), pch = ".",gp = gpar(cex=3,col=alpha("black",0.5)))
          })


# Final fit - diagnostics -------------------------------------------------
# Now finer diagnostics comparing moments and quantiles of the real data against data generated 
# by the fitted models. I am including the Gaussian and SHASH models for comparison. 

# Simulate data from fitted SHASH model
MLmu = U%*%coefs[1:ncol(U)] 

n_sim <- 500
cactus_sim<-cactus_sim_norm<-matrix(NA,nrow=nrow(CYIM),ncol=n_sim)
for(i in 1:n_sim){
  cactus_sim[,i] <- rSHASH(n = nrow(CYIM), 
                    mu = MLmu, 
					   sigma = exp(coefs["sigma_b0"] + coefs["sigma_b1"]*MLmu  + coefs["sigma_b2"]*MLmu^2), 
             nu = exp(coefs["nu_b0"] + coefs["nu_b1"]*MLmu + coefs["nu_b2"]*MLmu^2), 
             tau = exp(coefs["tau_b0"] + coefs["tau_b1"]*MLmu))
  cactus_sim_norm[,i] <- rnorm(n = nrow(CYIM),
                               mean = predict(CYIM_lmer_best),
                               sd = exp(pars[[best_model]][1] + pars[[best_model]][2]*predict(CYIM_lmer_best) + pars[[best_model]][3]*predict(CYIM_lmer_best)^2))
}

n_bins = 10
alpha_scale = 0.7
cactus_moments <- CYIM %>% 
  arrange(log(vol_t)) %>% 
  mutate(size_bin = cut_number(log(vol_t),n=n_bins)) %>% 
  group_by(size_bin) %>% 
  summarise(mean_t1 = mean(log(vol_t1)),
            sd_t1 = sd(log(vol_t1)),
            skew_t1 = NPskewness(log(vol_t1)),
            kurt_t1 = NPkurtosis(log(vol_t1)),
            bin_mean = mean(log(vol_t)),
            bin_n = n()) 

par(mfrow=c(2,2),mar=c(4,4,2,1),cex.axis=1.3,cex.lab=1.3,mgp=c(2,1,0),bty="l"); 
sim_bin_means=sim_moment_means=sim_moment_means_norm = matrix(NA,n_bins,n_sim); 
for(i in 1:n_sim){
    sim_moments <- bind_cols(CYIM,data.frame(sim=cactus_sim[,i],
                                             sim_norm=cactus_sim_norm[,i])) %>% 
    arrange(log(vol_t)) %>% 
    mutate(size_bin = cut_number(log(vol_t),n=n_bins)) %>% 
    group_by(size_bin) %>% 
    summarise(mean_t1 = mean(sim),
              mean_t1_norm = mean(sim_norm),
              bin_mean = mean(log(vol_t)))
	sim_bin_means[,i]=sim_moments$bin_mean; 
	sim_moment_means[,i]=sim_moments$mean_t1; sim_moment_means_norm[,i]=sim_moments$mean_t1_norm;		  
}

matplot(cactus_moments$bin_mean, sim_moment_means,col=alpha("gray",0.5),pch=16,xlab="Mean size t0",ylab="mean(Size t1)",cex=1.4)
matplot(cactus_moments$bin_mean+0.4, sim_moment_means_norm,col=alpha("cornflowerblue",0.5),pch=16,add=T)
points(cactus_moments$bin_mean+0.2, cactus_moments$mean_t1,pch=1,lwd=2,col=alpha("red",alpha_scale),cex=1.4)
points(cactus_moments$bin_mean, apply(sim_moment_means,1,median),pch=1,lwd=2,col=alpha("black",alpha_scale),cex=1.4)
points(cactus_moments$bin_mean+0.4, apply(sim_moment_means_norm,1,median),pch=1,lwd=2,col=alpha("black",alpha_scale),cex=1.4)
legend("topleft",legend=c("SHASH","Gaussian","Median","Data"),
col=c(alpha("gray",0.5),alpha("cornflowerblue",0.5),alpha("black",alpha_scale), alpha("red",alpha_scale)),pch=1,lwd=2,bty="n"); 
add_panel_label("a")

for(i in 1:n_sim){
    sim_moments <- bind_cols(CYIM,data.frame(sim=cactus_sim[,i],
                                             sim_norm=cactus_sim_norm[,i])) %>% 
    arrange(log(vol_t)) %>% 
    mutate(size_bin = cut_number(log(vol_t),n=n_bins)) %>% 
    group_by(size_bin) %>% 
    summarise(mean_t1 = sd(sim),
              mean_t1_norm = sd(sim_norm),
              bin_mean = mean(log(vol_t)))
	sim_bin_means[,i]=sim_moments$bin_mean; 
	sim_moment_means[,i]=sim_moments$mean_t1; sim_moment_means_norm[,i]=sim_moments$mean_t1_norm;		  
}
matplot(cactus_moments$bin_mean, sim_moment_means,col=alpha("gray",0.5),pch=16,xlab="Mean size t0",ylab="SD(Size t1)",cex=1.4) 
matplot(cactus_moments$bin_mean+0.4, sim_moment_means_norm,col=alpha("cornflowerblue",0.5),pch=16,add=T)
points(cactus_moments$bin_mean+0.2, cactus_moments$sd_t1,pch=1,lwd=2,col=alpha("red",alpha_scale),cex=1.4)
points(cactus_moments$bin_mean, apply(sim_moment_means,1,median),pch=1,lwd=2,col=alpha("black",alpha_scale),cex=1.4)
points(cactus_moments$bin_mean+0.4, apply(sim_moment_means_norm,1,median),pch=1,lwd=2,col=alpha("black",alpha_scale),cex=1.4)
add_panel_label("b")

for(i in 1:n_sim){
    sim_moments <- bind_cols(CYIM,data.frame(sim=cactus_sim[,i],
                                             sim_norm=cactus_sim_norm[,i])) %>% 
    arrange(log(vol_t)) %>% 
    mutate(size_bin = cut_number(log(vol_t),n=n_bins)) %>% 
    group_by(size_bin) %>% 
    summarise(mean_t1 = NPskewness(sim),
              mean_t1_norm = NPskewness(sim_norm),
              bin_mean = mean(log(vol_t)))
	sim_bin_means[,i]=sim_moments$bin_mean; 
	sim_moment_means[,i]=sim_moments$mean_t1; sim_moment_means_norm[,i]=sim_moments$mean_t1_norm;	  
}
matplot(cactus_moments$bin_mean, sim_moment_means,col=alpha("gray",0.5),pch=16,xlab="Mean size t0",ylab="Skew(Size t1)",cex=1.4)
matplot(cactus_moments$bin_mean+0.4, sim_moment_means_norm,col=alpha("cornflowerblue",0.5),pch=16,add=T)
points(cactus_moments$bin_mean+0.2, cactus_moments$skew_t1,pch=1,lwd=2,col=alpha("red",alpha_scale),cex=1.4)
points(cactus_moments$bin_mean, apply(sim_moment_means,1,median),pch=1,lwd=2,col=alpha("black",alpha_scale),cex=1.4)
points(cactus_moments$bin_mean+0.4, apply(sim_moment_means_norm,1,median),pch=1,lwd=2,col=alpha("black",alpha_scale),cex=1.4)
add_panel_label("c")

for(i in 1:n_sim){
    sim_moments <- bind_cols(CYIM,data.frame(sim=cactus_sim[,i],
                                             sim_norm=cactus_sim_norm[,i])) %>% 
    arrange(log(vol_t)) %>% 
    mutate(size_bin = cut_number(log(vol_t),n=n_bins)) %>% 
    group_by(size_bin) %>% 
    summarise(mean_t1 = NPkurtosis(sim),
              mean_t1_norm = NPkurtosis(sim_norm),
              bin_mean = mean(log(vol_t)))
	sim_bin_means[,i]=sim_moments$bin_mean; 
	sim_moment_means[,i]=sim_moments$mean_t1; sim_moment_means_norm[,i]=sim_moments$mean_t1_norm;	  
}
matplot(cactus_moments$bin_mean, sim_moment_means,col=alpha("gray",0.5),pch=16,xlab="Mean size t0",ylab="Log kurtosis(Size t1)",cex=1.4)
matplot(cactus_moments$bin_mean+0.4, sim_moment_means_norm,col=alpha("cornflowerblue",0.5),pch=16,add=T)
points(cactus_moments$bin_mean+0.2, cactus_moments$kurt_t1,pch=1,lwd=2,col=alpha("red",alpha_scale),cex=1.4)
points(cactus_moments$bin_mean, apply(sim_moment_means,1,median),pch=1,lwd=2,col=alpha("black",alpha_scale),cex=1.4)
points(cactus_moments$bin_mean+0.4, apply(sim_moment_means_norm,1,median),pch=1,lwd=2,col=alpha("black",alpha_scale),cex=1.4)
add_panel_label("d")

# Steve wrote this nifty function to do something similar but with quantiles. 
# I have not added the Gaussian comparison but I can print this out separately. Here is the skewet $t$. 
quantileComparePlot(log(CYIM$vol_t),log(CYIM$vol_t1),cactus_sim,n_bins)

#And here is the Gaussian. These look _really_ similar, except perhaps a mismatch at the smallest sizes. 
quantileComparePlot(log(CYIM$vol_t),log(CYIM$vol_t1),cactus_sim_norm,n_bins)

# IPM comparisons ---------------------------------------------------------
##Finally, we can compare IPM predictions between the best Gaussian model and the SHASH model. 

# Here are the size-dependent functions for survival, flowering, and flowerbud production, 
# fit in lme4 with year and plot random effects.

surv_mod <- glmer(Survival_t1 ~ log(vol_t) + (1|year_t) + (1|plot), family="binomial", data=CYIM_full)
flow_mod <- glmer(Goodbuds_t1>0 ~ log(vol_t1) + (1|year_t) + (1|plot), family="binomial", data=CYIM_full)
fert_mod <- glmer(Goodbuds_t1 ~ log(vol_t1) + (1|year_t) + (1|plot), family="poisson", data=subset(CYIM_full,Goodbuds_t1>0))

par(mfrow=c(2,2),mar=c(4,4,2,1))
plot(log(CYIM_full$vol_t),CYIM_full$Survival_t1,xlab="log Size_t",ylab="Survival",col=alpha("gray",0.5))
lines(size_dum,invlogit(fixef(surv_mod)[1]+fixef(surv_mod)[2]*size_dum),lwd=3)
plot(log(CYIM_full$vol_t),log(CYIM_full$vol_t1),xlab="log Size_t",ylab="log Size_t+1",col=alpha("gray",0.5))
lines(size_dum,mean(coefs[c(years,plots)]) + coefs["log(vol_t)"] * size_dum,lwd=3)
lines(size_dum,fixef(CYIM_lmer_best)[1]+fixef(CYIM_lmer_best)[2]*size_dum,col="red")
legend("topleft",legend=c("SHASH location","lmer mean"),lwd=c(2,1),col=c("black","red"))
plot(log(CYIM_full$vol_t),CYIM_full$Goodbuds_t1>0,xlab="log Size_t",ylab="Flowering",col=alpha("gray",0.5))
lines(size_dum,invlogit(fixef(flow_mod)[1]+fixef(flow_mod)[2]*size_dum),lwd=3)
plot(log(CYIM_full$vol_t[CYIM_full$Goodbuds_t1>0]),CYIM_full$Goodbuds_t1[CYIM_full$Goodbuds_t1>0],xlab="log Size_t",ylab="Fertility",col=alpha("gray",0.5))
lines(size_dum,exp(fixef(fert_mod)[1]+fixef(fert_mod)[2]*size_dum),lwd=3)

source("cactus_IPM_source_fns.R")

seeds_per_fruit<-read.csv("JO_fruit_data_final_dropplant0.csv",T)  %>% drop_na() %>% 
  summarise(seeds_per_fruit = mean(seed_count))
seed_survival <- read.csv("FruitSurvival.csv",T) %>% drop_na() %>% mutate(fruit_frac = Fr.on.grnd.not.chewed/Fr.on.plant) %>% 
  summarise(seed_survival = mean(fruit_frac))
germination <- read.csv("Germination.csv") %>% drop_na() %>% 
  mutate(germ1_frac = Seedlings04/Input,
         germ2_frac = Seedlings05/(Input-Seedlings04)) %>% 
  summarise(germ1 = mean(germ1_frac), germ2 = mean(germ2_frac))
precensus_survival <- read.csv("PrecensusSurvival.csv") %>% select(survive0405) %>% drop_na() %>% 
  summarise(precensus_survival = mean(survive0405))
seedling_size <- read_csv("cholla_demography_20042018_EDI.csv") %>% 
  ## subset seed germination plots
  filter(str_sub(Plot,1,1)=="H") %>% 
  mutate(vol_t = volume(Height_t,Width_t,Perp_t)) %>% 
  summarise(mean_size = mean(log(vol_t),na.rm=T),
            sd_size = sd(log(vol_t),na.rm=T))

cactus_params<-c()
## skewed t growth params
## intecept is the mean across years and plots, corresponds to lme4 coefs for other vital rates
cactus_params$grow.mu <- mean(coefs[c(years,plots)]) 
cactus_params$grow.bsize <- coefs["log(vol_t)"]
cactus_params$grow.bsize2 <- coefs["I(log(vol_t)^2)"]
cactus_params$sigma_b0 <- coefs["sigma_b0"]
cactus_params$sigma_b1 <- coefs["sigma_b1"]
cactus_params$sigma_b2 <- coefs["sigma_b2"]
cactus_params$nu_b0 <- coefs["nu_b0"]
cactus_params$nu_b1 <- coefs["nu_b1"]
cactus_params$nu_b2 <- coefs["nu_b2"]
cactus_params$tau_b0 <- coefs["tau_b0"]
cactus_params$tau_b1 <- coefs["tau_b1"]
## Gaussian growth from best lme4 model
cactus_params$grow.mu.norm <- fixef(CYIM_lmer_best)[1]
cactus_params$grow.bsize.norm <- fixef(CYIM_lmer_best)[2]
cactus_params$grow.bsize2.norm <- fixef(CYIM_lmer_best)[3]
cactus_params$grow.sd.b0 <- pars[[best_model]][1]
cactus_params$grow.sd.b1 <- pars[[best_model]][2]
cactus_params$grow.sd.b2 <- pars[[best_model]][3]
## survival
cactus_params$surv.mu <- fixef(surv_mod)[1]
cactus_params$surv.bsize <- fixef(surv_mod)[2]
## flowering
cactus_params$flow.mu <- fixef(flow_mod)[1]
cactus_params$flow.bsize <- fixef(flow_mod)[2]
## fruit production
cactus_params$fert.mu <- fixef(fert_mod)[1]
cactus_params$fert.bsize <- fixef(fert_mod)[2]
## seeds per fruit
cactus_params$mu_spf <- seeds_per_fruit$seeds_per_fruit
## seed survival
cactus_params$seedsurv <- seed_survival$seed_survival
## germination rates
cactus_params$germ1 <- germination$germ1
cactus_params$germ2 <- germination$germ2
## precensus survival
cactus_params$precenus_surv <- precensus_survival$precensus_survival
## seedling size
cactus_params$mu_sdlgsize <- seedling_size$mean_size
cactus_params$sigma_sdlgsize <- seedling_size$sd_size
## size bounds
cactus_params$min.size <- log(min(CYIM$vol_t)) 
cactus_params$max.size <- log(max(CYIM$vol_t1))

mat.size = 200
lower.extension = -1
upper.extension = 1.5
kernel_SHASH <- bigmatrix(params = cactus_params,
          lower.extension = lower.extension, 
          upper.extension = upper.extension,
          mat.size = mat.size,
          dist="SHASH")
kernel_norm <- bigmatrix(params = cactus_params,
          lower.extension = lower.extension, 
          upper.extension = upper.extension,
          mat.size = mat.size,
          dist="norm")


# And finally for IPM results. The mean model (averaging across years and plots) predicts very different growth rates:
 tibble(Growth_Dist = c("SHASH","Gaussian"),
       lambda = c(round(lambda(kernel_SHASH$IPMmat),4),round(lambda(kernel_norm$IPMmat),4)))

#Also, the SSD's are very different and the observed distribution is maybe in between 
# the two of them (in the data plot, colored bars are years). The Gaussian ssd is bi-modal but the 
# lower mode disappears with the skewed $t$ growth kernel, probably because the big reproductive plants 
# are rare in the ssd so you lose the signal of recruitment. 
ssd_SHASH <- stable.stage(kernel_SHASH$IPMmat)[3:(mat.size+2)] / sum(stable.stage(kernel_SHASH$IPMmat)[3:(mat.size+2)])
ssd_norm <- stable.stage(kernel_norm$IPMmat)[3:(mat.size+2)] / sum(stable.stage(kernel_norm$IPMmat)[3:(mat.size+2)])
empirical_sd <- density(log(CYIM$vol_t1),n=mat.size)

plot(kernel_norm$meshpts,ssd_norm,type="l",lty=2,lwd=3,xlab="log volume",ylab="Density")
lines(kernel_SHASH$meshpts,ssd_SHASH,type="l",lwd=3)
lines(empirical_sd$x,empirical_sd$y/sum(empirical_sd$y),col="red",lwd=3)
legend("topleft",c("Gaussian SSD","SHASH SSD", "Empirical SD"),lty=c(2,1,1),col=c("black","black","red"),lwd=3,bty="n")

# The difference in ssd's is pretty striking, so I just wanted to have a closer look at the two 
# growth kernels. The skewed $t$ is "peak-ier" than the Gaussian and this causes the Gaussian 
# to allow for large increases in size that don't happen with the skewed $t$. I think that is 
# where the big difference in SSD comes from. 
x = c(-4,0,4,8,12)
par(mfrow=c(2,3))
for(i in 1:length(x)){
plot(size_dum,gxy_SHASH(x=x[i],y=size_dum,params=cactus_params),type="l",main=paste("Size_t0 = ", x[i]),
     xlab="Future size",ylab="Pr density")
lines(size_dum,gxy_norm(x=x[i],y=size_dum,params=cactus_params),lty=2)
}
legend("topleft",legend=c("SHASH","Gaussian"),lty=c(1,2),bty="n")

