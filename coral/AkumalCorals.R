#############################################################
# The corals case study, based on Bruno et al. (2011)
#
# Original code by Steve Ellner (using Tom's diagnostic plots)
#
# Last modified by Steve Ellner August 20, 2020 
#############################################################

rm(list=ls(all=TRUE); 
setwd("c:/repos/IPM_size_transitions/coral"); 

require(car); require(zoo); require(moments); require(mgcv); 
require(gamlss); require(AICcmodavg); 
require(tidyverse); require(maxLik); 

source("../Diagnostics.R"); 
source("../fitChosenDists.R"); 
source("AkumalCoralsSetup.R"); # load the data frame on healthy corals 

# names(XH); 
#  "Site"  "Fan.number" "Year"  "Area1"  "Area2"  "logarea.t0" "logarea.t1"

#####################################################################
# Re-fit the Bruno et al. (2011) model that used a cube-root transform.  
# We fit all healthy fans, and look for effects of fate, site, year
#####################################################################
nlog=function(x) x^(1/3)
fitH1=lm(nlog(Area2)~nlog(Area1), data=XH); 
fitH2=lm(nlog(Area2)~nlog(Area1)+factor(Site):factor(Year),data=XH); 
fitH3=lm(nlog(Area2)~nlog(Area1)+factor(Site)*factor(Year)+nlog(Area1):factor(Site)+nlog(Area1):factor(Year),data=XH); 

################### Interrogate the residuals: not Gaussian 
jarque.test(fitH1$residuals) # normality test: FAILS, P < 0.001 
anscombe.test(fitH1$residuals) # kurtosis: FAILS, P < 0.01 
agostino.test(fitH1$residuals) # skewness: FAILS, P<0.001 

jarque.test(fitH3$residuals) # normality test: FAILS, P < 0.001 
anscombe.test(fitH3$residuals) # kurtosis: FAILS, P < 0.01 
agostino.test(fitH3$residuals) # skewness: FAILS, P<0.001 

######################################################################### 
# Since cube-root transformation isn't the cure, use log 
# First step is to fit a pilot Gaussian model. 
#########################################################################
fitGAU <- gam(list(logarea.t1~s(logarea.t0),~s(logarea.t0)), data=XH, gamma=1.4,family=gaulss())
summary(fitGAU); plot(fitGAU); 

## the mean looks almost linear; is there evidence against this? 
fitGAU0 <- gam(list(logarea.t1~logarea.t0,~s(logarea.t0)), data=XH, gamma=1.4, family=gaulss())
AIC(fitGAU); AIC(fitGAU0); # yes, Delta AIC of about 9 in favor of the spline 

## the log(sigma) fit looks almost linear; is there evidence against this? 
fitGAU00 <- gam(list(logarea.t1~s(logarea.t0),~logarea.t0), data=XH, gamma=1.4,family=gaulss())
AIC(fitGAU); AIC(fitGAU00);  # Delta AIC < 2, so perhaps weak evidence 

####  Extract values of the fitted splines to get polynomial approximations
####  to use as start values in fitting the final non-Gaussian model 
z_vals = XH$logarea.t0; 
fitted_vals = predict(fitGAU,type="response"); 

##### Mean is fitted almost exactly by linear, exactly by quadratic (spline has df just below 3). 
mean_fit1 = lm(fitted_vals[,1]~z_vals); 
mean_fit2 = lm(fitted_vals[,1]~z_vals+I(z_vals^2)); 

#### log(sigma) is fitted well by linear, perfectly by quadratic (spline has df just above 2) 
sigma_hat = 1/fitted_vals[,2]; 
sd_fit1 = lm(log(sigma_hat)~z_vals); # R^2 = 0.97 
sd_fit2 = lm(log(sigma_hat)~z_vals+I(z_vals^2)); # R^2 = 0.999 

###### Interrogate the scaled residuals: not Gaussian!  
fitted_all = predict(fitGAU,type="response"); 
fitted_mean = fitted_all[,1];
fitted_sd = 1/fitted_all[,2]; 
scaledResids=residuals(fitGAU,type="response")/fitted_sd;  
 
jarque.test(scaledResids) # normality test: FAILS, P =.03 
agostino.test(scaledResids) # skewness: passes, p=0.6 
anscombe.test(scaledResids) # kurtosis: FAILS, P < 0.03 

#####################################################################
#  Graph: data and pilot model 
#####################################################################
graphics.off(); dev.new(width=9,height=7); 
par(mfrow=c(2,2),mar=c(4,4,2,1), bty="l",cex.axis=1.3,cex.lab=1.3,mgp=c(2.2,1,0)); 

plot(Area2~Area1,data=XH,xlab="Initial size",ylab="Subsequent size") 
add_panel_label("a"); 

plot(I(Area2^0.3333)~I(Area1^0.3333),data=XH,xlab="(Initial size)^1/3",ylab="(Subsequent size)^(1/3)") 
add_panel_label("b"); 

e = order(z_vals); 
plot(log(Area2)~log(Area1),data=XH,xlab="log(Initial size)",ylab="log(Subsequent size)",ylim=c(1,8.3)) 
points(z_vals[e],fitted_vals[e,1],type="l",lty=1,col="red",lwd=2); 
points(log(XH$Area1), log(XH$Area2)) 
points(z_vals[e],fitted_vals[e,1]+2/fitted_vals[e,2],type="l",lty=2,col="blue",lwd=2); 
points(z_vals[e],fitted_vals[e,1]-2/fitted_vals[e,2],type="l",lty=2,col="blue",lwd=2); 
add_panel_label("c");

qqPlot(scaledResids,xlab="Normal quantiles",ylab="Scaled residuals"); 
add_panel_label("d");

dev.copy2pdf(file="../manuscript/figures/AkumalPilot.pdf"); 

#####################################################################
# Graph: rolling moments diagnostic using NP skew, kurtosis
#####################################################################
z = rollMomentsNP(XH$logarea.t0,scaledResids,windows=10,smooth=TRUE,scaled=TRUE,xlab="Initial log area") 
dev.copy2pdf(file="../manuscript/figures/AkumalRollingResiduals.pdf");
## mean and SD look good, skew is variable and small except at small sizes, kurtosis on both sides of Gaussian! 

###########################################################################
# Fit suitable distributions to binned data 
###########################################################################
logResids <- data.frame(init=XH$logarea.t0,resids=scaledResids); 
logResids <- logResids %>% mutate(size_bin = cut_number(init,n=8))

source("../fitChosenDists.R"); 

tryDists=c("GT","JSU", "SHASHo","SEP1","SEP3","SEP4"); 

bins = levels(logResids$size_bin); maxVals = matrix(NA,length(bins),length(tryDists)); 
for(j in 1:length(bins)){
for(k in 1:length(tryDists)) {
	Xj=subset(logResids,size_bin==bins[j])
	fitj = gamlssMaxlik(y=Xj$resids,DIST=tryDists[k]); 
	maxVals[j,k] = fitj$out[[1]]$maximum;
	cat("Finished ", tryDists[k]," ",j,k, fitj$maximum,"\n"); 
}
}

## best two for each bin 
for(j in 1:length(bins)){
	e = order(-maxVals[j,]); 
	cat(j, tryDists[e][1:8],"\n"); 
}	

# overall ranking 
e = order(-colSums(maxVals)); 
rbind(tryDists[e],round(colSums(maxVals)[e],digits=3)); 

# And the winner by a hair is: SEP2, just ahead of SEP1 and SEP3. 
# SEP2 has poor inferential properties (DiCiccio TJ and Monti AC 2004) 
# parameters are easily confounded - so we go with SEP1 for the sake of the
# diagnostic. In the paper, we exclude SEP2 for that reason.   

##########################################################################
# Graphical diagnostics to choose the form of the fitted model.
# Estimate the parameters of SEP1 for a set of initial size bins, and plot.
#
# ALERT: parameter estimates returned by our gamlssMaxlik() function are on 
# the family's link-transformed scale, e.g. log(sigma) rather than sigma. 
# Before usine the estimates in a distribution function (e.g., to simulate values 
# or compute a likelihood) the inverse link function needs to be
# applied to the parameter estimates. 
###########################################################################
DIST = "SEP1"; simfun=rSEP1; fam = as.gamlss.family(DIST);

XH <- XH %>% mutate(size_bin = cut_number(logarea.t0,n=8))
bins = levels(XH$size_bin); 
mus = sigmas = nus = taus = bin_means = numeric(length(bins)-1); 
for(j in 1:length(mus)){
	Xj=subset(XH,size_bin%in%c(bins[j-1],bins[j]))
	fitj = gamlssMaxlik(Xj$logarea.t1, DIST=DIST) 
	pars = fitj$out[[1]]$estimate; 
	mus[j]=pars[1]; sigmas[j]=pars[2]; nus[j]=pars[3]; taus[j]=pars[4]; 
	bin_means[j]=mean(Xj$logarea.t0); 
}	

graphics.off(); dev.new(width=8,height=6); 
par(mfrow=c(2,2),bty="l",mar=c(4,4,2,1),mgp=c(2.2,1,0),cex.axis=1.4,cex.lab=1.4);   
spline.scatter.smooth(bin_means,mus,xlab="Initial size",ylab=expression(paste("Location parameter  ", mu ))); #OK
add_panel_label("a"); 
spline.scatter.smooth(bin_means,sigmas,xlab="Initial size",
			ylab=expression(paste("log Scale parameter  ", sigma)));  
add_panel_label("b"); 
spline.scatter.smooth(bin_means,nus,xlab="Initial size",
			ylab=expression(paste("Skewness parameter  ", nu ))); 
add_panel_label("c"); 
spline.scatter.smooth(bin_means,taus,xlab="Initial size",
			ylab=expression(paste("log Kurtosis parameter  ", tau)));  
add_panel_label("d"); 

savePlot(file="../manuscript/figures/AkumalRollingSEP1pars.png", type="png"); 

## Based on the pilot Gaussian fit we let the mu be quadratic 
## and we try letting log sigma be linear (but will consider 
## quadratic if this model does not pass muster). Based on the binned data 
# parameter estimates, we let nu and log(tau) be linear 

# sigma.link = "log", nu.link = "identity", tau.link = "log"

KernelLogLik=function(pars,y,x){
	mu = pars[1]+ pars[2]*x + pars[3]*x^2  
	sigma = exp(pars[4] + pars[5]*x)
	nu = pars[6] + pars[7]*x
	tau = exp(pars[8]+pars[9]*x)
	val = dSEP1(y, mu=mu,sigma=sigma,nu=nu,tau=tau,log=TRUE)
	return(val); 
}

## The start values should be on the link-transformed scale (e.g., log sigma)
## because the inverse-link is applied in KernelLoglik(), so as to 
## guarantee that the parameters are valid in the distribution. 
start=numeric(9);
start[1:3]=coef(mean_fit2); # this is the pilot fit to the mean (link=identity) 
start[4:5]=coef(sd_fit1); # pilot fit to log (sigma)
start[6:7]=c(median(nus),0); # from binned data diagnostic
start[8:9]=c(median(taus),0) # from binned data diagnostic 

fit = maxLik(logLik=KernelLogLik,start=start*exp(0.1*rnorm(8)), y=XH$logarea.t1, x=XH$logarea.t0,method="BHHH",
		control=list(iterlim=5000,printLevel=2),finalHessian=FALSE);  
for(k in 1:5) {		
fit = maxLik(logLik=KernelLogLik,start=fit$estimate, y=XH$logarea.t1, x=XH$logarea.t0,method="NM",
		control=list(iterlim=25000,printLevel=2),finalHessian=FALSE);  
		
fit = maxLik(logLik=KernelLogLik,start=fit$estimate, y=XH$logarea.t1, x=XH$logarea.t0,method="BHHH",
		control=list(iterlim=5000,printLevel=2),finalHessian=FALSE);
}
fit = maxLik(logLik=KernelLogLik,start=fit$estimate, y=XH$logarea.t1, x=XH$logarea.t0,method="BHHH",
		control=list(iterlim=5000,printLevel=2),finalHessian=TRUE);
fit = maxLik(logLik=KernelLogLik,start=fit$estimate, y=XH$logarea.t1, x=XH$logarea.t0,method="NM",
		control=list(iterlim=25000,printLevel=2),finalHessian=TRUE);	
fit = maxLik(logLik=KernelLogLik,start=fit$estimate, y=XH$logarea.t1, x=XH$logarea.t0,method="BHHH",
		control=list(iterlim=5000,printLevel=2),finalHessian=TRUE);	
		
#######################################################################################
#  Plot the fits 
#######################################################################################
graphics.off(); dev.new(width=8,height=6); 
par(mfrow=c(2,2),bty="l",mar=c(4,4,2,1),mgp=c(2.2,1,0),cex.axis=1.4,cex.lab=1.4);   

x = seq(min(XH$logarea.t0),max(XH$logarea.t0),length=100); pars=fit$estimate; 
plot(x,pars[1]+ pars[2]*x+ pars[3]*x^2 ,xlab="Initial size", ylab="Location parameter mu",type="l"); 
plot(x,exp(pars[4] + pars[5]*x), xlab="Initial size", ylab="Scale parameter sigma",type="l"); 
plot(x,pars[6] + pars[7]*x, xlab="Initial size",ylab="Shape parameter nu",type="l");
plot(x,	exp(pars[8]+pars[9]*x), xlab="Initial size",ylab="Shape parameter tau",type="l");
		
		
############## Simulate data from fitted model

fitted_vals = predict(fitGAU,type="response"); 
sigma_hat = 1/fitted_vals[,2];

pars = fit$estimate; x = XH$logarea.t0; 
n_sim <- 500
coral_sim = normal_sim = matrix(NA,nrow=nrow(XH),ncol=n_sim)
for(i in 1:n_sim){
  if(i%%10==0) cat(i,"\n"); 
  coral_sim[,i] <- simfun(n = nrow(XH), 
                    mu =  pars[1]+ pars[2]*x + pars[3]*x^2,  
					sigma = exp(pars[4] + pars[5]*x),
					nu=pars[6] + pars[7]*x,
					tau=exp(pars[8]+pars[9]*x) )
  normal_sim[,i] =  rnorm(n=nrow(XH), mean=fitted_vals[,1], sd=sigma_hat)                 
}

save.image(file="AkumalCoralsModeling.Rdata");


out = quantileComparePlot(sortVariable=XH$logarea.t0,trueData=XH$logarea.t1,simData=coral_sim,
                        nBins=8,alpha_scale = 0.7) 		
                        
dev.copy2pdf(file="../manuscript/figures/CoralQuantileComparePlot.pdf")

out = quantileComparePlot(sortVariable=XH$logarea.t0,trueData=XH$logarea.t1,simData=normal_sim,
                        nBins=8,alpha_scale = 0.7)                      
   

source("../Diagnostics.R"); 
out = momentsComparePlot(sortVariable=XH$logarea.t0,trueData=XH$logarea.t1,simData=coral_sim,normData=normal_sim,
            nBins=10,alpha_scale = 0.7) 	
dev.copy2pdf(file="../manuscript/figures/CoralMomentsComparePlot.pdf")
		

		