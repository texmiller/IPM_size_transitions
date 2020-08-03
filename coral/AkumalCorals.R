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
# Following the original model, use a cube-root transform.   
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
fitGAU <- gam(list(logarea.t1~s(logarea.t0),~s(logarea.t0)), data=XH,family=gaulss())
summary(fitGAU); 

####  Extract values of the fitted splines to explore their properties 
z_vals = seq(min(XH$logarea.t0),max(XH$logarea.t0),length=250); 
fitted_vals = predict(fitGAU,type="response",newdata=data.frame(logarea.t0=z_vals)); 

##### Mean is fitted very well by a quadratic, slightly better by cubic (spline has df just above 3). 
mean_fit1 = lm(fitted_vals[,1]~z_vals); 
mean_fit2 = lm(fitted_vals[,1]~z_vals+I(z_vals^2)); 
mean_fit3 = lm(fitted_vals[,1]~z_vals+I(z_vals^2) + I(z_vals^3));  

#### log(sigma) is fitted well by linear, perfectly by quadratic (spline has df just above 2) 
sigma_hat = 1/fitted_vals[,2]; 
sd_fit1 = lm(log(sigma_hat)~z_vals); # R^2 = 0.97 
sd_fit2 = lm(log(sigma_hat)~z_vals+I(z_vals^2)); # R^2 = 0.999 

###### Interrogate the scaled residuals: not Gaussian!  
fitted_all = predict(fitGAU,type="response"); 
fitted_mean = fitted_all[,1];
fitted_sd = 1/fitted_all[,2]; 
scaledResids=residuals(fitGAU,type="response")/fitted_sd;  
 
jarque.test(scaledResids) # normality test: FAILS, P < 0.001 
agostino.test(scaledResids) # skewness: marginal, p=0.044 
anscombe.test(scaledResids) # kurtosis: FAILS, P < 0.01 

#####################################################################
#  Graph: data and pilot model 
#####################################################################
graphics.off(); dev.new(width=9,height=7); 
par(mfrow=c(2,2),mar=c(4,4,2,1), bty="l",cex.axis=1.3,cex.lab=1.3,mgp=c(2.2,1,0)); 

plot(Area2~Area1,data=XH,xlab="Initial size",ylab="Subsequent size") 
add_panel_label("a"); 

plot(I(Area2^0.3333)~I(Area1^0.3333),data=XH,xlab="(Initial size)^1/3",ylab="(Subsequent size)^(1/3)") 
add_panel_label("b"); 

plot(log(Area2)~log(Area1),data=XH,xlab="log(Initial size)",ylab="log(Subsequent size)",ylim=c(1,8.3)) 
points(z_vals,fitted_vals[,1],type="l",lty=3,col="blue",lwd=3); 
points(z_vals,fitted_vals[,1]+2/fitted_vals[,2],type="l",lty=2,col="blue",lwd=2); 
points(z_vals,fitted_vals[,1]-2/fitted_vals[,2],type="l",lty=2,col="blue",lwd=2); 
add_panel_label("c");

qqPlot(scaledResids,xlab="Normal quantiles",ylab="Scaled residuals"); 
add_panel_label("d");

# dev.copy2pdf(file="../manuscript/figures/AkumalPilot.pdf"); 

#####################################################################
# Graph: rolling moments diagnostic using NP skew, kurtosis
#####################################################################
z = rollMomentsNP(XH$logarea.t0,scaledResids,windows=10,smooth=TRUE,scaled=TRUE,xlab="Initial log area") 
# dev.copy2pdf(file="../manuscript/figures/AkumalRollingResiduals.pdf");
## mean and SD look good, skew is variable and small except at small sizes, kurtosis on both sides of Gaussian! 

###########################################################################
# Fit suitable distributions to binned data 
###########################################################################
logResids <- data.frame(init=XH$logarea.t0,resids=scaledResids); 
logResids <- logResids %>% mutate(size_bin = cut_number(init,n=5))

source("../fitChosenDists.R"); 

tryDists=c("EGB2","GT","JSU", "SHASHo","SEP1","SEP2","SEP3","SEP4"); 

bins = levels(logResids$size_bin); maxVals = matrix(NA,length(bins),length(tryDists)); 
for(j in 1:length(bins)){
for(k in 1:length(tryDists)) {
	Xj=subset(logResids,size_bin==bins[j])
	fitj = gamlssMaxlik(y=Xj$resids,DIST=tryDists[k]); 
	maxVals[j,k] = fitj$maximum;
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

# And the winner by a hair is: SEP2, just ahead of SEP1
# SEP2 has poor inferential properties (DiCiccio TJ and Monti AC 2004) 
# parameters are easily confounded - so we go with SEP1 for the sake of the
# diagnostic. In the paper, maybe exclude SEP2 for that reason?  

##########################################################################
# Graphical diagnostics to choose the form of the fitted model.
# Estimate the parameters of SEP1 for a set of initial size bins, and plot.
# Values returned by gamlssMaxlik are on the family's link-transformed scale, 
# e.g. log(sigma) rather than sigma. 
##########################################################################
DIST = "SHASHo"; density=dSHASHo; simfun=rSHASHo; 

XH <- XH %>% mutate(size_bin = cut_number(logarea.t0,n=8))
bins = levels(XH$size_bin); 
mus = sigmas = nus = taus = bin_means = numeric(length(bins)-1); 
for(j in 1:length(mus)){
	Xj=subset(XH,size_bin%in%c(bins[j-1],bins[j]))
	fitj = gamlssMaxlik(Xj$logarea.t1, DIST=DIST) 
	pars = fitj$estimate; 
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

## Based on the pilot Gaussian fit, we let the mean be cubic,
## even though we don't see it in the binned data diagnostics 
## Based on the binned data parameter estimates, we let 
## log(sigma), nu and log(tau) be linear 

# sigma.link = "log", nu.link = "identity", tau.link = "log"

KernelLogLik=function(pars,y,x){
	mu = pars[1]+ pars[2]*x + pars[3]*(x^2) + pars[4]*(x^3);   
	sigma = exp(pars[5] + pars[6]*x)
	nu = pars[7] + pars[8]*x
	tau = exp(pars[9]+pars[10]*x)
	val = density(y, mu=mu,sigma=sigma,nu=nu,tau=tau,log=TRUE)
	return(val); 
}

## These should be on the link-transformed scale (e.g., log sigma)
## because the inverse-link is applied in KernelLoglik(), so as to 
## guarantee that the parameters are valid in the distribution family. 
start=numeric(10);
start[1:4]=coef(mean_fit3); # this is the pilot fit to the mean (link=identity) 
start[5:6]=coef(sd_fit1); # pilot fit to log (sigma)
start[7:8]=c(median(nus),0); # from binned data diagnostic
start[9:10]=c(median(taus),0) # from binned data diagnostic 

fit = maxLik(logLik=KernelLogLik,start=start*exp(0.1*rnorm(10)), y=XH$logarea.t1, x=XH$logarea.t0,method="BHHH",
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
		
#######################################################################################
#  Plot the fits 
#######################################################################################
graphics.off(); dev.new(width=8,height=6); 
par(mfrow=c(2,2),bty="l",mar=c(4,4,2,1),mgp=c(2.2,1,0),cex.axis=1.4,cex.lab=1.4);   

x = seq(min(XH$logarea.t0),max(XH$logarea.t0),length=100); pars=fit$estimate; 
plot(x,pars[1]+ pars[2]*x + pars[3]*(x^2) + pars[4]*(x^3),xlab="Initial size", ylab="Location parameter mu",type="l"); 
plot(x,exp(pars[5] + pars[6]*x), xlab="Initial size", ylab="Scale parameter sigma",type="l"); 
plot(x,pars[7] + pars[8]*x, xlab="Initial size",ylab="Shape parameter nu",type="l");
plot(x,	exp(pars[9]+pars[10]*x), xlab="Initial size",ylab="Shape parameter tau",type="l");
		
		
############## Simulate data from fitted model
pars = fit$estimate; x = XH$logarea.t0; 
n_sim <- 500
coral_sim<-matrix(NA,nrow=nrow(XH),ncol=n_sim)
for(i in 1:n_sim){
  if(i%%10==0) cat(i,"\n"); 
  coral_sim[,i] <- simfun(n = nrow(XH), 
                    mu =  pars[1]+ pars[2]*x + pars[3]*(x^2) + pars[4]*(x^3),
					sigma = exp(pars[5] + pars[6]*x),
					nu=pars[7] + pars[8]*x,
					tau=exp(pars[9]+pars[10]*x) )
}

out = quantileComparePlot(sortVariable=XH$logarea.t0,trueData=XH$logarea.t1,simData=coral_sim,nBins=10,alpha_scale = 0.7) 		

# dev.copy2pdf(file="../manuscript/figures/CoralQuantileComparePlot.pdf")
		
		
		
		
		
		