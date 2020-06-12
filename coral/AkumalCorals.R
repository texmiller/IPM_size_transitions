setwd("c:/repos/IPM_size_transitions/coral"); 

require(car); require(lme4); require(zoo); require(moments); require(mgcv); 
require(gamlss); require(gamlss.tr); require(AICcmodavg); 
require(lmerTest); require(tidyverse); require(maxLik); 

source("../Diagnostics.R"); 

##########################################################
# Setup the main data frame: do this before any analysis
##########################################################
X=read.csv("e:\\pubs\\mardis\\bruno\\data\\lesion_final_condensed2.csv",as.is=TRUE); 
X$Area.T5=as.numeric(X$Area.T5)
#### NOTE: the l_f_c2 file has correct NA's for state=I and areas not measured. 

# Combine into one list of transitions with Year (=first year) as factor. 
x23=subset(X,select=c(Site,Fan.number,State.T2,Area.T2,Infected.T2,Lesion.T2,State.T3,Area.T3,Infected.T3,Lesion.T3)); x23$Year="2002"; 
x34=subset(X,select=c(Site,Fan.number,State.T3,Area.T3,Infected.T3,Lesion.T3,State.T4,Area.T4,Infected.T4,Lesion.T4)); x34$Year="2003";
x45=subset(X,select=c(Site,Fan.number,State.T4,Area.T4,Infected.T4,Lesion.T4,State.T5,Area.T5,Infected.T5,Lesion.T5)); x45$Year="2004";

newnames=c("Site","Fan.number","State1","Area1","Infected1","Lesion1","State2","Area2","Infected2","Lesion2","Year"); 
names(x23)=newnames; names(x34)=newnames; names(x45)=newnames; 
XC=rbind(x23,x34,x45);
e=!is.na(XC$State1); XC=XC[e,];  

XC$H1=XC$Area1-XC$Infected1; XC$H1[XC$H1<0]=0; 
XC$H2=XC$Area2-XC$Infected2; XC$H2[XC$H2<0]=0; 
XC$State1[XC$H1==0]="D"; 
XC$State2[XC$H2==0]="D";

XH = subset(XC, (State1=="H")&(State2=="H") ); 
XH = subset(XH, select=c(Site,Fan.number,Year,Area1,Area2))
XH = na.omit(XH); 

XH$Site=factor(XH$Site); XH$Year=factor(XH$Year); 

plot(log(Area2)~log(Area1),data=XH);


#####################################################################
# Follow the original model, using cube-root transform  
# Fit all healthy fans, look for effects of fate, site, year
#####################################################################
nlog=function(x) x^(1/3)
fitH1=lm(nlog(Area2)~nlog(Area1), data=XH); 
fitH2=lm(nlog(Area2)~nlog(Area1)+factor(Site):factor(Year),data=XH); 
fitH3=lm(nlog(Area2)~nlog(Area1)+factor(Site)*factor(Year)+nlog(Area1):factor(Site)+nlog(Area1):factor(Year),data=XH); 

################### Interrogate the residuals: not Gaussian 
jarque.test(fitH3$residuals) # normality test: FAILS, P < 0.001 
anscombe.test(fitH3$residuals) # kurtosis: FAILS, P < 0.01 
agostino.test(fitH3$residuals) # skewness: FAILS, P<0.001 

######################################################################### 
# Since 1/3 power transformation isn't the cure, use log 
#########################################################################
XH$logarea.t0 = log(XH$Area1); XH$logarea.t1 = log(XH$Area2); 
fitH2=lm(logarea.t1~logarea.t0 + I(logarea.t0^2), data=XH); # quadratic is significant 
fitH3=lm(logarea.t1~logarea.t0 + I(logarea.t0^2) + I(logarea.t0^3), data=XH); # also cubic
fitH4=lm(logarea.t1~logarea.t0 + I(logarea.t0^2) + I(logarea.t0^3)+ I(logarea.t0^4), data=XH); # nope
# cubic it is 


########################################################################## 
## Use iterative re-weighting to fit with nonconstant variance,  
##########################################################################

## NegLogLik function to fit variance model for residuals 
varPars = function(pars) {
    return(-sum(dnorm(resids, mean=0, sd=exp(pars[1] + pars[2]*fitted_vals + pars[3]*fitted_vals^2),log=TRUE)))
}	

err = 1; rep=0; log_model=lm(logarea.t1 ~ logarea.t0 + I(logarea.t0^2) + I(logarea.t0^3), weights=rep(1,length(XH$logarea.t0)),data=XH); 
while(err > 0.000001) {
	rep=rep+1; 
	fitted_vals = fitted(log_model); resids = residuals(log_model); 
	out=optim(c(log(sd(resids)),0,0),varPars,control=list(maxit=5000)); 
	pars=out$par; 
	new_sigma = exp(pars[1] + pars[2]*fitted_vals + pars[3]*fitted_vals^2); 
	new_weights = 1/((new_sigma)^2)
	new_weights = 0.5*(weights(log_model) + new_weights); # cautious update 
	new_model <- update(log_model,weights=new_weights); 
	err = weights(log_model)-weights(new_model); err=sqrt(mean(err^2)); 
	cat(rep,err,"\n") # check on convergence of estimated weights 
	log_model <- new_model; 
}
## refit and compare, using the estimated weights 
fitH2=lm(logarea.t1~logarea.t0 + I(logarea.t0^2), data=XH,weights=new_weights); 
fitH3=lm(logarea.t1~logarea.t0 + I(logarea.t0^2) + I(logarea.t0^3), data=XH,weights=new_weights); 
fitH4=lm(logarea.t1~logarea.t0 + I(logarea.t0^2) + I(logarea.t0^3)+ I(logarea.t0^4), data=XH,weights=new_weights); 
AIC(fitH2,fitH3,fitH4); # cubic wins (by a hair), so stick with it. 

scaledResids = residuals(log_model)*sqrt(weights(log_model))
fitted_vals = fitted(log_model); 

jarque.test(scaledResids) # normality test: FAILS, P < 0.02 
agostino.test(scaledResids) # skewness: OK, P=0.78 
anscombe.test(scaledResids) # kurtosis: FAILS, P=0.02, kurtosis=3.7 

z = rollMoments(fitted_vals,scaledResids,windows=8,smooth=TRUE,scaled=TRUE) 
## mean and SD look good, skew is variable, kurtosis on both sides of Gaussian! 

