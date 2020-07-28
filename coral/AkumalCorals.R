setwd("c:/repos/IPM_size_transitions/coral"); 

require(car); require(lme4); require(zoo); require(moments); require(mgcv); 
require(gamlss); require(gamlss.tr); require(AICcmodavg); 
require(lmerTest); require(tidyverse); require(maxLik); 

source("../Diagnostics.R"); 

#################################################################
# Recruit size data
#################################################################
recruitSizes=c(16.7, 6.4, 10, 4.1, 4.7, 11.2, 6.3, 9.7, 5); #these are areas of known recruits  

#xmin=min(log(recruitSizes)); xmax=max(log(recruitSizes));
#mins=maxs=numeric(500); 
#for(k in 1:500) {
#	z=rnorm(length(recruitSizes),mean=mean(log(recruitSizes)),sd=sd(log(recruitSizes)));
#	mins[k]=min(z); maxs[k]=max(z); 
#}	
# hist(mins); abline(v=xmin); hist(maxs); abline(v=xmax); 

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
XH$logarea.t0 = log(XH$Area1); 
XH$logarea.t1 = log(XH$Area2); 

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
#########################################################################

fitGAU <- gam(list(logarea.t1~s(logarea.t0),~s(logarea.t0)), data=XH,family=gaulss())
summary(fitGAU); 

####  get values of the fitted splines to explore their properties 
z_vals = seq(min(XH$logarea.t0),max(XH$logarea.t0),length=250); 
fitted_vals = predict(fitGAU,type="response",newdata=data.frame(logarea.t0=z_vals)); 

##### mean is fitted very well by a quadratic, slightly better by cubic (spline has df just above 3). 
mean_fit1 = lm(fitted_vals[,1]~z_vals); 
mean_fit2 = lm(fitted_vals[,1]~z_vals+I(z_vals^2)); 
mean_fit3 = lm(fitted_vals[,1]~z_vals+I(z_vals^2) + I(z_vals^3));  

#### log(sigma) is fitted well be a quadratic (spline has df just above 2) 
sigma_hat = 1/fitted_vals[,2]; 
sd_fit1 = lm(log(sigma_hat)~z_vals); # R^2 = 0.97 
sd_fit2 = lm(log(sigma_hat)~z_vals+I(z_vals^2)); # R^2 = 0.999 

###### Interrogate the scaled residuals 
fitted_all = predict(fitGAU,type="response"); 
fitted_mean = fitted_all[,1];
fitted_sd = 1/fitted_all[,2]; 
scaledResids=residuals(fitGAU,type="response")/fitted_sd;  
 
jarque.test(scaledResids) # normality test: FAILS, P < 0.001 
agostino.test(scaledResids) # skewness: marginal, p=0.044 
anscombe.test(scaledResids) # kurtosis: FAILS, P < 0.01 

#####################################################################
#  Data display figure and pilot model 
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

z = rollMomentsNP(XH$logarea.t0,scaledResids,windows=10,smooth=TRUE,scaled=TRUE,xlab="Initial log area") 
# dev.copy2pdf(file="../manuscript/figures/AkumalRollingResiduals.pdf");
## mean and SD look good, skew is variable and small except at small sizes, kurtosis on both sides of Gaussian! 


###########################################################################
# Try fitDist on binned data
###########################################################################
logResids <- data.frame(init=XH$logarea.t0,resids=scaledResids); 
logResids <- logResids %>% mutate(size_bin = cut_number(init,n=5))

bins = levels(logResids$size_bin); dists=list(length(bins)); 
for(j in 1:length(bins)){
	Xj=subset(logResids,size_bin==bins[j])
	dists[[j]]=fitDist(resids,data=Xj,type="realline"); 
	cat(j,"\n"); 
	print(dists[[j]]$fits); 
	cat("   ","\n"); 
}



