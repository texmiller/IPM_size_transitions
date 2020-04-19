# Original: PBA March 2016
# Modified for IPM_size_transitions by SPE, April 2020 
# NOTE: this is ARTR-specific; some shields against bringing in removals
# plots data for other species have been removed from the original code. 
rm(list=ls(all=TRUE));
setwd("c:/repos/IPM_size_transitions/idaho"); 

require(car); require(lme4); require(zoo); require(moments); require(mgcv); 

source("fetchGrowthData.R");
source("trimQuadrats.R"); 

## Steve's diagnostics code
### Scatterplot smoothing functions
spline.scatter.smooth=function(x,y,gamma=10,show.quadratic=FALSE,...) {
  fit=gam(y~s(x),gamma=gamma,method="REML")
  plot(x,y,type="p",...);
  out=predict(fit,type="response"); 
  points(x,out,type="l",lwd=1)
  if(show.quadratic){
    fit2 = lm(y~x+I(x^2));
    points(x,fit2$fitted,type="l",col="red",lwd=2,lty=2);
  }
}

##############################################################
#  1. Import data and calculate W's  - ARTR
##############################################################
sppList <- c("ARTR","HECO","POSE","PSSP","allcov","allpts")
doSpp <- "ARTR"   # !!!!!!!!!!!!!  Don't change this

dataDir1 <- "c:/repos/ExperimentTests/data/idaho"
dataDir2 <- "c:/repos/ExperimentTests/data/idaho_modern"
nonCompLength.s=5 #Number of columns in SppData that are not measures of competitors 

# set up distance weights------------------------------------------------
dists <- read.csv("IdahoModDistanceWeights_noExptl.csv");
dists$allcov <- rowMeans(dists[,1:4])  # for "other" polygons use average of big 4
dists$allpts <- dists$POSE  # set forb dist wts = smallest grass (POSE)

#########################################
# import old data------------------------
#########################################
D1 <- fetchGdat(doSpp=doSpp,speciesList=sppList,datadir=dataDir1,distWts=dists)
D1$Treatment <- "Control"

#########################################
# import modern data---------------------
#########################################
source("fetchGrowthData.R"); 
D2 <- fetchGdat(doSpp=doSpp,speciesList=sppList,datadir=dataDir2,distWts=dists)

ii <- which(D2$year>=2011 & D2$Treatment=="No_grass")
D2$W.HECO[ii] <- 0 ; D2$W.POSE[ii] <- 0 ; D2$W.PSSP[ii] <- 0

# merge in treatment identity for each quadrat 
tmp <- read.csv(paste(dataDir2,"/quad_info.csv",sep=""))
tmp <- tmp[,c("quad","Treatment")]
D2 <- merge(D2,tmp, all.x=T)

# limit to control and removals plots 
keep <- which(is.element(D2$Treatment,c("Control","No_grass")))
D2 <- D2[keep,]; 

#####################################
# merge and clean up 
#####################################
allD = rbind(D1,D2); 
allD$year[allD$year<2000] <- allD$year[allD$year<2000] + 1900

tmp <- which((allD$area.t0>100)&(allD$area.t1 < 0.26)); 
allD <- allD[-tmp,]; 

# get rid of seedlings 
allD <- trimQuadrats(allD)$data;
allD <- subset(allD,age>1); 

plot(logarea.t1~logarea.t0,data=allD); 

allD$Treatment[allD$Treatment=="Control" & allD$year>2000] <- "ControlModern"
allD$year <- as.factor(allD$year)

cols <- c("area.t0","area.t1","logarea.t0","logarea.t1","year","Group","W.ARTR","W.HECO","W.POSE","W.PSSP","W.allcov","W.allpts","Treatment")     
allD <- allD[,cols]; 
e = order(allD$area.t0); allD <- allD[e,]; 

#write.csv(allD,file="ARTR_growth_data.csv"); 

##############################################################
# 2. Diagnostic plots: first log-transform scale fit 
##############################################################
e = order(allD$area.t0); 
rollmean=rollapply(allD$area.t0[e]^0.5,50,mean,by=25);
rollvar=rollapply(allD$logarea.t1[e],50,sd,by=25); max(rollvar)/min(rollvar); 
rollkurt=rollapply(allD$logarea.t1[e],50,kurtosis,by=25);
rollskew=rollapply(allD$logarea.t1[e],50,skewness,by=25);

graphics.off(); dev.new(); 
par(mfcol=c(3,2),mar=c(5,5,2,1),cex.axis=1.3,cex.lab=1.3); 
spline.scatter.smooth(rollmean,rollvar,gamma=2,xlab="",ylab="Std Dev",ylim=c(0,max(rollvar))); 
title(main="Log-transformed size"); 
spline.scatter.smooth(rollmean,rollskew,gamma=2,xlab="",ylab="Skewness"); 
abline(h=0,col="blue",lty=2) 
spline.scatter.smooth(rollmean,rollkurt/3-1,gamma=2,xlab="Sqrt area t0",ylab="Excess kurtosis"); 
abline(h=0,col="blue",lty=2)

#####################################################
##  sqrt-transform scale fit 
#####################################################
e = order(allD$area.t0); 
rollmean=rollapply(allD$area.t0[e]^0.5,50,mean,by=25);
rollvar=rollapply(allD$area.t1[e]^0.5,50,sd,by=25); max(rollvar)/min(rollvar);
rollkurt=rollapply(allD$area.t1[e]^0.5,50,kurtosis,by=25);
rollskew=rollapply(allD$area.t1[e]^0.5,50,skewness,by=25);

spline.scatter.smooth(rollmean,rollvar,xlab="",ylab="Std Dev",ylim=c(0,max(rollvar))); 
title(main="Sqrt-transformed size"); 
spline.scatter.smooth(rollmean,rollskew,xlab="",ylab="Skewness"); 
abline(h=0,col="blue",lty=2) 
spline.scatter.smooth(rollmean,rollkurt/3-1,xlab="Sqrt area t0",ylab="Excess kurtosis"); 
abline(h=0,col="blue",lty=2)

###############################################################################
#  3. Fit models lmer and lm, Gaussian residuals. 
#     Note: size by year term is important (Tredennick et al. 2018)
###############################################################################
require(lme4); 
m0 <- lmer(sqrt(area.t1)~sqrt(area.t0)+W.ARTR + W.HECO + W.POSE + W.PSSP+  W.allcov + W.allpts + Treatment + 
              (1|Group)+(logarea.t0|year),control=lmerControl(optimizer="bobyqa"),data=allD) 
			  
m1 <- lmer(sqrt(area.t1)~sqrt(area.t0) + W.ARTR + W.POSE + W.PSSP + Treatment + 
              (1|Group)+(logarea.t0|year),control=lmerControl(optimizer="bobyqa"),data=allD) 

m2 <- lmer(sqrt(area.t1)~sqrt(area.t0) + Treatment + 
              (1|Group)+(logarea.t0|year),control=lmerControl(optimizer="bobyqa"),data=allD) 
			  
m3 <- lmer(sqrt(area.t1)~sqrt(area.t0) + Treatment + 
              (1|Group)+(1|year),control=lmerControl(optimizer="bobyqa"),data=allD) 			  
			  
m4 <- lmer(sqrt(area.t1)~sqrt(area.t0) + Treatment + 
              (1|Group)+(0+logarea.t0|year),control=lmerControl(optimizer="bobyqa"),data=allD)

m5 <- lmer(sqrt(area.t1)~sqrt(area.t0) + Treatment + 
              (0+logarea.t0|year),control=lmerControl(optimizer="bobyqa"),data=allD)

m6 <- lm(sqrt(area.t1) ~ sqrt(area.t0) + Treatment + logarea.t0:year,data=allD)
 			  
# yrVals=c(0,as.numeric(coef(m6)[-(1:4)])); sd(yrVals); 

#########################################################################
# 4. Residual diagnostics on scaled residuals: are they Gaussian?   
#########################################################################

spline.scatter.smooth(sqrt(allD$area.t0),abs(residuals(m5))); 
fit = lm(abs(residuals(m5))~sqrt(area.t0),data=allD); 
abline(fit,col="blue"); 

scaledResids = residuals(m5)/fit$fitted; 
qqPlot(scaledResids); # VERY bad in both tails

require(moments); 
jarque.test(scaledResids); # normality test: FAILS, P < 0.001 
anscombe.test(scaledResids); # kurtosis: FAILS, P < 0.001 
agostino.test(scaledResids); # skewness: FAILS, P<0.001 

e = order(allD$area.t0); 
rollmean=rollapply(sqrt(allD$area.t0[e]),50,mean,by=25);
rollkurt=rollapply(scaledResids[e],50,kurtosis,by=25);
rollskew=rollapply(scaledResids[e],50,skewness,by=25);
par(mfrow=c(2,1)); 
spline.scatter.smooth(rollmean,rollskew,gamma=2); 
spline.scatter.smooth(rollmean,rollkurt/3-1,gamma=2); 

######################################################################
# Based on diagnostics, try fitting a zero-truncated SHASH with
# constant skew and kurtosis parameters. The optimizers fail. 
######################################################################
require(gamlss); require(gamlss.tr); 

## truncated distribution; 
trun.SHASH = trun(0,family="SHASHo2",par=0,type="left",local=FALSE); 
 
# pilot fit on scaled residuals to get starting values for the shape parameters 
NegLogLik1 = function(par) {
	sigma=par[1]; nu=par[2]; tau=par[3]; 
	out=dSHASHo2(scaledResids,mu=0,sigma,nu,tau,log=TRUE);
	return(-sum(out))
}	

out=optim(par=c(1,1,1),fn=NegLogLik1,control=list(trace=4,maxit=5000));
sigma.h=out$par[1]; nu.h=out$par[2]; tau.h=out$par[3]
 
fit=gamlss(sqrt(area.t1) ~ sqrt(area.t0) + Treatment + logarea.t0:year,data=allD,family=trun.SHASH,
sigma.formula = ~sqrt(area.t0), nu.formula = ~1, tau.formula = ~1,nu.start=nu.h,tau.start=tau.h,method=RS(250));


# Re-fit: the optimizer takes one step, finds a better fit, and then crashes. 
fit2=gamlss(sqrt(area.t1) ~ sqrt(area.t0) + Treatment + logarea.t0:year,data=allD,family=trun.SHASH,
sigma.formula = ~sqrt(area.t0), nu.formula = ~1, tau.formula = ~1, start.from=fit,method=CG(250)); 

######################################################################
# BASEMENT: left over from a fit to only the modern Idaho data 
######################################################################
 #LL = function(b0,b1.area,b1.ARTR,b1.POSE,b1.Treatment,b0.sigma,b1.sigma,b0.nu,b1.nu) {
 LL=function(par) {
		b0=par[1]; b1.area=par[2]; b1.ARTR=par[3];b1.POSE=par[4];
		b1.Treatment=par[5]; b0.sigma=par[6]; b1.sigma=par[7];
		b0.nu=par[8];b1.nu=par[9]; 
		yhat = b0 + b1.area*sqrt(allD$area.t0)+ b1.ARTR*allD$W.ARTR+ b1.POSE*allD$W.POSE + b1.Treatment*allD$removal; 
		sigma.hat = exp(b0.sigma + b1.sigma*sqrt(allD$area.t0));
		nu.hat = exp(b0.nu + b1.nu*sqrt(allD$area.t0)); 
		out=trun.dT2(sqrt(allD$area.t1),mu=yhat,sigma=sigma.hat,nu=nu.hat,log=TRUE);
		return(sum(out))
}

m1 <- lm(sqrt(area.t1) ~ sqrt(area.t0) + W.ARTR + W.POSE + Treatment,data=allD)

start=c(b0=as.numeric(coef(m1)[1]),
	   b1.area=as.numeric(coef(m1)[2]),
	   b1.ARTR=0,b1.POSE=0,b1.Treatment=0,
	   b0.sigma=log(sigma.h),b1.sigma=0,
	   b0.nu=log(nu.h),b1.nu=0);

out2 = maxLik(logLik=LL, start=start, method="NM",control=list(iterlim=5000,printLevel=4));
for(j in 1:10){
out2 = maxLik(logLik=LL, start=out2$estimate, method="NM",control=list(iterlim=5000,printLevel=4));
}
 
 