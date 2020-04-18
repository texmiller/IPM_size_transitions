# Original: PBA March 2016
# Modified for IPM_size_transitions by SPE, April 2020 
# NOTE: this is ARTR-specific; some shields against bringing in removals
# plots data for other species have been removed. 
rm(list=ls(all=TRUE));
setwd("c:/repos/IPM_size_transitions/idaho"); 

source("fetchGrowthData.R");
source("trimQuadrats.R"); 

##############################################################
#  1. Import data and calculate W's  - ARTR, Idaho modern ONLY 
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

# import modern data--------------------------------------------------------
source("fetchGrowthData.R"); 
D2 <- fetchGdat(doSpp=doSpp,speciesList=sppList,datadir=dataDir2,distWts=dists)

# merge in treatment identity for each quadrat 
tmp <- read.csv(paste(dataDir2,"/quad_info.csv",sep=""))
tmp <- tmp[,c("quad","Treatment")]
allD <- merge(D2,tmp, all.x=T)

# remove outliers (large plants that obviously do not turn into tiny plants)
tmp <- which((allD$logarea.t0>2)&(allD$logarea.t1< -1.25)); 
allD <- allD[-tmp,]; 

# get rid of seedlings 
allD <- trimQuadrats(allD)$data;
allD <- subset(allD,age>1); 

# limit to control and removals plots 
keep <- which(is.element(allD$Treatment,c("Control","No_grass")))
allD <- allD[keep,]; 

cols <- c("area.t0","area.t1","logarea.t0","logarea.t1","year","Group","W.ARTR","W.HECO","W.POSE","W.PSSP","W.allcov","W.allpts","Treatment")     
allD <- allD[,cols]; 
e = order(allD$area.t0); allD <- allD[e,]; 

write.csv(allD,file="ARTR_modern.csv"); 

plot(logarea.t1~logarea.t0,data=allD); 

#########################################
#  2. Fit models
#########################################

m1 <- lm(sqrt(area.t1) ~ sqrt(area.t0) + Treatment,data=allD)

plot(sqrt(allD$area.t0),abs(m1$residuals)); 
fit = lm(abs(m1$residuals)~sqrt(area.t0),data=allD); 
abline(fit); 

scaledResids = m1$residuals/fit$fitted; 
qqPlot(scaledResids); # bad in both tails

require(moments); 
jarque.test(scaledResids); # normality test: FAILS, P < 0.001 
anscombe.test(scaledResids); # kurtosis: FAILS, P < 0.001 
agostino.test(scaledResids); # skewness: OK, P>0.3 

require(zoo);
e = order(m1$fitted); 
rollmean=rollapply(m1$fitted[e],30,mean,by=10);
rollkurt=rollapply(scaledResids[e],30,kurtosis,by=10);
scatter.smooth(rollmean,rollkurt); 
 
 ## truncated distribution; 
 trun.dT2 = trun.d(0,family="TF",par=0,type="left"); 
 
NegLogLik1 = function(par) {
	sigma=par[1]; nu=par[2];
	out=dTF(scaledResids,mu=0,sigma,nu,log=TRUE);
	return(-sum(out))
}	

out=optim(par=c(1,1),fn=NegLogLik1,control=list(trace=4,maxit=5000));
sigma.h=out$par[1]; nu.h=out$par[2];  
 
 allD$removal = 0*allD$area.t0;
 allD$removal[allD$Treatment=="No_grass"]<-1; 
 
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
 
 