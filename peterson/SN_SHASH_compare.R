#### Use red gorgonian data to compare SN and SHASH fits
#### for modeling growth data

rm(list=ls(all=TRUE))
setwd("c:/repos/IPM_size_transitions/peterson"); 

require(moments); require(dplyr); 
require(zoo); require(sn); require(bbmle); require(gamlss); 

### Load, manipulate, and sort data
dat <- read.csv("Gorgonian_raw_data.csv")
dat.size <- as.data.frame(dat)
colnames(dat.size)[6:8] <- c("Size"="t0", "Sizenext"="t1", "Survnext"= "survival") # rename column header
size <- arrange(dat.size, t0) # arrange initial column in ascending orders
size$t1[which(size$survival==0)] <- NA 
minsize=0 

size <- size[-which(is.na(size$t0)),]
size <- size[-which(is.na(size$t1)),]
size<-size[size$Site=='Portcros',] # For this example, use Portcros site 
plot(size$t0,size$t1); 

#### Diagnostics 

### Scatterplot smoothing functions
spline.scatter.smooth=function(x,y,...) {
	fit=gam(y~s(x),gamma=2,method="REML")
	plot(x,y,type="p",...);
    out=predict(fit,type="response"); 
	points(x,out,type="l",lwd=2)
    fit2 = lm(y~x+I(x^2));
    points(x,fit2$fitted,type="l",col="red",lwd=2,lty=2);
 
}

graphics.off(); 
size$z0 = size$t0; size$z1=size$t1; 
rollmean<-rollapply(size$z0,width=10,mean,na.rm=T)
rollmean1<-rollapply(size$z1,width=10,mean,na.rm=T)
rollsd <-rollapply(size$z1,width=10,sd,na.rm=T); 
rollskew<-rollapply(size$z1,width=10,skewness,na.rm=T)
rollkurt<-rollapply(size$z1,width=10,kurtosis,na.rm=T)/3 - 1
par(mfrow=c(2,2),bty="l",mar=c(4,4,1,1),mgp=c(2,1,0),cex.axis=1.3,cex.lab=1.3);
spline.scatter.smooth(rollmean,rollmean1,col="grey50",xlab="z0",ylab="z1"); 
spline.scatter.smooth(rollmean,rollsd,col="grey50",xlab="z0",ylab="SD"); 
spline.scatter.smooth(rollmean,rollskew,col="grey50",xlab="z0",ylab="Skew"); 
spline.scatter.smooth(rollmean,rollkurt,col="grey50",xlab="z0",ylab="Excess Kurtosis"); 

############################################################### 
# Skewed Normal by ML using bbmle() 
###############################################################

NLLsn <- function(a,b,c,av,bv,cv,aa,ba,ca) {
    mean = a + b*z0 + c*z0^2;
    sd = av + bv*z0 + cv*z0^2
    alpha = aa + ba*z0 + ca*z0^2; 
    -sum(dsn(z1,xi=mean,omega=sd,alpha=alpha,tau=0,log=TRUE))
}    
z0 = size$z0; z1=size$z1; 

fit0=lm(z1~z0+I(z0^2)); sigma0=sqrt(mean(fit0$residuals^2));
c0=as.vector(coef(fit0)); 
fitSN=mle2(NLLsn,start=list(a=c0[1],b=c0[2],c=c0[3],av=sigma0^2,bv=0.001,cv=0.001,aa=0,ba=0,ca=0),
method="Nelder-Mead",skip.hessian=TRUE, control=list(trace=4,maxit=10000)); 
for(k in 1:5) {
    c1 = as.list(coef(fitSN)); 
    fitSN=mle2(NLLsn,start=c1, method="Nelder-Mead",skip.hessian=TRUE, control=list(trace=4,maxit=10000)); 
    c1 = as.list(coef(fitSN)); 
    fitSN=mle2(NLLsn,start=c1, method="BFGS",skip.hessian=TRUE, control=list(trace=4,maxit=10000)); 
}
AIC(fitSN); 

### A qsn function that actually works, unlike the one in the sn package. 
### Vectorized in p, but not the distribution parameters 
my.qsn = function(p,xi,omega,alpha) {
    px = seq(-10,120,length=250); 
    py = psn(px,xi=xi,omega=omega,alpha=alpha,tau=0);
    F = splinefun(py,px,method="monoH.FC"); 
    return(F(p))
}    

c1=coef(fitSN); 
xiFun = function(x) c1[1]+c1[2]*x + c1[3]*x^2;
omegaFun = function(x) c1[4]+c1[5]*x + c1[6]*x^2 
alphaFun = function(x) c1[7]+ c1[8]*x + c1[9]*x^2;
px=seq(min(size$z0),max(size$z0),length=50);  
qSN = matrix(NA,length(px),3); 
for(j in 1:length(px)) {
    qSN[j,]=my.qsn(c(0.025,0.5,0.975),xiFun(px[j]),omegaFun(px[j]),alphaFun(px[j]))

}

########################################################### 
# Now do a SHASH fit using gamlss 
###########################################################
gamSHASH<-gamlss(z1~z0+I(z0^2),sigma.formula=~cs(z0,4), nu.formula=~z0+I(z0^2), 
        tau.formula=~z0+I(z0^2),data=size,family=SHASH,method=RS(250))

dev.new(); 
plot(size$z0,size$z1); 
matpoints(px,qSN,lty=c(2,1,2),col="forestgreen",type="l",lwd=2); 

c.mu = coef(gamSHASH,what="mu"); 
mu.hat = c.mu[1] + c.mu[2]*px + c.mu[3]*px^2; 

sig = predict(gamSHASH,what="sigma",type="response"); 
sigfun = approxfun(size$z0,sig);
sig.hat= sigfun(px); 

c.tau = coef(gamSHASH,what="tau"); 
tau.hat = exp(c.tau[1] + c.tau[2]*px + c.tau[3]*px^2);

c.nu = coef(gamSHASH,what="nu"); 
nu.hat = exp(c.nu[1] + c.nu[2]*px + c.nu[3]*px^2);

py1 = qSHASH(0.025,mu.hat,sig.hat,nu.hat,tau.hat); 
py2 = qSHASH(0.5,mu.hat,sig.hat,nu.hat,tau.hat); 
py3 = qSHASH(0.975,mu.hat,sig.hat,nu.hat,tau.hat); 

matpoints(px,cbind(py1,py2,py3),lty=c(2,1,2),col="blue",type="l",lwd=2); 

zvals=seq(10,80,length=200); 
out = dSHASH(zvals,mu.hat[32],sig.hat[32],nu.hat[32],tau.hat[32]); 
par(yaxs="i")
plot(zvals,out,ylim=c(0,max(out)),type="l"); 

############################################################
### ALL SITES: load, manipulate, and sort data
############################################################
dat <- read.csv("Gorgonian_raw_data.csv")
dat.size <- as.data.frame(dat)
colnames(dat.size)[6:8] <- c("Size"="t0", "Sizenext"="t1", "Survnext"= "survival") # rename column header
size <- arrange(dat.size, t0) # arrange initial column in ascending orders
size$t1[which(size$survival==0)] <- NA 
minsize=0 

size <- size[-which(is.na(size$t0)),]
size <- size[-which(is.na(size$t1)),]
sizeM<-size[size$Site=='Medas',] 
plot(sizeM$t0,sizeM$t1,xlim=range(size$t0),ylim=range(size$t1),col="black",pch=1); 
sizeP<-size[size$Site=='Portcros',] 
points(sizeP$t0,sizeP$t1,col="red",pch=2); 
sizeC<-size[size$Site=='CapdeCreus',] 
points(sizeC$t0,sizeC$t1,col="blue",pch=3); 

size$z0 = size$t0; size$z1=size$t1; 
fit_lm <- lm( z1 ~ z0 + I(z0^2) + Site + z0:Site,data=size);
size$r1 <- fit_lm$residuals; 

### Scatterplot smoothing functions
my.scatter.smooth=function(x,y,...) {
    scatter.smooth(x,y,...); 
    fit = lm(y ~ x + I(x^2) );
    points(x,fit$fitted,type="l",col="red",lwd=2,lty=2);
}

rollmean<-rollapply(size$z0,width=20,mean,na.rm=T)
rollmean1<-rollapply(size$r1,width=20,mean,na.rm=T)
rollsd <-rollapply(size$r1,width=20,sd,na.rm=T); 
rollskew<-rollapply(size$r1,width=20,skewness,na.rm=T)
rollkurt<-rollapply(size$r1,width=20,kurtosis,na.rm=T)/3 - 1
par(mfrow=c(2,2),bty="l",mar=c(4,4,1,1),mgp=c(2,1,0),cex.axis=1.3,cex.lab=1.3);
my.scatter.smooth(rollmean,rollmean1,col="grey50",xlab="z0",ylab="z1",degree=2); 
my.scatter.smooth(rollmean,rollsd,col="grey50",xlab="z0",ylab="SD",degree=2); 
my.scatter.smooth(rollmean,rollskew,col="grey50",xlab="z0",ylab="Skew",degree=2); 
my.scatter.smooth(rollmean,rollkurt,col="grey50",xlab="z0",ylab="Excess Kurtosis",degree=2); 

size$z0 = size$t0; size$z1=size$t1; 
gamSHASH2<-gamlss(z1~z0 + I(z0^2) + z0:Site, sigma.formula=~z0 + I(z0^2), nu.formula=~z0+I(z0^2), 
        tau.formula=~z0+I(z0^2),data=size,family=SHASHo2,method=RS(250)) 

################# Fit to the residuals from a linear fit 
fit <- gam(list(z1 ~ z0 + I(z0^2) + Site:z0, ~s(z0)),data=size,family=gaulss())
size$r1 <- residuals(fit); 
gamSHASH3 = gamlss(r1~z0, sigma.formula=~z0 + I(z0^2), nu.formula=~z0+I(z0^2), 
        tau.formula=~z0+I(z0^2),data=size,family=SHASHo2,method=RS(250)) 

############## Compare distribution estimates, direct with two-stage 
coef(gamSHASH2,what="sigma"); coef(gamSHASH3,what="sigma");
coef(gamSHASH2,what="nu"); coef(gamSHASH3,what="nu");
coef(gamSHASH2,what="tau"); coef(gamSHASH3,what="tau");