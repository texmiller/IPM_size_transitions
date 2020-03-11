#### Appendix 1: This script uses red gorgonian data to compare the beta vs. normal approach
#### to modeling growth data

rm(list=ls(all=TRUE))
require(quantreg)
require(betareg)
require(moments)
require(zoo)
require(dplyr)
require(MuMIn)
require(sn); 

setwd("e:\\admin\\reviews\\BES"); 

### Functions related to beta approach

# This function transforms continuous size data to a (0,1) interval
# based on size-dependent min and max values
betaFn<-function(x,min,max){
  y=(x-min)/(max-min)
}

# This function back-transforms data from a (0,1) interval
# to the original data scale based on size-dependent min and max values
backbeta<-function(x,min,max){
  (x*(max-min))+min
}

# This function back-transforms variance from a (0,1) interval
# to the original data scale based on size-dependent min and max values
backvar<-function(sigma,min,max){
  sigma*((max-min)^2)
}

# This function calculates the alpha and beta shape parameters of a 
# beta distribution based on the mean and variance
shapebeta<-function(mu,sigma){
  alpha<-((1-mu)/sigma-1/mu)*mu^2
  beta<-alpha*(1/mu-1)
  return(params=list(alpha=alpha,beta=beta))
}

# This function calculates the mean and variance of a beta distribution
# based on the alpha and beta shape parameters
revshapebeta<-function(alpha,beta){
  mu<-alpha/(alpha+beta)
  sigma<-(alpha*beta)/(((alpha+beta)^2)*(alpha+beta+1))
  return(params=list(mu=mu,sigma=sigma))
}

# This function calculates the skewness of a beta distribution
# based on the alpha and beta shape parameters
skewbeta<-function(alpha,beta){
  (2*(beta-alpha)*sqrt(alpha+beta+1))/((alpha+beta+2)*sqrt(alpha*beta))
}

### Load, manipulate, and sort data
dat <- read.csv("Gorgonian_raw_data.csv")
dat.size <- as.data.frame(dat)
colnames(dat.size)[6:8] <- c("Size"="t0", "Sizenext"="t1", "Survnext"= "survival") # rename column header
size <- arrange(dat.size, t0) # arrange initial column in ascending orders
size$t1[which(size$survival==0)] <- NA 
minsize=0 

#maxsize=75 # this is set manually to not have a large category at end that has no plants.
#size$t0[which(size$t0>maxsize)]=maxsize-0.1
#size$t1[which(size$t1>maxsize)]=maxsize-0.1

size <- size[-which(is.na(size$t0)),]
size <- size[-which(is.na(size$t1)),]
size<-size[size$Site=='Portcros',] # For this example, use Portcros site 

# Step 1: Look for departures from normality for a given starting size
par(mfrow=c(2,2)); 
plot(size$t0,size$t1); plot(size$t0^0.25, size$t1^0.25); 
rollmean<-rollapply(size$t0^0.25,width=30,mean,na.rm=T)
rollskew<-rollapply(size$t1^0.25,width=30,skewness,na.rm=T)
rollkurt<-rollapply(size$t1^0.25,width=30,kurtosis,na.rm=T)
plot(rollmean,rollskew,xlab='Size at time t',ylab='Skewness in size t+1',type='l'); 
abline(h=0,col='grey') # skew shifts from positive to negative with size
plot(rollmean,rollkurt/3-1,xlab='Size at time t',ylab='Excess Kurtosis in size t+1',type='l'); 
abline(h=0,col='grey') 

#### Start with the beta approach 

# Step 2: Fit size-dependent minimum and maximum values
quant1 <- rq(t1~ t0, data=size, tau=c(0.005,0.995)) # fit quantile regression
quant2 <- rq(t1~ t0 + I(t0^2), data=size,tau=c(0.005,0.995)) # fit quantile regression
quant <- rq(t1~ t0 + I(t0^2) + I(t0^3),data=size,tau=c(0.005,0.995)) # fit quantile regression
summary(quant)
minsizes=pmax(0,predict(quant,size)[,1]) # size dependent minimum
maxsizes=predict(quant,size)[,2] # size dependent maximum

par(bty="l",cex.axis=1.3,cex.lab=1.3,mar=c(5,5,1,1)); 
plot(size$t0,size$t1,xlab='Size at time t',ylab='Size at time t+1',col='grey20')
lines(0:100,pmax(predict(quant,data.frame(t0=0:100))[,1],0),col='black',lwd=2)
lines(0:100,predict(quant,data.frame(t0=0:100))[,2],col='black',lwd=2)
dev.copy2pdf(file="Fig2Arev.pdf"); ###############################################################################

# Step 3: Transform size at time t+1 to (0,1) interval
size$t1b<-betaFn(x=size$t1,min=minsizes,max=maxsizes)
range(size$t1b,na.rm=T) # produces some values outside (0,1) interval
size$t1b[size$t1b>0.99]=0.99
size$t1b[size$t1b<0.01]=0.01

plot(size$t0,size$t1b,xlab='Size at time t',ylab='Beta-transformed size t+1',col='grey')

# Step 4: Fit and compare beta regression models
m1<-betareg(t1b~t0+I(t0^2)|t0+I(t0^2),link='loglog',data=size) # size and size2 effects on mean and precision
m2<-betareg(t1b~t0|t0+I(t0^2),link='loglog',data=size) # size effects on mean and precision, size2 effect on precision
m3<-betareg(t1b~t0+I(t0^2)|t0,link='loglog',data=size) # size effects on mean and precision, size2 effect on mean
m4<-betareg(t1b~t0|t0,link='loglog',data=size) # size effects on mean and precision
m5<-betareg(t1b~t0+I(t0^2),link='loglog',data=size) # size and size2 effects on mean, constant precision
m6<-betareg(t1b~t0,link='loglog',data=size) # size effects on mean, constant precision
AICc(m1,m2,m3,m4,m5,m6) # m1 is best supported
bestgrowth=m2  # SPE change, resulting from cubic quantile regression  

plot(size$t0,size$t1b,xlab='Size at time t',ylab='Beta-transformed size t+1',col='grey')
lines(0:100,predict(bestgrowth,data.frame(t0=0:100),type='response'),col='darkred',lty=1,lwd=2)
lines(0:100,predict(bestgrowth,data.frame(t0=0:100),type='quantile',at=c(0.005)),col='darkred',lty=2,lwd=2) # 95th prediction interval
lines(0:100,predict(bestgrowth,data.frame(t0=0:100),type='quantile',at=c(0.995)),col='darkred',lty=2,lwd=2) # 95th prediction interval

# Step 5: Get shape parameters and pdf for a given starting size
newSize<-predict(quant,data.frame(t0=0:100))
newSize[,1]<-newSize[,1]
newSize[,2]<-newSize[,2]
betaMean<-predict(bestgrowth,data.frame(t0=0:100),type='response') # mean on (0,1) scale
betaVar<-predict(bestgrowth,data.frame(t0=0:100),type='variance') # variance on (0,1) scale
params<-shapebeta(mu=betaMean,sigma=betaVar) # alpha and beta parameters for a given starting size
zb<-seq(0,1,length.out=100)

db<-dbeta(zb,shape1=params$alpha[46],shape2=params$beta[46]) # get pdf for starting size = 45
zn<-backbeta(x=zb,min=newSize[46,1],max=newSize[46,2]) # back-transform to original data scale
db<-db/(zn[100]-zn[1]) # back-transform to original data scale
sub<-size[which(size$t0>42&size$t0<48),] # get actual size at t+1 for starting size between 44 and 48

hist(sub$t1,breaks=10,freq=F,xlim=c(newSize[46,1]*0.8,newSize[46,2]*1.1),main='',xlab="Size at time t+1")
lines(zn,db,type='l',lwd=2,col='darkred')

# Step 5b: Get shape parameters and pdf for a given starting size
newSize<-predict(quant,data.frame(t0=0:100))
newSize[,1]<-newSize[,1]
newSize[,2]<-newSize[,2]
betaMean<-predict(bestgrowth,data.frame(t0=0:100),type='response') # mean on (0,1) scale
betaVar<-predict(bestgrowth,data.frame(t0=0:100),type='variance') # variance on (0,1) scale
params<-shapebeta(mu=betaMean,sigma=betaVar) # alpha and beta parameters for a given starting size
zb<-seq(0,1,length.out=100)

db<-dbeta(zb,shape1=params$alpha[60],shape2=params$beta[60]) # get pdf for starting size = 60
zn<-backbeta(x=zb,min=newSize[60,1],max=newSize[60,2]) # back-transform to original data scale
db<-db/(zn[100]-zn[1]) # back-transform to original data scale
sub<-size[which(size$t0>53&size$t0<67),] # get actual size at t+1 for starting size between 53 and 67

hist(sub$t1,breaks=10,freq=F,xlim=c(newSize[60,1]*0.8,newSize[60,2]*1.1),main='',xlab="Size at time t+1")
lines(zn,db,type='l',lwd=2,col='darkred')

dev.copy2pdf(file="FigS2rev.pdf"); ###############################################################################

# Step 6: Back-transform model predictions to original data scale
betaLow<-predict(bestgrowth,data.frame(t0=0:100),type='quantile',at=c(0.005)) # 99th prediction interval on (0,1) scale
betaHi<-predict(bestgrowth,data.frame(t0=0:100),type='quantile',at=c(0.995)) # 99th prediction interval on (0,1) scale

plot(size$t0,size$t1,xlab='Size time t',ylab='Size time t+1',col='grey')
lines(0:100,backbeta(betaMean,min=newSize[,1],max=newSize[,2]),lty=1,lwd=2,col='darkred')
lines(0:100,backbeta(betaLow,min=newSize[,1],max=newSize[,2]),lty=2,lwd=2,col='darkred')
lines(0:100,backbeta(betaHi,min=newSize[,1],max=newSize[,2]),lty=2,lwd=2,col='darkred')

Params=shapebeta(mu=predict(bestgrowth,size,type='response'),sigma=predict(bestgrowth,size,type='variance'))
LikBeta<-dbeta(size$t1b,shape1=Params$alpha,shape2=Params$beta)
LLBeta<-sum(log(LikBeta))-sum(log(maxsizes-minsizes)) # transform to original datascale
KBeta<-length(c(coef(bestgrowth),coef(quant)))
AICBeta<- -2*LLBeta + 2*KBeta

#### Now compare to the Normal approach ################################################################

# Step 2N: Fit and compare linear regressions for mean growth
m1<-lm(t1~t0+I(t0^2),data=size) # size and size2 on mean 
m2<-lm(t1~t0,data=size) # size on mean
m3<-lm(t1~1,data=size) # no size effect
m4<-nls(t1 ~ a + b*(t0^c), data= size, start=list(a= 1, b=1,c=1)) # power function of size on mean

AICc(m1,m2,m3,m4) # m1 is best-supported
bestmean<-m1

plot(size$t0,size$t1,xlab='Size at time t',ylab='Size at time t+1',col='grey')
lines(0:100,predict(bestmean,data.frame(t0=0:100),type='response'),col='darkblue',lty=1,lwd=2)

# Step 3N: Get residuals
size$resid<-size$t1-predict(bestmean,size)
size$resid2<-size$resid^2

plot(size$t0,size$resid,xlab='Size at time t',ylab='Residuals in size t+1',col='grey')

# Step 4N: Fit and compare linear regressions for variance in growth
m1<-lm(resid2~-1+t0+I(t0^2),data=size) # size and size2 on variance 
m2<-lm(resid2~-1+t0,data=size) # size on variance
m3<-lm(resid2~1,data=size) # no size effect
m4<-lm(resid2~t0+I(t0^2),data=size) # size and size2 on variance 
m5<-lm(resid2~t0,data=size) # size on variance

AICc(m1,m2,m3,m4,m5) # m4 is best-supported
bestvar<-m4

plot(size$t0,size$resid2,xlab='Size at time t',ylab='Squared residuals in size t+1',col='grey')
lines(0:100,predict(bestvar,data.frame(t0=0:100),type='response'),col='darkblue',lty=1,lwd=2)

### do the Normal fit by ML ##########################################################
require(bbmle); 
z0 = size$t0; z1=size$t1; 
NLL <- function(a,b,c,av,bv,cv) {
    mean=a+b*z0+c*z0^2;
    sd = sqrt(av + bv*z0 + cv*z0^2)
    -sum(dnorm(z1,mean=mean,sd=sd,log=TRUE))
}    
fit0=lm(z1~z0+I(z0^2)); sigma0=sqrt(mean(fit0$residuals^2));
c0=as.vector(coef(fit0)); 
fitNorm=mle2(NLL,start=list(a=c0[1],b=c0[2],c=c0[3],av=sigma0^2,bv=0.001,cv=0.001),
method="Nelder-Mead",skip.hessian=TRUE,control=list(trace=4,maxit=10000)); 
c1=as.list(coef(fitNorm)); 
fitNorm=mle2(NLL,start=c1,method="Nelder-Mead",skip.hessian=TRUE,control=list(trace=4,maxit=10000)); 
c1=as.list(coef(fitNorm)); 
fitNorm=mle2(NLL,start=c1,method="BFGS",skip.hessian=TRUE,control=list(trace=4,maxit=10000)); 
AIC(fitNorm); 

### Skewed Normal ######################################################################
require(bbmle); 
NLLsn <- function(a,b,c,av,bv,cv,aa,ba,ca) {
    mean=a+b*z0+c*z0^2;
    sd = sqrt(av + bv*z0 + cv*z0^2)
    alpha = aa + ba*z0 + ca*z0^2; 
    -sum(dsn(z1,xi=mean,omega=sd,alpha=alpha,tau=0,log=TRUE))
}    
fit0=lm(z1~z0+I(z0^2)); sigma0=sqrt(mean(fit0$residuals^2));
c0=as.vector(coef(fit0)); 
fitSN=mle2(NLLsn,start=list(a=c0[1],b=c0[2],c=c0[3],av=sigma0^2,bv=0.001,cv=0.001,aa=0,ba=0,ca=0),
method="Nelder-Mead",skip.hessian=TRUE, control=list(trace=4,maxit=10000)); 
for(k in 1:200) {
    c1 = as.list(coef(fitSN)); 
    fitSN=mle2(NLLsn,start=c1, method="Nelder-Mead",skip.hessian=TRUE, control=list(trace=4,maxit=10000)); 
}
AIC(fitSN); 
# confint(fitSN); 

# Step 5N: Get parameters and pdf for a given starting size
normMean<-predict(bestmean,data.frame(t0=0:100)) 
normVar<-predict(bestvar,data.frame(t0=0:100)) 
# normVar[normVar<0.05]<-0.05 # keep variance positive

dn<-dnorm(0:100,mean=normMean[46],sd=sqrt(normVar[46])) # get pdf for starting size = 45
sub<-size[which(size$t0>42&size$t0<48),] # get actual size at t+1 for starting size between 42 and 48

hist(sub$t1,breaks=10,freq=F,xlim=c(newSize[46,1]*0.8,newSize[46,2]*1.1),main='',xlab="Size at time t+1")
lines(0:100,dn,type='l',lwd=2,col='darkblue')

size$t1n<-backbeta(size$t1b,min=minsizes,max=maxsizes) 
# Compare AIC to the exact same data by back-transforming to the original data scale
mn<-predict(bestmean,size)
mv<-predict(bestvar,size)
LikNorm<-dnorm(size$t1n,mean=mn,sd=sqrt(mv))
LLNorm<-sum(log(LikNorm))
KNorm<-length(c(coef(bestmean),coef(bestvar)))
AICNorm<- -2*LLNorm + 2*KNorm
