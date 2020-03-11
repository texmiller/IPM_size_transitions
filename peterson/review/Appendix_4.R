#### Appendix 4: This script uses red gorgonian data to illustrate using the 
#### gamlss package to fit the parameters of the beta distribution as cubic splines

rm(list=ls())

setwd("e:\\admin\\reviews\\BES"); 

require(quantreg)
require(betareg)
require(moments)
require(zoo)
require(dplyr)
require(MuMIn)
require(gamlss)

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

### Load, manipulate, and sort data
dat <- read.csv("Gorgonian raw data.csv")
dat.size <- as.data.frame(dat)
colnames(dat.size)[6:8] <- c("Size"="t0", "Sizenext"="t1", "Survnext"= "survival") # rename column header
size <- arrange(dat.size, t0) # arrange initial column in ascending orders
size$t1[which(size$survival==0)] <- NA 
minsize=0 
maxsize=75 # this is set manually to not have a large category at end that has no plants.
size$t0[which(size$t0>maxsize)]=maxsize-0.1
size$t1[which(size$t1>maxsize)]=maxsize-0.1
size <- size[-which(is.na(size$t0)),]
size <- size[-which(is.na(size$t1)),]
size<-size[size$Site=='Portcros',] # For this example, use Portcros site 


# Step 1: Fit size-dependent minimum and maximum values
quant<-rq(t1~t0+I(t0^2)+I(t0^3),data=size,tau=c(0.005,0.995)) # fit quantile regression
summary(quant)
minsizes=predict(quant,size)[,1] # size dependent minimum
maxsizes=predict(quant,size)[,2] # size dependent maximum

plot(size$t0,size$t1,xlab='Size at time t',ylab='Size at time t+1',col='grey')
lines(0:100,predict(quant,data.frame(t0=0:100))[,1],col='black')
lines(0:100,predict(quant,data.frame(t0=0:100))[,2],col='black')

# Step 2: Transform size at time t+1 to (0,1) interval
size$t1b<-betaFn(x=size$t1,min=minsizes,max=maxsizes)
range(size$t1b,na.rm=T) # produces some values outside (0,1) interval
size$t1b[size$t1b>0.99]=0.99
size$t1b[size$t1b<0.01]=0.01

plot(size$t0,size$t1b,xlab='Size at time t',ylab='Beta-transformed size t+1',col='grey')

# Step 3: Fit beta distribution parameters with cubic splines and compare
# to parametric models fit with betareg
gam3<-gamlss(t1b~cs(t0,3),sigma.formula=~cs(t0,3), data=size, family=BE) # mean and standard deviation fit as cubic splines of size
gam4<-gamlss(t1b~cs(t0,3), data=size, family=BE) # mean fit as cubic spline of size, constant standard deviation

gam1<-gamlss(t1b~pb(t0), sigma.formula=~pb(t0), data=size, family=BE) # mean and standard deviation fit as cubic splines of size
gam2<-gamlss(t1b~pb(t0), data=size, family=BE) # mean fit as cubic spline of size, constant standard deviation
bet1<-betareg(t1b~t0+I(t0^2)|t0+I(t0^2),data=size)
bet2<-betareg(t1b~t0|t0+I(t0^2),data=size)
bet3<-betareg(t1b~t0+I(t0^2)|t0,data=size)
bet4<-betareg(t1b~t0|t0,data=size)

AIC(gam1,gam2,gam3,gam4) # gam3 is best-supported
AIC(bet1,bet2,bet3,bet4) 

# Step 4: Compare model predictions to the data
gam1mean<-predict(gam1,newdata=data.frame(t0=0:100),type='response')
gam1lo<-centiles.pred(gam1,xname='t0',xvalues=0:100, type='centiles',cent=c(5)) # 5th prediction interval on (0,1) scale
gam1hi<-centiles.pred(gam1,xname='t0',xvalues=0:100, type='centiles',cent=c(95)) # 95th prediction interval on (0,1) scale
bet1mean<-predict(bet1,data.frame(t0=0:100),type='response')
bet1lo<-predict(bet1,data.frame(t0=0:100),type='quantile',at=c(0.05))
bet1hi<-predict(bet1,data.frame(t0=0:100),type='quantile',at=c(0.95))
mins<-predict(quant,data.frame(t0=0:100))[,1]
maxs<-predict(quant,data.frame(t0=0:100))[,2]

plot(size$t0,size$t1b,xlab='Size at time t',ylab='Beta-transformed size t+1',col='grey')
lines(0:100,gam1mean,col='black',lwd=2)
lines(0:100,gam1lo[,2],col='black',lwd=2,lty=2)
lines(0:100,gam1hi[,2],col='black',lwd=2,lty=2)
lines(0:100,bet1mean,col='blue',lwd=2)
lines(0:100,bet1lo,col='blue',lwd=2,lty=2)
lines(0:100,bet1hi,col='blue',lwd=2,lty=2)

plot(size$t0,size$t1,xlab='Size at time t',ylab='Size at time t+1',col='grey')
lines(0:100,backbeta(gam1mean,mins,maxs),col='black',lwd=2)
lines(0:100,backbeta(gam1lo[,2],mins,maxs),col='black',lwd=2,lty=2)
lines(0:100,backbeta(gam1hi[,2],mins,maxs),col='black',lwd=2,lty=2)
lines(0:100,backbeta(bet1mean,mins,maxs),col='blue',lwd=2)
lines(0:100,backbeta(bet1lo,mins,maxs),col='blue',lwd=2,lty=2)
lines(0:100,backbeta(bet1hi,mins,maxs),col='blue',lwd=2,lty=2)
legend('topleft',col=c('black','blue'),c('cubic spline','quadratic'),pch=16)

