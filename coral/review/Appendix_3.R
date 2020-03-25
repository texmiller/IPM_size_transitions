#### Appendix 3: This script uses red gorgonian data to illustrate one method
#### to simultaneously fit the minimum, maximum, mean, and precision of the beta model

rm(list=ls())
require(quantreg)
require(betareg)
require(dplyr)
require(MuMIn)
require(MASS)
require(MBESS)

setwd("e:\\admin\\reviews\\BES"); 

# This function transforms continuous size data to a (0,1) interval
# based on size-dependent min and max values
betaFn<-function(x,min,max){
  y=(x-min)/(max-min)
}

# Load, manipulate, and sort data
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

## Step 1: get covariance matrix to use for sampling parameters 
# for the minimum and maximum sizes
quant<-rq(t1~t0,data=size,tau=c(0.005,0.995))
quantcovmin<-summary.rq(rq(t1~t0,data=size,tau=c(0.005)),covariance=T)$cov 
quantcovmax<-summary.rq(rq(t1~t0,data=size,tau=c(0.995)),covariance=T)$cov
minsizes=predict(quant,size)[,1] # size dependent minimum
maxsizes=predict(quant,size)[,2] # size dependent maximum

plot(size$t0,size$t1)
lines(0:100,predict(quant,data.frame(t0=0:100))[,1],col='red')
lines(0:100,predict(quant,data.frame(t0=0:100))[,2],col='red')

## Step 2: transform to (0,1) interval and 
# keep track of the number of data points outside these bounds
size$t1b<-betaFn(size$t1,min=minsizes,max=maxsizes)
startnout<-length(c(size$t1b[size$t1b>=1],size$t1b[size$t1b<=0])) 

## Step 3: sample parameter values for max and min sizes
# from their covariance matrices 
ndraw<-10000 
plo<-mvrnorm(ndraw,mu=coef(quant)[,1],Sigma=quantcovmin)
phi<-mvrnorm(ndraw,mu=coef(quant)[,2],Sigma=quantcovmax)
LL<-numeric(ndraw)
nout<-numeric(ndraw)

## Step 4: for each set of parameter values, transform data to (0,1) interval
# and, if not excluding any more data points, fit the beta-regression
# then get the log likelihood of the entire model to the data
for(i in 1:ndraw){ cat(i,"\n"); 
  mins<-plo[i,1]+plo[i,2]*size$t0 # new min sizes
  maxs<-phi[i,1]+phi[i,2]*size$t0 # new max sizes
  size$t1bfit<-(size$t1-mins)/(maxs-mins) # new (0,1) interval
  nout[i]<-length(c(size$t1bfit[size$t1bfit>=1],size$t1bfit[size$t1bfit<=0])) # number of data points outside new interval
  if(nout[i]<=startnout){ # only consider parameters that don't exclude additional data points
    size$t1bfit[size$t1bfit>0.99]<-0.99
    size$t1bfit[size$t1bfit<0.01]<-0.01
    b1<-betareg(t1bfit~t0+I(t0^2)|t0+I(t0^2),data=size)
    b2<-betareg(t1bfit~t0|t0+I(t0^2),data=size)
    b3<-betareg(t1bfit~t0+I(t0^2)|t0,data=size)
    b4<-betareg(t1bfit~t0|t0,data=size)
    betamod<-list(b1,b2,b3,b4)[[which.min(AIC(b1,b2,b3,b4)$AIC)]]
    LL[i]<-logLik(betamod)-sum(log(maxs-mins))
  } else {
    LL[i]<-NA
  }
}

# Step 5: find the parameters that maximize the log likelihood
best<-which.max(LL)
plo[best,];phi[best,]
mins<-plo[best,1]+plo[best,2]*size$t0
maxs<-phi[best,1]+phi[best,2]*size$t0

plot(size$t0,size$t1)
lines(0:100,predict(quant,data.frame(t0=0:100))[,1],col='red')
lines(0:100,predict(quant,data.frame(t0=0:100))[,2],col='red')
lines(0:100,plo[best,1]+plo[best,2]*0:100,col='red',lty=2)
lines(0:100,phi[best,1]+phi[best,2]*0:100,col='red',lty=2)
legend('topleft',lty=c(1,2),c('Original','Updated'))

# Step 6: Transform to updated (0,1) interval
size$t1bfit<-betaFn(size$t1,min=mins,max=maxs)
size$t1bfit[size$t1bfit>0.99]<-0.99
size$t1bfit[size$t1bfit<0.01]<-0.01
b1<-betareg(t1bfit~t0+I(t0^2)|t0+I(t0^2),data=size)
b2<-betareg(t1bfit~t0|t0+I(t0^2),data=size)
b3<-betareg(t1bfit~t0+I(t0^2)|t0,data=size)
b4<-betareg(t1bfit~t0|t0,data=size)
betamod<-list(b1,b2,b3,b4)[[which.min(AIC(b1,b2,b3,b4)$AIC)]]
LLNew<-logLik(betamod)-sum(log(maxs-mins)) # log likelihood on the original data scale


