##### Appendix 2: This script simulates data representing different biological scenarios
##### and compares the fit of the beta vs. normal approach

rm(list=ls(all=T))

setwd("e:\\admin\\reviews\\BES"); 
require(quantreg)
require(betareg)
require(moments)
require(zoo)
require(MuMIn)
require(beanplot)
require(dplyr)
require(reshape2)
require(MuMIn)
require(reldist)

### Functions related to beta approach

# This function transforms continuous size data to a (0,1) interval
# based on size-dependent min and max values
betaFn<-function(x,min,max){
  y=(x-min)/(max-min)  
  n=length(x)
  if(any(y<0.001|y>0.999,na.rm=TRUE)){ # correct for 0s and 1s following Smithson & Verkuilen 2006
    y=(y*(n-1)+0.5)/n
  }
  return(y)
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

# This function simulates n values for each starting size z from a beta distribution
# defined by mean mu and variance sigma and backtransforms to the size scale based on min and max
simulatebeta<-function(mu,sigma,z,n,min,max){
  alpha<-((1-mu)/sigma-1/mu)*mu^2
  beta<-alpha*(1/mu-1)
  outs<-c()
  for(i in 1:length(alpha)){
    datb<-rbeta(n,shape1=alpha[i],shape2=beta[i])
    datb[datb==1]<-0.999
    datb[datb==0]<-0.001
    dat<-backbeta(datb,min=min[i],max=max[i])
    outs<-rbind(outs,cbind(dat,z=z[i]))
  }
  as.data.frame(outs)
}


#### Simulate 4 data scenarios
zb<-rep(seq(0,1,length.out=100))
z0<-rep(seq(1,100,length.out=100))
n=length(z0)
nsim=100

data<-list()
skew<-list()
variance<-list()
alpha<-list()
beta<-list()
mins<-list()
maxs<-list()

# 1) Constant and symmetric variance
mu=rep(0.5,n)
sigma=rep(0.03,n)
params=shapebeta(mu=mu,sigma=sigma)
minSize=z0-20
maxSize=z0+40
hist(rbeta(1000,shape1=params$alpha[1],shape2=params$beta[1]),
     main='Constant and symmetric variance',xlab='Beta-transformed size')
df1<-simulatebeta(mu=mu,sigma=sigma,z=z0,n=nsim,min=minSize,max=maxSize)
data[[1]]<-df1
plot(df1$z,df1$dat,main='Constant and symmetric variance')
skew[[1]]<-c(by(df1$dat,df1$z,skewness))
variance[[1]]<-c(by(df1$dat,df1$z,var))
alpha[[1]]<-params$alpha
beta[[1]]<-params$beta
mins[[1]]<-minSize
maxs[[1]]<-maxSize

# 2) Increasing and symmetric variance
mu=rep(0.5,n)
sigma=rep(0.03,n)
params=shapebeta(mu=mu,sigma=sigma)
minSize=z0*0.8 - 5
maxSize=z0*1.2 + 5
hist(rbeta(1000,shape1=params$alpha[1],shape2=params$beta[1]),
     main='Constant and symmetric variance',xlab='Beta-transformed size')
hist(rbeta(1000,shape1=params$alpha[100],shape2=params$beta[100]),
     main='Constant and symmetric variance',xlab='Beta-transformed size')
df2<-simulatebeta(mu=mu,sigma=sigma,z=z0,n=nsim,min=minSize,max=maxSize)
data[[2]]<-df2
plot(df2$z,df2$dat,main='Increasing and symmetric variance')
skew[[2]]<-c(by(df2$dat,df2$z,skewness))
variance[[2]]<-c(by(df2$dat,df2$z,var))
alpha[[2]]<-params$alpha
beta[[2]]<-params$beta
mins[[2]]<-minSize
maxs[[2]]<-maxSize

# 3) Constant and asymmetric variance
mu=rep(0.8,n)
sigma=rep(0.03,n)
params=shapebeta(mu=mu,sigma=sigma)
hist(rbeta(1000,shape1=params$alpha[1],shape2=params$beta[1]),
     main='Constant and asymmetric variance',xlab='Beta-transformed size')
hist(rbeta(1000,shape1=params$alpha[100],shape2=params$beta[100]),
     main='Constant and asymmetric variance',xlab='Beta-transformed size')
minSize=z0-20
maxSize=z0+40
df3<-simulatebeta(mu=mu,sigma=sigma,z=z0,n=nsim,min=minSize,max=maxSize)
data[[3]]<-df3
plot(df3$z,df3$dat,main='Constant and asymmetric variance')
skew[[3]]<-c(by(df3$dat,df3$z,skewness))
variance[[3]]<-c(by(df3$dat,df3$z,var))
alpha[[3]]<-params$alpha
beta[[3]]<-params$beta
mins[[3]]<-minSize
maxs[[3]]<-maxSize

# 4) Increasing and asymmetric variance
mu=0.5+z0*0.004
sigma=rep(0.03,n)
params=shapebeta(mu=mu,sigma=sigma)
hist(rbeta(1000,shape1=params$alpha[1],shape2=params$beta[1]))
hist(rbeta(1000,shape1=params$alpha[100],shape2=params$beta[100]))
minSize=z0*0.8 - 5
maxSize=z0*1.2 + 5
df4<-simulatebeta(mu=mu,sigma=sigma,z=z0,n=nsim,min=minSize,max=maxSize)
data[[4]]<-df4
plot(df4$z,df4$dat,main='Increasing and asymmetric variance')
skew[[4]]<-c(by(df4$dat,df4$z,skewness))
variance[[4]]<-c(by(df4$dat,df4$z,var))
alpha[[4]]<-params$alpha
beta[[4]]<-params$beta
mins[[4]]<-minSize
maxs[[4]]<-maxSize

### Fit beta and linear regressions and visualize fits
# store fitted beta and linear models
mb<-list()
mq<-list()
mn<-list()
mv<-list()
# store AIC table for model selection
mbAIC<-list()
mnAIC<-list()
mvAIC<-list()
# store AIC for each approach
LLBeta<-list()
LLNorm<-list()
KBeta<-list()
KNorm<-list()
AICBeta<-list()
AICNorm<-list()


for (i in 1:4){
  df<-data[[i]]
  # quantile regression for min and max size limits
  mq[[i]]<-rq(dat~z+I(z^2),data=df,tau=c(0,1))
  # transform to (0,1) interval
  SizeEst<-predict(mq[[i]],df)
  df$datb<-betaFn(df$dat,min=SizeEst[,1],max=SizeEst[,2])
  # beta regression model selection
  mba<-betareg(datb~z+I(z^2)|z+I(z^2),data=df)
  mbb<-betareg(datb~z|z+I(z^2),data=df)
  mbc<-betareg(datb~z+I(z^2)|z,data=df)
  mbd<-betareg(datb~z|z,data=df)
  mbe<-betareg(datb~z+I(z^2),data=df)
  mbf<-betareg(datb~z,data=df)
  mbAIC[[i]]<-AIC(mba,mbb,mbc,mbd,mbe,mbf)
  mb[[i]]<-list(mba,mbb,mbc,mbd,mbe,mbf)[[which.min(mbAIC[[i]]$AIC)]]
  # linear regression model selection
  mna<-lm(dat~z+I(z^2),data=df)
  mnb<-lm(dat~z,data=df)
  mnAIC[[i]]<-AIC(mna,mnb)
  mn[[i]]<-list(mna,mnb)[[which.min(mnAIC[[i]]$AIC)]]
  # variance regression model selection
  df$resid2<-(df$dat-predict(mn[[i]],df))^2
  mva<-lm(resid2~z+I(z^2),data=df)
  mvb<-lm(resid2~z,data=df)
  mvc<-lm(resid2~1,data=df)
  mvAIC[[i]]<-AIC(mva,mvb,mvc)
  mv[[i]]<-list(mva,mvb,mvc)[[which.min(mvAIC[[i]]$AIC)]]
  # min and max size predictions
  Size<-predict(mq[[i]],data.frame(z=z0))
  # best-fit beta predictions
  pb<-predict(mb[[i]],data.frame(z=z0),type='response')
  lb<-predict(mb[[i]],data.frame(z=z0),type='quantile',at=c(0.005)) # 99th prediction interval
  hb<-predict(mb[[i]],data.frame(z=z0),type='quantile',at=c(0.995)) # 99th prediction interval
  # best-fit linear predictions
  pn<-predict(mn[[i]],data.frame(z=z0))
  pv<-predict(mv[[i]],data.frame(z=z0))
  ln<-qnorm(0.005,mean=pn,sd=sqrt(pv)) # 99th prediction interval
  hn<-qnorm(0.995,mean=pn,sd=sqrt(pv)) # 99th prediction interval
  # plot model fits to data
  plot(df$z,df$dat,ylab='Size t+1',xlab='Size',main=i,col='grey')
  lines(z0,pn,lwd=1,lty=1,col='darkblue')
  lines(z0,ln,lwd=1,lty=2,col='darkblue')
  lines(z0,hn,lwd=1,lty=2,col='darkblue')
  lines(z0,backbeta(pb,min=Size[,1],max=Size[,2]),lwd=1,lty=1,col='darkred')
  lines(z0,backbeta(lb,min=Size[,1],max=Size[,2]),lwd=1,lty=2,col='darkred')
  lines(z0,backbeta(hb,min=Size[,1],max=Size[,2]),lwd=1,lty=2,col='darkred')
  # calculate AIC for each approach
  Params=shapebeta(mu=predict(mb[[i]],df,type='response'),sigma=predict(mb[[i]],df,type='variance'))
  LikBeta<-dbeta(df$datb,shape1=Params$alpha,shape2=Params$beta)
  LLBeta[[i]]<-sum(log(LikBeta))-sum(log(SizeEst[,2]-SizeEst[,1]))
  KBeta[[i]]<-length(c(coef(mb[[i]]),coef(mq[[i]])))
  AICBeta[[i]]<- -2*LLBeta[[i]] + 2*KBeta[[i]]
  df$datn<-backbeta(df$datb,min=SizeEst[,1],max=SizeEst[,2]) # compare to same dataset
  LikNorm<-dnorm(df$datn,mean=predict(mn[[i]],df),sd=sqrt(predict(mv[[i]],df)))
  LLNorm[[i]]<-sum(log(LikNorm))
  KNorm[[i]]<-length(c(coef(mn[[i]]),coef(mv[[i]])))
  AICNorm[[i]]<- -2*LLNorm[[i]] + 2*KNorm[[i]]
  # Plot normal approach after re-normalization to prevent eviction
  minsize<-min(df$dat); maxsize<-max(df$dat)
  vec.bin<-c(minsize, minsize+1:100*(maxsize-minsize)*(1/100))
  binmids <- 0.5*(vec.bin[2:length(vec.bin)] + vec.bin[1:(length(vec.bin)-1)])
  n.bin = length(binmids)
  pn.bin<-predict(mn[[i]],data.frame(z=binmids))
  pv.bin<-predict(mv[[i]],data.frame(z=binmids))
  pv.bin[pv.bin<=0.05]<-0.05
  pn.renorm=ln.renorm=hn.renorm=numeric(n.bin)
  for (ss in 1:(n.bin)) {
    growcdf <- pnorm(vec.bin,pn.bin[ss],sqrt(pv.bin[ss]))
    grows <- growcdf[2:length(vec.bin)]-growcdf[1:(length(vec.bin)-1)]
    if(sum(grows)>0){grows <- grows/sum(grows)}
    vecnew<-sample(binmids,prob=grows,replace=T,size=100000)
    pn.renorm[ss]<-mean(vecnew)
    ln.renorm[ss]<-quantile(vecnew,c(0.005))
    hn.renorm[ss]<-quantile(vecnew,c(0.995))
  }
  plot(df$z,df$dat,ylab='Size t+1',xlab='Size',main=i,col='grey')
  lines(z0,pn,lwd=1,lty=1,col='darkblue')
  lines(z0,ln,lwd=1,lty=2,col='darkblue')
  lines(z0,hn,lwd=1,lty=2,col='darkblue')
  lines(binmids,pn.renorm,lwd=1,lty=1,col='blue')
  lines(binmids,ln.renorm,lwd=1,lty=2,col='blue')
  lines(binmids,hn.renorm,lwd=1,lty=2,col='blue')
  lines(z0,backbeta(pb,min=Size[,1],max=Size[,2]),lwd=1,lty=1,col='darkred')
  lines(z0,backbeta(lb,min=Size[,1],max=Size[,2]),lwd=1,lty=2,col='darkred')
  lines(z0,backbeta(hb,min=Size[,1],max=Size[,2]),lwd=1,lty=2,col='darkred')
  # Plot predicted vs. observed distributions at a given starting size
  xsize<-c(10,30,50,60,70,80,90)
  znew<-seq(1,150,length.out=300)
  zSize<-predict(mq[[i]],data.frame(z=xsize))
  zpb<-predict(mb[[i]],data.frame(z=xsize),type='response')
  zvb<-predict(mb[[i]],data.frame(z=xsize),type='variance')
  params<-shapebeta(mu=zpb,sigma=zvb)
  zpn<-predict(mn[[i]],data.frame(z=xsize))
  zvn<-predict(mv[[i]],data.frame(z=xsize))
  for(j in 1:length(xsize)){
    x<-xsize[j]
    zn<-backbeta(x=zb,min=zSize[j,1],max=zSize[j,2])
    db<-dbeta(zb,shape1=params$alpha[j],shape2=params$beta[j])
    dn<-dnorm(znew,mean=zpn[j],sd=sqrt(zvn[j]))
    d.renorm<-pnorm(vec.bin,zpn[j],sqrt(zvn[j]))
    d.renorm<-d.renorm[2:length(vec.bin)]-d.renorm[1:(length(vec.bin)-1)]
    if(sum(d.renorm)>0){d.renorm <- d.renorm/sum(d.renorm)}
    hist(df$dat[df$z==x&df$z==x],xlim=c(zSize[j,1]*0.9,zSize[j,2]*1.1),ylim=c(0,max(dn)*2),freq=F,main=paste(i,x),xlab="Size t+1")
    lines(zn,db/(zn[100]-zn[1]),type='l',lwd=1,col='darkred')
    lines(znew,dn,type='l',lwd=1,col='darkblue')
    lines(binmids,d.renorm,type='l',lwd=1,lty=2,col='darkblue')
  }
  # Plot predicted vs. observed distributions at a given starting size after re-normalization
  xsize<-c(10,30,50,60,70,80,90,98)
  znew<-seq(1,150,length.out=300)
  zSize<-predict(mq[[i]],data.frame(z=xsize))
  zpb<-predict(mb[[i]],data.frame(z=xsize),type='response')
  zvb<-predict(mb[[i]],data.frame(z=xsize),type='variance')
  params<-shapebeta(mu=zpb,sigma=zvb)
  zpn<-predict(mn[[i]],data.frame(z=xsize))
  zvn<-predict(mv[[i]],data.frame(z=xsize))
  for(j in 1:length(xsize)){
    x<-xsize[j]
    vec.beta<-betaFn(x=vec.bin,min=zSize[j,1],max=zSize[j,2])
    db<-pbeta(vec.beta,shape1=params$alpha[j],shape2=params$beta[j])
    db<-db[2:length(vec.bin)]-db[1:(length(vec.bin)-1)]
    dn<-pnorm(vec.bin,zpn[j],sqrt(zvn[j]))
    dn<-dn[2:length(vec.bin)]-dn[1:(length(vec.bin)-1)]
    d.renorm<-dn
    if(sum(d.renorm)>0){d.renorm <- d.renorm/sum(d.renorm)}
    hist(df$dat[df$z==x&df$z==x],xlim=c(zSize[j,1]*0.9,zSize[j,2]*1.1),ylim=c(0,max(dn)*2),freq=F,main=paste(i,x),xlab="Size t+1")
    lines(binmids,db,type='l',lwd=1,col='darkred')
    lines(binmids,dn,type='l',lwd=1,col='darkblue')
    lines(binmids,d.renorm,type='l',lwd=1,lty=2,col='darkblue')
  }
}


##### Many simulations of 200 data points for scenario #4
n=10
z0<-rep(seq(1,100,length.out=n))
mu=0.5+z0*0.004
sigma=rep(0.03,n)
params=shapebeta(mu=mu,sigma=sigma)
minSize=z0*0.8 - 5
maxSize=z0*1.2 + 5

truemean<-backbeta(mu,min=minSize,max=maxSize)
truevar<-backvar(sigma,min=minSize,max=maxSize)
trueskew<-skewbeta(alpha=params$alpha,beta=params$beta)
bmean<-matrix(NA,nrow=500,ncol=n)
bvar<-matrix(NA,nrow=500,ncol=n)
bskew<-matrix(NA,nrow=500,ncol=n)
nmean<-matrix(NA,nrow=500,ncol=n)
nvar<-matrix(NA,nrow=500,ncol=n)

for (i in 1:500){
  df<-simulatebeta(mu=mu,sigma=sigma,z=z0,n=20,min=minSize,max=maxSize)
  mquant<-rq(dat~z+I(z^2),data=df,tau=c(0.05,0.995))
  # transform to (0,1) interval
  SizeEst<-predict(mquant,df)
  df$datb<-betaFn(df$dat,min=SizeEst[,1],max=SizeEst[,2])
  # beta regression model selection
  mba<-betareg(datb~z+I(z^2)|z+I(z^2),data=df)
  mbb<-betareg(datb~z|z+I(z^2),data=df)
  mbc<-betareg(datb~z+I(z^2)|z,data=df)
  mbd<-betareg(datb~z|z,data=df)
  mbe<-betareg(datb~z+I(z^2),data=df)
  mbf<-betareg(datb~z,data=df)
  mbeta<-list(mba,mbb,mbc,mbd,mbe,mbf)[[which.min(AIC(mba,mbb,mbc,mbd,mbe,mbf)$AIC)]]
  # linear regression model selection
  mna<-lm(dat~z+I(z^2),data=df)
  mnb<-lm(dat~z,data=df)
  mnorm<-list(mna,mnb)[[which.min(AIC(mna,mnb)$AIC)]]
  # variance regression model selection
  df$resid2<-(df$dat-predict(mnorm,df))^2
  mva<-lm(resid2~z+I(z^2),data=df)
  mvb<-lm(resid2~z,data=df)
  mvc<-lm(resid2~1,data=df)
  mvar<-list(mva,mvb,mvc)[[which.min(AIC(mva,mvb,mvc)$AIC)]]
  # get parameter estimates
  Size<-predict(mquant,data.frame(z=z0))
  nmean[i,]<-predict(mnorm,data.frame(z=z0))
  nvar[i,]<-predict(mvar,data.frame(z=z0))
  bmean[i,]<-backbeta(predict(mbeta,data.frame(z=z0),type='response'),min=Size[,1],max=Size[,2])
  bvar[i,]<-backvar(predict(mbeta,data.frame(z=z0),type='variance'),min=Size[,1],max=Size[,2])
  EstParams<-shapebeta(mu=predict(mbeta,data.frame(z=z0),type='response'),sigma=predict(mbeta,data.frame(z=z0),type='variance'))
  bskew[i,]<-skewbeta(alpha=EstParams$alpha,beta=EstParams$beta)
}

# compare estimated moments to the true moments
means<-bind_rows(list(norm=melt(nmean),beta=melt(bmean)),.id='source')
vars<-bind_rows(list(norm=melt(nvar),beta=melt(bvar)),.id='source')
skews<-bind_rows(list(beta=melt(bskew)),.id='source')

beanplot(value~source*Var2,main='',ylab="Mean",xlab='Size at time t',what=c(0,1,1,0),xaxt='n',
         side="both",border=NA,col=list("darkred", c("darkblue", "white")),data=means)
legend("topleft", fill = c("darkred", "darkblue"), legend = c("Beta", "Normal"))
axis(side=1,at=1:n,labels=z0)
lines(truemean,lwd=2,col='grey')

beanplot(value~source*Var2,main='',ylab="Variance",xlab='Size at time t',what=c(0,1,1,0),xaxt='n',
         side="both",border=NA,col=list("darkred", c("darkblue", "white")),data=vars,ylim=c(-10,150))
axis(side=1,at=1:n,labels=z0)
lines(truevar,lwd=2,col='grey')

beanplot(value~Var2,main='',ylab="Skewness",xlab='Size at time t',what=c(0,1,1,0),xaxt='n',
        side='first',border=NA,col=list("darkred", c("darkred", "white")),data=skews)
axis(side=1,at=1:n,labels=z0)
lines(trueskew,lwd=2,col='grey')
abline(h=0,lwd=2,col='darkblue')

