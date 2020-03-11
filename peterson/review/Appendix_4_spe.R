#### Appendix 4: This script uses red gorgonian data to illustrate using the 
#### gamlss package to fit the parameters of the beta distribution as cubic splines

rm(list=ls()); 

graphics.off(); setwd("e:\\admin\\reviews\\BES"); 

require(quantreg)
require(betareg)
require(moments)
require(zoo)
require(dplyr)
require(MuMIn)
require(gamlss)
require(stats4)
require(bbmle); 
require(sn); 

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
# quant<-rq(t1~t0,data=size,tau=c(0.005,0.995)) # fit quantile regression
summary(quant)
minsizes=pmax(0,predict(quant,size)[,1]) # size dependent minimum
maxsizes=predict(quant,size)[,2] # size dependent maximum

dev.new(); 
plot(size$t0,size$t1,xlab='Size at time t',ylab='Size at time t+1',col='grey')
lines(0:100,pmax(0,predict(quant,data.frame(t0=0:100))[,1]),col='black')
lines(0:100,predict(quant,data.frame(t0=0:100))[,2],col='black')

# Step 2: Transform size at time t+1 to (0,1) interval
size$t1b<-betaFn(x=size$t1,min=minsizes,max=maxsizes)
range(size$t1b,na.rm=T) # produces some values outside (0,1) interval
size$t1b[size$t1b>0.99]=0.99
size$t1b[size$t1b<0.01]=0.01

plot(size$t0,size$t1b,xlab='Size at time t',ylab='Beta-transformed size t+1',col='grey')

# Step 3: Fit beta distribution parameters with cubic splines and compare
# to parametric models fit with betareg
gam1<-gamlss(t1b~cs(t0,3),sigma.formula=~cs(t0,3), data=size, family=BE) # mean and standard deviation fit as cubic splines of size
gam2<-gamlss(t1b~cs(t0,3), data=size, family=BE) # mean fit as cubic spline of size, constant standard deviation

bet1<-betareg(t1b~t0+I(t0^2)|t0+I(t0^2),data=size)
bet2<-betareg(t1b~t0|t0+I(t0^2),data=size)
bet3<-betareg(t1b~t0+I(t0^2)|t0,data=size)
bet4<-betareg(t1b~t0|t0,data=size)

AIC(gam1,gam2) # gam1 is best-supported
AIC(bet1,bet2,bet3,bet4) #bet2 is best supported, and used below 

LikBeta<-dBE(size$t1b,mu=gam1$mu.fv,sigma=gam1$sigma.fv)
LLBeta<-sum(log(LikBeta))-sum(log(maxsizes-minsizes)) # transform to original datascale
KBeta<-length(coef(quant)) + gam1$df.fit; 
AICgamBeta<- -2*LLBeta + 2*KBeta;  

# Step 4: Compare model predictions to the data
gam1mean<-predict(gam1,newdata=data.frame(t0=0:100),type='response')
gam1lo<-centiles.pred(gam1,xname='t0',xvalues=0:100, type='centiles',cent=c(1)) # 5th prediction interval on (0,1) scale
gam1hi<-centiles.pred(gam1,xname='t0',xvalues=0:100, type='centiles',cent=c(99)) # 95th prediction interval on (0,1) scale
bet1mean<-predict(bet2,data.frame(t0=0:100),type='response')
bet1lo<-predict(bet2,data.frame(t0=0:100),type='quantile',at=c(0.01))
bet1hi<-predict(bet2,data.frame(t0=0:100),type='quantile',at=c(0.99))
mins<-pmax(0,predict(quant,data.frame(t0=0:100))[,1])
maxs<-predict(quant,data.frame(t0=0:100))[,2]

par(bty="l",cex.axis=1.4,cex.lab=1.4,mar=c(5,5,1,1)); 
plot(size$t0,size$t1b,xlab='Size at time t',ylab='Beta-transformed size t+1',col='grey')
lines(0:100,gam1mean,col='black',lwd=2)
lines(0:100,gam1lo[,2],col='black',lwd=2,lty=2)
lines(0:100,gam1hi[,2],col='black',lwd=2,lty=2)
lines(0:100,bet1mean,col='darkred',lwd=2)
lines(0:100,bet1lo,col='darkred',lwd=2,lty=2)
lines(0:100,bet1hi,col='darkred',lwd=2,lty=2)
legend('right',col=c('black','darkred'),c('Cubic spline Beta','Parametric Beta'),pch=16,cex=1.5)
dev.copy2pdf(file="Fig2Brev.pdf"); 

############### Last plot ############################################################################
par(bty="l",cex.axis=1.4,cex.lab=1.4,mar=c(5,5,1,1)); 
plot(size$t0,size$t1,xlab='Size at time t',ylab='Size at time t+1',col='grey')
lines(0:100,backbeta(gam1mean,mins,maxs),col='black',lwd=2)
lines(0:100,backbeta(gam1lo[,2],mins,maxs),col='black',lwd=2,lty=2)
lines(0:100,backbeta(gam1hi[,2],mins,maxs),col='black',lwd=2,lty=2)
lines(0:100,backbeta(bet1mean,mins,maxs),col='red',lwd=2)
lines(0:100,backbeta(bet1lo,mins,maxs),col='red',lwd=2,lty=2)
lines(0:100,backbeta(bet1hi,mins,maxs),col='red',lwd=2,lty=2)
legend('topleft',col=c('black','red','blue','forestgreen','purple'),
c('Spline Beta','Parametric Beta','Normal','Parametric SkewNormal','Spline Skew-Normal'),pch=16,cex=1.5)

#### Now compare to the normal approach
# Step 2N: Fit and compare linear regressions for mean growth
m1<-lm(t1~t0+I(t0^2),data=size) # size and size2 on mean 
m2<-lm(t1~t0,data=size) # size on mean
m3<-lm(t1~1,data=size) # no size effect
m4<-nls(t1 ~ a + b*(t0^c), data= size, start=list(a= 1, b=1,c=1)) # power function of size on mean

AICc(m1,m2,m3,m4) # m1 is best-supported
bestmean<-m1

# Step 3N: Get residuals
size$resid<-size$t1-predict(bestmean,size)
size$resid2<-size$resid^2

# Step 4N: Fit and compare linear regressions for variance in growth
m1<-lm(resid2~t0+I(t0^2),data=size) # size and size2 on variance 
m2<-lm(resid2~t0,data=size) # size on variance
m3<-lm(resid2~1,data=size) # no size effect

AICc(m1,m2,m3) # m1 is best-supported
bestvar<-m1

# Step 5N: Get parameters and pdf for a given starting size
normMean<-predict(bestmean,data.frame(t0=(0:100))) 
normVar<-predict(bestvar,data.frame(t0=(0:100))) 
normVar[normVar<0.001]<-0.001 # keep variance positive

# Step 6N: Visualize model predictions
normLow<-qnorm(0.01,mean=normMean,sd=sqrt(normVar))
normHi<-qnorm(0.99,mean=normMean,sd=sqrt(normVar))

lines(0:100,normMean,col="blue",lwd=2,lty=1)
lines(0:100,normLow,col="blue",lwd=2,lty=2)
lines(0:100,normHi,col="blue",lwd=2,lty=2)

################################ New Step: fit a Skewed Normal ###########################
NLLsn <- function(a,b,c,av,bv,cv,dv,aa,ba,ca) {
    mean=a+b*z0+c*z0^2;
    sd = sqrt(av + bv*z0 + cv*z0^2 + dv*z0^3) # let the variance be cubic; AIC says 'yes' to this. 
    alpha = aa + ba*z0 + ca*z0^2; 
    -sum(dsn(z1,xi=mean,omega=sd,alpha=alpha,tau=0,log=TRUE))
}    

if(!exists("fitSN")) {
z0 = size$t0; z1=size$t1; 
fit0=lm(z1~z0+I(z0^2)); sigma0=sqrt(mean(fit0$residuals^2));
c0=as.vector(coef(fit0)); 
fitSN=mle2(NLLsn,start=list(a=c0[1],b=c0[2],c=c0[3],av=sigma0^2,bv=0.001,cv=0.001,dv=0.001,aa=0,ba=0,ca=0),
method="Nelder-Mead",control=list(trace=4,maxit=10000)); 
for(k in 1:200) { # it really is that bad 
    c1 = as.list(coef(fitSN)); 
    fitSN=mle2(NLLsn,start=c1, method="Nelder-Mead",skip.hessian=TRUE,control=list(trace=4,maxit=10000)); 
    cat("===============================================",k,"\n"); 
}
} 
AIC(fitSN); summary(fitSN); 

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
omegaFun = function(x) sqrt(c1[4]+c1[5]*x + c1[6]*x^2 + c1[7]*x^3);  
alphaFun = function(x) c1[8]+ c1[9]*x + c1[10]*x^2;
px=seq(min(size$t0),max(size$t0),length=50);  
qSN = matrix(NA,length(px),3); 
for(j in 1:length(px)) {
    qSN[j,]=my.qsn(c(0.01,0.5,0.99),xiFun(px[j]),omegaFun(px[j]),alphaFun(px[j]))

}
matpoints(px,qSN,lty=c(2,1,2),col="forestgreen",type="l",lwd=2); 

if(!exists("gamSN")) {
gamSN<-gamlss(t1~t0+I(t0^2),sigma.formula=~cs(t0,3), nu.formula=~cs(t0,3), data=size, family=SN1,method=RS(250))
} 
   
muFun = approxfun(size$t0,predict(gamSN,what="mu"),rule=2);
sigmaFun = approxfun(size$t0,gamSN$sigma.fv,rule=2);  
nuFun = approxfun(size$t0,gamSN$nu.fv,rule=2); 

points(size$t0,size$t1,xlab='Size at time t',ylab='Size at time t+1',col='grey')
px=seq(min(size$t0),max(size$t0),length=150);  
SNlow=qSN1(0.01,mu=muFun(px),sigma=sigmaFun(px),nu=nuFun(px)); 
SNmid=qSN1(0.5,mu=muFun(px),sigma=sigmaFun(px),nu=nuFun(px)); 
SNhigh=qSN1(0.99,mu=muFun(px),sigma=sigmaFun(px),nu=nuFun(px));
matpoints(px,cbind(SNlow,SNmid,SNhigh),lty=c(2,1,2),col="purple",type="l",lwd=2);  
dev.copy2pdf(file="Fig2Erev.pdf"); 

graphics.off(); 
dev.new(width=11,height=5); par(mfrow=c(1,2),bty="l",cex.axis=1.3,cex.lab=1.3,mgp=c(2.2,1,0))
require(viridis); 
colors=plasma(10); 
px=seq(0,85,length=500);
jx=10; 
plot(px, dsn(px,xi=xiFun(jx),omega=omegaFun(jx),alpha=alphaFun(jx)),type="l",lty=1,col=colors[1],
    xlim=c(0,85),xlab="Size at time t",ylab="Size distribution at time t+1",lwd=2) 
for(j in 2:7) {
    jx = j*10;  
    points(px, dsn(px,xi=xiFun(jx),omega=omegaFun(jx),alpha=alphaFun(jx)),type="l",lty=1,col=colors[j],lwd=2) 
}

py <- dsn(px, xi = xiFun(60),omega=omegaFun(60),alpha=alphaFun(60))
sub<-size[which(size$t0>53&size$t0<67),] # get actual size at t+1 for starting size between 53 and 67
hist(sub$t1,breaks=10,freq=F,main='',xlab="Size at time t+1")
lines(px,py,type='l',lwd=2,col='darkred')
out=density(sub$t1,bw = "ucv"); points(out$x,out$y,col="darkred",lty=2,lwd=2,type="l"); 
legend("topleft",legend=c("Skew-Normal","Kernel density estimate"),col="darkred",lty=c(1,2),bty="n",inset=0.05); 
dev.copy2pdf(file="FigS2skewNormal.pdf"); 
