graphics.off(); 
rm(list=ls(all=TRUE)); 
library(quantreg); library(fda); library(sn); 
library(moments); library(dplyr); library(zoo); library(bbmle); 

############## Pull in gorgonian data, to fit an SN distribution and then 
############## generate synthetic data with known properties 
setwd("e:\\admin\\reviews\\BES"); 
dat <- read.csv("Gorgonian_raw_data.csv")
setwd("c:\\repos\\IPM_size_transitions\\gam practice"); 

size <- as.data.frame(dat)
colnames(size)[6:8] <- c("Size"="t0", "Sizenext"="t1", "Survnext"= "survival") # rename column header
size <- arrange(size, t0) # arrange initial column in ascending orders
size$t1[which(size$survival==0)] <- NA 
size <- size[-which(is.na(size$t0)),]
size <- size[-which(is.na(size$t1)),]
size <- size[size$Site=='Portcros',] # For this example, use Portcros site 
nx = nrow(size); 
size$z0 = size$t0^0.5; size$z1 = size$t1^0.5; 
plot(z1~z0,data=size); 

### Fit a skewed Normal ##################################################

par(mfrow=c(2,2)); 
rollmean <- rollapply(size$z0,width=50,mean,na.rm=T)
rollmu<-rollapply(size$z1,width=50,mean,na.rm=T)
rollsd <- rollapply(size$z1,width=50,sd,na.rm=T)
rollskew<-rollapply(size$z1,width=50,skewness,na.rm=T)
rollkurt<-rollapply(size$z1,width=50,kurtosis,na.rm=T)
plot(rollmean,rollmu,xlab='Size at time t',ylab='Mean size at size t+1',type='l'); 
plot(rollmu,log(rollsd),xlab='Size at time t',ylab='log(SD) size at size t+1',type='l'); 
plot(rollmu,rollskew,xlab='Size at time t',ylab='Skewness size at size t+1',type='l'); 
plot(rollmu,rollkurt/3 -1 ,xlab='Size at time t',ylab='Excess kurtosis size at size t+1',type='l'); 


NLLsn <- function(a,b,av,bv,aa,ba,ca) {
    mu = a+b*z0;
    sd = exp(av + bv*mu)
    alpha = aa + ba*mu + ca*mu^2; 
    -sum(dSHASHo(z1,mu=mu,sigma=sd,nu=alpha,tau=1,log=TRUE))
}    
z0 = size$z0; z1 = size$z1; 
fit0=lm(z1~z0); sigma0=sqrt(mean(fit0$residuals^2));
c0=as.vector(coef(fit0)); 
fitSN=mle2(NLLsn,start=list(a=c0[1],b=c0[2],av=log(sigma0),bv=0.001, aa=0,ba=-0.5, ca=0),
method="Nelder-Mead",skip.hessian=TRUE, control=list(trace=4,maxit=10000)); 
for(k in 1:25) {
    c1 = as.list(coef(fitSN)); 
    fitSN=mle2(NLLsn,start=c1, method="Nelder-Mead",skip.hessian=TRUE, control=list(trace=4,maxit=10000)); 
    c1 = as.list(coef(fitSN)); 
    cat("BFGS","\n"); 
    fitSN=mle2(NLLsn,start=c1, method="BFGS",skip.hessian=TRUE, control=list(REPORT=1,maxit=10000)); 
 }
c1 = as.list(coef(fitSN));
fitSN=mle2(NLLsn,start=c1, method="BFGS",skip.hessian=FALSE, control=list(REPORT=1,maxit=10000)); 
summary(fitSN); 
pars = round(coef(fitSN),digits=3); 
names(pars) <- c("a","b","av","bv","aa","ba","ca"); 
a = pars["a"]; b = pars["b"];  av = pars["av"]; bv = pars["bv"];  
aa = pars["aa"]; ba = pars["ba"];  ca = pars["ca"];  


###########################################################
#  Generate artificial data from fitted distribution 
########################################################### 
x = sample(size$z0,nrow(size),replace=FALSE); x=sort(x); nx = length(x); 
means = a + b*x; 
sds =  exp(av + bv*means); 
alphas = aa + ba*means  + ca*means^2;  
par(mfrow=c(3,1)); 
plot(x,means,type="l");  plot(x,log(sds),type="l");  plot(x,alphas,type="l");  

y = rSHASHo(length(x),means,sds,alphas,tau=1)
par(mfrow=c(2,1)); 
plot(x,y,type="p");

fsize = data.frame(z0=x,z1=y); 
par(mfrow=c(2,2)); 
rollmean <- rollapply(fsize$z0,width=50,mean,na.rm=T)
rollmu<-rollapply(fsize$z1,width=50,mean,na.rm=T)
rollsd <- rollapply(fsize$z1,width=50,sd,na.rm=T)
rollskew<-rollapply(fsize$z1,width=50,skewness,na.rm=T)
rollkurt<-rollapply(fsize$z1,width=50,kurtosis,na.rm=T)
plot(rollmean,rollmu,xlab='Size at time t',ylab='Mean size at size t+1',type='l'); 
plot(rollmu,log(rollsd),xlab='Size at time t',ylab='log(SD) size at size t+1',type='l'); 
plot(rollmu,rollskew,xlab='Size at time t',ylab='Skewness size at size t+1',type='l'); 
plot(rollmu,rollkurt/3 -1 ,xlab='Size at time t',ylab='Excess kurtosis size at size t+1',type='l'); 

breakpts = c(min(x),quantile(x,c(1/3,2/3)), max(x))
B=create.bspline.basis(rangeval=range(x), breaks=breakpts,  norder=3);
X = eval.basis(B,x);

f.50 <- rq(y ~ X-1, tau=0.5)
f.90 <- rq(y ~ X-1, tau=0.9)
f.10 <- rq(y ~ X-1, tau=0.1)
q.50 = fitted(f.50); q.90 = fitted(f.90); q.10 = fitted(f.10);
matpoints(x,cbind(q.10, q.50, q.90), type="l",lty=2, col="blue"); 
NPS = (q.10 + q.90 - 2*q.50)/(q.90 - q.10); 

boot.10=boot.rq(X,y,tau=0.1, R=501);
boot.50=boot.rq(X,y,tau=0.5,  R=501, U=boot.10$U);
boot.90=boot.rq(X,y,tau=0.9, R=501, U=boot.10$U);

bootq.10=X%*%t(boot.10$B);
bootq.50=X%*%t(boot.50$B);
bootq.90=X%*%t(boot.90$B);

bootNPS = (bootq.10 + bootq.90 - 2*bootq.50)/(bootq.90 - bootq.10); 
LS = apply(bootNPS,1,function(x) quantile(x,0.05)); 
MS = apply(bootNPS,1,function(x) quantile(x,0.5));
US = apply(bootNPS,1,function(x) quantile(x,0.95)); 
matplot(x,cbind(LS,MS,US),type="l",col="red",lty=2,lwd=2); 
points(x,NPS,type="l",col="black",lwd=2); 
abline(0,0,col="blue",lty=2); 

b0 = apply(bootNPS,2,mean); 
quantile(b0,c(0.01, 0.05,0.95,0.99)); 

b1 = apply(bootNPS,2,function(s) lm(s~x)$coef[2]);
quantile(b1,c(0.01, 0.05,0.95,0.99)); 

b2 = apply(bootNPS,2,function(s) lm(s~x + I(x^2))$coef[3]);
quantile(b2,c(0.01, 0.05,0.95,0.99)); 

require(mgcv); 
gam_fit = gam(list(y~x,~s(x), ~s(x), ~x), family="shash"); 
summary(gam_fit); 