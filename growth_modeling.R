## Model size transitions using data list prepared in dryad_data_prep.R
library(moments)
library(zoo)
library(mgcv)

home = "c:/repos" ## edit as needed 
setwd(home); setwd("IPM_size_transitions"); 

## load data
growth_dat <- readRDS("growth_dat.rds")

## Steve's diagnostics code
### Scatterplot smoothing functions
spline.scatter.smooth=function(x,y,...) {
  fit=gam(y~s(x),gamma=10,method="REML")
  plot(x,y,type="p",...);
  out=predict(fit,type="response"); 
  points(x,out,type="l",lwd=2)
  #fit2 = lm(y~x+I(x^2));
  #points(x,fit2$fitted,type="l",col="red",lwd=2,lty=2);
}

par(mfrow=c(4,4),bty="l",mar=c(4,4,1,1),mgp=c(2,1,0),cex.axis=1.3,cex.lab=1.3);
for(i in 1:4){
z0 = growth_dat[[i]]$t0; z1=growth_dat[[i]]$t1; 
e = order(z0); z0=z0[e]; z1 = z1[e]; 

rollmean<-rollapply(z0,width=20,mean,na.rm=T)
rollmean1<-rollapply(z1,width=20,mean,na.rm=T)
rollsd <-rollapply(z1,width=20,sd,na.rm=T); 
rollskew<-rollapply(z1,width=20,skewness,na.rm=T)
rollkurt<-rollapply(z1,width=20,kurtosis,na.rm=T)/3 - 1
spline.scatter.smooth(rollmean,rollmean1,col="grey50",xlab="z0",ylab="z1");title(main=names(growth_dat[i])) 
spline.scatter.smooth(rollmean,rollsd,col="grey50",xlab="z0",ylab="SD"); 
spline.scatter.smooth(rollmean,rollskew,col="grey50",xlab="z0",ylab="Skew"); 
spline.scatter.smooth(rollmean,rollkurt,col="grey50",xlab="z0",ylab="Excess Kurtosis"); 

#scatter.smooth(rollmean,rollmean1,col="grey50",xlab="z0",ylab="z1",degree=2);title(main=names(growth_dat[i])) 
#scatter.smooth(rollmean,rollsd,col="grey50",xlab="z0",ylab="SD",degree=2); 
#scatter.smooth(rollmean,rollskew,col="grey50",xlab="z0",ylab="Skew",degree=2); 
#scatter.smooth(rollmean,rollkurt,col="grey50",xlab="z0",ylab="Excess Kurtosis",degree=2); 
}
