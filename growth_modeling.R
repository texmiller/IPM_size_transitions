## Model size transitions using data list prepared in dryad_data_prep.R
library(moments)
library(zoo)
library(mgcv)
library(scales)

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


# coral data --------------------------------------------------------------
coral <- growth_dat[["coral"]]
## there are sites, plots, and years here - but mostly sites are visited in separate years
table(coral$Plot,coral$Year,coral$Site)

## a simple linear regression
coral_lm <- lm(t1 ~ t0, data = coral)
plot(t1 ~ t0, data = coral)
abline(coef(coral_lm))

## simluate data with fitted model
layout.matrix <- matrix(c(1, 2, 3, 1, 4, 5), nrow = 3, ncol = 2)
layout(mat = layout.matrix, heights = c(2, 1, 1))
layout.show(5)
plot(density(coral$t1),type="n")
n_sim <- 500
coral_sim <- matrix(NA,dim(coral)[1],n_sim); coral_mean <- coral_sd <- coral_skew <- coral_kurt <- vector("numeric",n_sim)
for(i in 1:n_sim){
  coral_sim[,i] <- rnorm(n = nrow(coral_sim), mean = coef(coral_lm)[1] + coef(coral_lm)[2] * coral$t0, sd = sigma(coral_lm))
  lines(density(coral_sim[,i]),col=alpha("gray",0.5))
  coral_mean[i] <- mean(coral_sim[,i])
  coral_sd[i] <- sd(coral_sim[,i])
  coral_skew[i] <- skewness(coral_sim[,i])
  coral_kurt[i] <- kurtosis(coral_sim[,i])
}
lines(density(coral$t1),lwd=4)
plot(density(coral_mean),lwd=2,col="gray",xlim=c(min(c(coral_mean,mean(coral$t1))),max(c(coral_mean,mean(coral$t1)))))
abline(v=mean(coral$t1),lwd=4)
plot(density(coral_sd),lwd=2,col="gray",xlim=c(min(c(coral_sd,sd(coral$t1))),max(c(coral_sd,sd(coral$t1)))))
abline(v=sd(coral$t1),lwd=4)
plot(density(coral_skew),lwd=2,col="gray",xlim=c(min(c(coral_skew,skewness(coral$t1))),max(c(coral_skew,skewness(coral$t1)))))
abline(v=skewness(coral$t1),lwd=4)
plot(density(coral_kurt),lwd=2,col="gray",xlim=c(min(c(coral_kurt,kurtosis(coral$t1))),max(c(coral_kurt,kurtosis(coral$t1)))))
abline(v=kurtosis(coral$t1),lwd=4)

## next see how well Peterson et al. did with beta regression