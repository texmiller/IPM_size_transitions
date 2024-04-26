# purpose: assemble data figures and tables using data and model outputs 
# combining multiple case studies
library(xtable)

rm(list=ls(all=TRUE))
require(scales); 

tom = "C:/Users/tm9/Dropbox/github/IPM_size_transitions"
steve = "c:/repos/IPM_size_transitions" 
home = ifelse(Sys.info()["user"] == "Ellner", steve, tom)
setwd(home)

source("code/variance_diagnostics.R"); 

## read in species-specific outputs
lichen<-readRDS("lichen/lichen_out.rds")
cactus<-readRDS("cactus/cactus_out.rds")
load("cactus/CYIM_grow.Rdata") ### <-------------WARNING! different pilot fit. 
orchid<-readRDS("orchid/orchid_out.rds")

pike<-readRDS("pike/pike_out.rds")
creosote<-readRDS("creosote/creosote_out.rds")

alpha_val <- 0.15
graphics.off(); dev.new(height=9,width=7.5); 
par(mfrow=c(3,2),mar = c(3.5, 5, 2, 4), mgp=c(2,1,0), cex.axis=1,cex.lab=1.4, oma=c(0,0,0,2)) 

################### Lichen 
pilot = lichen$lichen_grow
plot(pilot$fitted,pilot$scaledResids,pch=1,col=alpha("black",alpha_val),xlab="Fitted values", ylab="Scaled residuals",
ylim = quantile(pilot$scaledResids,c(0.025,0.975))); 
mfit = rsq.smooth.spline(pilot$fitted,pilot$scaledResids);
points(mfit$x,mfit$yhat,type="l",col="red",lty=1,lwd=2);  
title("A   lichen",font=3,adj=0); title(round(sd(mfit$yhat),3)); 
plot(pilot$fitted,abs(pilot$scaledResids),pch=1,col=alpha("black",alpha_val),xlab="Fitted values", ylab="|Scaled residuals|",
ylim = quantile(abs(pilot$scaledResids),c(0,0.975))); 
vfit = rsq.smooth.spline(pilot$fitted,abs(pilot$scaledResids));
points(vfit$x,vfit$yhat,type="l",col="blue",lty=1,lwd=2);  
title("B",font=3,adj=0); title(round(sd(vfit$yhat),3));  


################### Cactus 
pilot = CYIM_grow; 
# pilot = cactus$cactus_grow; 
pilot$fitted = pilot$fitted_norfx; 
plot(pilot$fitted,pilot$scaledResids,pch=1,col=alpha("black",alpha_val),xlab="Fitted values", ylab="Scaled residuals",
ylim = quantile(pilot$scaledResids,c(0.025,0.975)));  
mfit = rsq.smooth.spline(pilot$fitted,pilot$scaledResids);
points(mfit$x,mfit$yhat,type="l",col="red",lty=1,lwd=2);  
title("C   cactus",font=3,adj=0); title(round(sd(mfit$yhat),3)); 
plot(pilot$fitted,abs(pilot$scaledResids)^0.5,pch=1,col=alpha("black",alpha_val),xlab="Fitted values", ylab="|Scaled residuals|",
ylim = quantile(abs(pilot$scaledResids),c(0,0.975)));  
vfit = rsq.smooth.spline(pilot$fitted,abs(pilot$scaledResids));
points(vfit$x,vfit$yhat^0.5,type="l",col="blue",lty=1,lwd=2);  
title("D",font=3,adj=0); title(round(sd(vfit$yhat),3));  

################### Orchid 
pilot = orchid$orchid_grow; 
pilot$fitted = pilot$GAU_fitted; pilot$scaledResids = pilot$GAU_scaled_resids;  
plot(pilot$fitted,pilot$scaledResids,pch=1,col=alpha("black",alpha_val),xlab="Fitted values", ylab="Scaled residuals",
ylim = quantile(pilot$scaledResids,c(0.025,0.975)));  
mfit = rsq.smooth.spline(pilot$fitted,pilot$scaledResids);
points(mfit$x,mfit$yhat,type="l",col="red",lty=1,lwd=2);  
title("E   orchid",font=3,adj=0); title(round(sd(mfit$yhat),3)); 
plot(pilot$fitted,abs(pilot$scaledResids)^0.5,pch=1,col=alpha("black",alpha_val),xlab="Fitted values", ylab="|Scaled residuals|",
ylim = quantile(abs(pilot$scaledResids),c(0,0.975)));  
vfit = rsq.smooth.spline(pilot$fitted,abs(pilot$scaledResids));
points(vfit$x,vfit$yhat^0.5,type="l",col="blue",lty=1,lwd=2);  
title("F",font=3,adj=0); title(round(sd(vfit$yhat),3)); 