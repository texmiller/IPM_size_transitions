# purpose: assemble data figures and tables using data and model outputs 
# combining multiple case studies


rm(list=ls(all=TRUE))
require(scales); library(xtable)
library(oce); library(lme4) 

tom = "C:/Users/tm9/Dropbox/github/IPM_size_transitions"
steve = "c:/repos/IPM_size_transitions" 
home = ifelse(Sys.info()["user"] == "Ellner", steve, tom)
setwd(home)

## read in species-specific outputs
lichen<-readRDS("lichen/lichen_out.rds")
cactus<-readRDS("cactus/cactus_out.rds")
orchid<-readRDS("orchid/orchid_out.rds")
pike<-readRDS("pike/pike_out.rds")
creosote<-readRDS("creosote/creosote_out.rds")

# residuals plot
alpha_val <- 0.15

pdf("manuscript/figures/combo_resid_diagnostics.pdf",height = 9, width = 7.5,useDingbats = F)
par(mfrow=c(3,2),mar = c(3.5, 5, 2, 4), oma=c(0,0,0,2)) 
## lichen
plot(lichen$lichen_grow$t0,lichen$lichen_grow$t1,pch=1,col=alpha("black",alpha_val),cex.axis=0.8,
     xlab=" ",ylab=" ")
mtext("Thallus area, time t", side = 1, line = 2,cex=0.7)
mtext("Thallus area, time t+1", side = 2, line = 2,cex=0.7)
lines(lichen$lichen_grow$t0,predict(lichen$lichen_GAU_best,type="response")[,1],col="red3",lwd=2)
par(new = TRUE)                           
plot(lichen$lichen_grow$t0,lichen$lichen_grow$fitted_sd,col="blue3",type="l",lwd=2,
     axes = FALSE, xlab = "", ylab = "",cex.axis=0.8)
axis(side = 4, at = pretty(range(lichen$lichen_grow$fitted_sd)),cex.axis=0.8)
mtext("Standard deviation", side = 4, line = 2,cex=0.7)
legend("topleft",legend=c("Fitted mean","Fitted sd"),bg="white",lwd=1.5,col=c("red3","blue3"),cex=0.8)
title("A - lichen",font=3,adj=0)

plot(lichen$lichen_grow$t0,lichen$lichen_grow$scaledResids,col=alpha("black",alpha_val),cex.axis=0.8,
     xlab=" ",ylab=" ")
mtext("Thallus area, time t", side = 1, line = 2,cex=0.7)
mtext("Scaled residuals of\nthallus area, time t+1", side = 2, line = 2,cex=0.7)
lines(lichen$lichen_grow$t0,lichen$q.05,col="black")
lines(lichen$lichen_grow$t0,lichen$q.10,col="black")
lines(lichen$lichen_grow$t0,lichen$q.25,col="black")
lines(lichen$lichen_grow$t0,lichen$q.50,col="black")
lines(lichen$lichen_grow$t0,lichen$q.75,col="black")
lines(lichen$lichen_grow$t0,lichen$q.90,col="black")
lines(lichen$lichen_grow$t0,lichen$q.95,col="black")
par(new = TRUE)                           
matplot(cbind(lichen$lichen_grow$t0,lichen$lichen_grow$t0),
        cbind(lichen$NPS_hat,lichen$NPK_hat), type="l",lwd=2,
        col=c("blue3","red3"), lty=1, axes = FALSE, xlab = "", ylab = "",cex.axis=0.8)
axis(side = 4, at = pretty(range(c(lichen$NPS_hat,lichen$NPK_hat))),cex.axis=0.8)
mtext("NP skewness or kurtosis", side = 4, line = 2,cex=0.7)
legend("topleft",legend=c("Skewness","Excess kurtosis"),bg="white",lwd=1.5,col=c("blue3","red3"),cex=0.8)
title("B",font=3,adj=0)

## cactus
plot(cactus$cactus_grow$logvol_t,cactus$cactus_grow$logvol_t1,pch=1,col=alpha("black",alpha_val),
     xlab=" ",ylab=" ")
mtext("log(volume), time t", side = 1, line = 2,cex=0.7)
mtext("log(volume), time t+1", side = 2, line = 2,cex=0.7)
lines(cactus$cactus_grow$logvol_t,cactus$cactus_grow$fitted_norfx,col="red3",lwd=2)
par(new = TRUE)                           
plot(cactus$cactus_grow$logvol_t,cactus$cactus_grow$fitted_sd,col="blue3",lwd=2,
     axes = FALSE, xlab = "", ylab = "",type="l")
axis(side = 4, at = pretty(range(cactus$cactus_grow$fitted_sd)))
mtext("Standard deviation", side = 4, line = 2,cex=0.7)
title("C - cactus",font=3,adj=0)

plot(cactus$cactus_grow$logvol_t,cactus$cactus_grow$scaledResids,col=alpha("black",alpha_val),
     xlab=" ",ylab=" ")
mtext("log(volume), time t", side = 1, line = 2,cex=0.7)
mtext("Scaled residuals of\nlog(volume), time t+1", side = 2, line = 2,cex=0.7)
lines(cactus$cactus_grow$logvol_t,cactus$q.05,col="black")
lines(cactus$cactus_grow$logvol_t,cactus$q.10,col="black")
lines(cactus$cactus_grow$logvol_t,cactus$q.25,col="black")
lines(cactus$cactus_grow$logvol_t,cactus$q.50,col="black")
lines(cactus$cactus_grow$logvol_t,cactus$q.75,col="black")
lines(cactus$cactus_grow$logvol_t,cactus$q.90,col="black")
lines(cactus$cactus_grow$logvol_t,cactus$q.95,col="black")
par(new = TRUE)                           
matplot(cbind(cactus$cactus_grow$logvol_t,cactus$cactus_grow$logvol_t),
     cbind(cactus$NPS_hat,cactus$NPK_hat),
     type="l",lwd=2,
     col=c("blue3","red3"), lty=1, axes = FALSE, xlab = "", ylab = "",cex.axis=0.8)
axis(side = 4, at = pretty(range(c(cactus$NPS_hat,cactus$NPK_hat))),cex.axis=0.8)
mtext("NP skewness or kurtosis", side = 4, line = 2,cex=0.7)
title("D",font=3,adj=0)

## orchid
plot(orchid$orchid_grow$log_area_t,orchid$orchid_grow$log_area_t1,type="n",
     xlab=" ",ylab=" ")
mtext("log(leaf area), time t", side = 1, line = 2,cex=0.7)
mtext("log(leaf area), time t+1", side = 2, line = 2,cex=0.7)
points(orchid$orchid_grow$log_area_t[orchid$orchid_grow$flowering==0],
       orchid$orchid_grow$log_area_t1[orchid$orchid_grow$flowering==0],col=alpha("black",alpha_val))
points(orchid$orchid_grow$log_area_t[orchid$orchid_grow$flowering==1],
       orchid$orchid_grow$log_area_t1[orchid$orchid_grow$flowering==1],col=alpha("pink",alpha_val))
lines(orchid$veg_size,fixef(orchid$orchid_GAU_best)[1]+fixef(orchid$orchid_GAU_best)[2]*orchid$veg_size,
      col=alpha("red3",0.75),lwd=2,lty=1)
lines(orchid$flow_size,fixef(orchid$orchid_GAU_best)[1]+fixef(orchid$orchid_GAU_best)[3]+(fixef(orchid$orchid_GAU_best)[2]+fixef(orchid$orchid_GAU_best)[4])*orchid$flow_size,
      col=alpha("red3",0.75),lwd=2,lty=2)
legend("topleft",legend=c("Vegetative, time t","Flowering, time t"),bg="white",pch=1,col=c("black","pink"),cex=0.8)
plotInset(2.5, 0.25, 7, 2,
          expr = plot(orchid$orchid_grow$GAU_fitted,orchid$orchid_grow$GAU_sd,
                      type="l", xlab = "E[log(leaf area), time t+1]",col="blue3",lwd=2,
                      ylab = "Std Dev",cex.lab=0.8,
                      cex.axis = 0.5, mgp = c(3/2, 1/2, 0)),
          mar = c(0, 3, 0, 0))
title("E - orchid",font=3,adj=0)

plot(orchid$orchid_grow$GAU_fitted,orchid$orchid_grow$GAU_scaled_resids,
     col=alpha("black",alpha_val),xlab=" ",ylab=" ")
mtext("E[log(leaf area), time t+1]", side = 1, line = 2,cex=0.7)
mtext("Scaled residuals of\nlog(leaf area), time t+1", side = 2, line = 2,cex=0.7)
lines(orchid$orchid_grow$GAU_fitted,orchid$q.05,col="black")
lines(orchid$orchid_grow$GAU_fitted,orchid$q.10,col="black")
lines(orchid$orchid_grow$GAU_fitted,orchid$q.25,col="black")
lines(orchid$orchid_grow$GAU_fitted,orchid$q.50,col="black")
lines(orchid$orchid_grow$GAU_fitted,orchid$q.75,col="black")
lines(orchid$orchid_grow$GAU_fitted,orchid$q.90,col="black")
lines(orchid$orchid_grow$GAU_fitted,orchid$q.95,col="black")
par(new = TRUE)      
matplot(cbind(orchid$orchid_grow$GAU_fitted,orchid$orchid_grow$GAU_fitted),
        cbind(orchid$NPS_hat,orchid$NPK_hat),
        type="l",lwd=2,
        col=c("blue3","red3"), lty=1, axes = FALSE, xlab = "", ylab = "",cex.axis=0.8)
axis(side = 4,cex.axis=0.8,at = pretty(range(c(orchid$NPS_hat,orchid$NPK_hat))))
mtext("NP skewness or kurtosis", side = 4, line = 2,cex=0.7)
title("F",font=3,adj=0)
dev.off()

