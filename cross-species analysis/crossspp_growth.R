# purpose: assemble data figures and tables using data and model outputs 
# combining multiple case studies
library(xtable)

rm(list=ls(all=TRUE))
require(scales); 

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
title("A",font=3,adj=0)

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
title("C",font=3,adj=0)

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
title("E",font=3,adj=0)

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


### Table of life history traits from matrices
source("./code/metaluck_fns_CMH.R")

##lichen
lichenGAU_mat <- lichen$mat_GAU$K; lichenJSU_mat <- lichen$mat_JSU$K
lichenGAU_matU <- lichen$mat_GAU$P; lichenJSU_matU <- lichen$mat_JSU$P
lichenGAU_matF <- lichen$mat_GAU$F; lichenJSU_matF <- lichen$mat_JSU$F
lichen_c0 = rep(0,nrow(lichenGAU_matU)); lichen_c0[1]=1
lichenGAU_traits<-c(
  Re(eigen(lichenGAU_mat)$values[1]),
  mean_lifespan(lichenGAU_matU, mixdist=lichen_c0),
  mean_LRO(lichenGAU_matU,lichenGAU_matF,mixdist=lichen_c0),
  #var_LRO_mcr(lichenGAU_matU,lichenGAU_matF,mixdist=lichen_c0)^0.5,
  mean_age_repro(lichenGAU_matU,lichenGAU_matF,mixdist=lichen_c0),
  gen_time_mu1_v(lichenGAU_matU,lichenGAU_matF)
)
lichenJSU_traits<-c(
  Re(eigen(lichenJSU_mat)$values[1]),
  mean_lifespan(lichenJSU_matU, mixdist=lichen_c0),
  mean_LRO(lichenJSU_matU,lichenJSU_matF,mixdist=lichen_c0),
  #var_LRO_mcr(lichenJSU_matU,lichenJSU_matF,mixdist=lichen_c0)^0.5,
  mean_age_repro(lichenJSU_matU,lichenJSU_matF,mixdist=lichen_c0),
  gen_time_mu1_v(lichenJSU_matU,lichenJSU_matF)
)

##cactus
cactusGAU_mat <- cactus$mat_GAU$IPMmat; cactusSHASH_mat <- cactus$mat_SHASH$IPMmat
cactusGAU_matU <- cactus$mat_GAU$Tmat; cactusSHASH_matU <- cactus$mat_SHASH$Tmat
cactusGAU_matF <- cactus$mat_GAU$Fmat; cactusSHASH_matF <- cactus$mat_SHASH$Fmat
cactus_c0 = rep(0,nrow(cactusGAU_matU)); cactus_c0[1]=1
cactusGAU_traits<-c(
  Re(eigen(cactusGAU_mat)$values[1]),
  mean_lifespan(cactusGAU_matU, mixdist=cactus_c0),
  mean_LRO(cactusGAU_matU,cactusGAU_matF,mixdist=cactus_c0),
  #var_LRO_mcr(cactusGAU_matU,cactusGAU_matF,mixdist=cactus_c0)^0.5,
  mean_age_repro(cactusGAU_matU,cactusGAU_matF,mixdist=cactus_c0),
  gen_time_mu1_v(cactusGAU_matU,cactusGAU_matF)
)
cactusSHASH_traits<-c(
  Re(eigen(cactusSHASH_mat)$values[1]),
  mean_lifespan(cactusSHASH_matU, mixdist=cactus_c0),
  mean_LRO(cactusSHASH_matU,cactusSHASH_matF,mixdist=cactus_c0),
  #var_LRO_mcr(cactusSHASH_matU,cactusSHASH_matF,mixdist=cactus_c0)^0.5,
  mean_age_repro(cactusSHASH_matU,cactusSHASH_matF,mixdist=cactus_c0),
  gen_time_mu1_v(cactusSHASH_matU,cactusSHASH_matF)
)

##orchid
orchidGAU_mat <- orchid$mat_GAU$matrix; orchidSST_mat <- orchid$mat_SST$matrix
orchidGAU_matU <- orchid$mat_GAU$Tmatrix; orchidSST_matU <- orchid$mat_SST$Tmatrix
orchidGAU_matF <- orchid$mat_GAU$Fmatrix; orchidSST_matF <- orchid$mat_SST$Fmatrix
orchid_c0 = rep(0,nrow(orchidGAU_matU)); orchid_c0[1]=1
orchidGAU_traits<-c(
  Re(eigen(orchidGAU_mat)$values[1]),
  mean_lifespan(orchidGAU_matU, mixdist=orchid_c0),
  mean_LRO(orchidGAU_matU,orchidGAU_matF,mixdist=orchid_c0),
  #var_LRO_mcr(orchidGAU_matU,orchidGAU_matF,mixdist=orchid_c0)^0.5,
  mean_age_repro(orchidGAU_matU,orchidGAU_matF,mixdist=orchid_c0),
  gen_time_mu1_v(orchidGAU_matU,orchidGAU_matF)
)
orchidSST_traits<-c(
  Re(eigen(orchidSST_mat)$values[1]),
  mean_lifespan(orchidSST_matU, mixdist=orchid_c0),
  mean_LRO(orchidSST_matU,orchidSST_matF,mixdist=orchid_c0),
  #var_LRO_mcr(orchidSST_matU,orchidSST_matF,mixdist=orchid_c0)^0.5,
  mean_age_repro(orchidSST_matU,orchidSST_matF,mixdist=orchid_c0),
  gen_time_mu1_v(orchidSST_matU,orchidSST_matF)
)

##pike
pikeGAU_mat <- pike$GAU_IPM$IPMmat; pikeSHASH_mat <- pike$SHASH_IPM$IPMmat
pikeGAU_matU <- pike$GAU_IPM$Pmat; pikeSHASH_matU <- pike$SHASH_IPM$Pmat
pikeGAU_matF <- pike$GAU_IPM$Fmat; pikeSHASH_matF <- pike$SHASH_IPM$Fmat
pike_c0 = rep(0,nrow(pikeGAU_matU)); pike_c0[1]=1
pikeGAU_traits<-c(
  Re(eigen(pikeGAU_mat)$values[1]),
  mean_lifespan(pikeGAU_matU, mixdist=pike_c0),
  mean_LRO(pikeGAU_matU,pikeGAU_matF,mixdist=pike_c0),
  #var_LRO_mcr(pikeGAU_matU,pikeGAU_matF,mixdist=pike_c0)^0.5,
  mean_age_repro(pikeGAU_matU,pikeGAU_matF,mixdist=pike_c0),
  gen_time_mu1_v(pikeGAU_matU,pikeGAU_matF)
)
pikeSHASH_traits<-c(
  Re(eigen(pikeSHASH_mat)$values[1]),
  mean_lifespan(pikeSHASH_matU, mixdist=pike_c0),
  mean_LRO(pikeSHASH_matU,pikeSHASH_matF,mixdist=pike_c0),
  #var_LRO_mcr(pikeSHASH_matU,pikeSHASH_matF,mixdist=pike_c0)^0.5,
  mean_age_repro(pikeSHASH_matU,pikeSHASH_matF,mixdist=pike_c0),
  gen_time_mu1_v(pikeSHASH_matU,pikeSHASH_matF)
)

##creosote
creosoteGAU_mat <- creosote$mat_GAU$IPMmat; creosoteJSU_mat <- creosote$mat_JSU$IPMmat
creosoteGAU_matU <- creosote$mat_GAU$Pmat; creosoteJSU_matU <- creosote$mat_JSU$Pmat
creosoteGAU_matF <- creosote$mat_GAU$Fmat; creosoteJSU_matF <- creosote$mat_JSU$Fmat
creosote_c0 = rep(0,nrow(creosoteGAU_matU)); creosote_c0[1]=1
creosoteGAU_traits<-c(
  Re(eigen(creosoteGAU_mat)$values[1]),
  mean_lifespan(creosoteGAU_matU, mixdist=creosote_c0),
  mean_LRO(creosoteGAU_matU,creosoteGAU_matF,mixdist=creosote_c0),
  #var_LRO_mcr(creosoteGAU_matU,creosoteGAU_matF,mixdist=creosote_c0)^0.5,
  mean_age_repro(creosoteGAU_matU,creosoteGAU_matF,mixdist=creosote_c0),
  gen_time_mu1_v(creosoteGAU_matU,creosoteGAU_matF)
)
creosoteJSU_traits<-c(
  Re(eigen(creosoteJSU_mat)$values[1]),
  mean_lifespan(creosoteJSU_matU, mixdist=creosote_c0),
  mean_LRO(creosoteJSU_matU,creosoteJSU_matF,mixdist=creosote_c0),
  #var_LRO_mcr(creosoteJSU_matU,creosoteJSU_matF,mixdist=creosote_c0)^0.5,
  mean_age_repro(creosoteJSU_matU,creosoteJSU_matF,mixdist=creosote_c0),
  gen_time_mu1_v(creosoteJSU_matU,creosoteJSU_matF)
)

## compile species into table
traits_table<-round(as.matrix(rbind(
                      lichenGAU_traits,lichenJSU_traits,
                      cactusGAU_traits,cactusSHASH_traits,
                      orchidGAU_traits,orchidSST_traits,
                      pikeGAU_traits,pikeSHASH_traits,
                      creosoteGAU_traits,creosoteJSU_traits)),3) 
spp_names<-c("Lichen","(Vulpicida pinastri)",
             "Cactus","(Cylindriopunia imbricata)",
             "Orchid","(Orchis purpurea)",
             "Pike","(Esox Lucius)",
             "Creosote","(Larrea tridentata)")
out_table<-data.frame(cbind(spp_names,rep(c("Gaussian","Improved"),times=length(spp_names)/2),
               traits_table),row.names=NULL)
names(out_table)<-c("Species","Growth model","lambda","Lifespan",
                    "Lifetime reproductive output","Age at reproduction","Generation time")

## print latex code then paste into ms
bold <- function(x) {paste('{\\textbf{',x,'}}', sep ='')}
print(xtable(out_table,align="rrp{1.5cm}|p{1.5cm}p{1.5cm}p{1.5cm}p{1.5cm}p{1.5cm}"),
      size="\\fontsize{9pt}{10pt}\\selectfont",include.rownames=FALSE,
      hline.after=c(0,2,4,6,8),
      sanitize.colnames.function=bold)

