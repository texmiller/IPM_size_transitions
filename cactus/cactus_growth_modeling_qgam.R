## Returning to cactus growth, now estimating skew and kurtosis by quantile regression

## setwd
setwd("C:/Users/tm9/Dropbox/github/IPM_size_transitions")

## load libraries
library(tidyverse)
library(mgcv)
library(scales)
library(qgam)
library(gamlss.dist)

# function for converting cactus size measurements to volume
volume <- function(h, w, p){
  (1/3)*pi*h*(((w + p)/2)/2)^2
}

## read in cactus demography data
## these data are published on EDI: https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-sev.323.1
CYIM_full<-read_csv("cactus/cholla_demography_20042018_EDI.csv")%>% 
  ## filter out transplants and drop seed addition plots (which start with letter H)
  ## also drop Year_t==2018 because I don't have 2019 size data (not entered yet). 2018 data still included in 2017-2018 transition.
  filter(Transplant == 0,
         str_sub(Plot,1,1)!="H",
         Year_t!=2018) %>% 
  ## convert height, max width, perp width to volume of cone, take natural log
  mutate(vol_t = volume(Height_t,Width_t,Perp_t),
         vol_t1 = volume(Height_t1,Width_t1,Perp_t1),
         plot = as.factor(Plot),
         year_t = as.factor(Year_t),
         ID = interaction(TagID,plot)) %>%
  #select(ID,year_t,plot,vol_t,vol_t1,Survival_t1,Goodbuds_t1) %>% 
  ## sort by initial size
  arrange(vol_t) 


## In prelim analysis I inspected several unrealistic size transitions
## this file identifies plants to drop
CYIM_outliers<-read_csv("cactus/CYIM_outliers.csv") %>% 
  filter(FLAG==1) %>% select(ID) %>% unique()

## pull out and na.omit size transitions for growth modeling
## NAs come from new recruits (missing year_t size) and mortality (missing year_t1 size)
CYIM_full %>% 
  filter(!ID%in%CYIM_outliers$ID) %>% 
  mutate(logvol_t=log(vol_t),
         logvol_t1=log(vol_t1)) %>% 
  select(ID,year_t,plot,logvol_t,logvol_t1) %>% 
  ## drop rows with NAs
  drop_na() -> CYIM_grow

## use gam to fit Gaussian growth model with non-constant variance
CYIM_grow_m1 <- gam(list(logvol_t1 ~ s(logvol_t,k=4) + s(plot,bs="re") + s(year_t,bs="re"), ~s(logvol_t,k=4)), 
                    data=CYIM_grow, family=gaulss())
CYIM_gam_pred <- predict(CYIM_grow_m1,type="response",exclude=c("s(plot)","s(year_t)"))

## visualize mean and variance fit
par(mar = c(5, 4, 4, 4) + 0.3) 
plot(CYIM_grow$logvol_t,CYIM_grow$logvol_t1,pch=1,col=alpha("black",0.25))
points(CYIM_grow$logvol_t,CYIM_gam_pred[,1],col=alpha("red",0.25),pch=16,cex=.5)
par(new = TRUE)                           
plot(CYIM_grow$logvol_t,1/CYIM_gam_pred[,2],col=alpha("blue",0.25),pch=16,cex=.5,
     axes = FALSE, xlab = "", ylab = "")
axis(side = 4, at = pretty(range(1/CYIM_gam_pred[,2])))
mtext("sigma", side = 4, line = 3,col="blue")

## inspect residuals, scaled by sd; re-run predict now w/RFX
fitted_sd<-1/predict(CYIM_grow_m1,type="response")[,2]
CYIM_grow$scaledResids=residuals(CYIM_grow_m1,type="response")/fitted_sd

## fit qgam -- we will need several quantiles for skewness and kurtosis
gamma_param<-2
S.05<-qgam(scaledResids~s(logvol_t,k=4), data=CYIM_grow,qu=0.05)#,argGam=list(gamma=gamma_param)) 
S.10<-qgam(scaledResids~s(logvol_t,k=4), data=CYIM_grow,qu=0.1)#,argGam=list(gamma=gamma_param)) 
S.25<-qgam(scaledResids~s(logvol_t,k=4), data=CYIM_grow,qu=0.25)#,argGam=list(gamma=gamma_param)) 
S.50<-qgam(scaledResids~s(logvol_t,k=4), data=CYIM_grow,qu=0.5)#,argGam=list(gamma=gamma_param)) 
S.75<-qgam(scaledResids~s(logvol_t,k=4), data=CYIM_grow,qu=0.75)#,argGam=list(gamma=gamma_param)) 
S.90<-qgam(scaledResids~s(logvol_t,k=4), data=CYIM_grow,qu=0.9)#,argGam=list(gamma=gamma_param)) 
S.95<-qgam(scaledResids~s(logvol_t,k=4), data=CYIM_grow,qu=0.95)#,argGam=list(gamma=gamma_param)) 

## NP skewness
q.10<-predict(S.10);q.50<-predict(S.50);q.90<-predict(S.90)
NPS_hat = (q.10 + q.90 - 2*q.50)/(q.90 - q.10)

## NP kurtosis (relative to Gaussian)
q.05<-predict(S.05);q.25<-predict(S.25);q.75<-predict(S.75);q.95<-predict(S.95)
qN = qnorm(c(0.05,0.25,0.75,0.95))
KG = (qN[4]-qN[1])/(qN[3]-qN[2])
NPK_hat = ((q.95-q.05)/(q.75-q.25))/KG - 1

par(oma=c(0,0,0,4))
plot(CYIM_grow$logvol_t,CYIM_grow$scaledResids,col=alpha("black",0.25),
     xlab="Size t",ylab="Scaled residuals")
par(new = TRUE)                           
plot(CYIM_grow$logvol_t,NPS_hat,col=alpha("blue",0.25),pch=16,cex=.5,
     axes = FALSE, xlab = "", ylab = "")
axis(side = 4, at = pretty(range(c(NPS_hat,NPK_hat))))
mtext("NP Skewness", side = 4, line = 2,col="blue")
par(new = TRUE)                           
plot(CYIM_grow$logvol_t,NPK_hat,col=alpha("red",0.25),pch=16,cex=.5,
     axes = FALSE, xlab = "", ylab = "")
mtext("NP Excess Kurtosis", side = 4, line =3,col="red")

##I would like to see the actual qgams
plot(CYIM_grow$logvol_t,CYIM_grow$scaledResids,col=alpha("black",0.25),
     xlab="Size t",ylab="Scaled residuals")
points(CYIM_grow$logvol_t,q.05,col="red",pch=".")
points(CYIM_grow$logvol_t,q.10,col="red",pch=".")
points(CYIM_grow$logvol_t,q.25,col="red",pch=".")
points(CYIM_grow$logvol_t,q.50,col="red",pch=".")
points(CYIM_grow$logvol_t,q.75,col="red",pch=".")
points(CYIM_grow$logvol_t,q.90,col="red",pch=".")
points(CYIM_grow$logvol_t,q.95,col="red",pch=".")

## put it all together:
par(mar = c(5, 4, 2, 3), mfrow=c(1,2), oma=c(0,0,0,4)) 
plot(CYIM_grow$logvol_t,CYIM_grow$logvol_t1,pch=1,col=alpha("black",0.25))
points(CYIM_grow$logvol_t,CYIM_gam_pred[,1],col=alpha("red",0.25),pch=16,cex=.5)
par(new = TRUE)                           
plot(CYIM_grow$logvol_t,1/CYIM_gam_pred[,2],col=alpha("blue",0.25),pch=16,cex=.5,
     axes = FALSE, xlab = "", ylab = "")
axis(side = 4, at = pretty(range(1/CYIM_gam_pred[,2])))
mtext("sigma", side = 4, line = 2,col="blue")

plot(CYIM_grow$logvol_t,CYIM_grow$scaledResids,col=alpha("black",0.25),
     xlab="Size t",ylab="Scaled residuals")
par(new = TRUE)                           
plot(CYIM_grow$logvol_t,NPS_hat,col=alpha("blue",0.25),pch=16,cex=.5,
     axes = FALSE, xlab = "", ylab = "")
axis(side = 4, at = pretty(range(c(NPS_hat,NPK_hat))))
mtext("NP Skewness", side = 4, line = 2,col="blue")
par(new = TRUE)                           
plot(CYIM_grow$logvol_t,NPK_hat,col=alpha("red",0.25),pch=16,cex=.5,
     axes = FALSE, xlab = "", ylab = "")
mtext("NP Excess Kurtosis", side = 4, line =3,col="red")

## now I need to fit a distribution with negative skew and positive and negative excess kurtosis
## both skewness and kurtosis should be non-monotonic wrt size
## turns out mgcv can do this!
CYIM_gam_shash <- gam(list(logvol_t1 ~ s(logvol_t,k=4) + s(plot,bs="re") + s(year_t,bs="re"), # <- model for location 
                           ~ s(logvol_t,k=4),   # <- model for log-scale
                           ~ s(logvol_t,k=4),   # <- model for skewness
                           ~ s(logvol_t,k=4)), # <- model for log-kurtosis
                      data = CYIM_grow, 
                      family = shash,  
                      optimizer = "efs")
CYIM_shash_pred <- predict(CYIM_gam_shash,type="response",exclude=c("s(plot)","s(year_t)"))

## view parameter estimates
par(mfrow=c(2,2))
plot(CYIM_grow$logvol_t,CYIM_grow$logvol_t1,pch=1,col=alpha("black",0.25))
points(CYIM_grow$logvol_t,CYIM_shash_pred[,1],col=alpha("red",0.25),pch=16,cex=.5)
plot(CYIM_grow$logvol_t,exp(CYIM_shash_pred[,2]),col=alpha("red",0.25),pch=16,cex=.5)
plot(CYIM_grow$logvol_t,CYIM_shash_pred[,3],col=alpha("red",0.25),pch=16,cex=.5)
plot(CYIM_grow$logvol_t,exp(CYIM_shash_pred[,4]),col=alpha("red",0.25),pch=16,cex=.5)

## simulate data from fitted model and compare to real data
sim_dat <- rSHASHo2(n=nrow(CYIM_grow),
         mu=CYIM_shash_pred[,1],
         sigma=exp(CYIM_shash_pred[,2]),
         nu=CYIM_shash_pred[,3],
         tau=exp(CYIM_shash_pred[,4]))


plot(CYIM_grow$logvol_t,sim_dat)
points(CYIM_grow$logvol_t,CYIM_grow$logvol_t1,col="red")

# the basement ------------------------------------------------------------

## test that I understand how to get sigma from the gaulss fit
x <- runif(500,-1,1)
y <- rnorm(500,mean=8.6+3*x,sd=exp(0.5+0.8*x))
plot(x,y)
testgam<-gam(list(y~s(x),~s(x)),family=gaulss())
pred.response <- predict(testgam,type="response")
plot(exp(0.5+0.8*x),1/pred.response[,2]);abline(0,1)

pred.lpmat <- predict(testgam,type="lpmatrix")
grow_sd_index <- which(as.factor(names(coef(testgam)))=="(Intercept).1") ## this is where the sd coefficients start
gam_coef_length <- length(coef(testgam))
plot(exp(0.5+0.8*x),
exp(pred.lpmat[, grow_sd_index:length(coef(testgam))] %*% coef(testgam)[grow_sd_index:length(coef(testgam))]));abline(0,1)


## now curious to look at the outliers
CYIM_grow %>% 
  filter(scaledResids < quantile(scaledResids,probs=0.025) |
           scaledResids >  quantile(scaledResids,probs=0.975) ) %>% 
  select(ID,year_t) %>% 
  left_join(.,CYIM_full,
            by=c("ID","year_t"))-> outliers
write_csv(outliers,"cactus/CYIM_outliers.csv")

separate(as.character(CYIM_grow$ID))

df <- data.frame(x = c(NA, "x.y", "x.z", "y.z"))
df %>% separate(x, c("A", "B"))

CYIM_test<-read_csv("cactus/cholla_demography_20042018_EDI.csv")
CYIM_test %>% 
  filter(Plot=="3",TagID=="45") %>% View()
