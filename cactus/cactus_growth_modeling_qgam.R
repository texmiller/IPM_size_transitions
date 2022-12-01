## Returning to cactus growth, now estimating skew and kurtosis by quantile regression

## setwd
setwd("C:/Users/tm9/Dropbox/github/IPM_size_transitions")

## load libraries
library(tidyverse)
library(mgcv)
library(scales)
library(qgam)

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
CYIM_grow_m1 <- gam(list(logvol_t1 ~ s(logvol_t) + s(plot,bs="re") + s(year_t,bs="re"), ~s(logvol_t,k=4)), 
                    data=CYIM_grow, gamma=2, family=gaulss())
CYIM_gam_pred <- predict(CYIM_grow_m1,type="response",exclude=c("s(plot)","s(year_t)"))

## visualize mean and variance fit
par(mar = c(5, 4, 4, 4) + 0.3) 
plot(CYIM_grow$logvol_t,CYIM_grow$logvol_t1,pch=1,col=alpha("black",0.25))
points(CYIM_grow$logvol_t,CYIM_gam_pred[,1],col=alpha("red",0.25),pch=16,cex=.5)
par(new = TRUE)                           
plot(CYIM_grow$logvol_t,1/CYIM_gam_pred[,2],col=alpha("blue",0.25),pch=16,cex=.5,
     axes = FALSE, xlab = "", ylab = "")
axis(side = 4, at = pretty(range(1/CYIM_gam_pred[,2])))
mtext("sigma", side = 4, line = 3)

## inspect residuals, scaled by sd; re-run predict now w/RFX
fitted_sd<-1/predict(CYIM_grow_m1,type="response")[,2]
CYIM_grow$scaledResids=residuals(CYIM_grow_m1,type="response")/fitted_sd

plot(CYIM_grow$logvol_t,CYIM_grow$scaledResids,col=alpha("black",0.25),
     xlab="Size t",ylab="Scaled residuals")

## fit qgam
S.10<-qgam(scaledResids~s(logvol_t), data=CYIM_grow,qu=0.1,argGam=list(gamma=2)) 
S.50<-qgam(scaledResids~s(logvol_t), data=CYIM_grow,qu=0.5,argGam=list(gamma=2)) 
S.90<-qgam(scaledResids~s(logvol_t), data=CYIM_grow,qu=0.9,argGam=list(gamma=2)) 
q.10<-predict(S.10);q.50<-predict(S.50);q.90<-predict(S.90)
NPS_hat = (q.10 + q.90 - 2*q.50)/(q.90 - q.10)

lines(CYIM_grow$logvol_t,q.10,col=alpha("red",0.5))
lines(CYIM_grow$logvol_t,q.50,col=alpha("red",0.5))
lines(CYIM_grow$logvol_t,q.90,col=alpha("red",0.5))

par(new = TRUE)                           
plot(CYIM_grow$logvol_t,NPS_hat,col=alpha("blue",0.25),pch=16,cex=.5,
     axes = FALSE, xlab = "", ylab = "")
axis(side = 4, at = pretty(range(NPS_hat)))
mtext("NP Skewness", side = 4, line = 3)



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
