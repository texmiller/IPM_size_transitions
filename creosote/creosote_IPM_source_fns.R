## Source functions for creosote IPM, using gams fit in creosote_growth_modeling_gammit.R
## These functions are size (x) and density (d), with no other spatial or temporal random effects
## d is the weighted density of the window (sum of log volume)
setwd("C:/Users/tm9/Desktop/git local/IPM_size_transitions/creosote")

library(mgcv)
library(popbio)
library(plot.matrix)

## Inv logit function for transformations
invlogit <- function(x){exp(x)/(1+exp(x))}

## Lower and upper integration limits (size bounds)
LATR_size_bounds <- readRDS("LATR_size_bounds.rds")
size_dummy <- seq(LATR_size_bounds$min_size,LATR_size_bounds$max_size,0.1)

## Growth -- Gaussian using best gam
LATR_grow_best <- readRDS("LATR_grow_best.rds")
growth_NO <- function(x,y,d){
  xb=pmin(pmax(x,LATR_size_bounds$min_size),LATR_size_bounds$max_size)
  lpmat <- predict.gam(LATR_grow_best,
                      newdata = data.frame(
                        weighted.dens = d,
                        log_volume_t = xb,
                        unique.transect="1.FPS"),
                      type="lpmatrix",
                      exclude = "s(unique.transect)")
  ## linear predictor for mean and log sigma -- need to update so these indices are not hard-coded but for now they work
  grow_mu <- lpmat[,1:31]%*%coef(LATR_grow_best)[1:31]
  grow_sigma <- exp(lpmat[,32:50]%*%coef(LATR_grow_best)[32:50])
  return(dnorm(y,mean=grow_mu,sd=grow_sigma))
}

growth_SGT <- function(x, y, d){
  xb = pmin(pmax(x, LATR_size_bounds$min_size), LATR_size_bounds$max_size)
  lpmat <- predict.gam(LATR_grow_best,
                       newdata = data.frame(weighted.dens = d, log_volume_t = xb, unique.transect = "1.FPS"),
                       type = "lpmatrix",
                       exclude = "s(unique.transect)")
  # Linear predictor for mean, sigma, and lambda
  grow_mu <- lpmat[, 1:31] %*% coef_grow_best[1:31]
  grow_sigma <- exp(lpmat[, 32:50] %*% coef_grow_best[32:50])
  grow_lambda <- -invlogit(coef_grow_best[51]+coef_grow_best[52]*xb)
  return(dsgt(x = y, 
              mu=grow_mu,
              sigma=grow_sigma,
              lambda=grow_lambda,
              p=exp(coef_grow_best[53]),
              q=exp(coef_grow_best[54]),
              mean.cent=T,
              var.adj=T))}

## Survival-- prediction from naturally occuring plants
LATR_surv_best <- readRDS("LATR_surv_best.rds")
survival_fn <- function(x,d){
  xb=pmin(pmax(x,LATR_size_bounds$min_size),LATR_size_bounds$max_size)
  lpmat <- predict.gam(LATR_surv_best,
              newdata = data.frame(
                weighted.dens = d,
                log_volume_t = xb,
                transplant=F,
                unique.transect="1.FPS"),
              type="lpmatrix",
              exclude = "s(unique.transect)")
  pred <- lpmat[,1:21]%*%coef(LATR_surv_best)[1:21]
  return(invlogit(pred))
}

## combined growth and survival
pxy <- function(x,y,d,dist="SGT"){
  if(dist=="SGT"){return(survival_fn(x,d) * growth_SGT(x,y,d))}
  if(dist=="NO"){return(survival_fn(x,d) * growth_NO(x,y,d))}
}

## Flowering
LATR_flower_best <- readRDS("LATR_flower_best.rds")
flower_fn <- function(x,d){
  xb=pmin(pmax(x,LATR_size_bounds$min_size),LATR_size_bounds$max_size)
  lpmat <- predict.gam(LATR_flower_best,
              newdata = data.frame(
                weighted.dens = d,
                log_volume_t = xb,
                unique.transect="1.FPS"),
              type="lpmatrix",
              exclude = "s(unique.transect)")
  pred <- lpmat[,1:20]%*%coef(LATR_flower_best)[1:20]
  return(invlogit(pred))
}

## seed production (fruits * seeds/fruit)
LATR_fruits_best <- readRDS("LATR_fruits_best.rds")
seeds_fn <- function(x,d,seeds.per.fruit=6){
  xb=pmin(pmax(x,LATR_size_bounds$min_size),LATR_size_bounds$max_size)
  lpmat <- predict.gam(LATR_fruits_best,
                      newdata = data.frame(
                        weighted.dens = d,
                        log_volume_t = xb,
                        unique.transect="1.FPS"),
                      type="lpmatrix",
                      exclude = "s(unique.transect)")
  pred <- lpmat[,1:19]%*%coef(LATR_fruits_best)[1:19]
  return(exp(pred)*seeds.per.fruit)
}

## Seed-to-Seedling recruitment probability
LATR_recruit_best <- readRDS("LATR_recruit_best.rds")
recruitment_fn <- function(d){
  lpmat <- predict.gam(LATR_recruit_best,
                      newdata = data.frame(
                        weighted.dens = d,
                        unique.transect="1.FPS"),
                      type="lpmatrix",
                      exclude = "s(unique.transect)")
  pred <- lpmat[,1]%*%coef(LATR_recruit_best)[1]
  return(invlogit(pred[[1]]))
}

## Recruit size distribution
LATR_recruit_size <- readRDS("LATR_recruit_size.rds")
recruit_size<-function(y){
  dnorm(x=y,mean=LATR_recruit_size$recruit_mean,sd=LATR_recruit_size$recruit_sd)
}

## combined flowering, fertility, and recruitment
fxy <- function(x,y,d){
  flower_fn(x,d) * seeds_fn(x,d) * recruitment_fn(d) * recruit_size(y)
}

#PUT IT ALL TOGETHER
## projection matrix is a function of weighted density (dens)
bigmatrix<-function(lower.extension = -8, ## needs a large lower extension because growth variance (gaussian) is greater for smaller plants 
                    upper.extension = 2,
                    min.size = LATR_size_bounds$min_size,
                    max.size = LATR_size_bounds$max_size,
                    mat.size = 200,
                    dens,
                    dist="SGT"){
  
  n<-mat.size
  L<-min.size + lower.extension
  U<-max.size + upper.extension
  #these are the upper and lower integration limits
  h<-(U-L)/n                   #Bin size
  b<-L+c(0:n)*h;               #Lower boundaries of bins 
  y<-0.5*(b[1:n]+b[2:(n+1)]);  #Bin midpoints
  
  # Growth/Survival matrix
  Pmat<-t(outer(y,y,pxy,d=dens,dist=dist)) * h 
  
  # Fertility/Recruiment matrix
  Fmat<-t(outer(y,y,fxy,d=dens)) * h 

  # Put it all together
  IPMmat<-Pmat+Fmat

    #and transition matrix
  return(list(IPMmat=IPMmat,Fmat=Fmat,Pmat=Pmat,meshpts=y))
}
