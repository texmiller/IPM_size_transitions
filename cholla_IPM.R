# SETUP -------------------------------------------------------------------
library(tidyverse)
library(popbio)
volume <- function(h, w, p){
  (1/3)*pi*h*(((w + p)/2)/2)^2
}
invlogit<-function(x){exp(x)/(1+exp(x))}

dir <- "C:/Users/tm634/"
dir <- "C:/Users/tm9/"

cholla <- read.csv(paste0(dir,"Dropbox/IPM size transitions/cholla_demography_20042018.csv")) %>% 
  ## drop seed addition plots (don't want them in plot RFX)
  ## change plots and years to 1,...,N integers
  #select(TagID,Plot,Year_t,Height_t,Width_t,Perp_t,Height_t1,Width_t1,Perp_t1) %>% 
  filter(str_sub(Plot,1,1)!="H") %>% 
  mutate(year_int = Year_t - (min(Year_t,na.rm = T)-1),
         plot_int = ifelse(Plot=="T1",9,ifelse(Plot=="T2",10,ifelse(Plot=="T3",11,as.integer(Plot)))),
         ind_int = as.integer(as.numeric(interaction(plot_int,TagID))),
         vol_t = volume(h = Height_t, w = Width_t, p = Perp_t),
         vol_t1 = volume(h = Height_t1, w = Width_t1, p = Perp_t1)) %>% 
  #filter(!is.na(vol_t),
  #       !is.na(vol_t1)) %>% 
  mutate(size_change = vol_t1 - vol_t)
## find outliers
(high_grow <- cholla %>% filter(size_change > quantile(size_change,probs=c(0.99))))
drop_high <- high_grow[c(22,25),]## these two stand out by an order of magnitude and are clearly wrong
(low_grow <- cholla %>% filter(size_change < quantile(size_change,probs=c(0.01))))
drop_low <- low_grow[c(8,9,22,27,31,36),]## these three are errors
drop <- bind_rows(drop_high,drop_low)
cholla <- anti_join(cholla, drop)

     
## the seedling data set minimum is a little smaller than the rest of the data set, so
## I will use seedling min size as the lower bound
min.size <- min(log(cholla$vol_t),na.rm=T) 
max.size <- cholla %>% mutate(vol_t = log(volume(h = Height_t, w = Width_t, p = Perp_t))) %>% 
  summarise(max(vol_t,na.rm=T)) 

## fit simple models
surv_model <- glm(Survival_t1 ~ log(vol_t),family="binomial",data=cholla)
flow_model <- glm(TotFlowerbuds_t1>0 ~ log(vol_t1),family="binomial",data=cholla)
fert_model <- glm(TotFlowerbuds_t1 ~ log(vol_t1),family="poisson",data=subset(cholla,TotFlowerbuds_t1>0))

plot(log(cholla$vol_t),cholla$Survival_t1)

## other misc params
params_post <- read.csv("C:/Users/tm9/Desktop/git local/cholla_climate_IPM/allrates.selected.posterior.csv")


## parameter vector
cholla_params<- c()
cholla_params$surv.mu <- coef(surv_model)[1]
cholla_params$surv.bsize <- coef(surv_model)[2]
cholla_params$flow.mu <- coef(flow_model)[1]
cholla_params$flow.bsize <- coef(flow_model)[2]
cholla_params$fert.mu <- coef(fert_model)[1]
cholla_params$fert.bsize <- coef(fert_model)[2]
cholla_params$mu_spf <- mean(params_post$mu_spf)
cholla_params$seedsurv <- mean(params_post$seedsurv)
cholla_params$seedsurv <- mean(params_post$seedsurv)
cholla_params$germ1 <- mean(params_post$germ1)
cholla_params$germ2 <- mean(params_post$germ2)
cholla_params$precenus_surv <- mean(params_post$precenus_surv)
cholla_params$mu_sdlgsize <- mean(params_post$mu_sdlgsize)
cholla_params$sigma_sdlgsize <- mean(params_post$sigma_sdlgsize)
cholla_params$min.size <- min.size
cholla_params$max.size <- max.size$`max(vol_t, na.rm = T)`
## gaussian
cholla_params$grow.mean.gaussian <- mean(cholla_pred_MoM_gaussian$mu)
cholla_params$grow.sd.gaussian <- mean(cholla_pred_MoM_gaussian$sigma)
## sgt
cholla_params$grow.mean.sgt <- mean(cholla_pred_MoM$mu)
cholla_params$grow.sd.sgt <- mean(cholla_pred_MoM$sigma)
cholla_params$grow.lambda.sgt <- mean(cholla_pred_MoM$l)
cholla_params$grow.p.sgt <- mean(cholla_pred_MoM$p)
cholla_params$grow.q.sgt <- mean(cholla_pred_MoM$q)

# VITAL RATE FUNCTIONS ----------------------------------------------------

pxy_gaussian <- function(x,y,params){
  xb=pmin(pmax(x,params$min.size),params$max.size)
  surv <- invlogit(params$surv.mu + params$surv.bsize*xb)
  return(surv * dnorm(y,mean=xb+params$grow.mean.gaussian,
                      sd=params$grow.sd.gaussian))
}

pxy_sgt <- function(x,y,params){
  xb=pmin(pmax(x,params$min.size),params$max.size)
  surv <- invlogit(params$surv.mu + params$surv.bsize*xb)
  return(surv * dsgt(y,mu=xb+params$grow.mean.sgt,
                      sigma=params$grow.sd.sgt,
                      lambda=params$grow.lambda.sgt,
                      p=params$grow.p.sgt,
                      q=params$grow.q.sgt))
}

fx<-function(x,params){
  xb=pmin(pmax(x,params$min.size),params$max.size)
  flow <- invlogit(params$flow.mu + params$flow.bsize*xb)
  nfruits <- exp(params$fert.mu + params$fert.bsize*xb)
  return(flow*nfruits*params$mu_spf*params$seedsurv)  
}

recruit.size<-function(y,params){
  dnorm(x=y,mean=params$mu_sdlgsize,sd=params$sigma_sdlgsize)
}

# BIGMATRIX ---------------------------------------------------------------
bigmatrix<-function(params,
                    lower.extension = 0, ## I'll need to extend lower and upper beyond true size limits
                    upper.extension = 0,
                    mat.size,
                    sgt=F){
  
  n<-mat.size
  L<-params$min.size + lower.extension
  U<-params$max.size + upper.extension
  #these are the upper and lower integration limits
  h<-(U-L)/n                   #Bin size
  b<-L+c(0:n)*h;               #Lower boundaries of bins 
  y<-0.5*(b[1:n]+b[2:(n+1)]);  #Bins' midpoints
  #these are the boundary points (b) and mesh points (y)
  
  # Fertility matrix
  Fmat<-matrix(0,(n+2),(n+2))
  
  # 1-yo banked seeds go in top row
  Fmat[1,3:(n+2)]<-fx(y,params)
  
  # Growth/survival transition matrix
  Tmat<-matrix(0,(n+2),(n+2))
  
  # Graduation to 2-yo seed bank = pr(not germinating as 1-yo)
  Tmat[2,1]<-(1-params$germ1)
  
  # Graduation from 1-yo bank to cts size = germination * size distn * pre-census survival
  Tmat[3:(n+2),1]<- params$germ1 * params$precenus_surv * recruit.size(y,params) * h   
  
  # Graduation from 2-yo bank to cts size = germination * size distn * pre-census survival
  Tmat[3:(n+2),2]<- params$germ2 * params$precenus_surv * recruit.size(y,params) * h  
  
  # Growth/survival transitions among cts sizes
  if(sgt==F){Tmat[3:(n+2),3:(n+2)]<-t(outer(y,y,pxy_gaussian,params=params)) * h} 
  if(sgt==T){Tmat[3:(n+2),3:(n+2)]<-t(outer(y,y,pxy_sgt,params=params)) * h} 
  
  # Put it all together
  IPMmat<-Fmat+Tmat     #Full Kernel is simply a summation ot fertility
  #and transition matrix
  return(list(IPMmat=IPMmat,Fmat=Fmat,Tmat=Tmat,meshpts=y))
}

gaussian_mat <- bigmatrix(params=cholla_params,
          lower.extension = -.2, ## I'll need to extend lower and upper beyond true size limits
          upper.extension = 1,
          mat.size=200,
          sgt=F)$IPMmat
Re(eigen(gaussian_mat)$values[1])

sgt_mat <- bigmatrix(params=cholla_params,
                          lower.extension = -.2, ## I'll need to extend lower and upper beyond true size limits
                          upper.extension = 1,
                          mat.size=200,
                          sgt=T)$IPMmat
Re(eigen(sgt_mat)$values[1])

# posterior samples of lambda ----------------------------------------
n_post <- pmin(100,nrow(params_post)) ## number of posterior draws
rand.indices <- sample.int(nrow(params_post), size=n_post)
mat.size = 200
lower.extension = -0.2
upper.extension = 1

## Gaussian growth model
cholla_fit_norm <- read_rds(paste0(dir,"Dropbox/IPM size transitions/cholla_fit.rds"))
params_norm <- rstan::extract(cholla_fit_norm, pars = c("b_0","b_size","d_0","d_size"))

lambda_post <-c()

for(i in 1:n_post){
  print(i)
  ## now convert params to list for the rest of it
  sample.params <- as.list(params_post[i,])
  sample.params$flow.bclim <- params_post[rand.indices[i],] %>% 
    select(flow.bclim.1.1.:flow.bclim.3.3.) %>% 
    matrix(nrow=3)
  sample.params$fert.bclim <- params_post[rand.indices[i],] %>% 
    select(fert.bclim.1.1.:fert.bclim.3.3.) %>% 
    matrix(nrow=3)  
  sample.params$surv.bclim <- params_post[rand.indices[i],] %>% 
    select(surv.bclim.1.1.:surv.bclim.3.3.) %>% 
    matrix(nrow=3) 
  sample.params$grow.bclim <- params_post[rand.indices[i],] %>% 
    select(grow.bclim.1.1.:grow.bclim.3.3.) %>% 
    matrix(nrow=3) 
  sample.params$min.size <- min.size
  sample.params$max.size <- max.size$`max(vol_t)`
  
  ## growth
  #sample.params$grow.mu <- params_norm$b_0[i]
  #sample.params$grow.bsize <- params_norm$b_size[i]
  #sample.params$d0 <- params_norm$d_0[i]
  #sample.params$d_size <- params_norm$d_size[i]  
  
  lambda_post[i]<-lambda(bigmatrix(params = sample.params,
                                           PC1 = rep(0,2), PC2 = rep(0,2), PC3 = rep(0,2),
                                           random = F, 
                                           lower.extension = lower.extension, 
                                           upper.extension = upper.extension,
                                           mat.size = mat.size)$IPMmat)


}

hist(lambda_post)

params_post[rand.indices[i],] %>% 
  select(grow.mu) %>% summarise(mean(grow.mu))
