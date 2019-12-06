library(tidyverse)
library(rstan)
rstan_options( auto_write = TRUE )
options( mc.cores = parallel::detectCores() )
library(bayesplot)
library(moments)

volume <- function(h, w, p){
  (1/3)*pi*h*(((w + p)/2)/2)^2
}

dir <- "C:/Users/tm634/"
dir <- "C:/Users/tm9/"


# Cholla ------------------------------------------------------------------

## read in cholla data
cholla <- read.csv(paste0(dir,"Dropbox/IPM size transitions/cholla_demography_20042018.csv")) %>% 
## drop seed addition plots (don't want them in plot RFX)
## change plots and years to 1,...,N integers
  select(TagID,Plot,Year_t,Height_t,Width_t,Perp_t,Height_t1,Width_t1,Perp_t1) %>% 
  filter(str_sub(Plot,1,1)!="H") %>% 
  mutate(year_int = Year_t - (min(Year_t,na.rm = T)-1),
         plot_int = ifelse(Plot=="T1",9,ifelse(Plot=="T2",10,ifelse(Plot=="T3",11,as.integer(Plot)))),
         ind_int = as.integer(as.numeric(interaction(plot_int,TagID))),
         vol_t = volume(h = Height_t, w = Width_t, p = Perp_t),
         vol_t1 = volume(h = Height_t1, w = Width_t1, p = Perp_t1)) %>% 
  filter(!is.na(vol_t),
         !is.na(vol_t1)) %>% 
  mutate(size_change = vol_t1 - vol_t)

plot((cholla$vol_t),(cholla$vol_t1))

## find outliers
(high_grow <- cholla %>% filter(size_change > quantile(size_change,probs=c(0.99))))
drop_high <- high_grow[c(22,25),]## these two stand out by an order of magnitude and are clearly wrong
(low_grow <- cholla %>% filter(size_change < quantile(size_change,probs=c(0.01))))
drop_low <- low_grow[c(8,9,22,27,31,36),]## these three are errors
drop <- bind_rows(drop_high,drop_low)

cholla <- anti_join(cholla, drop)
plot((cholla$vol_t),(cholla$vol_t1))
plot(log(cholla$vol_t),log(cholla$vol_t1))

## prep model for Stan
cholla_dat <- list(cholla_N = nrow(cholla),
                   cholla_sizet = log(cholla$vol_t),
                   cholla_delta_size = log(cholla$vol_t1) - log(cholla$vol_t),
                   cholla_Nplots = max(cholla$plot_int),
                   cholla_Nyears = max(cholla$year_int),
                   cholla_plot = cholla$plot_int,
                   cholla_year = cholla$year_int)

sim_pars <- list(
  warmup = 1000, 
  iter = 10000, 
  thin = 3, 
  chains = 3
)

cholla_fit <- stan(
  file = 'cholla_growth.stan',
  data = cholla_dat,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains )
#write_rds(cholla_fit,paste0(dir,"Dropbox/IPM size transitions/cholla_fit.rds"))
cholla_fit <- read_rds(paste0(dir,"Dropbox/IPM size transitions/cholla_fit.rds"))

# Posterior predictive checks ---------------------------------------------
## need to generate simulated data, doing this in Stan gave me errors (problems with log_neg_binom_2_rng)
cholla_pred <- rstan::extract(cholla_fit, pars = c("cholla_pred","cholla_sd","b_0","b_size","d_0","d_size"))

plot(cholla_dat$cholla_sizet,1 + cholla_dat$cholla_sizet*0,ylim=c(0,5),type="n")
for(i in 1:n_post_draws){
lines(cholla_dat$cholla_sizet , cholla_pred$d_0[i] + cholla_pred$d_size[i] * cholla_dat$cholla_sizet)
  }

plot(cholla_dat$cholla_sizet, cholla_pred$cholla_sd[100,])

plot(exp(cholla_dat$cholla_sizet),exp(cholla_dat$cholla_sizet1))
mcmc_dens_overlay(cholla_fit,par=c("b_0","b_size","d_0","d_size"))

n_post_draws <- 500
post_draws <- sample.int(dim(cholla_pred$cholla_pred)[1], n_post_draws)
y_cholla_sim <- matrix(NA,n_post_draws,cholla_dat$cholla_N)
for(i in 1:n_post_draws){
  y_cholla_sim[i,] <- rnorm(n=cholla_dat$cholla_N, mean = cholla_pred$cholla_pred[i,],sd = cholla_pred$cholla_sd[i,])
}
ppc_dens_overlay(cholla_dat$cholla_delta_size, y_cholla_sim)
ppc_stat(cholla_dat$cholla_delta_size, y_cholla_sim,stat="mean")+theme(legend.position = "none")
ppc_stat(cholla_dat$cholla_delta_size, y_cholla_sim,stat="sd")+theme(legend.position = "none")
ppc_stat(cholla_dat$cholla_delta_size, y_cholla_sim,stat="skewness")+theme(legend.position = "none")
ppc_stat(cholla_dat$cholla_delta_size, y_cholla_sim,stat="kurtosis")+theme(legend.position = "none")



for(i in 1:n_post_draws){
  steve <- ((cholla_dat$cholla_delta_size - cholla_pred$cholla_pred[i,]) / cholla_pred$cholla_sd[i,])
  #steve <- steve / sd(steve)
  y_cholla_sim[i,] <- cholla_pred$cholla_pred[i,] + cholla_pred$cholla_sd[i,] * sample(steve,replace = T)#rnorm(n=cholla_dat$cholla_N, mean = cholla_pred$cholla_pred[i,],sd = cholla_pred$cholla_sd[i,])
}
ppc_dens_overlay(cholla_dat$cholla_delta_size, y_cholla_sim)
ppc_stat(cholla_dat$cholla_delta_size, y_cholla_sim,stat="skewness")+theme(legend.position = "none")


## what if I added an individual random effect?
## prep model for Stan
cholla_ind_dat <- list(cholla_N = nrow(cholla),
                   cholla_sizet = cholla$vol_t,
                   cholla_delta_size = (cholla$vol_t1 - cholla$vol_t),
                   cholla_Nplots = max(cholla$plot_int),
                   cholla_Nyears = max(cholla$year_int),
                   cholla_Nind = max(cholla$ind_int),
                   cholla_plot = cholla$plot_int,
                   cholla_year = cholla$year_int,
                   cholla_ind = cholla$ind_int)
cholla_ind_fit <- stan(
  file = 'cholla_growth_individual.stan',
  data = cholla_ind_dat,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains )

cholla_ind_pred <- rstan::extract(cholla_ind_fit, pars = c("cholla_pred","cholla_sd"))
y_cholla_ind_sim <- matrix(NA,n_post_draws,cholla_ind_dat$cholla_N)
for(i in 1:n_post_draws){
  y_cholla_ind_sim[i,] <- rnorm(n=cholla_ind_dat$cholla_N, mean = cholla_ind_pred$cholla_pred[i,],sd = cholla_ind_pred$cholla_sd[i,])
}
ppc_dens_overlay(cholla_ind_dat$cholla_delta_size, y_cholla_ind_sim)
ppc_stat(cholla_ind_dat$cholla_delta_size, y_cholla_ind_sim,stat="mean")+theme(legend.position = "none")
ppc_stat(cholla_ind_dat$cholla_delta_size, y_cholla_ind_sim,stat="sd")+theme(legend.position = "none")
ppc_stat(cholla_ind_dat$cholla_delta_size, y_cholla_ind_sim,stat="skewness")+theme(legend.position = "none")
ppc_stat(cholla_ind_dat$cholla_delta_size, y_cholla_ind_sim,stat="kurtosis")+theme(legend.position = "none")


# Orchis ------------------------------------------------------------------
## read in cholla data
orchid <- read.csv(paste0(dir,"Dropbox/IPM size transitions/Orchis_IPM_data.csv")) %>% 
  select(light,begin.year,total.leaf.area,flowering,end.total.leaf.area) %>% 
  mutate(size_t = (total.leaf.area),
         size_t1 = (end.total.leaf.area),
         year_int = begin.year - (min(begin.year,na.rm = T)-1),
         light = as.integer(ifelse(light=="L"|light=="#L",1,0))) %>% 
  filter(!is.na(size_t),
         !is.na(size_t1))  
  

plot(log(orchid$total.leaf.area),log(orchid$end.total.leaf.area))

orchid_dat <- list(N = nrow(orchid),
                   y = (orchid$size_t1 - orchid$size_t),
                   sizet = orchid$size_t,
                   flower = orchid$flowering,
                   light = orchid$light,
                   Nyears = max(orchid$year_int),
                   year = orchid$year_int)

orchid_fit <- stan(
  file = 'orchid_growth.stan',
  data = orchid_dat,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains )
#write_rds(orchid_fit,paste0(dir,"Dropbox/IPM size transitions/orchid_fit.rds"))
orchid_fit <- read_rds(paste0(dir,"Dropbox/IPM size transitions/orchid_fit.rds"))

orchid_pred <- rstan::extract(orchid_fit, pars = c("pred","std"))

y_orchid_sim <- matrix(NA,n_post_draws,orchid_dat$N)
for(i in 1:n_post_draws){
  y_orchid_sim[i,] <- rnorm(n=orchid_dat$N, mean = orchid_pred$pred[i,],sd = orchid_pred$std[i,])
}
ppc_dens_overlay(orchid_dat$y, y_orchid_sim)
ppc_stat(orchid_dat$y, y_orchid_sim,stat="mean")+theme(legend.position = "none")
ppc_stat(orchid_dat$y, y_orchid_sim,stat="sd")+theme(legend.position = "none")
ppc_stat(orchid_dat$y, y_orchid_sim,stat="skewness")+theme(legend.position = "none")
ppc_stat(orchid_dat$y, y_orchid_sim,stat="kurtosis")+theme(legend.position = "none")


for(i in 1:n_post_draws){
  steve <- ((orchid_dat$delta_size - orchid_pred$pred[i,]) / orchid_pred$std[i,])
  #steve <- steve / sd(steve)
  y_orchid_sim[i,] <- orchid_pred$pred[i,] + orchid_pred$std[i,] * sample(steve,replace = T)#rnorm(n=cholla_dat$cholla_N, mean = cholla_pred$cholla_pred[i,],sd = cholla_pred$cholla_sd[i,])
}
ppc_dens_overlay(cholla_dat$cholla_delta_size, y_cholla_sim)
ppc_stat(cholla_dat$cholla_delta_size, y_cholla_sim,stat="skewness")+theme(legend.position = "none")


# Creosote ----------------------------------------------------------------
creosote <- read.csv(paste0(dir,"Dropbox/IPM size transitions/creosote_dat.csv")) %>% 
  mutate(year_int = year_t - (min(year_t,na.rm = T)-1),
         site_int = as.integer(as.numeric(site)),
         unique_transect = as.integer(as.numeric(interaction(site,transect)))) %>% 
  filter(!is.na(volume_t),
         !is.na(volume_t1),
         !is.na(weighted.dens))   

creosote_dat <- list(N = nrow(creosote),
                   y = (creosote$volume_t1 - creosote$volume_t),
                   sizet = creosote$volume_t,
                   density = (creosote$weighted.dens - mean(creosote$weighted.dens)),
                   Nsites = max(creosote$site_int),
                   site = creosote$site_int,
                   Ntransects = max(creosote$unique_transect),
                   transect = creosote$unique_transect,
                   Nyears = max(creosote$year_int),
                   year = creosote$year_int)

creosote_fit <- stan(
  file = 'creosote_growth_rfx.stan',
  data = creosote_dat,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains )

#write_rds(creosote_fit,paste0(dir,"Dropbox/IPM size transitions/creosote_fit.rds"))
creosote_fit <- read_rds(paste0(dir,"Dropbox/IPM size transitions/creosote_fit.rds"))

creosote_pred <- rstan::extract(creosote_fit, pars = c("pred","std"))

y_creosote_sim <- matrix(NA,n_post_draws,creosote_dat$N)
for(i in 1:n_post_draws){
  y_creosote_sim[i,] <- rnorm(n=creosote_dat$N, mean = creosote_pred$pred[i,],sd = creosote_pred$std[i,])
}
ppc_dens_overlay(creosote_dat$y, y_creosote_sim)
ppc_stat(creosote_dat$y, y_creosote_sim,stat="mean")+theme(legend.position = "none")
ppc_stat(creosote_dat$y, y_creosote_sim,stat="sd")+theme(legend.position = "none")
ppc_stat(creosote_dat$y, y_creosote_sim,stat="skewness")+theme(legend.position = "none")
ppc_stat(creosote_dat$y, y_creosote_sim,stat="kurtosis")+theme(legend.position = "none")


# ML model selection ------------------------------------------------------
library(sn)
library(sgt)

## normal
cholla_normal <- function(params,log_ratio){
 -sum(dnorm(x=log_ratio,mean=params[1],sd=params[2],log=T))
}
MLE_cholla_normal <- optim(par=c(0,1),fn=cholla_normal,log_ratio=log(cholla$vol_t1/cholla$vol_t))
AIC_cholla_normal <- 2*MLE_cholla_normal $value+2*length(MLE_cholla_normal$par)

## skewed normal
cholla_sn <- function(params,log_ratio){
  -sum(dsn(x=log_ratio,xi=params[1],omega=params[2],alpha=params[3],log=T))
}
MLE_cholla_sn <- optim(par=c(0,1,0),fn=cholla_sn,log_ratio=log(cholla$vol_t1/cholla$vol_t))
AIC_cholla_sn <- 2*MLE_cholla_sn $value+2*length(MLE_cholla_sn$par)

## skewed t
cholla_st <- function(params,log_ratio){
  -sum(dst(x=log_ratio,xi=params[1],omega=params[2],alpha=params[3],nu=params[4],log=T))
}
MLE_cholla_st <- optim(par=c(0,1,0,100),fn=cholla_st,log_ratio=log(cholla$vol_t1/cholla$vol_t))
AIC_cholla_st <- 2*MLE_cholla_st $value+2*length(MLE_cholla_st$par)

## skewed generalized t
cholla_sgt <- function(params,log_ratio){
  -sum(dsgt(x=log_ratio, mu=params[1], sigma=params[2],
            lambda=params[3], p = params[4], q = params[5], 
            mean.cent = TRUE, var.adj = TRUE, log=T))
}
MLE_cholla_sgt <- optim(par=c(0,1,0,2,100),fn=cholla_sgt,log_ratio=log(cholla$vol_t1/cholla$vol_t))
AIC_cholla_sgt <- 2*MLE_cholla_sgt $value+2*length(MLE_cholla_sgt$par)

AIC_cholla_normal
AIC_cholla_sn
AIC_cholla_st
AIC_cholla_sgt


# Stan SGT ----------------------------------------------------------------
## run and write all of the following over a few days
sim_pars <- list(
  warmup = 2000, 
  iter = 15000, 
  thin = 3, 
  chains = 3)

## cholla MoM and linear predictor fits
cholla_dat_MoM <- list(N = nrow(cholla),
                   delta_size = log(cholla$vol_t1) - log(cholla$vol_t))
cholla_sgt_MoM <- stan(
  file = 'skewgent_MoM.stan',
  data = cholla_dat_MoM,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains )
write_rds(cholla_sgt_MoM,paste0(dir,"Dropbox/IPM size transitions/cholla_sgt_MoM.rds"))

## add linear predictor with random effects
cholla_dat <- list(cholla_N = nrow(cholla),
                   cholla_sizet = log(cholla$vol_t),
                   cholla_delta_size = log(cholla$vol_t1) - log(cholla$vol_t),
                   cholla_Nplots = max(cholla$plot_int),
                   cholla_Nyears = max(cholla$year_int),
                   cholla_plot = cholla$plot_int,
                   cholla_year = cholla$year_int)

cholla_sgt_linpred <- stan(
  file = 'skewgent_linpred_cholla.stan',
  data = cholla_dat,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains )
write_rds(cholla_sgt_linpred,paste0(dir,"Dropbox/IPM size transitions/cholla_sgt_linpred.rds"))

## orchids MoM and linear predictor fits
orchid_dat_MoM <- list(N = nrow(orchid),
                       delta_size = log(orchid$size_t1) - log(orchid$size_t))
orchid_sgt_MoM <- stan(
  file = 'skewgent_MoM.stan',
  data = orchid_dat_MoM,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains )
write_rds(orchid_sgt_MoM,paste0(dir,"Dropbox/IPM size transitions/orchid_sgt_MoM.rds"))

orchid_dat <- list(N = nrow(orchid),
                   y = log(orchid$size_t1) - log(orchid$size_t),
                   sizet = log(orchid$size_t),
                   flower = orchid$flowering,
                   light = orchid$light,
                   Nyears = max(orchid$year_int),
                   year = orchid$year_int)
orchid_sgt_linpred <- stan(
  file = 'skewgent_linpred_orchid.stan',
  data = orchid_dat,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains )
write_rds(orchid_sgt_linpred,paste0(dir,"Dropbox/IPM size transitions/orchid_sgt_linpred.rds"))

## creosote
creosote_dat_MoM <- list(N = nrow(creosote),
                         delta_size = creosote$volume_t1 - creosote$volume_t)
creosote_sgt_MoM <- stan(
  file = 'skewgent_MoM.stan',
  data = creosote_dat_MoM,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains )
write_rds(creosote_sgt_MoM,paste0(dir,"Dropbox/IPM size transitions/creosote_sgt_MoM.rds"))

creosote_dat <- list(N = nrow(creosote),
                     y = (creosote$volume_t1 - creosote$volume_t),
                     sizet = creosote$volume_t,
                     density = (creosote$weighted.dens - mean(creosote$weighted.dens)),
                     Nsites = max(creosote$site_int),
                     site = creosote$site_int,
                     Ntransects = max(creosote$unique_transect),
                     transect = creosote$unique_transect,
                     Nyears = max(creosote$year_int),
                     year = creosote$year_int)
creosote_sgt_linpred <- stan(
  file = 'skewgent_linpred_creosote.stan',
  data = creosote_dat,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains )
write_rds(creosote_sgt_linpred,paste0(dir,"Dropbox/IPM size transitions/creosote_sgt_linpred.rds"))



# SGT rstan fits ----------------------------------------------------------
n_post_draws <- 500

## cholla MOM
mcmc_dens_overlay(cholla_sgt_MoM,par=c("mu", "sigma","l","p","q"))
mcmc_trace(cholla_sgt_MoM,par=c("mu", "sigma","l","p","q"))
cholla_pred_MoM <- rstan::extract(cholla_sgt_MoM, pars = c("mu", "sigma", "l", "p", "q"))
post_draws <- sample.int(dim(cholla_pred_MoM$mu)[1], n_post_draws)
y_cholla_MoM <- matrix(NA,n_post_draws,cholla_dat_MoM$N)
for(i in 1:n_post_draws){
  y_cholla_MoM[i,] <- rsgt(n=cholla_dat_MoM$N, 
                           mu = cholla_pred_MoM$mu[post_draws[i]],
                           sigma = cholla_pred_MoM$sigma[post_draws[i]],
                           lambda = cholla_pred_MoM$l[post_draws[i]],
                           p = cholla_pred_MoM$p[post_draws[i]],
                           q = cholla_pred_MoM$q[post_draws[i]])}

ppc_dens_overlay(cholla_dat_MoM$delta_size, y_cholla_MoM) +xlim(-10, 10)
ppc_stat(cholla_dat_MoM$delta_size, y_cholla_MoM,stat="mean")+theme(legend.position = "none")
ppc_stat(cholla_dat_MoM$delta_size, y_cholla_MoM,stat="sd")+theme(legend.position = "none")
ppc_stat(cholla_dat_MoM$delta_size, y_cholla_MoM,stat="skewness")+theme(legend.position = "none")
ppc_stat(cholla_dat_MoM$delta_size, y_cholla_MoM,stat="kurtosis")+theme(legend.position = "none")

## cholla lin pred
mcmc_dens_overlay(cholla_sgt_linpred,par=c("b_0","b_size","d_0","d_size","l","p","q"))
mcmc_trace(cholla_sgt_linpred,par=c("b_0","b_size","d_0","d_size","l","p","q"))
## not sure why but chain 1 is F'd up in this run -- drop it? -- this is annoying but it works
cholla_linpred <- rstan::extract(cholla_sgt_linpred, pars = c("cholla_pred","cholla_sd","l","p","q"),permuted=F)

cholla_linpred_cholla_pred <- rstan::extract(cholla_sgt_linpred, pars = c("cholla_pred"),permuted=F)
cholla_linpred_cholla_sd <- rstan::extract(cholla_sgt_linpred, pars = c("cholla_sd"),permuted=F)
cholla_linpred_l <- rstan::extract(cholla_sgt_linpred, pars = c("l"),permuted=F)
cholla_linpred_p <- rstan::extract(cholla_sgt_linpred, pars = c("p"),permuted=F)
cholla_linpred_q <- rstan::extract(cholla_sgt_linpred, pars = c("q"),permuted=F)

post_draws <- sample.int(dim(cholla_sgt_linpred)[1], n_post_draws)
y_cholla_linpred <- matrix(NA,n_post_draws,cholla_dat$cholla_N)
for(i in 1:n_post_draws){
  y_cholla_linpred[i,] <- rsgt(n=cholla_dat$cholla_N, 
                           mu = cholla_linpred_cholla_pred[post_draws[i],2:3,],
                           sigma = cholla_linpred_cholla_sd[post_draws[i],2:3,],
                           lambda = cholla_linpred_l[post_draws[i],2:3,1],
                           p = cholla_linpred_p[post_draws[i],2:3,1],
                           q = cholla_linpred_q[post_draws[i],2:3,1])
}
ppc_dens_overlay(cholla_dat$cholla_delta_size, y_cholla_linpred) +xlim(-10, 10)
ppc_stat(cholla_dat$cholla_delta_size, y_cholla_linpred,stat="mean")+theme(legend.position = "none")
ppc_stat(cholla_dat$cholla_delta_size, y_cholla_linpred,stat="sd")+theme(legend.position = "none")
ppc_stat(cholla_dat$cholla_delta_size, y_cholla_linpred,stat="skewness")+theme(legend.position = "none")
ppc_stat(cholla_dat$cholla_delta_size, y_cholla_linpred,stat="kurtosis")+theme(legend.position = "none")



## orchid MoM
mcmc_dens_overlay(orchid_sgt_MoM,par=c("mu", "sigma","l","p","q"))
mcmc_trace(orchid_sgt_MoM,par=c("mu", "sigma","l","p","q"))
orchid_pred_MoM <- rstan::extract(orchid_sgt_MoM, pars = c("mu", "sigma", "l", "p", "q"))
post_draws <- sample.int(dim(orchid_pred_MoM$mu)[1], n_post_draws)
y_orchid_MoM <- matrix(NA,n_post_draws,orchid_dat_MoM$N)
for(i in 1:n_post_draws){
  y_orchid_MoM[i,] <- rsgt(n=orchid_dat_MoM$N, 
                           mu = orchid_pred_MoM$mu[post_draws[i]],
                           sigma = orchid_pred_MoM$sigma[post_draws[i]],
                           lambda = orchid_pred_MoM$l[post_draws[i]],
                           p = orchid_pred_MoM$p[post_draws[i]],
                           q = orchid_pred_MoM$q[post_draws[i]])}

ppc_dens_overlay(orchid_dat_MoM$delta_size, y_orchid_MoM) +xlim(-10, 10)
ppc_stat(orchid_dat_MoM$delta_size, y_orchid_MoM,stat="mean")+theme(legend.position = "none")
ppc_stat(orchid_dat_MoM$delta_size, y_orchid_MoM,stat="sd")+theme(legend.position = "none")
ppc_stat(orchid_dat_MoM$delta_size, y_orchid_MoM,stat="skewness")+theme(legend.position = "none")
ppc_stat(orchid_dat_MoM$delta_size, y_orchid_MoM,stat="kurtosis")+theme(legend.position = "none")

## orchid lin pred -- THIS DID NOT RUN (???)
mcmc_dens_overlay(orchid_sgt_linpred,par=c("b_0","b_size","d_0","d_size","l","p","q"))
mcmc_trace(orchid_sgt_linpred,par=c("b_0","b_size","d_0","d_size","l","p","q"))



## creosote MoM
mcmc_dens_overlay(creosote_sgt_MoM,par=c("mu", "sigma","l","p","q"))
mcmc_trace(creosote_sgt_MoM,par=c("mu", "sigma","l","p","q"))
creosote_pred_MoM <- rstan::extract(creosote_sgt_MoM, pars = c("mu", "sigma", "l", "p", "q"))
post_draws <- sample.int(dim(creosote_pred_MoM$mu)[1], n_post_draws)
y_creosote_MoM <- matrix(NA,n_post_draws,creosote_dat_MoM$N)
for(i in 1:n_post_draws){
  y_creosote_MoM[i,] <- rsgt(n=creosote_dat_MoM$N, 
                           mu = creosote_pred_MoM$mu[post_draws[i]],
                           sigma = creosote_pred_MoM$sigma[post_draws[i]],
                           lambda = creosote_pred_MoM$l[post_draws[i]],
                           p = creosote_pred_MoM$p[post_draws[i]],
                           q = creosote_pred_MoM$q[post_draws[i]])}

ppc_dens_overlay(creosote_dat_MoM$delta_size, y_creosote_MoM) +xlim(-10, 10)
ppc_stat(creosote_dat_MoM$delta_size, y_creosote_MoM,stat="mean")+theme(legend.position = "none")
ppc_stat(creosote_dat_MoM$delta_size, y_creosote_MoM,stat="sd")+theme(legend.position = "none")
ppc_stat(creosote_dat_MoM$delta_size, y_creosote_MoM,stat="skewness")+theme(legend.position = "none")
ppc_stat(creosote_dat_MoM$delta_size, y_creosote_MoM,stat="kurtosis")+theme(legend.position = "none")

## creosote lin pred
mcmc_dens_overlay(creosote_sgt_linpred,par=c("b_0","b_size","d_0","d_size","l","p","q"))
mcmc_trace(creosote_sgt_linpred,par=c("b_0","b_size","d_0","d_size","l","p","q"))
## maybe chain 3 was the 'right chain'?? Have a look
creosote_linpred <- rstan::extract(creosote_sgt_linpred, pars = c("creosote_pred","creosote_sd","l","p","q"),permuted=F)

creosote_linpred_creosote_pred <- rstan::extract(creosote_sgt_linpred, pars = c("pred"),permuted=F)
creosote_linpred_creosote_sd <- rstan::extract(creosote_sgt_linpred, pars = c("std"),permuted=F)
creosote_linpred_l <- rstan::extract(creosote_sgt_linpred, pars = c("l"),permuted=F)
creosote_linpred_p <- rstan::extract(creosote_sgt_linpred, pars = c("p"),permuted=F)
creosote_linpred_q <- rstan::extract(creosote_sgt_linpred, pars = c("q"),permuted=F)

post_draws <- sample.int(dim(creosote_sgt_linpred)[1], n_post_draws)
y_creosote_linpred <- matrix(NA,n_post_draws,creosote_dat$N)
for(i in 1:n_post_draws){
  y_creosote_linpred[i,] <- rsgt(n=creosote_dat$N, 
                               mu = creosote_linpred_creosote_pred[post_draws[i],3,],
                               sigma = creosote_linpred_creosote_sd[post_draws[i],3,],
                               lambda = creosote_linpred_l[post_draws[i],3,1],
                               p = creosote_linpred_p[post_draws[i],3,1],
                               q = creosote_linpred_q[post_draws[i],3,1])
}
ppc_dens_overlay(creosote_dat$y, y_creosote_linpred) +xlim(-10, 10)
ppc_stat(creosote_dat$y, y_creosote_linpred,stat="mean")+theme(legend.position = "none")
ppc_stat(creosote_dat$y, y_creosote_linpred,stat="sd")+theme(legend.position = "none")
ppc_stat(creosote_dat$y, y_creosote_linpred,stat="skewness")+theme(legend.position = "none")
ppc_stat(creosote_dat$y, y_creosote_linpred,stat="kurtosis")+theme(legend.position = "none")

