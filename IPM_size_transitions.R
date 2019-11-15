library(tidyverse)
library(rstan)
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
  select(Plot,Year_t,Height_t,Width_t,Perp_t,Height_t1,Width_t1,Perp_t1) %>% 
  filter(str_sub(Plot,1,1)!="H") %>% 
  mutate(year_int = Year_t - (min(Year_t,na.rm = T)-1),
         plot_int = ifelse(Plot=="T1",9,ifelse(Plot=="T2",10,ifelse(Plot=="T3",11,as.integer(Plot)))),
         vol_t = log(volume(h = Height_t, w = Width_t, p = Perp_t)),
         vol_t1 = log(volume(h = Height_t1, w = Width_t1, p = Perp_t1))) %>% 
  filter(!is.na(vol_t),
         !is.na(vol_t1))

## prep model for Stan
cholla_dat <- list(cholla_N = nrow(cholla),
                   cholla_sizet = cholla$vol_t,
                   cholla_delta_size = (cholla$vol_t1 - cholla$vol_t),
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

# Posterior predictive checks ---------------------------------------------
## need to generate simulated data, doing this in Stan gave me errors (problems with log_neg_binom_2_rng)
cholla_pred <- rstan::extract(cholla_fit, pars = c("cholla_pred","cholla_sd","b_0","b_size","d_0","d_size"))

mcmc_dens_overlay(cholla_fit,par=c("b_0","b_size","d_0","d_size"))

n_post_draws <- 500
post_draws <- sample.int(dim(cholla_pred$cholla_pred)[1], n_post_draws)
y_cholla_sim <- matrix(NA,n_post_draws,cholla_dat$cholla_N)
for(i in 1:n_post_draws){
  y_cholla_sim[i,] <- rnorm(n=cholla_dat$cholla_N, mean = cholla_pred$cholla_pred[i,],sd = cholla_pred$cholla_sd[i,])
}
ppc_dens_overlay(cholla_dat$cholla_delta_size, y_cholla_sim)
ppc_stat(cholla_dat$cholla_delta_size, y_cholla_sim,stat="mean")
ppc_stat(cholla_dat$cholla_delta_size, y_cholla_sim,stat="sd")
ppc_stat(cholla_dat$cholla_delta_size, y_cholla_sim,stat="skewness")
ppc_stat(cholla_dat$cholla_delta_size, y_cholla_sim,stat="kurtosis")


# Orchis ------------------------------------------------------------------
## read in cholla data
orchid <- read.csv(paste0(dir,"Dropbox/IPM size transitions/Orchis_IPM_data.csv")) %>% 
  select(light,begin.year,total.leaf.area,flowering,end.total.leaf.area) %>% 
  mutate(size_t = log(total.leaf.area),
         size_t1 = log(end.total.leaf.area),
         year_int = begin.year - (min(begin.year,na.rm = T)-1),
         light = as.integer(ifelse(light=="L"|light=="#L",1,0))) %>% 
  filter(!is.na(size_t),
         !is.na(size_t1))  
  
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

orchid_pred <- rstan::extract(orchid_fit, pars = c("pred","std"))

y_orchid_sim <- matrix(NA,n_post_draws,orchid_dat$N)
for(i in 1:n_post_draws){
  y_orchid_sim[i,] <- rnorm(n=orchid_dat$N, mean = orchid_pred$pred[i,],sd = orchid_pred$std[i,])
}
ppc_dens_overlay(orchid_dat$y, y_orchid_sim)
ppc_stat(orchid_dat$y, y_orchid_sim,stat="mean")
ppc_stat(orchid_dat$y, y_orchid_sim,stat="sd")
ppc_stat(orchid_dat$y, y_orchid_sim,stat="skewness")
ppc_stat(orchid_dat$y, y_orchid_sim,stat="kurtosis")


# Creosote ----------------------------------------------------------------
creosote <- read.csv(paste0(dir,"Dropbox/IPM size transitions/creosote_dat.csv"))

