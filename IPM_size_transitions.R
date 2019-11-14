library(tidyverse)
library(rstan)
library(bayesplot)
library(moments)

volume <- function(h, w, p){
  (1/3)*pi*h*(((w + p)/2)/2)^2
}

## read in cholla data
cholla <- read.csv("C:/Users/tm634/Dropbox/IPM size transitions/cholla_demography_20042018.csv") %>% 
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

plot(cholla$vol_t,cholla$vol_t1)

## prep model for Stan
cholla_dat <- list(cholla_N = nrow(cholla),
                   cholla_sizet = cholla$vol_t,
                   cholla_sizet1 = cholla$vol_t1,
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
  file = 'size_transition.stan',
  data = cholla_dat,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains )

# Posterior predictive checks ---------------------------------------------
## need to generate simulated data, doing this in Stan gave me errors (problems with log_neg_binom_2_rng)
cholla_pred <- rstan::extract(cholla_fit, pars = c("cholla_pred","cholla_sd","b_0","b_size","d_0","d_size"))

mcmc_dens_overlay(cholla_fit,par=c("b_0","b_size","d_0","d_size"))

n_post_draws <- 1000
post_draws <- sample.int(dim(cholla_pred$cholla_pred)[1], n_post_draws)
y_cholla_sim <- matrix(NA,n_post_draws,cholla_dat$cholla_N)
for(i in 1:n_post_draws){
  y_cholla_sim[i,] <- rnorm(n=cholla_dat$cholla_N, mean = cholla_pred$cholla_pred[i,],sd = cholla_pred$cholla_sd[i])
}
ppc_dens_overlay(cholla_dat$cholla_sizet1, y_cholla_sim)
ppc_stat(cholla_dat$cholla_sizet1, y_cholla_sim,stat="mean")
ppc_stat(cholla_dat$cholla_sizet1, y_cholla_sim,stat="sd")
ppc_stat(cholla_dat$cholla_sizet1, y_cholla_sim,stat="skewness")
ppc_stat(cholla_dat$cholla_sizet1, y_cholla_sim,stat="kurtosis")
