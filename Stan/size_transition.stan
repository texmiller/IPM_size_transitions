
data {
  int<lower=1> cholla_N;
  real cholla_sizet[cholla_N];
  real cholla_sizet1[cholla_N];
  int<lower=1> cholla_Nplots;
  int<lower=1> cholla_Nyears;
  int<lower=1, upper=cholla_Nplots> cholla_plot[cholla_N];
  int<lower=1, upper=cholla_Nyears> cholla_year[cholla_N];
}

parameters {
  real b_0;
  real b_size;
  //real<lower=0> sigma_res;
  real<lower=0> sigma_plot;
  real<lower=0> sigma_year;
  real plot_rfx[cholla_Nplots];
  real year_rfx[cholla_Nyears];
  real d_0;
  real d_size;
}

transformed parameters{
  real cholla_pred[cholla_N];
  real<lower=0> cholla_sd[cholla_N];
  for(i in 1:cholla_N){
    cholla_pred[i] = b_0 + b_size * cholla_sizet[i] + plot_rfx[cholla_plot[i]] + year_rfx[cholla_year[i]];
    cholla_sd[i] = exp(d_0 + d_size * cholla_sizet[i]);
  }
}

model {
  b_0 ~ normal(0, 100);    
  b_size ~ normal(0, 100);    
  //sigma_res ~ inv_gamma(0.001, 0.001);
  sigma_plot ~ inv_gamma(0.001, 0.001);
  for(i in 1:cholla_Nplots){
    plot_rfx[i]~normal(0,sigma_plot);
  }
  sigma_year ~ inv_gamma(0.001, 0.001);
  for(i in 1:cholla_Nyears){
    year_rfx[i]~normal(0,sigma_year);
  }
  
  //cholla_sizet1 ~ normal(cholla_pred,sigma_res);
  cholla_sizet1 ~ normal(cholla_pred,cholla_sd);
}

