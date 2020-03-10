
data {
  int<lower=1> N;
  real y[N];
  real sizet[N];
  real density[N];
  int<lower=1> Ntransects;
  int<lower=1, upper=Ntransects> transect[N];
  int<lower=1> Nsites;
  int<lower=1, upper=Nsites> site[N];  
  int<lower=1> Nyears;
  int<lower=1, upper=Nyears> year[N];  
}

parameters {
  real b_0;
  real b_size;
  real b_density;
  real b_size_density;
  real<lower=0> sigma_site;
  real site_rfx[Nsites];
  real<lower=0> sigma_transect;
  real transect_rfx[Ntransects];
  real<lower=0> sigma_year;
  real year_rfx[Nyears];
  real d_0;
  real d_size;
}

transformed parameters{
  real pred[N];
  real<lower=0> std[N];
  for(i in 1:N){
    pred[i] = b_0 + b_size*sizet[i] + b_density*density[i] + b_size_density*sizet[i]*density[i] + 
    site_rfx[site[i]] + transect_rfx[transect[i]] + year_rfx[year[i]];
    
    std[i] = exp(d_0 + d_size * sizet[i]);
  }
}

model {

  b_0~normal(0,100);
  b_size~normal(0,100);
  b_density~normal(0,100);
  b_size_density~normal(0,100);
  sigma_site ~ inv_gamma(0.001, 0.001);
  for(i in 1:Nsites){
    site_rfx[i]~normal(0,sigma_site);
  }  
  sigma_transect ~ inv_gamma(0.001, 0.001);
  for(i in 1:Ntransects){
    transect_rfx[i]~normal(0,sigma_transect);
  }
  sigma_year ~ inv_gamma(0.001, 0.001);
  for(i in 1:Nyears){
    year_rfx[i]~normal(0,sigma_year);
  }
  
  y ~ normal(pred,std);
}

