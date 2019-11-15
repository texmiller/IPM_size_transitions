
data {
  int<lower=1> N;
  real y[N];
  real sizet[N];
  int<lower=0, upper=1> density[N];
  int<lower=1> Ntransects;
  int<lower=1, upper=Ntransects> transect[N];
  int<lower=1> Nyears;
  int<lower=1, upper=Nyears> year[N];  
}

parameters {
  real b_0[Nyears];
  real b_size[Nyears];
  real b_density[Nyears];
  real b_size_density[Nyears];
  real<lower=0> sigma_transect;
  real transect_rfx[Ntransects];
  real d_0;
  real d_size;
}

transformed parameters{
  real pred[N];
  real<lower=0> std[N];
  for(i in 1:N){
    pred[i] = b_0[year[i]] + b_size[year[i]]*sizet[i] + b_density[year[i]]*density[i] + 
    b_size_density[year[i]]*sizet[i]*density[i] + transect_rfx[transect[i]];
    
    std[i] = exp(d_0 + d_size * sizet[i]);
  }
}

model {

  for(i in 1:Nyears){
    b_0[i]~normal(0,100);
    b_size[i]~normal(0,100);
    b_density[i]~normal(0,100);
    b_size_density[i]~normal(0,100);
  }
  sigma_transect ~ inv_gamma(0.001, 0.001);
  for(i in 1:Ntransects){
    transect_rfx[i]~normal(0,sigma_transect);
  }
  
  y ~ normal(pred,std);
}

