
data {
  int<lower=1> N;
  real y[N];
  real sizet[N];
  int<lower=0, upper=1> light[N];
  int<lower=0, upper=1> flower[N];
  int<lower=1> Nyears;
  int<lower=1, upper=Nyears> year[N];
}

parameters {
  real b_0;
  real b_size;
  real b_light;
  real b_flower;
  real b_size_light;
  real b_size_flower;
  real b_light_flower;
  real b_size_light_flower;
  
  //real<lower=0> sigma_res;
  real<lower=0> sigma_year;
  real year_rfx[Nyears];
  real d_0;
  real d_size;
}

transformed parameters{
  real pred[N];
  real<lower=0> std[N];
  for(i in 1:N){
    pred[i] = b_0 + b_size*sizet[i] + b_light*light[i] + b_flower*flower[i] +
    b_size_light*sizet[i]*light[i] + b_size_flower*sizet[i]*flower[i] + 
    b_light_flower*light[i]*flower[i] + b_size_light_flower*sizet[i]*light[i]*flower[i] +
    year_rfx[year[i]];
    
    std[i] = exp(d_0 + d_size * sizet[i]);
  }
}

model {
  b_0 ~ normal(0, 100);    
  b_size ~ normal(0, 100);    
  b_light ~ normal(0, 100);
  b_flower ~ normal(0, 100);
  b_size_light ~ normal(0, 100);
  b_size_flower ~ normal(0, 100);
  b_light_flower ~ normal(0, 100);
  b_size_light_flower ~ normal(0, 100);

  sigma_year ~ inv_gamma(0.001, 0.001);
  for(i in 1:Nyears){
    year_rfx[i]~normal(0,sigma_year);
  }
  
  y ~ normal(pred,std);
}

