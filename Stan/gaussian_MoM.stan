data {
  int<lower=0> N;
  vector[N] delta_size;
}
parameters {
  real mu; 
  real<lower=0> sigma;
}
model {
  mu ~ normal(0,100);
  sigma ~ inv_gamma(0.001, 0.001);
  delta_size ~ normal(mu, sigma);
}