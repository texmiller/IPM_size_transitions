functions {
  real sgt_log(real x, real mu, real s, real l, real p, real q) {
    // Skewed generalised t
    //int N;
    real lz1;
    real lz2;
    real v;
    real m;
    real r;
    real out;
    //N = dims(x)[1];
    lz1=lbeta(1.0/p,q);
    lz2=lbeta(2.0/p,q-1.0/p);
    v=q^(-1.0/p)*((3*l^2+1)*exp(lbeta(3.0/p,q-2.0/p)-lz1)-4*l^2*exp(lz2-lz1)^2)^(-0.5);
    m=2*v*s*l*q^(1.0/p)*exp(lz2-lz1);
    out=0;
    //for (n in 1:N) {
      r=x-mu+m;
      if (r<0)
      	     //out=out+log(p)-log(2*v*s*q^(1.0/p)*exp(lz1)*(fabs(r)^p /(q*(v*s)^p*(l*(-1)+1)^p)+1)^(1.0/p+q));
      	     out=log(p)-log(2*v*s*q^(1.0/p)*exp(lz1)*(fabs(r)^p /(q*(v*s)^p*(l*(-1)+1)^p)+1)^(1.0/p+q));
      else
      	     //out=out+log(p)-log(2*v*s*q^(1.0/p)*exp(lz1)*(fabs(r)^p /(q*(v*s)^p*(l*(1)+1)^p)+1)^(1.0/p+q));
      	     out=log(p)-log(2*v*s*q^(1.0/p)*exp(lz1)*(fabs(r)^p /(q*(v*s)^p*(l*(1)+1)^p)+1)^(1.0/p+q));
    //}
    return out;
  }
}

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
  real<lower=-0.99,upper=0.99> l; 
  real<lower=0.01> p; 
  real<lower=2/p> q; 
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
  l ~ uniform(-0.99,0.99);
  p ~ inv_gamma(0.001, 0.001);
  q ~ inv_gamma(0.001, 0.001);

  for(i in 1:N){  
  y ~ sgt(pred[i], std[i], l, p, q);
  }
}

