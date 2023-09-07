require(gamlss); 

f1 = function(x) dSHASHo2(x,mu=0,sigma=1,nu=0,tau=1); 
f2 = function(x) dSHASHo2(x,mu=0,sigma=1,nu=-1,tau=1); 

mu1 = integrate(function(x) x*f1(x), -20, 20)$value; mu1; 
mu2 = integrate(function(x) x*f2(x), -20, 20)$value; mu2; 
sigma = sqrt(integrate(function(x) (x^2)*f2(x), -20, 20)$value - mu2^2); sigma; 


f3 = function(x) dSST(x,mu=0,sigma=2,nu=10,tau=15); 
f4 = function(x) dSST(x,mu=0,sigma=2,nu=.1,tau=15); 

mu3 = integrate(function(x) x*f3(x), -20, 20)$value; mu3; 
mu4 = integrate(function(x) x*f4(x), -20, 20)$value; mu4; 
sigma = sqrt(integrate(function(x) (x^2)*f4(x), -20, 20)$value - mu4^2); sigma; 

