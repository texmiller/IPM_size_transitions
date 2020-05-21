####################################################
#  SHASH functions, from Jones and Pewsey 
#  epsilon is the skewness parameter (real)
#  delta controls tail weight (positive) 
#  delta < 1  <--> fatter than Gaussian 
####################################################

#####################################################################
# Moments of base SHASH distribution (without location/scale) 
# Jones and Pewsey, p. 764 
#####################################################################
require(Bessel); 

Pq = function(q) { # Jones and Pewsey, p. 764 
	a = exp(0.24)/sqrt(8*pi);
	B1 = BesselK(z = 0.25, nu = (q+1)/2)
	B2 = BesselK(z = 0.25, nu = (q-1)/2)
	return(a*(B1 + B2))
}	

m1B = function(epsilon,delta) {
	sinh(epsilon/delta)*Pq(1/delta)
}

m2B = function(epsilon,delta) {
		0.5*(cosh(2*epsilon/delta)*Pq(2/delta) - 1)
}		

m3B = function(epsilon,delta) {
	0.25*(sinh(3*epsilon/delta)*Pq(3/delta) - 3*sinh(epsilon/delta)*Pq(1/delta))
}	

m4B  = function(epsilon,delta) {
	0.125*(cosh(4*epsilon/delta)*Pq(4/delta) - 4*cosh(2*epsilon/delta)*Pq(2/delta) + 3) 
}	

varB  = function(epsilon,delta) {
	m2B(epsilon,delta) - m1B(epsilon,delta)^2
}

skewB = function(epsilon,delta) {
	M1 = m1B(epsilon,delta); 
	M2 = m2B(epsilon,delta); 
	M3 = m3B(epsilon,delta); 
	M3central = M3 - 3*M1*M2 + 2*M1^3;
	return( M3central/(varB(epsilon,delta)^(3/2)) ); 
}	

kurtB = function(epsilon,delta) {
	M1 = m1B(epsilon,delta); 
	M2 = m2B(epsilon,delta); 
	M3 = m3B(epsilon,delta); 
	M4 = m4B(epsilon,delta) 
	M4central = M4 - 4*M1*M3 + 6*M2*(M1^2) - 3*(M1^4); 
	return( M4central/(varB(epsilon,delta)^2) )
}	

nus = seq(-2,2,length=31); 
invdeltas = seq(0.5,2,length=30);
skews = kurts = matrix(0,31,30); 
for(j in 1:31){
for(k in 1:30){
    delta=1/invdeltas[k];
	epsilon=nus[j]; 
	skews[j,k]=skewB(epsilon,delta) 
	kurts[j,k]=kurtB(epsilon,delta)/3-1; 
}}	


graphics.off();
dev.new(); contour(nus,invdeltas,skews,xlab="nu",ylab="theta",main="Skew",labcex=1);

levels=round( seq(min(kurts),max(kurts),length=12), digits=1); 
dev.new(); contour(nus,invdeltas,kurts,xlab="nu",ylab="theta",main="Excess Kurtosis",labcex=1,
	levels=sort(c(0,levels))); 

#####################################################################
# Base density function 
#####################################################################
dJPB=function(x,epsilon,delta) {
		S = sinh(delta*asinh(x)-epsilon);
		C = cosh(delta*asinh(x)-epsilon);
		fac = sqrt(2*pi*(1+x^2));
		return(delta*C*exp(-0.5*S^2)/fac)
}

######## this should be the same as dSHASHo, and it is. 
eps=runif(1,-1,1); delta=1 + runif(1,-0.3,0.3);
plot(function(x) dJPB(x,eps,delta),-8,8); 
x=seq(-8,8,length=500); 
points(x, dSHASHo(x,mu=0,sigma=1,eps,delta),type="l",lty=2); 

integrate(function(x) x^2*dJP(x,0,2),-10,10,subdivisions=4000); 
m2B(0,2); 

#####################################################################
# Density with mu=0 parameterized by sigma, epsilon, and tau=1/delta 
#####################################################################
dJPBs = function(x,sigma,nu,tau,log=FALSE) {
		sigmaB=sqrt(varB(nu,tau))
		(sigmaB/sigma)*dSHASHo(sigmaB*x/sigma,nu,tau,log=log)
}		
plot(function(x) dJPBs(x,1,0.2,2),-3,3); 

#### test -- seems to be OK, relative to likely integration errors
sigma=sqrt(2.5); nu=runif(1,-.3,.3); tau=0.5+runif(1); 
mu = sigma*m1B(nu,tau); 
endpts = mu+seq(-50*sigma,50*sigma,length=50001); h=endpts[2]-endpts[1];
meshpts = endpts[-1]-h/2; 
h*sum(((meshpts-mu)^2)*dJPBs(meshpts,sigma,nu,tau));

#####################################################################
# Density parameterized by mu, sigma, nu= epsilon, and tau=delta 
#####################################################################
dSHASHsc = function(x,mu,sigma,nu,tau,log=FALSE) {
    if (any(sigma <= 0)) 
        stop(paste("sigma must be positive", "\n", ""))
    if (any(tau <= 0)) 
        stop(paste("tau must be positive", "\n", ""))

	JPsmean = sigma* m1B(nu,tau); 
	
	dJPBs(x - mu + JPsmean,sigma,nu,tau,log=log)
} 

mu = 0; sigma=sqrt(1.5); nu=runif(1,-.3,.3); tau=0.75+0.5*runif(1); 
endpts = mu+seq(-50*sigma,50*sigma,length=125000); h=endpts[2]-endpts[1];
meshpts = endpts[-1]-h/2; 
h*sum(meshpts*dSHASHsc(meshpts,mu,sigma,nu,tau));
mu; 

	

