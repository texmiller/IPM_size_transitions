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
dJP09B=function(x,epsilon,delta) {
		if (any(delta <= 0)) stop("need tau > 0");
		S = sinh(delta*asinh(x)-epsilon);
		C = cosh(delta*asinh(x)-epsilon);
		fac = sqrt(2*pi*(1+x^2));
		return(delta*C*exp(-0.5*S^2)/fac)
}

######## this should be the same as dSHASHo with mu=0, sigma=1, and it is. 
eps=runif(1,-1,1); delta=1 + runif(1,-0.3,0.3);
plot(function(x) dJP09B(x,eps,delta),-8,8); 
x=seq(-8,8,length=500); 
points(x, dSHASHo(x,mu=0,sigma=1,eps,delta),type="l",lty=2); 

integrate(function(x) (x^2)*dJP09B(x,0,2),-10,10,subdivisions=4000); 
m2B(0,2); 

#####################################################################
# Full SHASH density function, using gamlss notation 
# This should be identical to SHASHo 
#####################################################################
dJP09=function(x,mu,sigma,nu,tau) { 
   if (any(sigma <= 0)) stop("need sigma > 0")
   if (any(tau <= 0)) stop("need tau > 0")
   return( dJP09B((x-mu)/sigma, epsilon=nu, delta=tau)/sigma )
}


xvals=seq(-6,6,length=1000);
gof = function(pars) {
	y1 = dSHASHo(xvals,pars[1],exp(pars[2]),pars[3],exp(pars[4]))
	err = (y1-y2)^2; 
	return(sum(err))
}	

df=12; 
#y2 = dlnorm(xvals,mu=0,sigma=1,nu=6,tau=df);
y2 = 0.5*exp(-abs(xvals)^1.2);
out=optim(par=c(0,0,0,0),fn=gof,method="Nelder-Mead",control=list(maxit=10000,trace=4))$par; 
y1 = dSHASHo(xvals,out[1],exp(out[2]),out[3],exp(out[4]));

matplot(xvals,cbind(y1,y2),col=c("black","blue"),type="l",lty=c(1,2)); 
legend("topright",legend=c("SHASH","True"),col=c("black","blue"),lty=c(1,2))
abline(h=0,col="black",lty=3); 


xvals=seq(-6,6,length=1000);
gof = function(pars) {
	y1 = dsgt(xvals,pars[1],exp(pars[2]),2*pnorm(pars[3])-1,exp(pars[4]),exp(pars[5]))
	err = (y1-y2)^2; 
	return(sum(err))
}	

df=12; 
y2 = 0.5*exp(-abs(xvals)^1.2);
pars=optim(par=c(0,0,0,1,2),fn=gof,method="Nelder-Mead",control=list(maxit=10000,trace=4))$par; 
y1 = dsgt(xvals,pars[1],exp(pars[2]),2*pnorm(pars[3])-1,exp(pars[4]),exp(pars[5]))
matplot(xvals,cbind(y1,y2),col=c("black","blue"),type="l",lty=c(1,2)); 
legend("topright",legend=c("SHASH","True"),col=c("black","blue"),lty=c(1,2))
abline(h=0,col="black",lty=3); 











