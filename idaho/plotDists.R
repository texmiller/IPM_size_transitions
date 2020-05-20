tau=seq(0.2,1,length=90); 
kurt = numeric(length(tau));
for(j in 1:length(tau)){
	z = rSHASHo2(120000,mu=0,sigma=4,nu=0,tau=tau[j])
	kurt[j]=skewness(z)
}
plot(tau,kurt); abline(h=3); 


nus = seq(-2,2,length=21); 
taus = seq(0.5,5,length=20);
skews = kurts = matrix(0,21,20); 
for(j in 1:21){
for(k in 1:20){
	zmu=integrate(function(x) x*dSHASHo(x,mu=0,sigma=1,nu=nus[j],tau=taus[k]), -Inf, Inf)$value; 
	zvar = integrate(function(x) ((x-zmu)^2)*dSHASHo2(x,mu=0,sigma=1,nu=nus[j],tau=taus[k]), -Inf, Inf)$value; 
	z3 = integrate(function(x) ((x-zmu)^3)*dSHASHo2(x,mu=0,sigma=1,nu=nus[j],tau=taus[k]), -Inf, Inf)$value; 
	z4 = integrate(function(x) ((x-zmu)^4)*dSHASHo2(x,mu=0,sigma=1,nu=nus[j],tau=taus[k]), -Inf, Inf)$value; 
	sigma=sqrt(zvar);
	skew = z3/(sigma^(3/2)); 
	kurt = z4/(sigma^4)
	# cat(nus[j],taus[k],zmu,sigma,skew,kurt,"\n");
	skews[j,k]=skew; kurts[j,k]=kurt;
}}	
contour(nus,taus,skews,xlab="nu",ylab="tau");
contour(nus,taus,kurts,xlab="nu",ylab="tau");

	
nus = seq(-2,2,length=20); 
taus = seq(1.25,4,length=20);
skews = kurts = sigmas = matrix(0,20,20); 
for(j in 1:20){
for(k in 1:20){
	zmu=integrate(function(x) x*dJSU(x,mu=0,sigma=1,nu=nus[j],tau=taus[k]), -12, 12)$value; 
	zvar = integrate(function(x) ((x-zmu)^2)*dJSU(x,mu=0,sigma=1,nu=nus[j],tau=taus[k]), -Inf, Inf)$value; 
	z3 = integrate(function(x) ((x-zmu)^3)*dJSU(x,mu=0,sigma=1,nu=nus[j],tau=taus[k]), -Inf, Inf)$value; 
	z4 = integrate(function(x) ((x-zmu)^4)*dJSU(x,mu=0,sigma=1,nu=nus[j],tau=taus[k]), -Inf, Inf)$value; 
	sigma=sqrt(zvar);
	skew = z3/(sigma^(3/2)); 
	kurt = z4/(sigma^4)
	# cat(nus[j],taus[k],zmu,sigma,skew,kurt,"\n");
	skews[j,k]=skew; kurts[j,k]=kurt;
	sigmas[j,k]=sigma; 
}}	
contour(nus,-1/taus,skews,xlab="nu",ylab="tau"); 

contour(nus,-1/taus,kurts,xlab="nu",ylab="tau"); 

	


