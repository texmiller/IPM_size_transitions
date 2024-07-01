#########################################################
#  Testing JP2(epsilon, delta)  
#########################################################

## Testing the pdf and mean 
for(j in 1:10) {
 epsilon=runif(1,-1,1); delta=exp(runif(1,-1,1)); 
 u = integrate(function(x) dJP2(x,epsilon,delta),-Inf,Inf)$value # should equal 1. 
 v = integrate(function(x) x*dJP2(x,epsilon,delta), -Inf, Inf)$value   # should equal JP2mean 
cat(round(u,digits=6), round(v,digits=6), JP2mean(epsilon,delta), "\n"); 
}

## Testing the RNG and variance 
for(j in 1:10) {
	epsilon=runif(1,-2,2); delta=exp(runif(1,-2,2)); 
	z1 = rJP2(200000, epsilon, delta)
	z2 = rJP2(200000, epsilon, delta)
	cat( sd(z1),sd(z2),JP2sd(epsilon,delta), "     ", mean(z1), mean(z2), JP2mean(epsilon,delta), "\n");
}

## Testing the cdf 
for(j in 1:10) {
 epsilon=runif(1,-1,1); delta=exp(runif(1,-1,1)); x = rnorm(1,0,2);  
 u = integrate(function(x) dJP2(x,epsilon,delta),-Inf,x)$value  
 v = pJP2(x,epsilon,delta); 
 cat(u, v, "       ", log10(2*abs(u-v)/(u+v)),"\n"); 
}

## Comparing RNG and cdf 
for(j in 1:10) {
	epsilon=runif(1,-2,2); delta=exp(runif(1,-2,2)); 
	z1 = rJP2(500000, epsilon, delta)
	q = runif(1)
	x = quantile(z1,q);
	v = pJP2(x,epsilon,delta); 	
	cat(q,v,"\n"); 
}

## Testing the quantile function 
for(j in 1:10) {
	epsilon=runif(1,-2,2); delta=exp(runif(1,-2,2)); 
	z1 = rJP2(500000, epsilon, delta);
	p = runif(1); 
	qhat = quantile(z1,p); 
	q = qJP2(p,epsilon,delta)
 	cat(qhat,q,"\n"); 
}
	
####################################################
#  Testing JPS(mean, sd, epsilon, delta)  
####################################################	

## Testing the pdf and mean 
for(j in 1:10) {
 mu = 2*rnorm(1); sigma = exp(runif(1,-3,3)); epsilon=runif(1,-2,2); delta=exp(runif(1,-2,2)); 
 u = integrate(function(x) dJPS(x,mean=mu,sd=sigma,epsilon,delta),-Inf,Inf)$value # should equal 1. 
 v = integrate(function(x) x*dJPS(x,mean=mu,sd=sigma,epsilon,delta), -Inf, Inf)$value   # should equal JP2mean 
cat(round(u,digits=6), "     ", round(v,digits=6), mu, "\n"); 
}	

## Testing the RNG and variance 
foo = matrix(NA,25,6); 
for(j in 1:25) {
	mu = 2*rnorm(1); sigma = exp(runif(1,-3,3)); epsilon=runif(1,-2,2); delta=exp(runif(1,-2,2)); 
	z1 = rJPS(500000, mu, sigma, epsilon, delta)
	z2 = rJPS(500000, mu, sigma, epsilon, delta)
	cat( sd(z1),sd(z2),sigma, "     ", mean(z1), mean(z2), mu, "\n");
	foo[j,] = c(sd(z1),sd(z2),sigma, mean(z1), mean(z2), mu);
}
par(mfrow=c(2,1));
matplot(foo[,3],foo[,1:2],xlab="sigma",ylab="s-hat"); abline(0,1); 
matplot(foo[,6],foo[,4:5],xlab = "mu",ylab = "x_bar"); abline(0,1); 


## Testing the cdf 
for(j in 1:10) {
 mu = rnorm(1); sigma = exp(runif(1,-1,1));
 epsilon=runif(1,-1,1); delta=exp(runif(1,-1,1)); x = rnorm(1,0,1);  
 u = integrate(function(u) dJPS(u,mu,sigma,epsilon,delta),-Inf,x)$value  
 v = pJPS(x,mu,sigma,epsilon,delta); 
 cat(u, v, "       ", log10(2*abs(u-v)/(u+v)),"\n"); 
}

## Comparing RNG and cdf 
for(j in 1:10) {
	mu = 2*rnorm(1); sigma = exp(runif(1,-1,1));
	epsilon=runif(1,-2,2); delta=exp(runif(1,-2,2)); 
	z1 = rJPS(500000, mu,sigma, epsilon, delta)
	q = runif(1); x = quantile(z1,q);   # fraction q of random values are <x. 
	v = pJPS(x,mu,sigma,epsilon,delta); # cdf of x (should be the same). 	
	cat(q,v,"\n"); 
}
	
## Testing the quantile function 
for(j in 1:10) {
	mu = 2*rnorm(1); sigma = exp(runif(1,-1,1));
	epsilon=runif(1,-2,2); delta=exp(runif(1,-2,2)); 
	z1 = rJPS(500000, mu,sigma, epsilon, delta);
	p = runif(1); 
	qhat = quantile(z1,p); 
	q = qJPS(p,mu,sigma,epsilon,delta)
 	cat(qhat,q,"\n"); 
}
	
