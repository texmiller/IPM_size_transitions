####################################################
#  Testing JP2 
####################################################

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
	cat( sd(z1),sd(z2),JP2sd(epsilon,delta), "     ", mean(z1) mean(z2), JP2mean(epsilon,delta), "\n");
}

## Testing the pdf 
for(j in 1:10) {
 epsilon=runif(1,-1,1); delta=exp(runif(1,-1,1)); x = rnorm(1,0,2);  
 u = integrate(function(x) dJP2(x,epsilon,delta),-Inf,x)$value  
 v = pJP2(x,epsilon,delta); 
 cat(u, v, "       ", log10(2*abs(u-v)/(u+v)),"\n"); 
}