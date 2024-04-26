x = sort(2*rbeta(1000,2,3)); u = rnorm(1000,mean=0,sd=1-0.2*abs(x-1));
plot(x,abs(u));

fit = smooth.spline(abs(u)~x,penalty=1.2);
points(fit$x,fit$y,type="l",col="red"); 
sd(fit$y)/mean(fit$y);
rsq.smooth.spline(x,abs(u)); 


