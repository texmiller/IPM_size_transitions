library(quantreg); library(fda); library(sn); library(moments); 

x=seq(0,2,length=500); 
alphas = -3 + 10*x
y = rsn(length(x),xi=0,omega=1,alpha=alphas);
par(mfrow=c(2,1)); 
plot(x,y,type="p"); 

y1 = rsn(100000,xi=0,omega=1,alpha=min(alphas));
y2 = rsn(100000,xi=0,omega=1,alpha=max(alphas));
skewness(y1); skewness(y2); skewness(c(y1,y2)); 

B=create.bspline.basis(rangeval=range(x), nbasis=4, norder=3);
X = eval.basis(B,x);

f.50 <- rq(y ~ X-1, tau=0.5)
f.90 <- rq(y ~ X-1, tau=0.9)
f.10 <- rq(y ~ X-1, tau=0.1)
q.50 = fitted(f.50); q.90 = fitted(f.90); q.10 = fitted(f.10);
matpoints(x,cbind(q.10, q.50, q.90), type="l",lty=2, col="blue"); 
NPS = (q.10 + q.90 - 2*q.50)/(q.90 - q.10); 

boot.10=boot.rq(X,y,tau=0.1, R=501);
boot.50=boot.rq(X,y,tau=0.5,  R=501, U=boot.10$U);
boot.90=boot.rq(X,y,tau=0.9, R=501, U=boot.10$U);

bootq.10=X%*%t(boot.10$B);
bootq.50=X%*%t(boot.50$B);
bootq.90=X%*%t(boot.90$B);

bootNPS = (bootq.10 + bootq.90 - 2*bootq.50)/(bootq.90 - bootq.10); 
LS = apply(bootNPS,1,function(x) quantile(x,0.05)); 
MS = apply(bootNPS,1,function(x) quantile(x,0.5));
US = apply(bootNPS,1,function(x) quantile(x,0.95)); 
matplot(x,cbind(LS,MS,US),type="l",col="red",lty=2,lwd=2); 
points(x,NPS,type="l",col="black",lwd=2); 
abline(0,0,col="blue",lty=2); 

b0 = apply(bootNPS,2,mean); 
quantile(b0,c(0.01, 0.05,0.95,0.99)); 

b1 = apply(bootNPS,2,function(s) lm(s~x)$coef[2]);
quantile(b1,c(0.01, 0.05,0.95,0.99)); 

b2 = apply(bootNPS,2,function(s) lm(s~x + I(x^2))$coef[3]);
quantile(b2,c(0.01, 0.05,0.95,0.99)); 

require(mgcv); 
gam_fit = gam(list(y~1,~1, ~s(x), ~1), family="shash"); 
summary(gam_fit); 