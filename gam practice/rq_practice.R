x=seq(0,2,length=500); 
y = rsn(length(x),xi=0,omega=1,alpha=-3+x);
par(mfrow=c(2,1)); 
plot(x,y,type="p"); 

B=create.bspline.basis(rangeval=range(x), nbasis=4, norder=4);
X = eval.basis(B,x);

f.50 <- rq(y ~ X-1, tau=0.5)
f.90 <- rq(y ~ X-1, tau=0.95)
f.10 <- rq(y ~ X-1, tau=0.05)
q.50 = fitted(f.50); q.90 = fitted(f.90); q.10 = fitted(f.10);
matpoints(x,cbind(q.10, q.50, q.90), type="l",lty=2, col="blue"); 
NPS = (q.10 + q.90 - 2*q.50)/(q.90 - q.10); 



boot.10=boot.rq(X,y,tau=0.05, R=501);
boot.50=boot.rq(X,y,tau=0.5,  R=501, U=boot.10$U);
boot.90=boot.rq(X,y,tau=0.95, R=501, U=boot.10$U);

bootq.10=X%*%t(boot.10$B);
bootq.50=X%*%t(boot.50$B);
bootq.90=X%*%t(boot.90$B);

bootNPS = (bootq.10 + bootq.90 - 2*bootq.50)/(bootq.90 - bootq.10); 
LS = apply(bootNPS,1,function(x) quantile(x,0.05)); 
US = apply(bootNPS,1,function(x) quantile(x,0.95)); 
matplot(x,cbind(LS,US),type="l",col="red",lty=2); points(x,NPS,type="l",col="black"); 
abline(0,0,col="blue",lty=2); 

