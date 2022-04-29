library(quantreg); library(fda); library(sn); library(moments); library(mgcv); 

### A qsn function that actually works, unlike the one in the sn package. 
### Vectorized in p, but not the distribution parameters 
my.qsn = function(p,xi,omega,alpha) {
    px = seq(-50,50,length=500); 
    py = psn(px,xi=xi,omega=omega,alpha=alpha,tau=0);
    F = splinefun(py,px,method="monoH.FC"); 
    return(F(p))
}  
 
x = sort(2*rbeta(500,3,3)); 
alphas = 5 - 5*x; 
Y = matrix(NA,length(alphas),3)
for(i in 1:length(alphas)) {
    out = my.qsn(c(0.2,0.5,0.8), xi=0, omega=1, alpha=alphas[i]);
    Y[i,]=out; 
    if(i%%20==0) cat(i, "\n"); 
} 

graphics.off(); 
y = rsn(length(x),xi=0,omega=1,alpha=alphas); 
xydata = data.frame(x = x, y = y); 



plot(x, y-qo.50,type="p");

y1 = rsn(100000,xi=0,omega=1,alpha=min(alphas));
y2 = rsn(100000,xi=0,omega=1,alpha=max(alphas));
skewness(y1); skewness(y2); skewness(c(y1,y2)); 

xycdata = data.frame(x = x, yc = y - qo.50); 
S.10 = qgam(yc~s(x),data=xycdata,qu=0.2,argGam=list(gamma=1.4)); q.10 = predict(S.10,newdata=xycdata); points(x,q.10,type="l"); 
S.50 = qgam(yc~s(x),data=xycdata,qu=0.5,argGam=list(gamma=1.4)); q.50 = predict(S.50,newdata=xycdata); points(x,q.50,type="l"); 
S.90 = qgam(yc~s(x),data=xycdata,qu=0.8,argGam=list(gamma=1.4)); q.90 = predict(S.90,newdata=xycdata); points(x,q.90,type="l"); 

matpoints(x, Y - cbind(qo.50,qo.50,qo.50), type="l",lty=2, col="blue"); 

NPS.hat = (q.10 + q.90 - 2*q.50)/(q.90 - q.10); 
NPS = (Y[,1]+Y[,3]-2*Y[,2])/(Y[,3]-Y[,1]); 
dev.new(); 
matplot(x,cbind(NPS,NPS.hat),type="l", lty=c(1,2)); 




X = ns(x,df=4, intercept=TRUE, Boundary.knots=c(0,2)); 

boot.10=boot.rq(X,y,tau=0.1, R=501);
boot.50=boot.rq(X,y,tau=0.5,  R=501, U=boot.10$U);
boot.90=boot.rq(X,y,tau=0.9, R=501, U=boot.10$U);

bootq.10=X%*%t(boot.10$B);
bootq.50=X%*%t(boot.50$B);
bootq.90=X%*%t(boot.90$B);

matplot(x, cbind(bootq.10,bootq.50, bootq.90), type="l"); 

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
gam_fit = gam(list(y~s(x),~s(x), ~s(x), ~s(x)), family="shash"); 
summary(gam_fit); 