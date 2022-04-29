library(quantreg); library(fda); library(sn); library(moments); 

### A qsn function that actually works, unlike the one in the sn package. 
### Vectorized in p, but not the distribution parameters 
my.qsn = function(p,xi,omega,alpha) {
    px = seq(-50,50,length=500); 
    py = psn(px,xi=xi,omega=omega,alpha=alpha,tau=0);
    F = splinefun(py,px,method="monoH.FC"); 
    return(F(p))
}  
 
## AIC model averaging of df=(1,2,3,4) quantile regression
## df=3 or 4 are natural B-spline bases including intercept.  
## Model df are over-weighted by default, as recommended for 
## spline smoothing in ordinary regression. Here it is just
## a weak hedge against overfitting.
rqAIC = function(x,y,tau, L=NULL, U=NULL, gamma=1.4) {
    e = order(x); x=x[e]; y=y[e]; 

  # fit constant (df=1)
    fit.1 = rq(y~1, tau = tau); 

  # fit linear (df=2)
     fit.2 = rq(y~x, tau = tau); 

  # fit with df=3
    if(is.null(U)) U = max(x);
    if(is.null(L)) L = min(x) 
    X3 = ns(x,df=3, intercept=TRUE, Boundary.knots=c(L,U)); 
    fit.3 = rq(y ~ X3-1, tau=tau) 

  # fit with df=4
    X4 = ns(x,df=4, intercept=TRUE, Boundary.knots=c(L,U)); 
    fit.4 = rq(y ~ X4-1, tau=tau) 
   
    AICs = c(AIC(fit.1)[1],AIC(fit.2)[1], AIC(fit.3)[1], AIC(fit.4)[1]); 
  # For AIC weights with over-weighted model df 
    z = AICs + 2*(gamma-1)*c(1,2,3,4)  
    w = exp(-0.5*(z-mean(z))); w = w/sum(w); 
    
    px = seq(L,U,length=100); 
    X3 = ns(px,df=3, intercept=TRUE, Boundary.knots=c(L,U)); 
    X4 = ns(px,df=4, intercept=TRUE, Boundary.knots=c(L,U)); 
    qhat = w[1]*coef(fit.1)[1] + w[2]*(coef(fit.2)[1]+px*coef(fit.2)[2])  + w[3]*X3%*%coef(fit.3) + w[4]*X4%*%coef(fit.4); 
    qfun = splinefun(px,qhat,method="natural"); 
    return(qfun); 
} 
 
x=seq(0,2,length=250); alphas = -3 + 3*x
y = rsn(length(x),xi=0,omega=1,alpha=alphas);
plot(x,y,type="p"); 

y1 = rsn(100000,xi=0,omega=1,alpha=min(alphas));
y2 = rsn(100000,xi=0,omega=1,alpha=max(alphas));
skewness(y1); skewness(y2); skewness(c(y1,y2)); 

qfun = rqAIC(x,y,0.1,L=0,U=2); q.10 = qfun(x); points(x,q.10,type="l"); 
qfun = rqAIC(x,y,0.5,L=0,U=2); q.50 = qfun(x); points(x,q.50,type="l"); 
qfun = rqAIC(x,y,0.9,L=0,U=2); q.90 = qfun(x); points(x,q.90,type="l"); 

Y = matrix(NA,length(alphas),3)
for(i in 1:length(alphas)) {
    out = my.qsn(c(0.1,0.5,0.9), xi=0, omega=1, alpha=alphas[i]);
    Y[i,]=out; 
    if(i%%20==0) cat(i, "\n"); 
} 
matpoints(x,Y,type="l",lty=2, col="blue"); 

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