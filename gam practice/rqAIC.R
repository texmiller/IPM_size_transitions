library(quantreg); library(fda); library(sn); library(moments); 

### A qsn function that actually works, unlike the one in the sn package. 
### Vectorized in p, but not the distribution parameters 
my.qsn = function(p,xi,omega,alpha) {
    px = seq(-10,120,length=250); 
    py = psn(px,xi=xi,omega=omega,alpha=alpha,tau=0);
    F = splinefun(py,px,method="monoH.FC"); 
    return(F(p))
}  

rqAIC = function(x,y,tau, L=NULL, U=NULL) {
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
    
    # Compute AIC weights 
    z = c(AIC(fit.1)[1],AIC(fit.2)[1], AIC(fit.3)[1], AIC(fit.4)[1]); 
    w = exp(-0.5*(z-mean(z))); w = w/sum(w); 
    
    px = seq(L,U,length=100); 
    X3 = ns(px,df=3, intercept=TRUE, Boundary.knots=c(L,U)); 
    X4 = ns(px,df=4, intercept=TRUE, Boundary.knots=c(L,U)); 
    qhat = w[1]*coef(fit.1)[1] + w[2]*(coef(fit.2)[1]+px*coef(fit.2)[2])  + w[3]*X3%*%coef(fit.3) + w[4]*X4%*%coef(fit.4); 
    qfun = splinefun(px,qhat,method="natural"); 
    return(qfun); 
} 
 
x=seq(0,2,length=250); alphas = -3 + 5*x
y = rsn(length(x),xi=0,omega=1,alpha=alphas);
plot(x,y,type="p"); 

qfun = rqAIC(x,y,0.25,L=0,U=2); points(x,qfun(x),type="l"); 
qfun = rqAIC(x,y,0.5,L=0,U=2); points(x,qfun(x),type="l"); 
qfun = rqAIC(x,y,0.75,L=0,U=2); points(x,qfun(x),type="l"); 

py = qsn(0.25,xi=0, omega=1, alpha=alphas,solver="RFB"); points(x,py,type="l",col="blue"); 



y1 = rsn(100000,xi=0,omega=1,alpha=min(alphas));
y2 = rsn(100000,xi=0,omega=1,alpha=max(alphas));
skewness(y1); skewness(y2); skewness(c(y1,y2)); 

# B=create.bspline.basis(rangeval=range(x), nbasis=4, norder=4);
# X = eval.basis(B,x);

X = ns(x,df=4, intercept=TRUE, Boundary.knots=c(0,2)); 
matplot(x,X,type="l"); 

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