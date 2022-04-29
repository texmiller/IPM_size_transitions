library(quantreg); library(fda); library(sn); library(moments); library(mgcv); 

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
rqAIC = function(x,y,tau, L=NULL, U=NULL, gamma=2) {
    e = order(x); x=x[e]; y=y[e]; 
    if(is.null(U)) U = max(x); if(is.null(L)) L = min(x) 

  # fit constant (df=1)
    fit.1 = rq(y~1, tau = tau); 

  # fit linear (df=2)b
     fit.2 = rq(y~x, tau = tau); 

  # fit with df=3
    B3=create.bspline.basis(rangeval=c(L,U), nbasis=3, norder=2);
    X3 = eval.basis(B3,x);
    # X3 = ns(x,df=3, intercept=TRUE, Boundary.knots=c(L,U)); 
    fit.3 = rq(y ~ X3-1, tau=tau) 

  # fit with df=4
    B4=create.bspline.basis(rangeval=c(L,U), nbasis=4, norder=3);
    X4 = eval.basis(B4,x);
    # X4 = ns(x,df=4, intercept=TRUE, Boundary.knots=c(L,U)); 
    fit.4 = rq(y ~ X4-1, tau=tau) 
    
  # fit with df=5
    B5=create.bspline.basis(rangeval=c(L,U), nbasis=5, norder=4);
    X5 = eval.basis(B5,x);
    # X4 = ns(x,df=4, intercept=TRUE, Boundary.knots=c(L,U)); 
    fit.5 = rq(y ~ X5-1, tau=tau)     
   
    AICs = c(AIC(fit.1)[1],AIC(fit.2)[1], AIC(fit.3)[1], AIC(fit.4)[1],AIC(fit.5)[1]); 
  # For AIC weights with over-weighted model df 
    z = AICs + 2*(gamma-1)*c(1:5)  
    w = exp(-0.5*(z-mean(z))); w = w/sum(w); 
    
    px = seq(L,U,length=200); 
    # X3 = ns(px,df=3, intercept=TRUE, Boundary.knots=c(L,U)); 
    # X4 = ns(px,df=4, intercept=TRUE, Boundary.knots=c(L,U)); 
    X3 = eval.basis(B3,px); 
    X4 = eval.basis(B4,px); 
    X5 = eval.basis(B5,px);
    
    qhat = w[1]*coef(fit.1)[1] + w[2]*(coef(fit.2)[1]+px*coef(fit.2)[2]) + w[3]*X3%*%coef(fit.3) + w[4]*X4%*%coef(fit.4) + w[5]*X5%*%coef(fit.5)
    qfun = approxfun(px,qhat,rule=1); 
    return(list(qfun=qfun,wts=w)); 
} 

x = sort(runif(200,0,2)); 
alphas = 5 - 5*x; 
Y = matrix(NA,length(alphas),3)
for(i in 1:length(alphas)) {
    out = my.qsn(c(0.1,0.5,0.9), xi=0, omega=1, alpha=alphas[i]);
    Y[i,]=out; 
    if(i%%20==0) cat(i, "\n"); 
} 

graphics.off(); 
y = rsn(length(x),xi=0,omega=1,alpha=alphas); 

fit_gam = gam(y~s(x), gamma=1.4); 
plot(x,y); points(x, fitted(fit_gam)); 

y = y - fitted(fit_gam); 
plot(x,y,type="p"); 

y1 = rsn(100000,xi=0,omega=1,alpha=min(alphas));
y2 = rsn(100000,xi=0,omega=1,alpha=max(alphas));
skewness(y1); skewness(y2); skewness(c(y1,y2)); 

S.10 = rqAIC(x,y,0.1,L=0,U=2); q.10 = S.10$qfun(x); points(x,q.10,type="l"); 
S.50 = rqAIC(x,y,0.5,L=0,U=2); q.50 = S.50$qfun(x); points(x,q.50,type="l"); 
S.90 = rqAIC(x,y,0.9,L=0,U=2); q.90 = S.90$qfun(x); points(x,q.90,type="l"); 


matpoints(x,Y-fitted(fit_gam),type="l",lty=2, col="blue"); 

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