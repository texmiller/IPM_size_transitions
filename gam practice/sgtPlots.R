require(sgt); 
######### Visualize effects of changing p and q 
par(mfrow=c(3,1)); 
x=seq(-5,5,length=500); 

py1 = dsgt(x,mu=0,sigma=1,lambda=0,p=1,q=2.2);
py2 = dsgt(x,mu=0,sigma=1,lambda=0,p=1,q=4);
py3 = dsgt(x,mu=0,sigma=1,lambda=0,p=1,q=20);
matplot(x,cbind(py1,py2,py3),type="l",col=c("black","red","blue")); 

py1 = dsgt(x,mu=0,sigma=1,lambda=0,p=2,q=1.2);
py2 = dsgt(x,mu=0,sigma=1,lambda=0,p=2,q=4);
py3 = dsgt(x,mu=0,sigma=1,lambda=0,p=2,q=20);
matplot(x,cbind(py1,py2,py3),type="l",col=c("black","red","blue")); 


x=seq(-5,5,length=500); 
py1 = dsgt(x,mu=0,sigma=1,lambda=0,p=5,q=.5);
py2 = dsgt(x,mu=0,sigma=1,lambda=0,p=5,q=4);
py3 = dsgt(x,mu=0,sigma=1,lambda=0,p=5,q=20);
matplot(x,cbind(py1,py2,py3),type="l",col=c("black","red","blue")); 


######### Create covariate for residuals 
require(fda); require(sgt); require(maxLik); 
z = rt(500,df=10); z=sort(z); hist(z); 

########### Create artificial "residuals" with known sgt parameters 
resids = rsgt(length(z),mu=0,sigma=1,lambda= 1-exp(z/10), p = 2 + 0.1*z, q=3 - 0.1*z); 
hist(resids); plot(z,resids); 

########## Create B-spline basis, and plot it 
B = create.bspline.basis(rangeval=range(z), norder=4, nbasis=6)
x = seq(min(z),max(z),length=200); 
out = eval.basis(x,B);   
matplot(x,out,type="l"); 


########### function to convert coefficients into lambda, p, q values 
X = eval.basis(z,B); P = bsplinepen(B); 

notExp2 = function (x, d, b=1/d) exp(d * sin(x * b))  # from mgcv

make_lpq = function(pars){
    u = X%*%pars[1:6]; lambda = -1 + 2*exp(u)/(1+exp(u));
    u = X%*%pars[7:12]; pz = notExp2(u,d=log(50)); 
    u = X%*%pars[13:18]; qz = 0.5 + notExp2(u,d=log(100));
    
    # evaluate penalties 
    pars=matrix(pars,ncol=1); 
    Plam = t(pars[1:6])%*%P%*%pars[1:6]
    Pp = t(pars[7:12])%*%P%*%pars[7:12]
    Pq = t(pars[13:18])%*%P%*%pars[13:18]
    
    df1 = lambda2df(pars[1:6],B)
    
    r = rep(1,length(pz));  
    e = which(pz*qz < 2.01); r[e]=sqrt(2.05/(pz[e]*qz[e]));
    pz=pz*r; qz = qz*r; 
    return(list(lambda=lambda,p=pz,q=qz,Plam=Plam,Pp=Pp,Pq=Pq)); 
 }
 
 NegLogLik = function(pars){
    out=make_lpq(pars);
    NLL = -sum(dsgt(z,mu=0,sigma=1,lambda=out$lambda,p=out$p,q=out$q,log=TRUE))  
    return(NLL)
}    
 
 PenNegLogLik = function(pars,pen) {
     out = make_lpq(pars); 
     NLL = NegLogLik(pars); 
     PNLL = NLL + pen[1]*out$Plam + pen[2]*out$Pp + pen[3]*out$Pq
     return(PNLL)
 } 

 evalAIC = function(logpen){
        pen=exp(logpen); 
        fit = optim(par=rep(0,18),fn=PenNegLogLik,control=list(maxit=10000,trace=0),pen=pen)
        pars=fit$par; 
        edf1 = lambda2df(argvals=z,basisobj=B,lambda=pen[1])
        edf2 = lambda2df(argvals=z,basisobj=B,lambda=pen[2])
        edf3 = lambda2df(argvals=z,basisobj=B,lambda=pen[3])
        AICval = NegLogLik(pars) + 2*(edf1+edf2+edf3);
        return(AICval); 
 }       
evalAIC(c(-1,0,-2));


 

## log likelihood function for maxLik 
sgt_logLik = function(pars) {
    lpq = make_lpq(pars); 
    lik = dsgt(resids, mu = 0, sigma = 1, lambda = lpq$lambda, p = lpq$p, q = lpq$q, log = TRUE)
    return(lik)    
}

## log likelihood function for optim 
sgt_logLik1 = function(pars) {
    val = sum(sgt_logLik(pars))
    if(!is.finite(val)) val = -(10^16)
    return(-val)
}  





out=optim(par=rep(1,9),sgt_logLik1,control=list(trace=4,maxit=25000)); 
out$estimate=out$par; 

out = maxLik(sgt_logLik, start = rep(-3,9), method="BHHH",control=list(printLevel=2)); 
for(k in 1:5) {
    #out = optim(par=out$estimate,sgt_logLik1,control=list(trace=4,maxit=25000)); 
    out = maxLik(sgt_logLik, start = out$estimate, method="BHHH",control=list(printLevel=2)); 
  
}    

fit = make_lpq(out$estimate); 
matplot(x, cbind(fit$lambda,fit$p,fit$q), type="l",lty=1, col=c("black","red","blue")); 
 