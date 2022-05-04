######### Visualize effects of changing nu and tau 
par(mfrow=c(3,1),bty="l",mar=c(4,4,1,1)); 
x=seq(-5,5,length=500); 
py1 = dSHASHo2(x,mu=0,sigma=1,nu=0,tau=1) 
py2 = dSHASHo2(x,mu=0,sigma=1,nu=0,tau=.25)
py3 = dSHASHo2(x,mu=0,sigma=1,nu=0,tau=4)
matplot(x,cbind(py1,py2,py3),type="l",col=c("black","red","blue"),main="Vary tau"); 


py1 = dSHASHo2(x,mu=0,sigma=1,nu=-1,tau=1) 
py2 = dSHASHo2(x,mu=0,sigma=1,nu=0,tau=1)
py3 = dSHASHo2(x,mu=0,sigma=1,nu=1,tau=1)
matplot(x,cbind(py1,py2,py3),type="l",col=c("black","red","blue"),main="Vary nu"); 

########### Create artificial "residuals" with known parameters 
x = rnorm(500); x = sort(x); 
resids = rSHASHo2(length(x),mu=0,sigma=1,nu= 0.1*x^3 + (x^2)/5,tau=0.7); 
hist(resids); plot(x,resids); 


pbc = cs.control(cv=FALSE,nknots=5,penalty=2); 
fit = gamlss(resids~1, sigma.formula=~1, 
            mu.start=0, mu.fix=TRUE, sigma.start=1,sigma.fix=TRUE,
            nu.formula=~cs(x,control=pbc), 
            tau.formula=~cs(x,control=pbc), family=SHASHo2,method=CG(20))

val = numeric(100); val[1]=fit$G.deviance; fits=list(100); fits[[1]]=fit;          
for(j in 2:100){
        cat("<-------- Starting next run -------->","\n");
fit2 = gamlss(resids~1, sigma.formula=~1, 
            mu.start=0, mu.fix=TRUE, sigma.start=1,sigma.fix=TRUE,
            nu.formula=~cs(x,control=pbc), 
            tau.formula=~cs(x,control=pbc), family=SHASHo2,method=CG(10),
            start.from=fit)
            fit=fit2; val[j]=fit$G.deviance; fits[[j]]=fit; 
}            

jmin = which(val==min(val))[1]; fit = fits[[val]]



Y = cbind(fit$mu.fv,fit$sigma.fv,fit$nu.fv,fit$tau.fv); 
matplot(x,Y,type="l",col=c("black","red","blue","green3"), lty=1,lwd=2,xlab="x values", ylab="fitted pars"); 
abline(h=c(0, 1,0.7), col=c("black","red","green3")) 
nu=0.1*x^3 + x^2/5; points(x,nu,type="l",lty=1,col="blue"); 


L=min(x); U=max(x); 
B = create.bspline.basis(rangeval=c(L,U), norder=2, nbasis=2)
X = eval.basis(B,x); 


make_args = function(pars){
    mu = X%*%pars[1:3]; 
    sigma = exp(X%*%pars[4:6]);  
    nu = X%*%pars[7:9]; 
    tau = exp(X%*%pars[10:12]); 
    return(list(mu=mu,sigma=sigma,nu=nu,tau=tau)); 
 }
 

## log likelihood function for maxLik 
SHASH_logLik = function(pars) {
    args = make_args(pars); 
    lik = dSHASHo2(resids, mu = args$mu, sigma = args$sigma, nu=args$nu, tau=args$tau, log = TRUE)
    return(lik)    
}

## log likelihood function for optim 
SHASH_logLik1 = function(pars) {
    val = sum(SHASH_logLik(pars))
    if(!is.finite(val)) val = -(10^16)
    return(-val)
}  

out=optim(par=rep(0,12),SHASH_logLik1,control=list(trace=4,maxit=25000)); 
out$estimate=out$par; 
out = maxLik(SHASH_logLik, start = out$par, method="BHHH",control=list(printLevel=2)); 


for(k in 1:5) {
    #out = optim(par=out$estimate,SHASH_logLik1,control=list(trace=4,maxit=25000)); 
    out = maxLik(SHASH_logLik, start = out$estimate, method="BHHH",control=list(printLevel=2)); 
  
}    
fit = make_args(out$estimate); 
matplot(x, cbind(fit$mu,fit$sigma,fit$nu,fit$tau), type="l",lty=1, col=c("black","red","blue","green3")); 
 