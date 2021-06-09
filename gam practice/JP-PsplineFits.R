rm(list=ls(all=TRUE)); 

setwd("c:/repos/IPM_size_transitions/gam practice"); 
require(minqa); require(fda);
source("JPfuns.R"); 


######### Create covariate for residuals 

z = rt(500,df=10); z=sort(z); hist(z); 

########### Create artificial "residuals" with known sgt parameters 
resids = rSJP(length(z), epsilon=-0.5, delta = 0.5 + 0.05*z^2); 
par(mfrow=c(2,1)); 
hist(resids); plot(z,resids); 

########## Create B-spline basis, and plot it 
B = create.bspline.basis(rangeval=range(z), norder=4, nbasis=6)
x = seq(min(z),max(z),length=200); 
out = eval.basis(x,B);   
matplot(x,out,type="l"); 

########### function to convert coefficients into epsilon and delta values 
X = eval.basis(z,B); P = bsplinepen(B,Lfdobj=2); 

notExp2 = function (x, d, b=1/d) exp(d * sin(x * b))  # from mgcv

make_eps_delta = function(pars){
    epsilon = X%*%pars[1:6]; 
    u = X%*%pars[7:12]; delta = notExp2(u,d=log(10));  # 1/10 to 10 
      
    # evaluate roughness 
    pars=matrix(pars,ncol=1); 
    P.e = t(pars[1:6])%*%P%*%pars[1:6]
    P.d = t(pars[7:12])%*%P%*%pars[7:12]

    return(list(epsilon=epsilon,delta=delta,P.epsilon=P.e,P.delta=P.d)); 
 }
 
 NegLogLik = function(pars){
    out=make_eps_delta(pars);
    NLL = -sum(log(dSJP(z,epsilon=out$epsilon,out$delta)))  
    return(NLL)
}

########################################################################
# Unconstrained fit 
########################################################################
fit = bobyqa(par=rep(0,12), fn=NegLogLik, control = list(iprint=1000,maxfun=250000,rhobeg=1))
fit = bobyqa(par=fit$par, fn=NegLogLik, control = list(iprint=1000,maxfun=250000))

bestSplineFits = make_eps_delta(fit$par); 
par(mfrow=c(2,1)); 
plot(z,bestSplineFits$epsilon,type="l"); points(z,rep(-0.5,length(z)), type="l",lty=2,col="blue") 
rug(z); 

plot(z,bestSplineFits$delta,type="l"); points(z,0.5+0.05*z^2, type="l",lty=2,col="blue"); 
rug(z); 
    
########################################################################
# Penalized fit 
########################################################################
 PenNegLogLik = function(pars,pen) {
     out = make_eps_delta(pars); 
     NLL = NegLogLik(pars); 
     PNLL = NLL + pen[1]*out$P.epsilon + pen[2]*out$P.delta
     return(PNLL)
 } 

Pfit = optim(par=rep(0,12), fn=PenNegLogLik, control = list(trace=4,maxit=250000), pen=10^c(0,0))
Pfit = optim(par=Pfit$par, fn=PenNegLogLik, control = list(trace=4,maxit=250000), pen=10^c(0,0))

Pfit = bobyqa(par=rep(0,12), fn=PenNegLogLik, control = list(iprint=100,maxfun=250000,rhobeg=1), pen=10^c(2,2))
Pfit = bobyqa(par=Pfit$par, fn=PenNegLogLik, control = list(iprint=100,maxfun=250000,rhobeg=1), pen=10^c(2,2))

bestSplineFits = make_eps_delta(Pfit$par); 
par(mfrow=c(2,1)); 
plot(z,bestSplineFits$epsilon,type="l"); points(z,rep(-0.5,length(z)), type="l",lty=2,col="blue") 
rug(z); 

plot(z,bestSplineFits$delta,type="l"); points(z,0.5+0.05*z^2, type="l",lty=2,col="blue"); 
rug(z); 


















 evalAIC = function(logpen){
    pen=10^(logpen); 
       
    fit = bobyqa(par=rep(0,8), fn=PenNegLogLik, control = list(iprint=1000,maxfun=250000,rhobeg=1), pen=pen)
        
    pars=fit$par; 
              
    edf1 = lambda2df(argvals=z,basisobj=B,lambda=pen[1])
    edf2 = lambda2df(argvals=z,basisobj=B,lambda=pen[2])
    AICval = NegLogLik(pars) + log(length(z))*(edf1+edf2);
    write.table(round(pars,digits=8), append=FALSE, row.names=FALSE,col.names=FALSE,file="splinePars.txt")
        
    return(AICval); 
 }      

evalAIC = function(theta) {
    pars=theta[1:8]; pen=10^theta[9:10]; 
 

Val = matrix(NA,5,5); 
lambdas = seq(-6,6,length=5); 
for(i in 1:5){
for(j in 1:5){
    x = c(lambdas[i],lambdas[j])
    Val[i,j]=evalAIC(x)
    cat(i,j,Val[i,j],"\n"); 
}}

bestPlace = which(Val==min(Val),arr.ind=TRUE);
bestVal = lambdas[bestPlace]; 

# optimize the (highly) approximate AIC as a function of penalty 
bestPen = optim(par=bestVal,fn=evalAIC, control=list(maxit=10000,trace=4)); 

# refit with the best penalty 
start = read.table("splinePars.txt")$V1; 
bestFit = bobyqa(par=start, fn=PenNegLogLik, control = list(iprint=100,maxfun=250000),pen=10^bestPen$par)

bestSplineFits = make_lpq(bestSplinePars$V1); 
par(mfrow=c(3,1)); 
plot(z,bestSplineFits$lambda);
plot(z,bestSplineFits$p);
plot(z,bestSplineFits$q);






 
 