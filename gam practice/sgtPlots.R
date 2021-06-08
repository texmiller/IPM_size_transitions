rm(list=ls(all=TRUE)); 

setwd("c:/repos/IPM_size_transitions/gam practice"); 
require(sgt); require(minqa); 

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
B = create.bspline.basis(rangeval=range(z), norder=4, nbasis=4)
x = seq(min(z),max(z),length=200); 
out = eval.basis(x,B);   
matplot(x,out,type="l"); 


########### function to convert coefficients into lambda, p, q values 
X = eval.basis(z,B); P = bsplinepen(B); 

notExp2 = function (x, d, b=1/d) exp(d * sin(x * b))  # from mgcv

make_lpq = function(pars){
    u = X%*%pars[1:4]; lambda = -1 + 2*exp(u)/(1+exp(u)); # -1 to 1 (open interval)
    u = X%*%pars[5:8]; pz = notExp2(u,d=log(50));  # 1/50 to 50 
    u = X%*%pars[9:12]; qz = 0.5 + notExp2(u,d=log(100)); #1/100 to 100 
    
    # evaluate penalties 
    pars=matrix(pars,ncol=1); 
    Plam = t(pars[1:4])%*%P%*%pars[1:4]
    Pp = t(pars[5:8])%*%P%*%pars[5:8]
    Pq = t(pars[9:12])%*%P%*%pars[9:12]
    

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

########## Unconstrained fit 
fit = bobyqa(par=rep(0,12), fn=NegLogLik, control = list(iprint=0,maxfun=5000,rhobeg=1))
    
 
 PenNegLogLik = function(pars,pen) {
     out = make_lpq(pars); 
     NLL = NegLogLik(pars); 
     PNLL = NLL + pen[1]*out$Plam + pen[2]*out$Pp + pen[3]*out$Pq
     return(PNLL)
 } 

 evalAIC = function(logpen){
        pen=10^(logpen); 
        
        #fit = optim(par=rep(0,12),fn=PenNegLogLik,control=list(maxit=10000,trace=0),pen=pen)

        fit = bobyqa(par=rep(0,12), fn=PenNegLogLik, control = list(iprint=0,maxfun=5000,rhobeg=1), pen=pen)
        
        pars=fit$par; 
               
        edf1 = lambda2df(argvals=z,basisobj=B,lambda=pen[1])
        edf2 = lambda2df(argvals=z,basisobj=B,lambda=pen[2])
        edf3 = lambda2df(argvals=z,basisobj=B,lambda=pen[3])
        AICval = NegLogLik(pars) + log(length(z))*(edf1+edf2+edf3);
        write.table(round(pars,digits=8), append=FALSE, row.names=FALSE,col.names=FALSE,file="splinePars.txt")
        
        return(AICval); 
 }       

Val = array(NA,c(6,6,6)); 
lambdas = seq(-6,6,length=6); 
for(i in 1:6){
for(j in 1:6){
for(k in 1:6){
    x = c(lambdas[i],lambdas[j],lambdas[k])
    Val[i,j,k]=evalAIC(x)
    cat(i,j,k,Val[i,j,k],"\n"); 
}}}

 
bestPlace = which(Val==min(Val),arr.ind=TRUE);
bestVal = lambdas[bestPlace]; 

bestPen = optim(par=bestVal,fn=evalAIC, control=list(maxit=10000,trace=4)); 

bestAIC = evalAIC(bestPen$par); # last thing done 

bestSplinePars=read.table("splinePars.txt"); 

bestSplineFits = make_lpq(bestSplinePars$V1); 
par(mfrow=c(3,1)); 
plot(z,bestSplineFits$lambda);
plot(z,bestSplineFits$p);
plot(z,bestSplineFits$q);








 
 