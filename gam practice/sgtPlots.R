rm(list=ls(all=TRUE)); 

setwd("c:/repos/IPM_size_transitions/gam practice"); 
require(sgt); require(minqa); 

######### Visualize effects of changing p and q 
dev.new(); 
par(mfrow=c(3,1),mar=c(4,4,1,1),mgp=c(2.2,1,0),bty="l"); 
x=seq(-5,5,length=500); 

py1 = dsgt(x,mu=0,sigma=1,lambda=0.5,p=1,q=2.2);
py2 = dsgt(x,mu=0,sigma=1,lambda=0.5,p=1,q=4);
py3 = dsgt(x,mu=0,sigma=1,lambda=0.5,p=1,q=20);
matplot(x,cbind(py1,py2,py3),type="l",col=c("black","red","blue")); 
title(main="p=1"); 
legend("left",legend=c("q=2.2", "q=4", "q=20"), col=c("black","red","blue"), lty=1, pch=NULL,cex=2,bty="n",inset=0.04); 

py1 = dsgt(x,mu=0,sigma=1,lambda=0.5,p=2,q=1.2);
py2 = dsgt(x,mu=0,sigma=1,lambda=0.5,p=2,q=4);
py3 = dsgt(x,mu=0,sigma=1,lambda=0.5,p=2,q=20);
matplot(x,cbind(py1,py2,py3),type="l",col=c("black","red","blue")); 
title(main="p=2"); 
legend("left",legend=c("q=1.2", "q=4", "q=20"), col=c("black","red","blue"), lty=1, pch=NULL,cex=2,bty="n",inset=0.04); 

x=seq(-5,5,length=500); 
py1 = dsgt(x,mu=0,sigma=1,lambda=0.5,p=5,q=.5);
py2 = dsgt(x,mu=0,sigma=1,lambda=0.5,p=5,q=4);
py3 = dsgt(x,mu=0,sigma=1,lambda=0.5,p=5,q=20);
matplot(x,cbind(py1,py2,py3),type="l",col=c("black","red","blue")); 
title(main="p=5"); 
legend("left",legend=c("q=0.5", "q=4", "q=20"), col=c("black","red","blue"), lty=1, pch=NULL,cex=2,bty="n",inset=0.04); 


sgt_NPskewness = function(lambda,p,q,tail.p=0.1) {
	qt = qsgt(c(tail.p,0.5,1-tail.p),lambda=lambda,p=p,q = q)
	u = (qt[3]+qt[1]-2*qt[2])/(qt[3]-qt[1]);
	return(as.numeric(u)); 
	
}	

sgt_NPkurtosis=function(lambda,p,q, tail.p=0.05) {
	qt = qsgt(c(tail.p,0.25,0.75,1-tail.p),lambda=lambda,p=p,q=q)
	qN = qnorm(c(tail.p,0.25,0.75,1-tail.p))
	u = (qt[4]-qt[1])/(qt[3]-qt[2]);
	uN = (qN[4]-qN[1])/(qN[3]-qN[2]);
	return (as.numeric(u/uN-1)) 
}

############# Vary lambda and p 
SkewMat = KurtMat = matrix(NA,50,51);
lambda = seq(-0.95,0.95,length=50)
p = seq(1.5,10,length=51);  
for(i in 1:50) {
for(j in 1:51){
    SkewMat[i,j]=sgt_NPskewness(lambda=lambda[i],p=p[j],q=2);
    KurtMat[i,j]=sgt_NPkurtosis(lambda=lambda[i],p=p[j],q=2);
}}

require(viridisLite); require(fields); 
graphics.off(); dev.new(); 
image.plot(lambda,p,SkewMat,col=plasma(64));  title(main="sgt NP Skewness"); 
contour(lambda,p,SkewMat,add=TRUE,labcex=1) 


dev.new(); 
image.plot(lambda,p,KurtMat,col=plasma(64));  title(main="sgt NP Kurtosis"); 
contour(lambda,p,KurtMat,add=TRUE,labcex=1) 

############ Vary lambda and q 
SkewMat = KurtMat = matrix(NA,50,51);
lambda = seq(-0.95,0.95,length=50)
q = seq(1.5,10,length=51);  
for(i in 1:50) {
for(j in 1:51){
    SkewMat[i,j]=sgt_NPskewness(lambda=lambda[i],p=2,q=q[j]);
    KurtMat[i,j]=sgt_NPkurtosis(lambda=lambda[i],p=2,q=q[j]);
}}

dev.new(); 
image.plot(lambda,q,SkewMat,col=plasma(64));  title(main="sgt NP Skewness"); 
contour(lambda,p,SkewMat,add=TRUE,labcex=1) 

dev.new(); 
image.plot(lambda,q,KurtMat,col=plasma(64));  title(main="sgt NP Kurtosis"); 
contour(lambda,p,KurtMat,add=TRUE,labcex=1) 


############ Vary p and q 
SkewMat = KurtMat = matrix(NA,60,61);
p = seq(1.5,10,length=60);  
q = seq(1.5,10,length=61);  
for(i in 1:60) {
for(j in 1:61){
    SkewMat[i,j]=sgt_NPskewness(lambda=0.35, p=p[i],q=q[j]);
    KurtMat[i,j]=sgt_NPkurtosis(lambda=0.35, p=p[i],q=q[j]);
}}
graphics.off(); 

dev.new(); 
image.plot(p,q,SkewMat);  title(main="sgt NP Skewness"); 
contour(p,q,SkewMat,add=TRUE,labcex=1) 

dev.new(); 
image.plot(p,q,KurtMat);  title(main="sgt NP Kurtosis"); 
contour(p,q,KurtMat,add=TRUE,labcex=1) 




#######################################################
#  Display how NP skew and kurtosis depend on parameters
#######################################################

SkewMat = KurtMat = matrix(NA,150,151);
nu = seq(-3,3,length=150)
tau = seq(-1,1,length=151);  
for(i in 1:150) {
for(j in 1:151){
    delta=exp(-0.5*tau[j]); epsilon=delta*nu[i]; 
    SkewMat[i,j]=JP_NPskewness(epsilon,delta);
    KurtMat[i,j]=JP_NPkurtosis(epsilon,delta);
}}

graphics.off(); 
dev.new(); 
image.plot(nu,tau,SkewMat,col=plasma(64));  title(main="NP Skewness"); 
contour(nu,tau,SkewMat,add=TRUE) 

dev.new(); 
image.plot(nu,tau,KurtMat,col=plasma(64)); title(main = "NP Kurtosis"); 
contour(nu,tau,KurtMat,add=TRUE) 












######### Create covariate z for artificial "residuals" 
require(fda); require(sgt); require(maxLik); 
z = rt(500,df=10); z=sort(z); hist(z); 

########### Create artificial "residuals" with sgt parameters p(z) and q(z), mu=0, sigma=1 
True.lambda= -0.5 + plogis(z); True.p = 3 + 0.5*z^2; True.q = 3 - 0.5*z

graphics.off(); dev.new(height=9,width=7); par(mfrow=c(3,1)); 
matplot(z,cbind(True.lambda, True.p, True.q), type="l", col=1:3); range(True.p*True.q); 
resids = rsgt(length(z),mu=0,sigma=1,lambda= True.lambda, p = True.p, q= True.q); 
hist(resids); plot(z,resids); 

########## Create B-spline basis for fitting p and q, and plot it 
B = create.bspline.basis(rangeval=range(z), norder=4, nbasis=4)
x = seq(min(z),max(z),length=200); 
out = eval.basis(x,B);   
matplot(x,out,type="l"); 

########### function to convert spline coefficients into lambda, p, q values 
X = eval.basis(z,B); P = bsplinepen(B); 

notExp2 = function (x, d, b=1/d) exp(d * sin(x * b))  # from mgcv

make_lpq = function(pars){
    u = X%*%pars[1:4]; lambda = -1 + 2*exp(u)/(1+exp(u)); # -1 to 1 (open interval)
    u = X%*%pars[5:8]; pz = notExp2(u,d=log(50));  # 1/50 to 50 
    u = X%*%pars[9:12]; qz = 0.5 + notExp2(u,d=log(100)); #1/100 to 100 
    
    # evaluate penalties: integrated squared 2nd derivative  
    pars=matrix(pars,ncol=1); 
    Plam = t(pars[1:4])%*%P%*%pars[1:4]
    Pp = t(pars[5:8])%*%P%*%pars[5:8]
    Pq = t(pars[9:12])%*%P%*%pars[9:12]
    

    # p*q must be bigger than 2 to have finite variance 
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
 
 
PenNegLogLik1 = function(pars,pen) {
     out = make_lpq(pars); 
     NLL = NegLogLik(pars); 
     PNLL = NLL + pen*(out$Plam + out$Pp + out$Pq); 
     return(PNLL)
 } 


fits=list(10); pens=10^seq(-3,0,length=10)
for(k in 1:10) {
    foo = optim(par=rep(0,12), PenNegLogLik1, control = list(maxit=10000, trace=0), pen=pens[k]) 
    for(rr in 1:10) {
    foo = bobyqa(par=foo$par, fn=PenNegLogLik1, control = list(iprint=100,maxfun=10000,rhobeg=1), pen=pens[k]) 
    foo = optim(foo$par, PenNegLogLik1, control = list(maxit=10000, trace=0), pen=pens[k]) 
    cat(k, rr, "<---------------------------------------->", "\n"); 
    }
    fits[[k]]= bobyqa(par=foo$par, fn=PenNegLogLik1, control = list(iprint=100,maxfun=10000,rhobeg=1), pen=pens[k]) 
    cat("Finished fit ", k, "<------------------------------->", "\n"); 
}

par(mfrow=c(3,1),mar=c(4,4,1,1),bty="l"); 
out = make_lpq(fits[[1]]$par); ps = out$p; lambs=out$lambda; qs=out$q; 
for(k in 2:10) {
    out = make_lpq(fits[[k]]$par); 
    ps = cbind(ps,out$p); lambs=cbind(lambs,out$lambda); qs=cbind(qs,out$q); 
}
require(viridisLite); 
matplot(z,lambs,lty=1,type="l",col=magma(10)); points(z,True.lambda,col="red",lty=1,lwd=2,type="l"); 
matplot(z,ps,lty=1,type="l",col=magma(10));  points(z,True.p,col="red",lty=1,lwd=2,type="l"); 
matplot(z,qs,lty=1,type="l",col=magma(10));   points(z,True.q,col="red",lty=1,lwd=2,type="l"); 


##################################################
# Basement: 3-penalty fitting attempt. 
##################################################

if(FALSE){
## This isn't really AIC because it's based on the 2nd derivative complexity measure
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

Val = array(NA,c(4,4,4)); 
lambdas = seq(-6,6,length=4); 
for(i in 1:4){
for(j in 1:4){
for(k in 1:4){
    x = c(lambdas[i],lambdas[j],lambdas[k])
    Val[i,j,k]=evalAIC(x)
    cat(i,j,k,Val[i,j,k],"\n"); 
}}}

 
bestPlace = which(Val==min(Val),arr.ind=TRUE);
bestVal = lambdas[bestPlace]; 

bestPen = optim(par=bestVal,fn=evalAIC, method="BFGS", control=list(maxit=10000,trace=4)); 

bestAIC = evalAIC(bestPen$par); # last thing done 

bestSplinePars=read.table("splinePars.txt"); 

bestSplineFits = make_lpq(bestSplinePars$V1); 
par(mfrow=c(3,1)); 
plot(z,bestSplineFits$lambda);
plot(z,bestSplineFits$p);
plot(z,bestSplineFits$q);
} ########################################################## End of the basement 







 
 