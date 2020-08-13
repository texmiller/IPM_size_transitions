################################################################################# 
# Explore how skew and kurtosis depend on nu and tau for gamlss distributions
#################################################################################
rm(list=ls(all=TRUE))
setwd("c:/repos/IPM_size_transitions/manuscript"); #edit as needed 
source("../Diagnostics.R"); 

require(gamlss); 

NPshapes = function(nu,tau,qfun) {
	q = qfun(p=c(.1,0.5,0.9),0,1,nu,tau)
	sq = (q[3]+q[1]-2*q[2])/(q[3]-q[1]);

	q = qfun(p=c(0.05,0.25,0.75,0.95),0,1,nu,tau)
	qN = qnorm(p=c(0.05,0.25,0.75,0.95))
	u = (q[4]-q[1])/(q[3]-q[2]);
	uN = (qN[4]-qN[1])/(qN[3]-qN[2]);
	kurto = as.numeric(u/uN-1);
 
    return(list(sq=sq,excess.kurtosis=kurto)) 
}

objfun.min=function(p,qfun,target.skew) {
    nu=exp(p[1]); tau=exp(p[2])
    out=NPshapes(nu,tau,qfun)
    val = out$excess.kurtosis
    err = abs(out$sq - target.skew)
    val = ifelse(err<0.01, val, (10^7)*err)
    return(val)
}    

objfun.max=function(p,qfun,target.skew) {
    nu=exp(p[1]); tau=exp(p[2])
    out=NPshapes(nu,tau,qfun)
    val = out$excess.kurtosis
    if(val>1) val = 1 - 0.5*(val-1)/(1 + (val-1)); 
    val = -val; 
    err = abs(out$sq - target.skew)
    val = ifelse(err<0.01, val, (10^7)*err)
    return(val)
}    

minKurtNP=function(qfun,target.skew) {
    out=optim(par=runif(2),fn=objfun.min,qfun=qfun,target.skew=target.skew,
            control=list(maxit=20000,trace=0)); 
    for(j in 1:3) {
        out=optim(par=out$par,fn=objfun.min,qfun=qfun,target.skew=target.skew, 
           control=list(maxit=20000,trace=0));             
    }
    return(out$val); 
}

maxKurtNP=function(qfun,target.skew) { 
    out=optim(par=runif(2),fn=objfun.max,qfun=qfun,target.skew=target.skew,
            control=list(maxit=20000,trace=0)); 
    for(j in 1:3) {
        out=optim(par=out$par,fn=objfun.max,qfun=qfun,target.skew=target.skew, 
           control=list(maxit=20000,trace=0));             
    }
    return(-out$val); 
}

funs=list(12);
funs[[1]]=qJSU; funs[[2]]=qSHASH; funs[[3]]=qSHASHo2;
funs[[4]]=qSEP1; funs[[5]]=qSEP2; funs[[6]]=qSEP3; funs[[7]]=qSEP4; 
funs[[8]]=qST1; funs[[9]]=qST2; funs[[10]]=qST3; funs[[11]]=qST4; funs[[12]]=qST5; 

theMins = theMaxs = matrix(NA,6,12);
target.skews=c(0.01,0.05,seq(0.1,0.7,by=0.2)); 

for(k in 1:12) {
  qfun = funs[[k]]
  for(j in 1:6) { 
    theMins[j,k] = tryCatch(minKurtNP(qfun,target.skews[j]), error=function(e) NA); 
    theMaxs[j,k] = tryCatch(maxKurtNP(qfun,target.skews[j]), error=function(e) NA); 
    cat(j,k,"\n"); 
  }
}    


graphics.off(); 
dev.new(width=9,height=7); 
par(mfrow=c(2,3),bty="l",mgp=c(2.2,1,0),mar=c(4,4,1,1),cex.axis=1.3,cex.lab=1.3); 

ylim1=c(min(theMins,na.rm=TRUE),1.5); 

matplot(target.skews,theMins[,1:3],ylim=ylim1, type="o",lty=1,col=1:4,lwd=2,pch=1:4,xlab="NP Skew",
ylab = "Minimum NP Kurtosis",cex=1.3);
legend("topleft",legend=c("JSU","SHASH","SHASHo"),lty=1,col=1:4,lwd=2,pch=1:4,bty="n");

matplot(target.skews,theMins[,4:7],ylim=ylim1, type="o",lty=1,col=1:4,lwd=2,pch=1:4,xlab="NP Skew",
ylab = "Minimum NP Kurtosis",cex=1.3);
legend("topleft",legend=c("SEP1","SEP2","SEP3","SEP4"),lty=1,col=1:4,lwd=2,pch=1:4,bty="n");

matplot(target.skews,theMins[,8:12],ylim=ylim1, type="o",lty=1,col=1:5,lwd=2,pch=1:5,xlab="NP Skew",
ylab = "Minimum NP Kurtosis",cex=1.3);
legend("topleft",legend=c("ST1","ST2","ST3","ST4","ST5"),lty=1,col=1:5,lwd=2,pch=1:5,bty="n");


ylim1=c(min(theMaxs,na.rm=TRUE),1); 

matplot(target.skews,theMaxs[,1:3],ylim=ylim1, type="o",lty=1,col=1:4,lwd=2,pch=1:4,xlab="NP Skew",
ylab = "Maximum NP Kurtosis",cex=1.3);

matplot(target.skews,theMaxs[,4:7],ylim=ylim1, type="o",lty=1,col=1:4,lwd=2,pch=1:4,xlab="NP Skew",
ylab = "Maximum NP Kurtosis",cex=1.3);

matplot(target.skews,theMaxs[,8:12],ylim=ylim1, type="o",lty=1,col=1:5,lwd=2,pch=1:5,xlab="NP Skew",
ylab = "Maximum NP Kurtosis",cex=1.3);



