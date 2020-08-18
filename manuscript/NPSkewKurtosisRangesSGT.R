################################################################################# 
# Explore how skew and kurtosis depend on nu and tau for gamlss distributions
#################################################################################
rm(list=ls(all=TRUE))
setwd("c:/repos/IPM_size_transitions/manuscript"); #edit as needed 
source("../Diagnostics.R"); 

require(sgt); 

NPshapes = function(lambda,p,q) {
	qv = qsgt(p=c(.1,0.5,0.9),0,1,lambda,p,q)
	sq = (qv[3]+qv[1]-2*qv[2])/(qv[3]-qv[1]);

	qv = qsgt(p=c(0.05,0.25,0.75,0.95),0,1,lambda,p,q)
	qN = qnorm(p=c(0.05,0.25,0.75,0.95))
	u = (qv[4]-qv[1])/(qv[3]-qv[2]);
	uN = (qN[4]-qN[1])/(qN[3]-qN[2]);
	kurto = as.numeric(u/uN-1);
 
    return(list(sq=sq,excess.kurtosis=kurto)) 
}

objfun.min=function(par,target.skew) {
    lambda=exp(par[1]); lambda=-1 + 2*lambda/(1+lambda); 
    p =exp(par[2]); q = exp(par[3]); 
    out=NPshapes(lambda,p,q)
    val = out$excess.kurtosis
    err = abs(out$sq - target.skew)
    val = ifelse(err<0.01, val, (10^7)*err)
    return(val)
}    

objfun.max=function(par,target.skew) {
    lambda=exp(par[1]); lambda=-1 + 2*lambda/(1+lambda); 
    p =exp(par[2]); q = exp(par[3]); 
    out=NPshapes(lambda,p,q)
    val = out$excess.kurtosis
    if(val>1) val = 1 - 0.5*(val-1)/(1 + (val-1)); 
    val = -val; 
    err = abs(out$sq - target.skew)
    val = ifelse(err<0.01, val, (10^7)*err)
    return(val)
}    


minKurtNP=function(target.skew) {
    vals=numeric(3); 
    for(k in 1:3){
    out=optim(par=rnorm(3),fn=objfun.min,target.skew=target.skew,
            control=list(maxit=20000,trace=0)); 
    for(j in 1:2) {
        out=optim(par=out$par,fn=objfun.min,target.skew=target.skew, 
           control=list(maxit=20000,trace=0));             
    }
    vals[k]=out$val; 
    }
    return(min(vals)); 
}

maxKurtNP=function(target.skew) {
    vals=numeric(3); 
    for(k in 1:3) {
    out=optim(par=rnorm(3),fn=objfun.max,target.skew=target.skew,
            control=list(maxit=20000,trace=0)); 
    for(j in 1:2) {
        out=optim(par=out$par,fn=objfun.max,qfun=qfun,target.skew=target.skew, 
           control=list(maxit=20000,trace=0));             
    }
    vals[k]=out$val; 
    }
    val = min(vals); 
    return(-val); 
}

sgtMins = sgtMaxs = matrix(NA,6,1);
target.skews=c(0.01,0.05,seq(0.1,0.7,by=0.2)); 


for(j in 1:6) { 
    sgtMins[j,k] = tryCatch(minKurtNP(qfun,target.skews[j]), error=function(e) NA); 
    sgtMaxs[j,k] = tryCatch(maxKurtNP(qfun,target.skews[j]), error=function(e) NA); 
    cat(j,"\n"); 
}   
