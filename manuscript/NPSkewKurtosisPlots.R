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

doPlotNP=function(qfun,nu,tau,ylim=c(-0.5,4)) { 
SkewMat = KurtMat = matrix(NA,length(nu),length(tau));
for(i in 1:length(nu)){
cat(i,"\n"); 
for(j in 1:length(tau)){
    out=NPshapes(nu[i],tau[j],qfun);
    SkewMat[i,j]=out$sq;
    KurtMat[i,j]=out$excess.kurtosis
}}    
    e = which(SkewMat>=0); 
    plot(SkewMat[e],KurtMat[e],xlab="NP Skewness",ylab="NP Excess Kurtosis",ylim=ylim); 
}

graphics.off(); 
par(mfrow=c(3,4),bty="l",cex.axis=1.3,cex.lab=1.3,mar=c(4,4,2,1),mgp=c(2.2,1,0),yaxs="i"); 

nu=exp(seq(log(0.3),log(35),length=31))
tau=exp(seq(log(0.3),log(35),length=32))

## EGB2 
#nu=seq(.5,20,length=31); tau=seq(.5,20,length=32); 
doPlotNP(qEGB2,nu,tau,ylim=c(0,0.5));  title(main="EGB2"); 

## JSU 
#nu=seq(.5,20,length=31); tau=seq(.5,20,length=32); 
doPlotNP(qJSU,nu,tau);  title(main="JSU"); 

### SHASH 
#nu=seq(.5,20,length=31); tau=seq(.5,20,length=32); 
doPlotNP(qSHASH,nu,tau);  title(main="SHASH"); 

### SHASHo 
#nu=seq(.5,20,length=31); tau=seq(.5,20,length=32); 
doPlotNP(qSHASHo,nu,tau);   title(main="SHASHo"); 

### SEP1 
#nu=seq(.5,20,length=31); tau=seq(.5,20,length=32); 
doPlotNP(qSEP1,nu,tau);  title(main="SEP1"); 

## SEP2 
#nu=seq(.5,20,length=31); tau=seq(.5,20,length=32); 
doPlotNP(qSEP2,nu,tau);  title(main="SEP2"); 

## SEP3
doPlotNP(qSEP3,nu,tau);  title(main="SEP3"); 

## SEP4 
#nu=seq(.5,20,length=31); tau=seq(.5,20,length=32); 
doPlotNP(qSEP4,nu,tau);  title(main="SEP4"); 

## ST1 
#nu=seq(.5,20,length=31); tau=seq(.5,20,length=32); 
doPlotNP(qST1,nu,tau);  title(main="ST1"); 

## ST2 
#nu=seq(.5,20,length=31); tau=seq(.5,20,length=32); 
doPlotNP(qST2,nu,tau);  title(main="ST2"); 

## ST3 
#nu=seq(.5,20,length=31); tau=seq(.5,20,length=32); 
doPlotNP(qST3C,nu,tau);  title(main="ST3"); 

## ST4 
#nu=seq(.5,20,length=31); tau=seq(.5,20,length=32); 
nu=exp(seq(log(.1),log(1.2),length=31))
tau=exp(seq(log(.1),log(25),length=32))
doPlotNP(qST4,nu,tau);  title(main="ST4");   

## ST5 
nu=exp(seq(log(.1),log(1.2),length=31))
tau=exp(seq(log(.1),log(25),length=32))
doPlotNP(qST5,nu,tau,ylim=c(0,5));  title(main="ST5"); 



