### demo of sampling properties for NP skew and kurtosis
rm(list=ls(all=TRUE))
setwd("c:/repos/IPM_size_transitions"); #edit as needed 

require(moments); require(car); 

NPskewness=function(x) (mean(x)-median(x))/sd(x) 

NPkurtosis=function(x) {
	q = quantile(x,c(0.05,0.25,0.75,0.95))
	qN = qnorm(c(0.05,0.25,0.75,0.95))
	u = (q[4]-q[1])/(q[3]-q[2]);
	uN = (qN[4]-qN[1])/(qN[3]-qN[2]);
	return (as.numeric(u - uN)) 
}

m2 = integrate(function(x) (x^2)*dt(x,df=10), -30, 30)$value;
m4 = integrate(function(x) (x^4)*dt(x,df=10), -30, 30)$value;
trueK = m4/(m2^2); 

q = qt(c(0.05,0.25,0.75,0.95),df=10)
qN = qnorm(c(0.05,0.25,0.75,0.95))
u = (q[4]-q[1])/(q[3]-q[2]);
uN = (qN[4]-qN[1])/(qN[3]-qN[2]);
trueNPK = u - uN; 

n = 200; 
sampK = sampNPK = sampSkew = sampNPSkew = numeric(5000); 
for(j in 1:5000) {
	z = rt(n,df=10); 
	sampK[j]=kurtosis(z) - 3;
	sampNPK[j]=NPkurtosis(z);
	sampSkew[j]=skewness(z);
	sampNPSkew[j]=NPskewness(z);
}

qB = function(z) quantile(z,c(0.05,0.95)); 

graphics.off(); 
dev.new(height=5,width=8); 
par(mfcol=c(2,2),mar=c(4,4,2,2),mgp=c(2.2,1,0),cex.axis=1.3,cex.lab=1.3);
hist(sampSkew,20,xlab="Skewness",main="Skewness");
abline(v=qB(sampSkew),col="blue",lty=2,lwd=2); 
abline(v=range(sampSkew),col="red",lty=3,lwd=2); 
points(0,0,col="black",pch=16,cex=1.5); 

hist(sampNPSkew,20,xlab="NP Skewness",main="NP Skewness");
abline(v=qB(sampNPSkew),col="blue",lty=2,lwd=2); 
abline(v=range(sampNPSkew),col="red",lty=3,lwd=2); 
points(0,0,col="black",pch=16,cex=1.5);

hist(sampK,20,xlab="Excess Kurtosis",main="Excess Kurtosis"); 
abline(v=qB(sampK),col="blue",lty=2,lwd=2); 
abline(v=range(sampK),col="red",lty=3,lwd=2); 
points(trueK-3,0,col="black",pch=16,cex=1.5);

hist(sampNPK,20,xlab="Excess NP Kurtosis",main="Excess NP Kurtosis"); 
abline(v=qB(sampNPK),col="blue",lty=2,lwd=2); 
abline(v=range(sampNPK),col="red",lty=3,lwd=2); 
points(trueNPK,0,col="black",pch=16,cex=1.5);

dev.copy2pdf(file="manuscript/figures/NPmoments.pdf"); 

kurtosis(sampSkew)-3; kurtosis(sampNPSkew)-3;
kurtosis(sampK)-3; kurtosis(sampNPK)-3; 

dev.new();  
par(mfcol=c(2,2));
qqPlot(sampSkew,main="Skewness");
qqPlot(sampNPSkew,main="NP Skewness");
qqPlot(sampK,main="Kurtosis"); 
qqPlot(sampNPK,main="NP Kurtosis"); 
	