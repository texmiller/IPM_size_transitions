### demo of sampling properties for NP skew and kurtosis
rm(list=ls(all=TRUE))
setwd("c:/repos/IPM_size_transitions/manuscript"); #edit as needed 

require(moments); require(car); source("../Diagnostics.R"); 

m2 = integrate(function(x) (x^2)*dt(x,df=8), -Inf, Inf)$value;
m4 = integrate(function(x) (x^4)*dt(x,df=8), -Inf, Inf)$value;
trueK = m4/(m2^2); 

q = qt(c(0.1,0.25,0.75,0.9),df=8)
qN = qnorm(c(0.1,0.25,0.75,0.9))
u = (q[4]-q[1])/(q[3]-q[2]);
uN = (qN[4]-qN[1])/(qN[3]-qN[2]);
trueNPK = u - uN; 

n = 200; 
sampK = sampNPK = sampSkew = sampNPSkew = numeric(5000); 
for(j in 1:5000) {
	z = rt(n,df=8); 
	sampK[j]=kurtosis(z) - 3;
	sampNPK[j]=NPkurtosis(z);
	sampSkew[j]=skewness(z);
	sampNPSkew[j]=NPskewness(z);
}

qB = function(z) quantile(z,c(0.05,0.95)); 
qB2 = function(z) range(z); 

graphics.off(); 
dev.new(height=5,width=8); 
par(mfcol=c(2,2),mar=c(4,4,2,2),mgp=c(2.2,1,0),cex.axis=1.3,cex.lab=1.3);
hist(sampSkew,50,xlab="Skewness",main="Skewness");
abline(v=qB(sampSkew),col="blue",lty=2,lwd=2); 
abline(v=qB2(sampSkew),col="red",lty=3,lwd=3); 
points(0,0,col="black",pch=16,cex=1.75); 

hist(sampNPSkew,50,xlab="NP Skewness",main="NP Skewness");
abline(v=qB(sampNPSkew),col="blue",lty=2,lwd=2); 
abline(v=qB2(sampNPSkew),col="red",lty=3,lwd=3); 
points(0,0,col="black",pch=16,cex=1.75);

hist(sampK,50,xlab="Excess Kurtosis",main="Excess Kurtosis"); 
abline(v=qB(sampK),col="blue",lty=2,lwd=2); 
abline(v=qB2(sampK),col="red",lty=3,lwd=3); 
points(trueK-3,0,col="black",pch=16,cex=1.75);

hist(sampNPK,50,xlab="NP Excess Kurtosis",main="Excess NP Kurtosis"); 
abline(v=qB(sampNPK),col="blue",lty=2,lwd=2); 
abline(v=qB2(sampNPK),col="red",lty=3,lwd=3); 
points(trueNPK,0,col="black",pch=16,cex=1.75);

dev.copy2pdf(file="figures/NPmoments.pdf"); 

dev.new();  
par(mfcol=c(2,2));
qqPlot(sampSkew,main="Skewness");
qqPlot(sampNPSkew,main="NP Skewness");
qqPlot(sampK,main="Kurtosis"); 
qqPlot(sampNPK,main="NP Kurtosis"); 
	