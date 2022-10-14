require(gamlss); require(moments); 

epsilons = seq(-0.5,0.5,length=20); 

zbar = zsd = zskew = zkurt = numeric(length(epsilons)); 
for(k in 1:20){
    z = rSHASHo(500000, mu = 0, sigma = 1, nu = epsilons[k], tau = 1)
    if(k==1) z1 = z; 
    if(k==20) z20 = z; 
    zbar[k] = mean(z); 
    zsd[k] = sd(z);
    zskew[k] = skewness(z);
    zkurt[k] = kurtosis(z); 
    cat(k,"\n"); 
} 

par(mfrow=c(3,2),mar=c(4,4,1,1),bty="l",mgp=c(2.1,1,0), cex.axis=1.3, cex.lab=1.3); 
hist(z1, main="Largest left skew"); hist(z20, main="Largest right skew"); 
plot(epsilons,zskew,xlab="Skewness parameter nu", ylab="Skewness",type="o"); abline(0,1,col="blue", lty=2); title(main="Skewness"); 
plot(zskew, zbar, xlab="Skewness",ylab="Mean",type="o"); abline(0,1,col="blue", lty=2); title(main="Mean"); 
plot(zskew, zsd, xlab="Skewness",ylab="Std Dev",type="o"); title(main="Standard Deviation") 
plot(zskew, zkurt, xlab="Skewness",ylab="Kurtosis",type="o"); title(main="Kurtosis"); 


epsilons = seq(-0.5,0.5,length=20); 
zbar = zsd = zskew = zkurt = numeric(length(epsilons)); 
for(k in 1:20){
    z = rRSJP(500000, lambda = epsilons[k], tau = .2)
    if(k==1) z1 = z; 
    if(k==20) z20 = z; 
    zbar[k] = mean(z); 
    zsd[k] = sd(z);
    zskew[k] = skewness(z);
    zkurt[k] = kurtosis(z); 
    cat(k,"\n"); 
} 

par(mfrow=c(3,2),mar=c(4,4,1,1),bty="l",mgp=c(2.1,1,0), cex.axis=1.3, cex.lab=1.3); 
hist(z1, main="Largest left skew"); hist(z20, main="Largest right skew"); 
plot(epsilons,zskew,xlab="Skewness parameter nu", ylab="Skewness",type="o"); abline(0,1,col="blue", lty=2); title(main="Skewness"); 
plot(zskew, zbar, xlab="Skewness",ylab="Mean",type="o"); title(main="Mean"); 
plot(zskew, zsd, xlab="Skewness",ylab="Std Dev",type="o"); title(main="Standard Deviation") 
plot(zskew, zkurt, xlab="Skewness",ylab="Kurtosis",type="o"); title(main="Kurtosis"); 

