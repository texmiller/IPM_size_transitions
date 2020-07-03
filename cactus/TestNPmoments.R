kurt = NPkurt = sq = NPsq = numeric(2500); 
for(k in 1:2500) {
	z=rt(200,df=10); 
	kurt[k]=kurtosis(z);
	NPkurt[k]=NPkurtosis(z);
	sq[k]=skewness(z); 
	NPsq[k]=NPskewness(z);
}
par(mfrow=c(2,1));
hist(kurt); hist(NPkurt); 

dev.new(); 
par(mfrow=c(2,1));
hist(sq); hist(NPsq); 