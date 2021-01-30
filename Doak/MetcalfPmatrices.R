require(statmod); # for Gaussian quadrature 

home=ifelse(.Platform$OS.type=="windows", "c:/repos/fanova","~/fanova")
setwd(home); 
setwd("envVar/MetcalfModels"); 

source("MatrixImage.R"); 
source("LightTransitions.R"); 
source("KernelFunctions.R"); 

gpars.CI = make_growthPars("CI"); spars.CI = make_survPars("CI"); 
gpars.DP = make_growthPars("DP"); spars.DP = make_survPars("DP"); 

Q.CI = make_Qmatrix(X,"CI"); # pioneer canopy species 
Q.DP = make_Qmatrix(X,"DP"); # emergent species in the main text 

sizes=c(30,50,100,300,500); L=-2; U=7.5; 
meanLife = varLife = all_meshpts = list(5); 

for(isize in 1:5) {
mz = m = sizes[isize]; h = (U - L)/m; meshpts  = L + ((1:m) - 1/2) * h
Plist.CI = Plist.DP = list(6);
for(q in 1:6) {
    Plist.CI[[q]] = mk_P(L,U,m,gpars.CI,spars.CI,q)$P 
    ## Plist.DP[[q]] = mk_P(L,U,m,gpars.DP,spars.DP,q)$P 
} 

bigmz = 6*m
Marray = matrix(NA, bigmz, bigmz)
for (i in 1:6) {
  for (j in 1:6) {
    Marray[(i-1)*mz + 1:mz, (j-1)*mz + 1:mz] = Plist.CI[[j]]*Q.CI[i,j]
  }
}

N = solve(diag(bigmz)-Marray); 
meanLife[[isize]] = apply(N,2,sum); all_meshpts[[isize]]=meshpts;

P=rbind(Marray,1-apply(Marray,2,sum)); s=nrow(P); P=cbind(P,rep(0,s)); P[s,s]=1; # this is what we call P+ 
fT = matrix(1,1,bigmz); fT=cbind(fT,0); 
vec1s = matrix(1,s,1); 
R1 = R2 = (vec1s%*%fT); 

# Do the calculations for total LRO variance ###########################
tau = s-1; alpha = 1; 
Z = cbind(diag(tau),matrix(0,tau,alpha)); 
tildeRho1 = t(N) %*% Z %*% t(P*R1) %*% vec1s; 
tildeR1 = Z %*% R1 %*% t(Z) 
tildeRho2 = t(N) %*% ( Z %*% t(P*R2) %*% vec1s + 2*t(Marray*tildeR1)%*%tildeRho1) 

rho1 = apply (N, 2, sum)
rbarPib = rho1 %*% Marray
sigsqb=0; pb=1; b=1; 
r2 = (sigsqb + (pb*b)^2 + 2*pb*b*rbarPib) %*% N
varLife[[isize]] = tildeRho2 - tildeRho1^2; 
varLife[[isize]] = r2 - rho1^2; 

# sanity check: these should be the same 
cat(isize, range(tildeRho1 - meanLife[[isize]]), "\n"); 
cat(isize, range(tildeRho2 - t(r2)), "\n"); 


}

## graphing for CI
graphics.off(); 
dev.new(width=8,height=4.5); 
par(mfrow=c(1,2),xaxs="i",bty="l",yaxs="i",cex.axis=1.3,cex.lab=1.3,mgp=c(2.2,1,0),mar=c(4,4,2,2)); 
for(isize in 1:5) {
    m = sizes[isize]
    meshpts=all_meshpts[[isize]]; 
    tau = meanLife[[isize]]
    if(isize==1) {plot(meshpts,tau[3*m+(1:m)],type="l",col=isize,xlim=c(1,7),lwd=2,ylim=c(10,90),
                    xlab="log size (mm dbh)",ylab="Mean remaining lifespan")}
    if(isize>1) points(meshpts,tau[3*m+(1:m)],type="l",col=isize,xlim=c(1,7),ylim=c(10,max(tau)),lwd=2);
}
legend("topleft",legend=as.character(sizes),col=1:5,lwd=2,lty=1,bty="n"); 


for(isize in 1:5) {
    m = sizes[isize]
    meshpts=all_meshpts[[isize]]; 
    tau = varLife[[isize]]^0.5; 
    if(isize==1) plot(meshpts,tau[3*m+(1:m)],type="l",col=isize,xlim=c(1,7),lwd=2,ylim=c(20,80), 
    xlab="log size (mm dbh)",ylab="Std Dev(remaining lifespan)");
    if(isize>1) points(meshpts,tau[3*m+(1:m)],type="l",col=isize,xlim=c(1,7),lwd=2,ylim=c(20,80));
}



