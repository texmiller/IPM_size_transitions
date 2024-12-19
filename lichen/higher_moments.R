 ### Not ready for prime time, uses scripts from metaLuck and testLuck 
 ### assumes that the Gaussian and JSU IPM's have been made already 
 
 Umat = IPM_G$P; Fmat = IPM_G$F; 
 lichen_c0 = rep(0,nrow(matU)); lichen_c0[1]=1
 
 a = mean_lifespan(Umat, mixdist=lichen_c0);
 b = var_lifespan(Umat, mixdist=lichen_c0)^0.5;
 c = skew_lifespan(Umat,mixdist=lichen_c0);
 d = mean_LRO(Umat,Fmat,mixdist=lichen_c0);
 e = var_LRO_mcr(Umat,Fmat,mixdist=lichen_c0)^0.5;
 f = skew_LRO(Umat,Fmat,mixdist=lichen_c0);

Plist = list(150); for(j in 1:150) Plist[[j]]=Umat; 
meanReward=matrix(Fmat[1,],150,ncol(Fmat),byrow=TRUE); 
birthpops=numeric(150); birthpops[1]=250000; 
sims_G = manyLivesAllYears(Plist, meanReward, birthpops)
 

 Umat = IPM_J$P; Fmat = IPM_J$F; 
 lichen_c0 = rep(0,nrow(matU)); lichen_c0[1]=1
 
 a = mean_lifespan(Umat, mixdist=lichen_c0);
 b = var_lifespan(Umat, mixdist=lichen_c0)^0.5;
 c = skew_lifespan(Umat,mixdist=lichen_c0);
 d = mean_LRO(Umat,Fmat,mixdist=lichen_c0);
 e = var_LRO_mcr(Umat,Fmat,mixdist=lichen_c0)^0.5;
 f = skew_LRO(Umat,Fmat,mixdist=lichen_c0);
 
Plist = list(150); for(j in 1:150) Plist[[j]]=Umat; 
meanReward=matrix(Fmat[1,],150,ncol(Fmat),byrow=TRUE); 
birthpops=numeric(150); birthpops[1]=250000; 
sims_S = manyLivesAllYears(Plist, meanReward, birthpops)
 
par(mfrow=c(2,1));
R_G = sims_G$R; R_S = sims_S$R; xmax = (max(c(R_G,R_S))); 
hist((R_G[R_G>0]),breaks=seq(0,129,by=5),xlim=c(0,xmax),ylim=c(0,0.16), xlab="Offspring", ylab="Frequency",main="Gaussian model",probability=TRUE); 
hist((R_S[R_S>0]),breaks=seq(0,129,by=5),xlim=c(0,xmax),ylim=c(0,0.16), xlab="Offspring", ylab="Frequency",main="JSU model",probability=TRUE); 