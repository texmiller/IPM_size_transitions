############################################################################
## Corals (Bruno et al., Akumal) case study. Data from Bruno et al. (2011) 
###########################################################################
rm(list=ls(all=TRUE))

### move to the right local directory 
tom = "C:/Users/tm9/Dropbox/github/IPM_size_transitions"
steve = "c:/repos/IPM_size_transitions" 
home = ifelse(Sys.info()["user"] == "Ellner", steve, tom)
setwd(home); setwd("coral"); 

require(mgcv); require(gamlss.dist); 

source("../code/metaluck_fns_CMH.R"); 
source("../code/bca.R"); 

## source("Akumal_corals_boot.R"); 

#######################################################################
## Create complete data frame XC. Bootstrapping will draw from this. 
###################################################################### 
source("AkumalCoralsSetup.R"); 
recruitSizes_true = recruitSizes 

######################################################################
# Demographic functions from Bruno et al. (2011) 
######################################################################
fan.height=function(area) {1.3*sqrt(area)}
fan.area=function(height) {height*height/1.69}

nlog=function(x) x^(1/3) 

### Size distribution of new recruits 
recruitDensity=density(log(recruitSizes),bw="SJ"); 
recruitFun=approxfun(recruitDensity$x,recruitDensity$y,yleft=0,yright=0);

### Additional mortality of large fans, based on 1997 Keys data 
### Note, here x is cube-root of area! 
### NOTE: estimating this pulled in multiple sources of data, in complicated ways, 
### and was based on making the fitted IPM match observed size distributions.  
### Here we DO NOT include uncertainty arising from uncertainty in bigmu. 
bigmu=function(x) {
	p=c(2.3,8.7,0.093); 
	u=exp(p[1]*(x-p[2]))
	return(p[3]*u/(1+u))
}

## Base survival function. Here z is log area.  
s_z <- function(z) {
    x=exp(z)^(1/3); ## cube root of area, as z = log(area)  
    return(1 - mH(x))
}  

### Mortality of healthy fans. Here x is cube root of area. 
mH=function(x) {
	p=surv_coefs; 
	u= p[1]+p[2]*x; 
	exp(u)/(1+exp(u)) + bigmu(x) 
}  

### Fecundity function for age at reproduction, so reproducing
### is guaranteed once you get above the threshold   
b_z <- function(z) {
      x=exp(z); ### area 
      xmin=fan.area(20); ### cutoff for reproduction 
      ifelse(x>xmin,40,0)
}

########## Growth functions and P kernel functions 
G_z1zPilot <- function(z1,z,pars=NULL){
  pred = predict(fitGAU, newdata = data.frame(t0=z))
  return(dnorm(z1,mean=pred[,1],sd=exp(pred[,2])))
}  

G_z1zSHASH <- function(z1,z,pars=NULL){
  pred = predict(fitSHASH, newdata = data.frame(t0=z))
  return(dSHASHo2(x=z1, mu=pred[,1],
                  sigma = exp(pred[,2]), 
                  nu = pred[,3], 
                  tau = exp(pred[,4])))
} 

P_z1zSHASH = function(z1,z) {s_z(z) * G_z1zSHASH(z1, z)}
P_z1zPilot = function(z1,z) {s_z(z) * G_z1zPilot(z1, z)}

############################################################
### Function to make the P kernel, with floor and ceiling 
## to avoid eviction
############################################################
mk_P_ceiling <- function(m, L, U, L1, U1, Pfun) {
	# mesh points 
	h <- (U - L)/m;
	meshpts <- L + ((1:m) - 1/2) * h;
    truncMesh = pmin(meshpts,U1); 
    truncMesh = pmax(meshpts,L1); 
	P <- h * (outer(meshpts, truncMesh, Pfun));
    return(list(meshpts = meshpts, h=h, P = P))
}

############ Sanity check: do we re-create the fitted survival function? YES! 
############ Create data frame of healthy fans 
XH=subset(XC,(State1=="H"));
e=(XH$State2=="H")|(XH$State2=="D")|(XH$State2=="I"); 
e[is.na(e)]=FALSE; XH=XH[e,]; 
XH=XH[!is.na(XH$Area1),]; 
XH$Died=as.numeric(XH$State2=="D")
XH$t0 = log(XH$Area1); XH$t1 = log(XH$Area2); 

fitD=glm(Died~nlog(Area1),data=XH,family="binomial");  
coef(fitD); ## same as in Bruno et al. 2011 

#################################################################
# Fit growth models to get smoothing parameters 
#################################################################
fitGAU <- gam(list(t1~s(t0),~s(t0)), data=XH, gamma=1.4,family=gaulss())
sp.GAU = fitGAU$sp; 

fitSHASH <- gam(list(t1 ~ s(t0), # <- location 
                     ~ s(t0),   # <- log-scale
                     ~ s(t0,k=4),   # <- skewness
                     ~ s(t0,k=4)), # <- log-kurtosis
                data = XH, gamma=1.4, 
                family = shash,  
                optimizer = "efs")
sp.SHASH = fitSHASH$sp; 

##################################################################
##  Doing the BOOTSTRAP 
##################################################################
XC_true = XC; ### save the complete data set! 
recruitDensity=density(log(recruitSizes_true),bw="SJ"); 
bw_true = recruitDensity$bw; 

#ncores = detectCores(logical=FALSE)-2; 
#c1 = makeCluster(ncores)
#registerDoParallel(c1); 
#clusterExport(c1,varlist=objects()); 

bootreps = 501; 
traits_G = traits_S = matrix(NA,bootreps,2); 

for(bootrep in 1:bootreps){
e = sample(1:nrow(XC),nrow(XC), replace=TRUE); 
XC = XC_true[e,]; 
e = order(XC$Area1); XC = XC[e,]; 
recruitSizes = sample(recruitSizes_true,replace=TRUE) 

### first time, use the real data! 
if(bootrep==1) {XC = XC_true; recruitSizes = recruitSizes_true} 

# Create data frame of healthy fans 
XH=subset(XC,(State1=="H"));
e=(XH$State2=="H")|(XH$State2=="D")|(XH$State2=="I"); 
e[is.na(e)]=FALSE; XH=XH[e,]; 
XH=XH[!is.na(XH$Area1),]; 
XH$Died=as.numeric(XH$State2=="D")
XH$t0 = log(XH$Area1); XH$t1 = log(XH$Area2); 

############### Survival modeling 
fitD=glm(Died~XH$t0,data=XH,family="binomial");  
surv_coefs = coef(fitD); 

############### size distribution of new recruits 
recruitDensity=density(log(recruitSizes),bw=bw_true); ## Don't refit bandwidth (standard advice) 
recruitFun=approxfun(recruitDensity$x,recruitDensity$y,yleft=0,yright=0);

################################################################# 
## Make the GAUSSIAN kernels, compute traits and store 
#################################################################
fitGAU <- gam(list(t1~s(t0),~s(t0)), data=XH, sp=sp.GAU,family=gaulss())
m=500; L=0; U=10; L1 = 1.35; U1 = 7.9; 
out = mk_P_ceiling(m, L, U, L1, U1,P_z1zPilot); 

matU = out$P; 
matF = 0*matU; matF[1,] = b_z(out$meshpts); 
coral_c0 = recruitFun(out$meshpts); 
coral_c0 = coral_c0/sum(coral_c0); 

traits_G[bootrep,] <-c(
  mean_lifespan(matU, mixdist=coral_c0),
  mean_age_repro(matU,matF,mixdist=coral_c0)
)

################################################################# 
## Make the SHASH kernels, compute traits and store 
#################################################################
fitSHASH <- gam(list(t1 ~ s(t0), # <- location 
                     ~ s(t0),   # <- log-scale
                     ~ s(t0,k=4),   # <- skewness
                     ~ s(t0,k=4)), # <- log-kurtosis
                data = XH, sp=sp.SHASH,
                family = shash,  
                optimizer = "efs")

m=500; L=0; U=10; L1 = 1.35; U1 = 7.9; 
out = mk_P_ceiling(m, L, U, L1, U1,P_z1zSHASH); 

matU = out$P; 
matF = 0*matU; matF[1,] = b_z(out$meshpts); 

traits_S[bootrep,] <-c(
  mean_lifespan(matU, mixdist=coral_c0),
  mean_age_repro(matU,matF,mixdist=coral_c0)
)

cat(bootrep, signif(traits_G[bootrep,],3), signif(traits_S[bootrep,],3), "\n"); 
if(bootrep%%100 == 0) save.image(file="corals_boot.Rdata"); 

}
traits_G_true = traits_G[1,]; 
traits_S_true = traits_S[1,]; 
traits_G_boot = traits_G[-1,]; 
traits_S_boot = traits_S[-1,]

save.image(file="corals_boot.Rdata"); 

################################################### 
## Output results: GAUSSIAN 
###################################################

xbar = apply(traits_G_boot,2,mean); 
xsd = apply(traits_G_boot,2,var)^0.5; 

### Compute BCA intervals 
CI_G = matrix(NA,2,2); 
for(j in 1:2) {
	CI_G[1:2,j]= bca(theta = traits_G_boot[,j], theta_hat = traits_G_true[j], a=0, conf.level = 0.95) 
}

cat("GAUSSIAN", "\n"); 
cat("point    ", signif(traits_G_true,3),"\n"); 
cat("boot mean", signif(xbar,3),"\n"); 
cat("boot sd  ", signif(xsd,3), "\n"); 
cat("BCA 95% confidence intervals", "\n") 
print(signif(CI_G,3)); 

graphics.off(); par(mfrow=c(2,1)); 
for(j in 1:2) hist(traits_G_boot[,j]); 

################################################### 
## Output results: SHASH
###################################################
xbar = apply(traits_S_boot,2,mean); 
xsd = apply(traits_S_boot,2,var)^0.5; 

### Compute BCA intervals 
CI_S = matrix(NA,2,2); 
for(j in 1:2) {
	CI_S[1:2,j]=bca(theta = traits_S_boot[,j], theta_hat = traits_S_true[j], a=0, conf.level = 0.95) 
}

cat("SHASH", "\n"); 
cat("point    ", signif(traits_S_true,3),"\n"); 
cat("boot mean", signif(xbar,3),"\n"); 
cat("boot sd  ", signif(xsd,3), "\n"); 
cat("BCA 95% confidence intervals", "\n") 
print(signif(CI_S,3)); 

graphics.off(); par(mfrow=c(2,1)); 
for(j in 1:2) hist(traits_S_boot[,j]); 
