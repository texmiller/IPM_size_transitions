################################################################
# Lichen case study, data from Peterson (de Marche) et al. MEE 
################################################################
rm(list=ls(all=TRUE))

### move to the right local directory 
tom = "C:/Users/tm9/Dropbox/github/IPM_size_transitions"
steve = "c:/repos/IPM_size_transitions" 
home = ifelse(Sys.info()["user"] == "Ellner", steve, tom)
setwd(home); setwd("lichen"); 

## source("Vulpicida_boot.R"); 

# require(car); require(zoo); require(moments); require(mgcv); 
# require(gamlss); require(AICcmodavg); 
# require(tidyverse); require(maxLik); require(qgam); require(exactLTRE)

require(mgcv); require(gamlss.dist); require(exactLTRE); 
require(maxLik); require(doParallel); require(moments); 

source("../code/metaluck_fns_CMH.R"); 
source("../code/bca.R"); 

##################################################################
##  Functions for making the IPMs 
##################################################################

############### Fecundity function (from Shriver et al. model) 
## Size (sqrt area) can be negative and fecundity evaluates at zero
## Recruitment is 0.047*circumference
## Note, bootstrap does not include uncertainty in this! 
fx = function(t0) {
	area = pmax(0,t0)^2
	r = sqrt(area/pi);
	u = 2*pi*r; 
	return(0.047*u)
}		

## Survival, best-fitting regression model on full data set 
## Size x (=sqrt area) can be negative and survival evaluates at zero
sx = function(x)  {
  a = pmax(0,x)
  u1 = surv_coefs[1] + surv_coefs[2]*sqrt(a) + surv_coefs[3]*a 
  p1 = exp(u1)/(1+exp(u1)); 	 
  p1[x < 0]=0; 
  return(p1);
}	

############### JSU log-likelihood function 
LogLikJSU=function(pars,response){
  dJSU(response, 
       mu=pars[1] + pars[2]*XH$t0 + pars[3]*(XH$t0^2),
       sigma=exp(pars[4] + pars[5]*XH$t0 + pars[6]*(XH$t0^2)),
       nu = pars[7] + pars[8]*XH$t0,
       tau = exp(pars[9]), log=TRUE)
}

############## Functions to make the Gaussian IPM 
mu_G = function(x) {coef(fitGAU)[1]  +  coef(fitGAU)[2]*x  + coef(fitGAU)[3]*x^2} 
sd_G = function(x) {exp(coef(fitGAU)[4]   + coef(fitGAU)[5]*x + coef(fitGAU)[6]*x^2)}  
G_z1z_G = function(z1,z) {sx(z)*dnorm(z1,mu_G(z),sd_G(z)) }

mk_K_G <- function(m, L, U, L1, U1) {
	# mesh points 
	h <- (U - L)/m
	meshpts <- L + ((1:m) - 1/2) * h
	K <- P <- h * (outer(meshpts, pmax(pmin(meshpts, U1),L1), G_z1z_G))
	K[1,] = K[1,] + matrix(fx(meshpts),nrow=1); F = K - P; 
	return(list(K = K, meshpts = meshpts, P = P, F = F))
}

############## Functions to make the JSU IPM 
muJSU_fun <- function(z){outJSU$estimate[1]+outJSU$estimate[2]*z+outJSU$estimate[3]*z^2}
sigJSU_fun <- function(z){exp(outJSU$estimate[4]+outJSU$estimate[5]*z+outJSU$estimate[6]*z^2)}
epsJSU_fun <- function(z){outJSU$estimate[7]+outJSU$estimate[8]*z}
delJSU_fun <- function(z){exp(outJSU$estimate[9])}

G_z1z_J = function(z1,z) {sx(z)*dJSU(z1,muJSU_fun(z),sigJSU_fun(z),epsJSU_fun(z),delJSU_fun(z))}

mk_K_J <- function(m, L, U, L1, U1) {
  # mesh points 
  h <- (U - L)/m
  meshpts <- L + ((1:m) - 1/2) * h
  K <- P <- h * (outer(meshpts, pmax(pmin(meshpts, U1),L1),G_z1z_J))
  K[1,] = K[1,] + matrix(fx(meshpts),nrow=1); F = K - P; 
  return(list(K = K, meshpts = meshpts, P = P, F = F))
}

##################################################################
##  Read in the data 
##################################################################
XH_true = read.csv("Vulpicida raw data.csv"); 
e = order(XH_true$t0); XH_true = XH_true[e,]; 

##################################################################
##  Doing the jacknife 
##################################################################
ncores = detectCores(logical=FALSE)-2; 
c1 = makeCluster(ncores)
registerDoParallel(c1); 
clusterExport(c1,varlist=objects()); 

jackreps = 1000 
e_jack = sample(1:nrow(XH_true),jackreps, replace=FALSE); 
clusterExport(c1,varlist=objects()); 

################### START jackknife loop !!!!!!!!!!!!
traits = foreach(jackrep = 1:jackreps,.packages = c("maxLik","gamlss.dist","mgcv"))%dopar%
{ 

omitted = e_jack[jackrep]; 
XH = XH_true[-omitted,]; 
e = order(XH$t0); XH = XH[e,]; 

############### Survival modeling 
fit2 = glm(survival~sqrt(t0) + t0,data=XH,family="binomial"); 
surv_coefs = coef(fit2); 

########################################################################
#  Gaussian growth model: final choice  
########################################################################
XH = XH[XH$survival==1,]; ## discard those dead at (t+1); 
fitGAU <- gam(list(t1~t0 + I(t0^2), ~t0 + I(t0^2)), data=XH, gamma=1.4,family=gaulss())

################################################################# 
## Make the GAUSSIAN IPM, compute traits and store 
#################################################################
L=-1; U=10; L1 = 0.2; U1 = 7; # limits of the data for eviction-prevention 

IPM_G = mk_K_G(200,L=L,U=U,L1 = L1, U1 = U1); lambda = Re(eigen(IPM_G$K)$values[1]);
matU = IPM_G$P; matF = IPM_G$F; 
lichen_c0 = rep(0,nrow(matU)); lichen_c0[1]=1

traits_G <-c(
  lambda, 
  mean_lifespan(matU, mixdist=lichen_c0),
  mean_LRO(matU,matF,mixdist=lichen_c0),
  mean_age_repro(matU,matF,mixdist=lichen_c0),
  gen_time_mu1_v(matU,matF)
)

################################################################# 
## Make the JSU IPM, compute traits and store 
#################################################################

p0<-c(coef(fitGAU)[1:6],0,0,0)
## fit with maxlik
outJSU=maxLik(logLik=LogLikJSU,start=p0*exp(0.2*rnorm(length(p0))), response=XH$t1,
              method="BHHH",control=list(iterlim=5000,printLevel=0),finalHessian=FALSE); 
outJSU=maxLik(logLik=LogLikJSU,start=outJSU$estimate,response=XH$t1,
              method="NM",control=list(iterlim=5000,printLevel=0),finalHessian=FALSE); 
outJSU=maxLik(logLik=LogLikJSU,start=outJSU$estimate,response=XH$t1,
              method="BHHH",control=list(iterlim=5000,printLevel=0),finalHessian=FALSE); 

IPM_J =mk_K_J(200,L=L,U=U,L1 = L1, U1 = U1); lambda = Re(eigen(IPM_J$K)$values[1]);

### JSU
matU = IPM_J$P; matF = IPM_J$F; c0 = rep(0,nrow(matU)); c0[1]=1; 
traits_JSU<-c(
  lambda, 
  mean_lifespan(matU, mixdist=lichen_c0),
  mean_LRO(matU,matF,mixdist=lichen_c0),
  mean_age_repro(matU,matF,mixdist=lichen_c0),
  gen_time_mu1_v(matU,matF)
)

c(traits_G,traits_JSU); 

} 
################### END jackknife loop !!!!!!!!!!!!

#################################################################
# Estimate acceleration parameter a from jackknifed estimates
#################################################################
theta = matrix(unlist(traits),ncol=10,byrow=TRUE); 
a_G = apply(theta[,1:5],2, skewness)/6; 
a_JSU = apply(theta[,6:10],2,skewness)/6; 


##################################################################
##  Doing the bootstrap
##################################################################

bootreps = 501; 
################### START bootstrap loop !!!!!!!!!!!!
traits = foreach(bootrep = 1:bootreps,.packages = c("maxLik","gamlss.dist","mgcv"))%dopar%
{ 
e = sample(1:nrow(XH_true),nrow(XH_true), replace=TRUE); 
XH = XH_true[e,]; 
e = order(XH$t0); XH = XH[e,]; 
if(bootrep==1) XH = XH_true; ### first time, use the real data! 

############### Survival modeling 
fit2 = glm(survival~sqrt(t0) + t0,data=XH,family="binomial"); 
surv_coefs = coef(fit2); 

########################################################################
#  Gaussian growth model: final choice  
########################################################################
XH = XH[XH$survival==1,]; ## discard those dead at (t+1); 
fitGAU <- gam(list(t1~t0 + I(t0^2), ~t0 + I(t0^2)), data=XH, gamma=1.4,family=gaulss())

################################################################# 
## Make the GAUSSIAN IPM, compute traits and store 
#################################################################
L=-1; U=10; L1 = 0.2; U1 = 7; # limits of the data for eviction-prevention 

IPM_G = mk_K_G(200,L=L,U=U,L1 = L1, U1 = U1); lambda = Re(eigen(IPM_G$K)$values[1]);
matU = IPM_G$P; matF = IPM_G$F; 
lichen_c0 = rep(0,nrow(matU)); lichen_c0[1]=1

traits_G <-c(
  lambda, 
  mean_lifespan(matU, mixdist=lichen_c0),
  mean_LRO(matU,matF,mixdist=lichen_c0),
  mean_age_repro(matU,matF,mixdist=lichen_c0),
  gen_time_mu1_v(matU,matF)
)

cat(bootrep, signif(lambda,3), "\n"); 

################################################################# 
## Make the JSU IPM, compute traits and store 
#################################################################

p0<-c(coef(fitGAU)[1:6],0,0,0)
## fit with maxlik
outJSU=maxLik(logLik=LogLikJSU,start=p0*exp(0.2*rnorm(length(p0))), response=XH$t1,
              method="BHHH",control=list(iterlim=5000,printLevel=0),finalHessian=FALSE); 
outJSU=maxLik(logLik=LogLikJSU,start=outJSU$estimate,response=XH$t1,
              method="NM",control=list(iterlim=5000,printLevel=0),finalHessian=FALSE); 
outJSU=maxLik(logLik=LogLikJSU,start=outJSU$estimate,response=XH$t1,
              method="BHHH",control=list(iterlim=5000,printLevel=0),finalHessian=FALSE); 

IPM_J =mk_K_J(200,L=L,U=U,L1 = L1, U1 = U1); lambda = Re(eigen(IPM_J$K)$values[1]);

### JSU
matU = IPM_J$P; matF = IPM_J$F; c0 = rep(0,nrow(matU)); c0[1]=1; 
traits_JSU<-c(
  lambda, 
  mean_lifespan(matU, mixdist=lichen_c0),
  mean_LRO(matU,matF,mixdist=lichen_c0),
  mean_age_repro(matU,matF,mixdist=lichen_c0),
  gen_time_mu1_v(matU,matF)
)

cat(bootrep, signif(lambda,3), "\n")
c(traits_G,traits_JSU); 

} 
stopCluster(c1); ################### END bootstrap loop !!!!!!!!!!!!

theta = matrix(unlist(traits),ncol=10,byrow=TRUE); 
traits_G = theta[,1:5]; traits_JSU = theta[,6:10]; 

traits_G_true = traits_G[1,]; ## should match the paper 
traits_JSU_true = traits_JSU[1,]; ## should match the paper 

################################################### 
## Output results: GAUSSIAN 
###################################################

traits_G_boot = traits_G[-1,]; 
xbar = apply(traits_G_boot,2,mean); 
xsd = apply(traits_G_boot,2,var)^0.5; 

### Compute BCA intervals 
CI_G = matrix(NA,2,5); 
for(j in 1:5) {
	CI_G[1:2,j]=bca(theta = traits_G_boot[,j], theta_hat = traits_G_true[j], a = a_G[j], conf.level = 0.95) 
}

cat("GAUSSIAN", "\n"); 
cat("point    ", signif(traits_G_true,3),"\n"); 
cat("boot mean", signif(xbar,3),"\n"); 
cat("boot sd  ", signif(xsd,3), "\n"); 
cat("BCA 95% confidence intervals", "\n") 
print(signif(CI_G,3)); 

graphics.off(); par(mfrow=c(3,2)); 
for(j in 1:5) hist(traits_G_boot[,j]); 


###################################################
##### Output results: JSU
###################################################

traits_JSU_boot = traits_JSU[-1,]; 
xbar = apply(traits_JSU_boot,2,mean); 
xsd = apply(traits_JSU_boot,2,var)^0.5; 

### Compute BCA intervals 
CI_J = matrix(NA,2,5); 
for(j in 1:5) {
	CI_J[1:2,j]=bca(theta  = traits_JSU_boot[,j], theta_hat = traits_JSU_true[j], a = a_JSU[j], conf.level = 0.95) 
}

cat("JSU", "\n");  
cat("point    ", signif(traits_JSU_true,3),"\n"); 
cat("boot mean", signif(xbar,3),"\n"); 
cat("boot sd  ", signif(xsd,3), "\n"); 
cat("BCA 95% confidence intervals", "\n") 
print(signif(CI_J,3)); 

dev.new(); par(mfrow=c(3,2)); 
for(j in 1:5) hist(traits_JSU_boot[,j]); 

save.image(file="Vulpicida_boot.Oct.16.Rdata"); 