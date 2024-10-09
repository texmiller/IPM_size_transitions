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

require(car); require(zoo); require(moments); require(mgcv); 
require(gamlss); require(AICcmodavg); 
require(tidyverse); require(maxLik); require(qgam); require(exactLTRE)
library(gamlss.dist)

source("../code/Diagnostics.R"); 
source("../code/fitChosenDists.R"); 
source("../code/variance_diagnostics.R"); 

################# Read in data 
XH = read.csv("Vulpicida raw data.csv"); 
e = order(XH$t0); XH = XH[e,]; 
XH_true = XH; 

############################# bootstrap! ################
bootreps = 101; 
XG = XS = matrix(NA,bootreps,12); ## to hold life history outputs 

for(bootrep in 1:bootreps){

e = sample(1:nrow(XH_true),nrow(XH_true), replace=TRUE); 
XH = XH_true[e,]; 
e = order(XH$t0); XH = XH[e,]; 
if(bootrep==1) XH = XH_true; ### first time, use the real data! 

############### Survival modeling 
fit2 = glm(survival~sqrt(t0) + t0,data=XH,family="binomial"); 
surv_coefs = coef(fit2); 

## survival, set up so that area can be negative and survival evaluates at zero
sx = function(x)  {
  a = pmax(0,x)
  u1 = surv_coefs[1] + surv_coefs[2]*sqrt(a) - surv_coefs[3]*a 
  p1 = exp(u1)/(1+exp(u1)); 	 
  p1[x < 0]=0; 
  return(p1);
}	

############### Fecundity function (from Shriver et al. model) 
## Area can be negative and fecundity evaluates at zero
fx = function(Area) {
    a = pmax(0,Area)
	r = sqrt(a/pi);
	u = 2*pi*r; 
	u[Area <0] = 0; 
	return(0.047*u)
}		

########################################################################
#  Gaussian growth model: final choice  
########################################################################
XH = XH[XH$survival==1,]; ## discard those dead at (t+1); 
fitGAU <- gam(list(t1~t0 + I(t0^2), ~t0 + I(t0^2)), data=XH, gamma=1.4,family=gaulss())

########################################################################
#  Non-Gaussian growth model: final choice  
########################################################################

fitSHASH2 <- gam(list(t1 ~ t0 + I(t0^2), # <- location 
                        ~ s(t0),  # <- log-scale
                        ~ s(t0),  # <- skewness
                        ~ 1), # <- log-kurtosis
                      data = XH, gamma=1.4, family = shash,  optimizer = "efs")
SHASH_pred2<-predict(fitSHASH2,type="response")

## Save the fitted size-dependent parameters of the best SHASH model. 
muE <- fitSHASH2$fitted[ , 1]
sigE <- exp(fitSHASH2$fitted[ , 2])
epsE <- fitSHASH2$fitted[ , 3]
delE <- exp(fitSHASH2$fitted[ , 4])
t0 <- XH$t0; 
# save(muE,sigE,epsE,delE,t0,file="SHASHfuns.Rdata"); 

################################################################# 
## Make the IPMs 
#################################################################


## GAUSSIAN
L=-1; U=10; L1 = 0.2; U1 = 7; # limits of the data for eviction-prevention 

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
IPM_G = mk_K_G(200,L=L,U=U,L1 = L1, U1 = U1); Re(eigen(IPM_G$K)$values[1]);
chop_meshpts = pmax(pmin(IPM_G$meshpts,U1),L1) 
  
########################################### 
## SHASH 
###########################################
muE_fun <- splinefun(t0,muE, method="natural")
log_sigE_fun <- splinefun(t0,log(sigE), method="natural")
epsE_fun <- splinefun(t0,epsE, method="natural")
delE_fun <- splinefun(t0,delE, method="natural")

G_z1z_S = function(z1,z) {sx(z)*dSHASHo2(z1,muE_fun(z),exp(log_sigE_fun(z)),epsE_fun(z),delE_fun(z))}

mk_K_S <- function(m, L, U, L1, U1) {
	# mesh points 
	h <- (U - L)/m
	meshpts <- L + ((1:m) - 1/2) * h
	K <- P <- h * (outer(meshpts, pmax(pmin(meshpts, U1),L1),G_z1z_S))
	K[1,] = K[1,] + matrix(fx(meshpts),nrow=1); F = K - P; 
	return(list(K = K, meshpts = meshpts, P = P, F = F))
}

IPM_S =mk_K_S(200,L=L,U=U,L1 = L1, U1 = U1); Re(eigen(IPM_S$K)$values[1]);

################################################################# 
## Compare life history attributes 
#################################################################

source("../code/metaluck_fns_CMH.R"); 

### Gaussian
matU = IPM_G$P; matF = IPM_G$F; c0 = rep(0,nrow(matU)); c0[1]=1; 
XG[bootrep,] = c(mean_lifespan(matU, mixdist=c0),
var_lifespan(matU, mixdist=c0)^0.5,
skew_lifespan(matU, mixdist=c0),
mean_LRO(matU,matF,mixdist=c0), 
var_LRO_mcr(matU,matF,mixdist=c0)^0.5,
skew_LRO(matU,matF,mixdist=c0), 
prob_repro(matU,matF)[1], 
mean_age_repro(matU,matF,mixdist=c0), 
lifespan_reproducers(matU,matF,mixdist=c0), 
gen_time_Ta(matU,matF), 
gen_time_mu1_v(matU,matF), 
gen_time_R0(matU,matF))  

### SHASH
matU = IPM_S$P; matF = IPM_S$F; c0 = rep(0,nrow(matU)); c0[1]=1; 
XS[bootrep,] = c(
mean_lifespan(matU, mixdist=c0),
var_lifespan(matU, mixdist=c0)^0.5,
skew_lifespan(matU, mixdist=c0),
mean_LRO(matU,matF,mixdist=c0), 
var_LRO_mcr(matU,matF,mixdist=c0)^0.5,
skew_LRO(matU,matF,mixdist=c0), 
prob_repro(matU,matF)[1], 
mean_age_repro(matU,matF,mixdist=c0), 
lifespan_reproducers(matU,matF,mixdist=c0), 
gen_time_Ta(matU,matF), 
gen_time_mu1_v(matU,matF), 
gen_time_R0(matU,matF))  


cat(bootrep,"\n"); 

}