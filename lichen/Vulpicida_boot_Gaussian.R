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

source("../code/metaluck_fns_CMH.R"); 
source("../code/bca.R"); 

############### Fecundity function (from Shriver et al. model) 
## Area can be negative and fecundity evaluates at zero
fx = function(Area) {
    a = pmax(0,Area)
	r = sqrt(a/pi);
	u = 2*pi*r; 
	u[Area <0] = 0; 
	return(0.047*u)
}	

################# Read in data 
XH_true = read.csv("Vulpicida raw data.csv"); 
e = order(XH_true$t0); XH_true = XH_true[e,]; 

bootreps = 501; 
traits_G = traits_JSU = matrix(NA,bootreps,5); ## to hold life history outputs 

for(bootrep in 1:bootreps){ ################### START bootstrap loop !!!!!!!!!!!!

e = sample(1:nrow(XH_true),nrow(XH_true), replace=TRUE); 
XH = XH_true[e,]; 
e = order(XH$t0); XH = XH[e,]; 
if(bootrep==1) XH = XH_true; ### first time, use the real data! 

# hist(XH$t0,25,xlim=range(XH_true$t0)); 

############### Survival modeling 
fit2 = glm(survival~sqrt(t0) + t0,data=XH,family="binomial"); 
surv_coefs = coef(fit2); 

## survival, set up so that area can be negative and survival evaluates at zero
sx = function(x)  {
  a = pmax(0,x)
  u1 = surv_coefs[1] + surv_coefs[2]*sqrt(a) + surv_coefs[3]*a 
  p1 = exp(u1)/(1+exp(u1)); 	 
  p1[x < 0]=0; 
  return(p1);
}	
	
########################################################################
#  Gaussian growth model: final choice  
########################################################################
XH = XH[XH$survival==1,]; ## discard those dead at (t+1); 
fitGAU <- gam(list(t1~t0 + I(t0^2), ~t0 + I(t0^2)), data=XH, gamma=1.4,family=gaulss())

################################################################# 
## Make the GAUSSIAN IPM, compute traits and store 
#################################################################
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

IPM_G = mk_K_G(200,L=L,U=U,L1 = L1, U1 = U1); lambda = Re(eigen(IPM_G$K)$values[1]);
matU = IPM_G$P; matF = IPM_G$F; 
lichen_c0 = rep(0,nrow(matU)); lichen_c0[1]=1

traits_G[bootrep,]<-c(
  lambda, 
  mean_lifespan(matU, mixdist=lichen_c0),
  mean_LRO(matU,matF,mixdist=lichen_c0),
  mean_age_repro(matU,matF,mixdist=lichen_c0),
  gen_time_mu1_v(matU,matF)
)

cat(bootrep, signif(lambda,3), "\n"); 

} ################### END bootstrap loop !!!!!!!!!!!!

traits_G_true = traits_G[1,]; ## should match the paper 

traits_G_boot = traits_G[-1,]; 
xbar = apply(traits_G_boot,2,mean); 
xsd = apply(traits_G_boot,2,var)^0.5; 

#### Compute a 'basic bootstrap interval' 
BB = matrix(NA,bootreps-1,5); 
for(j in 1:5) BB[,j] = 2*traits_G_true[j] - traits_G_boot[,j] 
LQ = apply(BB,2,function(x) quantile(x,0.025)); 
UQ = apply(BB,2,function(x) quantile(x,0.975));  


### Compute BCA intervals 
CI = matrix(NA,2,5); 
for(j in 1:5) {
	CI[1:2,j]=bca(traits_G_boot[,j], conf.level = 0.95) 
}

cat("point    ", signif(traits_G_true,3),"\n"); 
cat("boot mean", signif(xbar,3),"\n"); 
cat("boot sd  ", signif(xsd,3), "\n"); 
cat("BCA 95% confidence intervals", "\n") 
print(signif(CI,3)); 

graphics.off(); par(mfrow=c(3,2)); 
for(j in 1:5) hist(traits_G_boot[,j]); 