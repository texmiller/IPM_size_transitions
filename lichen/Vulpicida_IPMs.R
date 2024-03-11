################################################################
# Lichen case study, data from Peterson (de Marche) et al. MEE 
################################################################

rm(list=ls(all=TRUE))

### move to the right local directory 
tom = "C:/Users/tm9/Dropbox/github/IPM_size_transitions"
steve = "c:/repos/IPM_size_transitions" 
home = ifelse(Sys.info()["user"] == "Ellner", steve, tom)
setwd(home); setwd("lichen"); 

require(car); require(zoo); require(moments); require(mgcv); 
require(gamlss); require(AICcmodavg); 
require(tidyverse); require(maxLik); require(qgam); require(exactLTRE)
library(gamlss.dist)

source("../code/Diagnostics.R"); 
source("../code/fitChosenDists.R"); 
source("../code/variance_diagnostics.R"); 

PLOTTING = TRUE; 

## quartile-based estimates of mean and sd
## see https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/1471-2288-14-135
Q.mean<-function(q.25,q.50,q.75){(q.25+q.50+q.75)/3}
Q.sd<-function(q.25,q.75){(q.75-q.25)/1.35}

## Steve's functions for NP skew and kurtosis
Q.skewness<-function(q.10,q.50,q.90){(q.10 + q.90 - 2*q.50)/(q.90 - q.10)}
Q.kurtosis<-function(q.05,q.25,q.75,q.95){
  qN = qnorm(c(0.05,0.25,0.75,0.95))
  KG = (qN[4]-qN[1])/(qN[3]-qN[2])
  return(((q.95-q.05)/(q.75-q.25))/KG - 1)
}

################# Read in data 
XH = read.csv("Vulpicida raw data.csv"); 
e = order(XH$t0); XH = XH[e,]; 


############### Survival modeling 
plot(XH$t0,XH$survival); 

fit1 = glm(survival~sqrt(t0),data=XH,family="binomial"); 
fit2 = glm(survival~sqrt(t0) + t0,data=XH,family="binomial"); 
fit3 = gam(survival~s(sqrt(t0)),data=XH,family="binomial"); 
## AIC goes for fit2, \Delta AIC of about 3. 
#(Intercept)    sqrt(t0)          t0 
#  -1.649455    5.363115   -1.249375 

## survival, set up so that area can be negative and survival evaluates at zero
sx = function(x)  {
  a = pmax(0,x)
  u1 = -1.649455  + 5.363115*sqrt(a) - 1.249375*a 
  p1 = exp(u1)/(1+exp(u1)); 	 
  p1[x < 0]=0; 
  return(p1);
}	

nbins = 12; 
indexx = seq_along(XH$t0)
u = nbins*indexx/(1 + max(indexx)); 
u = floor(u); pbar = tbar = numeric(nbins); 
for(j in 1:nbins) {pbar[j] = mean(XH$survival[u==j]); tbar[j]=mean(XH$t0[u==j])}
plot(tbar,pbar,pch=16,xlim=range(XH$t0),ylim=c(0.9,1)); 
points(XH$t0,sx(XH$t0),type="l",lty=1); 

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
#  Growth modeling 
########################################################################

XH = XH[XH$survival==1,]; ## discard those dead at (t+1); 

fitGAU <- gam(list(t1~s(t0),~s(t0)), data=XH, gamma=1.4,family=gaulss())
summary(fitGAU); # plot(fitGAU); 

## The mean looks almost linear; is there evidence against this? 
fitGAU0 <- gam(list(t1~t0,~s(t0)), data=XH, gamma=1.4, family=gaulss())
AIC(fitGAU); AIC(fitGAU0); # Somewhat, Delta AIC of about 4 in favor of the spline 

## the log(sigma) fit looks almost linear; is there evidence against this? 
fitGAU00 <- gam(list(t1~s(t0),~t0), data=XH, gamma=1.4,family=gaulss())
AIC(fitGAU); AIC(fitGAU00);  # Delta AIC about 37, so very strong evidence

## what about quadratic? 
fitGAU01 <- gam(list(t1~s(t0,k=12), ~t0 + I(t0^2)), data=XH, gamma=1.4,family=gaulss())
AIC(fitGAU); AIC(fitGAU01);  # Delta AIC about 1

## what about quadratic for the mean, as well? 
fitGAU22 <- gam(list(t1~t0 + I(t0^2), ~t0 + I(t0^2)), data=XH, gamma=1.4,family=gaulss())
AIC(fitGAU); AIC(fitGAU22);  # the quadratic-quadratic wins by a hair, Delta AIC about 1.  

## proceeding with fitGAU22 as "best" Gaussian model
XH$fitted_sd <- 1/predict(fitGAU22,type="response")[,2]
XH$fitted = predict(fitGAU22,type="response")[,1]
XH$scaledResids=residuals(fitGAU22,type="pearson")

## checking that I know how to use the residuals function
plot(residuals(fitGAU22,type="response"),XH$t1-XH$fitted)
plot(XH$scaledResids,residuals(fitGAU22,type="response")/XH$fitted_sd)

mean(XH$scaledResids); sd(XH$scaledResids); ## all good 

### Parametric SD function is nearly the same as the spline 
plot(XH$t0, log(XH$fitted_sd)); 
points(XH$t0, log( 1/predict(fitGAU,type="response")[,2]), col="red" )

####################### Save the best-fitting Gaussian as fitGAU. 
fitGAU = fitGAU22; rm(fitGAU22); rm(fitGAU0); rm(fitGAU00);  

######### Diagnostics on fitted parametric SD function: no problems! 
c1<- makeCluster(8); 
registerDoParallel(c1);
out = multiple_levene_test(XH$fitted, XH$scaledResids, 3, 8, 2000);
out$p_value; ## > 0.89

out = multiple_bs_test(XH$fitted, XH$scaledResids, 4, 8, 2000) 
out$p_value; ## > 0.78; 
stopCluster(c1); 

## quantile regressions on standardized  residuals 
S.05<-qgam(scaledResids~s(t0,k=6), data=XH,qu=0.05)
S.10<-qgam(scaledResids~s(t0,k=6), data=XH,qu=0.1)
S.25<-qgam(scaledResids~s(t0,k=6), data=XH,qu=0.25)
S.50<-qgam(scaledResids~s(t0,k=6), data=XH,qu=0.5) 
S.75<-qgam(scaledResids~s(t0,k=6), data=XH,qu=0.75)
S.90<-qgam(scaledResids~s(t0,k=6), data=XH,qu=0.9) 
S.95<-qgam(scaledResids~s(t0,k=6), data=XH,qu=0.95)

## NP skewness
NPS_hat = Q.skewness(q.10=predict(S.10),
                     q.50=predict(S.50),
                     q.90=predict(S.90))

## NP kurtosis (relative to Gaussian)
NPK_hat = Q.kurtosis(q.05=predict(S.05),
                     q.25=predict(S.25),
                     q.75=predict(S.75),
                     q.95=predict(S.95))

## view quantile diagnostics of scaled residuals
pdf("../manuscript/figures/lichen_qgam_diagnostics.pdf",height = 5, width = 11,useDingbats = F)
par(mfrow=c(1,2),mar = c(5, 5, 2, 3), oma=c(0,0,0,2)) 
plot(XH$t0,XH$t1,pch=1,col=alpha("black",0.25),cex.axis=0.8,
     xlab="area, time t",ylab="area, time t+1")
lines(XH$t0,predict(fitGAU,type="response")[,1],col="red",lwd=2)
par(new = TRUE)                           
plot(XH$t0,XH$fitted_sd,col="blue",type="l",lwd=2,
     axes = FALSE, xlab = "", ylab = "",cex.axis=0.8)
axis(side = 4, at = pretty(range(XH$fitted_sd)),cex.axis=0.8)
mtext("std dev", side = 4, line = 2)
legend("topleft",legend=c("Fitted mean","Fitted sd"),bg="white",pch=1,col=c("red","blue"),cex=0.8)
title("A",font=3,adj=0)

plot(XH$t0,XH$scaledResids,col=alpha("black",0.25),cex.axis=0.8,
     xlab="area, time t",ylab="Scaled residuals of size at time t+1")
points(XH$t0,predict(S.05),col="black",pch=".")
points(XH$t0,predict(S.10),col="black",pch=".")
points(XH$t0,predict(S.25),col="black",pch=".")
points(XH$t0,predict(S.50),col="black",pch=".")
points(XH$t0,predict(S.75),col="black",pch=".")
points(XH$t0,predict(S.90),col="black",pch=".")
points(XH$t0,predict(S.95),col="black",pch=".")
par(new = TRUE)                           
matplot(cbind(XH$t0,XH$t0),cbind(NPS_hat,NPK_hat), type="l",lwd=2,
     col=c("blue","red"), lty=1, axes = FALSE, xlab = "", ylab = "",cex.axis=0.8)
abline(h=0,col="lightgray",lty=3)
axis(side = 4, at = pretty(range(c(NPS_hat,NPK_hat))),cex.axis=0.8)
mtext("Skewness and kurtosis", side = 4, line = 2)
legend("topleft",legend=c("Quantiles","NP skewness","NP excess kurtosis"),bg="white",pch=c(20,1,1),col=c("black","blue","red"),cex=0.8)
title("B",font=3,adj=0)
dev.off()

## collect outputs to write
lichen_out <- list(
  lichen_grow = XH[,c("t0","t1","fitted","fitted_sd","scaledResids")],
  lichen_GAU_best = fitGAU,
  q.05 = predict(S.05),
  q.10 = predict(S.10),
  q.25 = predict(S.25),
  q.50 = predict(S.50),
  q.75 = predict(S.75),
  q.90 = predict(S.90),
  q.95 = predict(S.95),
  NPS_hat = NPS_hat,
  NPK_hat = NPK_hat
)

###### Candidate for an improved model: gam SHASH
###### We begin with some model selection on skewness and kurtosis vs. size 
fitSHASH <- gam(list(t1 ~ s(t0), # <- location 
                        ~ s(t0),  # <- log-scale
                        ~ s(t0),  # <- skewness
                        ~ 1), # <- log-kurtosis
                      data = XH, gamma=1.4, family = shash,  optimizer = "efs")
SHASH_pred<-predict(fitSHASH,type="response")
# plot(fitSHASH,scale=FALSE); 
AIC(fitGAU,fitSHASH); 

fitSHASH2 <- gam(list(t1 ~ t0 + I(t0^2), # <- location 
                        ~ s(t0),  # <- log-scale
                        ~ s(t0),  # <- skewness
                        ~ 1), # <- log-kurtosis
                      data = XH, gamma=1.4, family = shash,  optimizer = "efs")
SHASH_pred2<-predict(fitSHASH,type="response")
# plot(fitSHASH2,scale=FALSE); 

fitSHASH3 <- gam(list(t1 ~ t0 + I(t0^2), # <- location 
                        ~ t0 + I(t0^2),  # <- log-scale
                        ~ t0+ I(t0^2),  # <- skewness
                        ~ 1), # <- log-kurtosis
                      data = XH, gamma=1.4, family = shash,  optimizer = "efs")
SHASH_pred3<-predict(fitSHASH3,type="response")

## add spline for log-kurtosis
fitSHASH4 <- gam(list(t1 ~ t0 + I(t0^2), # <- location 
                      ~ s(t0),  # <- log-scale
                      ~ s(t0),  # <- skewness
                      ~ s(t0)), # <- log-kurtosis
                 data = XH, gamma=1.4, family = shash,  optimizer = "efs")
SHASH_pred4<-predict(fitSHASH4,type="response")

## or linear effect of size
fitSHASH5 <- gam(list(t1 ~ t0 + I(t0^2), # <- location 
                      ~ s(t0),  # <- log-scale
                      ~ s(t0),  # <- skewness
                      ~ t0), # <- log-kurtosis
                 data = XH, gamma=1.4, family = shash,  optimizer = "efs")
SHASH_pred5<-predict(fitSHASH5,type="response")

fitSHASH6 <- gam(list(t1 ~ t0 + I(t0^2), # <- location 
                      ~ t0 + I(t0^2),  # <- log-scale
                      ~ s(t0),  # <- skewness
                      ~ 1), # <- log-kurtosis
                 data = XH, gamma=1.4, family = shash,  optimizer = "efs")
SHASH_pred6<-predict(fitSHASH6,type="response")

AIC(fitSHASH,fitSHASH2,fitSHASH3,fitSHASH4,fitSHASH5,fitSHASH6); # Number 2 is still the winner. 

## Save the fitted size-dependent parameters of the best SHASH model. 
muE <- fitSHASH2$fitted[ , 1]
sigE <- exp(fitSHASH2$fitted[ , 2])
epsE <- fitSHASH2$fitted[ , 3]
delE <- exp(fitSHASH2$fitted[ , 4])
t0 <- XH$t0; 
save(muE,sigE,epsE,delE,t0,file="SHASHfuns.Rdata"); 

##########################################################
# Compare quantiles of the two growth models 
##########################################################

par(mfrow=c(2,1))
Y = cbind(qSHASHo2(0.01,muE,sigE,epsE,delE), qSHASHo2(0.25,muE,sigE,epsE,delE), qSHASHo2(0.5,muE,sigE,epsE,delE), 
			qSHASHo2(0.75,muE,sigE,epsE,delE), qSHASHo2(0.99,muE,sigE,epsE,delE)) 
matplot(XH$t0,Y,type="l",lty=1,col="black",xlab="Area t", ylab = "Area t+1"); 
points(XH$t0,XH$t1); title(main="SHASH percentiles 1, 25, 50, 75, 99"); 

mE <- fitGAU$fitted[ , 1]
sE <- sqrt(1/(fitGAU$fitted[ , 2]))
Y = cbind(qnorm(0.01,mE,sE), qnorm(0.25,mE,sE), qnorm(0.5,mE,sE), 
			qnorm(0.75,mE,sE), qnorm(0.99,mE,sE)) 
matplot(XH$t0,Y,type="l",lty=1,col="black", xlab="Area t", ylab = "Area t+1"); 
points(XH$t0,XH$t1); title(main="Gaussian  percentiles 1, 25, 50, 75, 99"); 
dev.copy2pdf(file="Compare_quantiles.pdf")


################################################################# 
## Make the IPMs 
#################################################################

### move to the right local directory 
tom = "C:/Users/tm9/Dropbox/github/IPM_size_transitions"
steve = "c:/repos/IPM_size_transitions" 
home = ifelse(Sys.info()["user"] == "Ellner", steve, tom)
setwd(home); setwd("lichen"); 

library(gamlss.dist); load(file="SHASHfuns.Rdata")

########################################### 
## GAUSSIAN
###########################################
graphics.off(); par(mfrow=c(2,1)); 
L=-1; U=10; L1 = 0.2; U1 = 7; # limits of the data for eviction-prevention 
## see below for where these limits come from

#coef(fitGAU)
# (Intercept)            t0       I(t0^2) (Intercept).1          t0.1     I(t0^2).1 
#   0.07053852    1.03383973   -0.01763714   -2.28977280    0.74370780   -0.05767208 
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
plot(IPM_G$meshpts,apply(IPM_G$P,2,sum)/sx(chop_meshpts), type="l",lty=1,col=c("black","red")); 
abline(v=c(L1,U1));  

########################################### 
## SHASH 
###########################################
muE_fun <- splinefun(t0,muE, method="natural")
sigE_fun <- splinefun(t0,sigE, method="natural")
epsE_fun <- splinefun(t0,epsE, method="natural")
delE_fun <- splinefun(t0,delE, method="natural")

## how well to these splines reproduce the fits? -- very well
plot(t0,muE);lines(t0,muE_fun(t0),col="red")
plot(t0,sigE);lines(t0,sigE_fun(t0),col="red")
plot(t0,epsE);lines(t0,epsE_fun(t0),col="red")
plot(t0,delE);lines(t0,delE_fun(t0),col="red")

G_z1z_S = function(z1,z) {sx(z)*dSHASHo2(z1,muE_fun(z),sigE_fun(z),epsE_fun(z),delE_fun(z))}

mk_K_S <- function(m, L, U, L1, U1) {
	# mesh points 
	h <- (U - L)/m
	meshpts <- L + ((1:m) - 1/2) * h
	K <- P <- h * (outer(meshpts, pmax(pmin(meshpts, U1),L1),G_z1z_S))
	K[1,] = K[1,] + matrix(fx(meshpts),nrow=1); F = K - P; 
	return(list(K = K, meshpts = meshpts, P = P, F = F))
}

IPM_S =mk_K_S(200,L=L,U=U,L1 = L1, U1 = U1); Re(eigen(IPM_S$K)$values[1]);
chop_meshpts = pmax(pmin(IPM_S$meshpts,U1),L1)
plot(IPM_S$meshpts,apply(IPM_S$P,2,sum)/sx(chop_meshpts), type="l",lty=1,col=c("black","red")); 
abline(v=c(L1,U1)); 

source("../code/matrixImage.R"); 
graphics.off(); 
dev.new(); matrix.image((IPM_G$K)^0.25, IPM_G$meshpts, IPM_G$meshpts, main="Gaussian"); 
dev.new(); matrix.image((IPM_S$K)^0.25, IPM_S$meshpts, IPM_S$meshpts,main ="SHASH"); 

################################################################# 
## Compare life history attributes 
#################################################################

source("../code/metaluck_fns_CMH.R"); 

X = matrix(NA, 3,12); 

### Gaussian
matU = IPM_G$P; matF = IPM_G$F; c0 = rep(0,nrow(matU)); c0[1]=1; 
X[1,] = c(mean_lifespan(matU, mixdist=c0),
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

X[2,] = c(
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

X = data.frame(round(X,digits=2)); 
names(X) = c("Mean lifespan", "SD lifespan", "Skew lifespan", "Mean LRO", "SD LRO", "Skew LRO", "Prob repro", "Mean age repro", "Conditional lifespan", 
"Gen time Ta", "Gen time mu1(v)", "Gen time R0"); 
row.names(X) = c("Gaussian", "SHASH"); 

View(X); 

##########################################################################
##  Simulate discrete population dynamics and compare extinction risk
##########################################################################
update_pop = function(sizes,SHASH=TRUE) {
  N = length(sizes); 
  if(N==0){
	return(sizes); 
  }else{
	new_sizes = numeric(0); 
	total_F = sum(fx(sizes)); 
	kids = rpois(1,total_F);
	if(kids>0) new_sizes = rep(L1,kids)

	survival_probs = sx(sizes); 
	live = which(runif(N)<survival_probs) 
	if(length(live)>0) {
		living = sizes[live];
		if(SHASH) {
			surv_sizes = rSHASHo2(length(living), muE_fun(living),sigE_fun(living),epsE_fun(living),delE_fun(living) ) 
		}else{
			surv_sizes = rnorm(length(living), mu_G(living),sd_G(living))
        }			
		surv_sizes = pmax(pmin(surv_sizes, U1),L1)
		new_sizes = c(new_sizes,surv_sizes)
	}
	return(new_sizes) 
  }
}	
	
	
popsize_S = popsize_G =matrix(NA,100,5000);
for(k in 1:5000) {
sizes_S = sizes_G = runif(12,L1,2*L1); 
for(j in 1:100) {
		sizes_G=update_pop(sizes_G,SHASH=FALSE);
		popsize_G[j,k]=length(sizes_G);
		sizes_S=update_pop(sizes_S,SHASH=TRUE);
		popsize_S[j,k]=length(sizes_S);
}		
if(k%%100==0) cat(k,"\n"); 
}
extinct_SHASH = apply(popsize_S,1,function(x) sum(x==0)); 
extinct_GAU =  apply(popsize_G,1,function(x) sum(x==0)); 

dev.new(width=8,height=6); par(bty="l",cex.axis=1.3,cex.lab=1.3,mgp=c(2.1,1,0)); 
matplot(1:100,cbind(extinct_SHASH,extinct_GAU)/5000,col=c("black","red"),type="l",lty=1,
	xlab="Years", ylab="Extinction probability",lwd=2); 
legend("topleft",legend = c("SHASH", "Gaussian"), col=c("black","red"),lty=1,lwd=2,inset=0.03)
dev.copy2pdf(file="exctinction_risk.pdf"); 

###################################################################
## Exploring the JSU as an alternative
###################################################################
## since kurtosis is constant, a JSU might be a more natural alternative to Gaussian. 
## We could hold over the mean and sd from pilot gaussian because in this parameterization mu and sigma are mean and sd. 
## However, we want to compare with SHASH, and we cannot make the same transfer with SHASH, 
## therefore re-fitting polynomials for mean and sd. 
LogLikJSU=function(pars,response){
  dJSU(response, 
       mu=pars[1] + pars[2]*XH$t0 + pars[3]*(XH$t0^2),
       sigma=exp(pars[4] + pars[5]*XH$t0 + pars[6]*(XH$t0^2)),
       nu = pars[7] + pars[8]*XH$t0,
       tau = exp(pars[9]), log=TRUE)
}
## compare to SHASH
LogLikSHASH=function(pars,response){
  dSHASHo2(response, 
           mu=pars[1]+pars[2]*XH$t0+pars[3]*(XH$t0^2),
           sigma=exp(pars[4]+pars[5]*XH$t0+pars[6]*(XH$t0^2)),
           nu = pars[7]+pars[8]*XH$t0,
           tau = exp(pars[9]), log=TRUE)
}

## starting parameters
p0<-c(coef(fitGAU)[1:6],0,0,0)
## fit with maxlik
outJSU=maxLik(logLik=LogLikJSU,start=p0*exp(0.2*rnorm(length(p0))), response=XH$t1,
              method="BHHH",control=list(iterlim=5000,printLevel=2),finalHessian=FALSE); 
outJSU=maxLik(logLik=LogLikJSU,start=outJSU$estimate,response=XH$t1,
              method="NM",control=list(iterlim=5000,printLevel=1),finalHessian=FALSE); 
outJSU=maxLik(logLik=LogLikJSU,start=outJSU$estimate,response=XH$t1,
              method="BHHH",control=list(iterlim=5000,printLevel=2),finalHessian=FALSE); 

outSHASH=maxLik(logLik=LogLikSHASH,start=p0*exp(0.2*rnorm(length(p0))), response=XH$t1,
                method="BHHH",control=list(iterlim=5000,printLevel=2),finalHessian=FALSE); 
outSHASH=maxLik(logLik=LogLikSHASH,start=outSHASH$estimate,response=XH$t1,
                method="NM",control=list(iterlim=5000,printLevel=1),finalHessian=FALSE); 
outSHASH=maxLik(logLik=LogLikSHASH,start=outSHASH$estimate,response=XH$t1,
                method="BHHH",control=list(iterlim=5000,printLevel=2),finalHessian=FALSE); 

AIC(outJSU);AIC(outSHASH) ##looks like JSU is the winner
AIC(fitSHASH2) ##it's also better than the best gam SHASH

##########################################################
# Compare quantiles
##########################################################
muE_SHASH <- outSHASH$estimate[1]+outSHASH$estimate[2]*XH$t0+outSHASH$estimate[3]*(XH$t0^2)
sigE_SHASH <- exp(outSHASH$estimate[4]+outSHASH$estimate[5]*XH$t0+outSHASH$estimate[6]*(XH$t0^2))
epsE_SHASH <- outSHASH$estimate[7]+outSHASH$estimate[8]*XH$t0
delE_SHASH <- exp(outSHASH$estimate[9])
Y_SHASH = cbind(qSHASHo2(0.01,muE_SHASH,sigE_SHASH,epsE_SHASH,delE_SHASH), qSHASHo2(0.25,muE_SHASH,sigE_SHASH,epsE_SHASH,delE_SHASH), qSHASHo2(0.5,muE_SHASH,sigE_SHASH,epsE_SHASH,delE_SHASH), 
                qSHASHo2(0.75,muE_SHASH,sigE_SHASH,epsE_SHASH,delE_SHASH), qSHASHo2(0.99,muE_SHASH,sigE_SHASH,epsE_SHASH,delE_SHASH)) 

muE_JSU <- outJSU$estimate[1]+outJSU$estimate[2]*XH$t0+outJSU$estimate[3]*(XH$t0^2)
sigE_JSU <- exp(outJSU$estimate[4]+outJSU$estimate[5]*XH$t0+outJSU$estimate[6]*(XH$t0^2))
epsE_JSU <- outJSU$estimate[7]+outJSU$estimate[8]*XH$t0
delE_JSU <- exp(outJSU$estimate[9])
Y_JSU = cbind(qJSU(0.01,muE_JSU,sigE_JSU,epsE_JSU,delE_JSU), qJSU(0.25,muE_JSU,sigE_JSU,epsE_JSU,delE_JSU), qJSU(0.5,muE_JSU,sigE_JSU,epsE_JSU,delE_JSU), 
              qJSU(0.75,muE_JSU,sigE_JSU,epsE_JSU,delE_JSU), qJSU(0.99,muE_JSU,sigE_JSU,epsE_JSU,delE_JSU)) 


mE <- fitGAU$fitted[ , 1]
sE <- sqrt(1/(fitGAU$fitted[ , 2]))
Y_GAU = cbind(qnorm(0.01,mE,sE), qnorm(0.25,mE,sE), qnorm(0.5,mE,sE), 
          qnorm(0.75,mE,sE), qnorm(0.99,mE,sE)) 


par(mfrow=c(2,1));  
matplot(XH$t0,Y_SHASH,type="l",lty=1,col="black",xlab="Area t", ylab = "Area t+1"); 
matplot(XH$t0,Y_JSU,type="l",lty=1,col="red",xlab="Area t", ylab = "Area t+1",add=T); 
points(XH$t0,XH$t1); title(main="SHASH percentiles 1, 25, 50, 75, 99"); 

matplot(XH$t0,Y_GAU,type="l",lty=1,col="black", xlab="Area t", ylab = "Area t+1"); 
points(XH$t0,XH$t1); title(main="Gaussian  percentiles 1, 25, 50, 75, 99")


## just out of curiosity, what would happen if we fix JSU mean and sd at the Gaussian fitted mean and sd?
LogLikJSU2=function(pars,response){
  dJSU(response, 
       mu=XH$fitted,
       sigma=XH$fitted_sd,
       nu = pars[1]+pars[2]*XH$t0,
       tau = exp(pars[3]), log=TRUE)
}
p0<-c(0,0,0)
## fit with maxlik
outJSU2=maxLik(logLik=LogLikJSU2,start=p0*exp(0.2*rnorm(length(p0))), response=XH$t1,
              method="BHHH",control=list(iterlim=5000,printLevel=2),finalHessian=FALSE) 
outJSU2=maxLik(logLik=LogLikJSU2,start=outJSU2$estimate,response=XH$t1,
              method="NM",control=list(iterlim=5000,printLevel=1),finalHessian=FALSE) 
outJSU2=maxLik(logLik=LogLikJSU2,start=outJSU2$estimate,response=XH$t1,
              method="BHHH",control=list(iterlim=5000,printLevel=2),finalHessian=FALSE) 

## compare nu and tau params
outJSU$estimate[7:9];outJSU2$estimate ##different!
muE_JSU2 <- XH$fitted
sigE_JSU2 <- XH$fitted_sd
epsE_JSU2 <- outJSU2$estimate[1]+outJSU2$estimate[2]*XH$t0
delE_JSU2 <- exp(outJSU2$estimate[3])
Y_JSU2 = cbind(qJSU(0.01,muE_JSU2,sigE_JSU2,epsE_JSU2,delE_JSU2), qJSU(0.25,muE_JSU2,sigE_JSU2,epsE_JSU2,delE_JSU2), qJSU(0.5,muE_JSU2,sigE_JSU2,epsE_JSU2,delE_JSU2), 
              qJSU(0.75,muE_JSU2,sigE_JSU2,epsE_JSU2,delE_JSU2), qJSU(0.99,muE_JSU2,sigE_JSU2,epsE_JSU2,delE_JSU2)) 

matplot(XH$t0,Y_JSU,type="l",lty=1,col="black",xlab="Area t", ylab = "Area t+1"); 
matplot(XH$t0,Y_JSU2,type="l",lty=1,col="red",xlab="Area t", ylab = "Area t+1",add=T); 
points(XH$t0,XH$t1) 

##########################################################
# Compare real and simulated data
##########################################################
n_sim<-100
GAUsim_mean<-GAUsim_sd<-GAUsim_skew<-GAUsim_kurt<-matrix(NA,nrow=nrow(XH),ncol=n_sim)
JSUsim_mean<-JSUsim_sd<-JSUsim_skew<-JSUsim_kurt<-matrix(NA,nrow=nrow(XH),ncol=n_sim)
for(i in 1:n_sim){
  ## add this iteration of sim data to real df
  XH$t1.simGAU <- rnorm(n=nrow(XH),mean=XH$fitted,sd=XH$fitted_sd)
  XH$t1.simJSU <- rJSU(n=nrow(XH),mu=muE_JSU,sigma=sigE_JSU,nu=epsE_JSU,tau=delE_JSU)
  ## Qreg on simGAU data
  q.05<-predict(qgam(t1.simGAU~s(t0,k=4),data=XH,qu=0.05)) 
  q.10<-predict(qgam(t1.simGAU~s(t0,k=4),data=XH,qu=0.10)) 
  q.25<-predict(qgam(t1.simGAU~s(t0,k=4),data=XH,qu=0.25))
  q.50<-predict(qgam(t1.simGAU~s(t0,k=4),data=XH,qu=0.5))
  q.75<-predict(qgam(t1.simGAU~s(t0,k=4),data=XH,qu=0.75)) 
  q.90<-predict(qgam(t1.simGAU~s(t0,k=4),data=XH,qu=0.90))
  q.95<-predict(qgam(t1.simGAU~s(t0,k=4),data=XH,qu=0.95))
  GAUsim_mean[,i]<-Q.mean(q.25,q.50,q.75)
  GAUsim_sd[,i]<-Q.sd(q.25,q.75)
  GAUsim_skew[,i]<-Q.skewness(q.10,q.50,q.90)
  GAUsim_kurt[,i]<-Q.kurtosis(q.05,q.25,q.75,q.95)
  
  ## same on simJSU data (writing over the q objects)
  q.05<-predict(qgam(t1.simJSU~s(t0,k=4),data=XH,qu=0.05)) 
  q.10<-predict(qgam(t1.simJSU~s(t0,k=4),data=XH,qu=0.10)) 
  q.25<-predict(qgam(t1.simJSU~s(t0,k=4),data=XH,qu=0.25))
  q.50<-predict(qgam(t1.simJSU~s(t0,k=4),data=XH,qu=0.5))
  q.75<-predict(qgam(t1.simJSU~s(t0,k=4),data=XH,qu=0.75)) 
  q.90<-predict(qgam(t1.simJSU~s(t0,k=4),data=XH,qu=0.90))
  q.95<-predict(qgam(t1.simJSU~s(t0,k=4),data=XH,qu=0.95))
  JSUsim_mean[,i]<-Q.mean(q.25,q.50,q.75)
  JSUsim_sd[,i]<-Q.sd(q.25,q.75)
  JSUsim_skew[,i]<-Q.skewness(q.10,q.50,q.90)
  JSUsim_kurt[,i]<-Q.kurtosis(q.05,q.25,q.75,q.95)
}

## and now the real data
q.05<-predict(qgam(t1~s(t0,k=4),data=XH,qu=0.05)) 
q.10<-predict(qgam(t1~s(t0,k=4),data=XH,qu=0.10)) 
q.25<-predict(qgam(t1~s(t0,k=4),data=XH,qu=0.25)) 
q.50<-predict(qgam(t1~s(t0,k=4),data=XH,qu=0.5))
q.75<-predict(qgam(t1~s(t0,k=4),data=XH,qu=0.75))
q.90<-predict(qgam(t1~s(t0,k=4),data=XH,qu=0.90))
q.95<-predict(qgam(t1~s(t0,k=4),data=XH,qu=0.95))

pdf("../manuscript/figures/lichen_JSU_fit.pdf",height = 6, width = 6,useDingbats = F)
par(mfrow=c(2,2),mar=c(4,4,1,1),mgp=c(2.1,1,0),cex.lab=1.2); 
plot(XH$t0,Q.mean(q.25,q.50,q.75),type="n",
     xlab="Size(t)",ylab="NP mean size(t+1)",
     ylim=c(min(c(GAUsim_mean,JSUsim_mean)),1 + max(c(GAUsim_mean,JSUsim_mean))))
matpoints(XH$t0,GAUsim_mean,col=alpha("tomato",0.25),type="l",lty=1)
matpoints(XH$t0,1 + JSUsim_mean,col=alpha("cornflowerblue",0.25),type="l",lty=1)
points(XH$t0,Q.mean(q.25,q.50,q.75),col="black",type="l",lwd=2)
points(XH$t0,1 + Q.mean(q.25,q.50,q.75),col="black",type="l",lwd=2)
legend("topleft",legend=c("Real data","Gaussian simulation","JSU simulation + offset"),
       lty=1,col=c("black","tomato","cornflowerblue"),cex=0.8,bty="n")

plot(XH$t0,Q.sd(q.25,q.75),type="n",
     xlab="Size(t)",ylab="NP SD Size(t+1)",
     ylim=c(min(c(GAUsim_sd,JSUsim_sd)),1 + max(c(GAUsim_sd,JSUsim_sd))))
matpoints(XH$t0,GAUsim_sd,col=alpha("tomato",0.25),type="l",lty=1)
matpoints(XH$t0,1 + JSUsim_sd,col=alpha("cornflowerblue",0.25),type="l",lty=1)
points(XH$t0,Q.sd(q.25,q.75),col="black",type="l",lwd=2)
points(XH$t0,1 + Q.sd(q.25,q.75),col="black",type="l",lwd=2)

plot(XH$t0,Q.skewness(q.10,q.50,q.90),type="n",
     xlab="Size(t)",ylab="NP skewness size(t+1)",
     ylim=c(min(c(GAUsim_skew,JSUsim_skew)),1 + max(c(GAUsim_skew,JSUsim_skew))))
matpoints(XH$t0,GAUsim_skew,col=alpha("tomato",0.25),type="l",lty=1)
matpoints(XH$t0,1+ JSUsim_skew,col=alpha("cornflowerblue",0.25),type="l",lty=1)
points(XH$t0,Q.skewness(q.10,q.50,q.90),col="black",type="l",lwd=2)
points(XH$t0,1 + Q.skewness(q.10,q.50,q.90),col="black",type="l",lwd=2)

plot(XH$t0,Q.kurtosis(q.05,q.25,q.75,q.95),type="n",
     xlab="Size(t)",ylab="NP kurtosis size(t+1)",
     ylim=c(min(GAUsim_kurt,JSUsim_kurt),1 + max(GAUsim_kurt,JSUsim_kurt)))
matpoints(XH$t0,GAUsim_kurt,col=alpha("tomato",0.25),type="l",lty=1)
matpoints(XH$t0,1 + JSUsim_kurt,col=alpha("cornflowerblue",0.25),type="l",lty=1)
points(XH$t0,Q.kurtosis(q.05,q.25,q.75,q.95),col="black",type="l",lwd=2)
points(XH$t0,1 + Q.kurtosis(q.05,q.25,q.75,q.95),col="black",type="l",lwd=2)
dev.off()

########################################### 
## JSU IPM
###########################################
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

IPM_J =mk_K_J(200,L=L,U=U,L1 = L1, U1 = U1); Re(eigen(IPM_J$K)$values[1]);
chop_meshpts = pmax(pmin(IPM_J$meshpts,U1),L1)
plot(IPM_J$meshpts,apply(IPM_J$P,2,sum)/sx(chop_meshpts), type="l",lty=1,col=c("black","red")); 
abline(v=c(L1,U1)); 

### JSU
matU = IPM_J$P; matF = IPM_J$F; c0 = rep(0,nrow(matU)); c0[1]=1; 

X[3,] = c(
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

X = data.frame(round(X,digits=2)); 
names(X) = c("Mean lifespan", "SD lifespan", "Skew lifespan", "Mean LRO", "SD LRO", "Skew LRO", "Prob repro", "Mean age repro", "Conditional lifespan", 
             "Gen time Ta", "Gen time mu1(v)", "Gen time R0"); 
row.names(X) = c("Gaussian", "SHASH", "JSU"); 

# proportional differences
(X$`Mean lifespan`[1]-X$`Mean lifespan`[3]) / X$`Mean lifespan`[1]
(X$`Mean LRO`[1]-X$`Mean LRO`[3]) / X$`Mean LRO`[1]
(X$`Mean age repro`[1]-X$`Mean age repro`[3]) / X$`Mean age repro`[1]

##########################################################################
##  Simulate discrete population dynamics and compare extinction risk
##########################################################################
update_pop = function(sizes,JSU=TRUE) {
  N = length(sizes); 
  if(N==0){
    return(sizes); 
  }else{
    new_sizes = numeric(0); 
    total_F = sum(fx(sizes)); 
    kids = rpois(1,total_F);
    if(kids>0) new_sizes = rep(L1,kids)
    
    survival_probs = sx(sizes); 
    live = which(runif(N)<survival_probs) 
    if(length(live)>0) {
      living = sizes[live];
      if(JSU) {
        surv_sizes = rJSU(length(living), muJSU_fun(living),sigJSU_fun(living),epsJSU_fun(living),delJSU_fun(living) ) 
      }else{
        surv_sizes = rnorm(length(living), mu_G(living),sd_G(living))
      }			
      surv_sizes = pmax(pmin(surv_sizes, U1),L1)
      new_sizes = c(new_sizes,surv_sizes)
    }
    return(new_sizes) 
  }
}	

popsize_J = popsize_G =matrix(NA,100,5000);
for(k in 1:5000) {
  sizes_J = sizes_G = runif(12,L1,2*L1); 
  for(j in 1:100) {
    sizes_G=update_pop(sizes_G,JSU=FALSE);
    popsize_G[j,k]=length(sizes_G);
    sizes_J=update_pop(sizes_J,JSU=TRUE);
    popsize_J[j,k]=length(sizes_J);
  }		
  if(k%%100==0) cat(k,"\n"); 
}
extinct_JSU = apply(popsize_J,1,function(x) sum(x==0)); 
extinct_GAU =  apply(popsize_G,1,function(x) sum(x==0)); 

dev.new(width=10,height=10); par(bty="l",cex.axis=1.3,cex.lab=1.3,mgp=c(2.1,1,0)); 
matplot(1:100,cbind(extinct_JSU,extinct_GAU)/5000,col=c("black","red"),type="l",lty=1,
        xlab="Years", ylab="Extinction probability",lwd=2); 
legend("topleft",legend = c("JSU", "Gaussian"), col=c("black","red"),lty=1,lwd=2,inset=0.03)
dev.copy2pdf(file="../manuscript/figures/lichen_extinction_risk.pdf"); 

## add GAU and JSU kernels to output and write to file
lichen_out$mat_GAU <- IPM_G
lichen_out$mat_JSU <- IPM_J
write_rds(lichen_out,"lichen_out.rds")


##### find eviction limits that work for GAU and JSU
par(mfrow=c(1,2))
IPM_G = mk_K_G(200,L=-1,U=10,L1 = L1, U1 = U1); Re(eigen(IPM_G$K)$values[1]);
chop_meshpts = pmax(pmin(IPM_G$meshpts,U1),L1) 
plot(IPM_G$meshpts,apply(IPM_G$P,2,sum)/sx(chop_meshpts), type="l",lty=1,col=c("black","red")); 
abline(v=c(L1,U1))

IPM_J =mk_K_J(200,L=-1,U=10,L1 = L1, U1 = U1); Re(eigen(IPM_J$K)$values[1]);
chop_meshpts = pmax(pmin(IPM_J$meshpts,U1),L1)
plot(IPM_J$meshpts,apply(IPM_J$P,2,sum)/sx(chop_meshpts), type="l",lty=1,col=c("black","red")); 
abline(v=c(L1,U1))
