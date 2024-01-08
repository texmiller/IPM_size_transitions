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
## AIC goes for fit2, by a hair. 
#(Intercept)    sqrt(t0)          t0 
#  -1.649455    5.363115   -1.249375 

sx = function(x)  {
	u1 = -1.649455  + 5.363115*sqrt(x) - 1.249375*x 
	u2 = 0.0954 + 2.252*sqrt(x); 
	p1 = exp(u1)/(1+exp(u1)); 	p2 = exp(u1)/(1+exp(u1)); ##<--what's going on here?
	return( 0.2*p1 + 0.8*p2 )
}	

nbins = 12; 
indexx = seq_along(XH$t0)
u = nbins*indexx/(1 + max(indexx)); 
u = floor(u); pbar = tbar = numeric(nbins); 
for(j in 1:nbins) {pbar[j] = mean(XH$survival[u==j]); tbar[j]=mean(XH$t0[u==j])}
plot(tbar,pbar,pch=16,xlim=range(XH$t0),ylim=c(0.9,1)); 
points(XH$t0,sx(XH$t0),type="l",lty=1); 

############### Fecundity function (from Shriver et al. model) 
fx = function(Area) {
		r = sqrt(Area/pi);
		u = 2*pi*r; 
		return(0.047*u)
}		

########################################################################
#  Growth modeling 
########################################################################

XH = XH[XH$survival==1,]; ## discard the dead 

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

mean(XH$scaledResids); sd(XH$scaledResids); ## all good 

### Parametric SD function is nearly the same as the spline 
plot(XH$t0, log(XH$fitted_sd)); 
points(XH$t0, log( 1/predict(fitGAU,type="response")[,2]), col="red" )

####################### Save the best-fitting Gaussian 
fitGAU = fitGAU22; rm(fitGAU22); rm(fitGAU0); rm(fitGAU00);  

######### Diagnostics on fitted parametric SD function: no problems! 
#stopCluster(c1); ##<--c1 not defined
c1<- makeCluster(8); 
registerDoParallel(c1);
out = multiple_levene_test(XH$fitted, XH$scaledResids, 3, 8, 2000);
out$p_value; ## > 0.89

out = multiple_bs_test(XH$fitted, XH$scaledResids, 4, 8, 2000) 
out$p_value; ## > 0.78; 
stopCluster(c1); 

## quantile regressions on stand resids
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
pdf("lichen_qgam_diagnostics.pdf",height = 5, width = 11,useDingbats = F)
par(mfrow=c(1,2),mar = c(5, 5, 2, 3), oma=c(0,0,0,2)) 
plot(XH$t0,XH$t1,pch=1,col=alpha("black",0.25),cex.axis=0.8,
     xlab="log area, time t",ylab="log area, time t+1")
points(XH$t0,predict(fitGAU,type="response")[,1],col=alpha("red",0.25),pch=16,cex=.5)
par(new = TRUE)                           
plot(XH$t0,XH$fitted_sd,col=alpha("blue",0.25),pch=16,cex=.5,
     axes = FALSE, xlab = "", ylab = "",cex.axis=0.8)
axis(side = 4, at = pretty(range(XH$fitted_sd)),cex.axis=0.8)
mtext("std dev", side = 4, line = 2)
legend("topleft",legend=c("Fitted mean","Fitted sd"),bg="white",pch=1,col=c("red","blue"),cex=0.8)
title("A",font=3,adj=0)

plot(XH$t0,XH$scaledResids,col=alpha("black",0.25),cex.axis=0.8,
     xlab="log area, time t",ylab="Scaled residuals of size at time t+1")
points(XH$t0,predict(S.05),col="black",pch=".")
points(XH$t0,predict(S.10),col="black",pch=".")
points(XH$t0,predict(S.25),col="black",pch=".")
points(XH$t0,predict(S.50),col="black",pch=".")
points(XH$t0,predict(S.75),col="black",pch=".")
points(XH$t0,predict(S.90),col="black",pch=".")
points(XH$t0,predict(S.95),col="black",pch=".")
par(new = TRUE)                           
matplot(cbind(XH$t0,XH$t0),cbind(NPS_hat,NPK_hat), type="l",
     col=c("blue","red"), lty=1, axes = FALSE, xlab = "", ylab = "",cex.axis=0.8)
abline(h=0,col="lightgray",lty=3)
axis(side = 4, at = pretty(range(c(NPS_hat,NPK_hat))),cex.axis=0.8)
mtext("Skewness and kurtosis", side = 4, line = 2)
legend("topleft",legend=c("Quantiles","NP skewness","NP excess kurtosis"),bg="white",pch=c(20,1,1),col=c("black","blue","red"),cex=0.8)
title("B",font=3,adj=0)
dev.off()

###### improved model: gam SHASH
## TM: why assume kurtosis is constant?
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


AIC(fitSHASH,fitSHASH2,fitSHASH3,fitSHASH4,fitSHASH5); # Number 2 is still the winner. 


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

dev.new(); par(mfrow=c(2,1));  
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
L=0; U=10; L1 = 0.2; U1 = 7; # limits of the data for eviction-prevention 

#coef(fitGAU)
# (Intercept)            t0       I(t0^2) (Intercept).1          t0.1     I(t0^2).1 
#   0.07053852    1.03383973   -0.01763714   -2.28977280    0.74370780   -0.05767208 
mu_G = function(x) {0.07053852  +  1.03383973*x  - 0.01763714*x^2} 
sd_G = function(x) {exp(-2.28977280   + 0.74370780*x - 0.05767208*x^2)}  
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

X = matrix(NA, 2,12); 

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


