rm(list=ls(all=TRUE)); 

setwd("c:/repos/IPM_size_transitions/coral"); 
require(zoo); require(moments); require(mgcv); 
require(gamlss); require(tidyverse); require(maxLik); 

source("../Diagnostics.R");
source("../fitChosenDists.R"); 

#############   load the data frame on healthy corals 
source("AkumalCoralsSetup.R"); 
# names(XH); 
#  "Site"  "Fan.number" "Year"  "Area1"  "Area2"  "logarea.t0" "logarea.t1"

######################################################################### 
# Fit the pilot Gaussian model. 
#########################################################################
fitGAU <- gam(list(logarea.t1~s(logarea.t0),~s(logarea.t0)), data=XH, gamma=1.4,family=gaulss())
summary(fitGAU); 

####  Extract values of the fitted splines to get polynomial approximations
z_vals = XH$logarea.t0; 
fitted_vals = predict(fitGAU,type="response"); 
##### Mean is fitted almost exactly by linear, exactly by quadratic (spline has df just below 3). 
mean_fit1 = lm(fitted_vals[,1]~z_vals); 
mean_fit2 = lm(fitted_vals[,1]~z_vals+I(z_vals^2)); 

#### log(sigma) is fitted well by linear, perfectly by quadratic (spline has df just above 2) 
sigma_hat = 1/fitted_vals[,2]; 
sd_fit1 = lm(log(sigma_hat)~z_vals); # R^2 = 0.97 
sd_fit2 = lm(log(sigma_hat)~z_vals+I(z_vals^2)); # R^2 = 0.999 

##########################################################################
# Estimate the parameters of SEP1 for a set of initial size bins.
##########################################################################
XH <- XH %>% mutate(size_bin = cut_number(logarea.t0,n=8))
bins = levels(XH$size_bin); 
mus = sigmas = nus = taus = numeric(length(bins)-1); 
for(j in 1:length(mus)){
	Xj=subset(XH,size_bin%in%c(bins[j-1],bins[j]))
	fitj = gamlssMaxlik(Xj$logarea.t1, DIST="SEP1") 
	pars = fitj$estimate; 
	mus[j]=pars[1]; sigmas[j]=pars[2]; nus[j]=pars[3]; taus[j]=pars[4]; 
}	

######################################################################
# Fit the SEP1 growth model, using those as start values 
######################################################################
KernelLogLik=function(pars,y,x){
	mu = pars[1]+ pars[2]*x + pars[3]*x^2  
	sigma = exp(pars[4] + pars[5]*x)
	nu = pars[6] + pars[7]*x
	tau = exp(pars[8]+pars[9]*x)
	val = dSEP1(y, mu=mu,sigma=sigma,nu=nu,tau=tau,log=TRUE)
	return(val); 
}

## The start values should be on the link-transformed scale (e.g., log sigma)
## because the inverse-link is applied in KernelLoglik(), so as to 
## guarantee that the parameters are valid in the distribution. 
start=numeric(9);
start[1:3]=coef(mean_fit2); # this is the pilot fit to the mean (link=identity) 
start[4:5]=coef(sd_fit1); # pilot fit to log (sigma)
start[6:7]=c(median(nus),0); # from binned data diagnostic
start[8:9]=c(median(taus),0) # from binned data diagnostic 

fit = maxLik(logLik=KernelLogLik,start=start*exp(0.1*rnorm(8)), y=XH$logarea.t1, x=XH$logarea.t0,method="BHHH",
		control=list(iterlim=5000,printLevel=2),finalHessian=FALSE);  
for(k in 1:5) {		
fit = maxLik(logLik=KernelLogLik,start=fit$estimate, y=XH$logarea.t1, x=XH$logarea.t0,method="NM",
		control=list(iterlim=25000,printLevel=2),finalHessian=FALSE);  
		
fit = maxLik(logLik=KernelLogLik,start=fit$estimate, y=XH$logarea.t1, x=XH$logarea.t0,method="BHHH",
		control=list(iterlim=5000,printLevel=2),finalHessian=FALSE);
}
fit = maxLik(logLik=KernelLogLik,start=fit$estimate, y=XH$logarea.t1, x=XH$logarea.t0,method="BHHH",
		control=list(iterlim=5000,printLevel=2),finalHessian=TRUE);
fit = maxLik(logLik=KernelLogLik,start=fit$estimate, y=XH$logarea.t1, x=XH$logarea.t0,method="NM",
		control=list(iterlim=25000,printLevel=2),finalHessian=TRUE);	
fit = maxLik(logLik=KernelLogLik,start=fit$estimate, y=XH$logarea.t1, x=XH$logarea.t0,method="BHHH",
		control=list(iterlim=5000,printLevel=2),finalHessian=TRUE);	

growPars = fit$estimate; 
save.image(file="AkumaCoralsIPMs.Rdata"); 

#################################################################
# Alternatively, fit SHASH gam. 
##################################################################
fitSHASH <- gam(list(logarea.t1 ~ s(logarea.t0,k=4), # <- location 
                     ~ s(logarea.t0,k=4),   # <- log-scale
                     ~ s(logarea.t0,k=4),   # <- skewness
                     ~ s(logarea.t0,k=4)), # <- log-kurtosis
                data = XH, 
                family = shash,  
                optimizer = "efs")

#################################################################
# Start here if the model is already fitted. 
##################################################################
setwd("c:/repos/IPM_size_transitions/coral"); 
require(mgcv); require(gamlss); 
load("AkumaCoralsIPMs.Rdata"); 


######################################################################
# Demographic functions from Bruno et al. (2011) 
######################################################################

fan.height=function(area) {1.3*sqrt(area)}
fan.area=function(height) {height*height/1.69}

### Additional mortality of large fans, based on 1997 Keys data 
bigmu=function(x) {
	p=c(2.3,8.7,0.093); 
	u=exp(p[1]*(x-p[2]))
	return(p[3]*u/(1+u))
}

### Mortality of healthy fans, based on 1997 Keys data 
mH=function(x) {u= -0.55 - 0.62*x; exp(u)/(1+exp(u)) + bigmu(x) }  

### Size distribution of new recruits 
recruitDensity=density(log(recruitSizes),bw="SJ"); 
recruitFun=approxfun(recruitDensity$x,recruitDensity$y,yleft=0,yright=0);

######################################################################
# Construct the IPM, using code from the Springer (2016) book 
######################################################################

####### SEP1 growth function 
G_z1z <- function(z1,z,pars=growPars){
	mu = pars[1]+ pars[2]*z + pars[3]*z^2  
	sigma = exp(pars[4] + pars[5]*z)
	nu = pars[6] + pars[7]*z
	tau = exp(pars[8]+pars[9]*z)
	return(dSEP1(z1, mu=mu,sigma=sigma,nu=nu,tau=tau))
}    

################### Gaussian pilot model growth function 
z_vals = XH$logarea.t0; 
fitted_vals = predict(fitGAU,type="response"); 
e = order(z_vals); 
meanFun = splinefun(z_vals[e],fitted_vals[e,1],method="natural"); 
sigmaFun = splinefun(z_vals[e],1/fitted_vals[e,2],method="natural"); 


####### growth functions 
#G_z1zPilot <- function(z1,z,pars=NULL){
#	mu = meanFun(z)
#	sigma = sigmaFun(z)
#	return(dnorm(z1, mean=mu,sd=sigma))
#}    

G_z1zPilot <- function(z1,z,pars=NULL){
  pred = predict(fitGAU,
                 newdata = data.frame(logarea.t0=z))
  return(dnorm(z1,mean=pred[,1],sd=exp(pred[,2])))
}  

G_z1zSHASH <- function(z1,z,pars=NULL){
  pred = predict(fitSHASH,
                 newdata = data.frame(logarea.t0=z))
  return(dSHASHo2(x=z1, 
                  mu=pred[,1],
                  sigma = exp(pred[,2]), 
                  nu = pred[,3], 
                  tau = exp(pred[,4])))
} 

## Survival function, logistic regression
s_z <- function(z) {
    x=exp(z)^(1/3); 
    return(1 - mH(x))
}    

P_z1zSHASH = function(z1,z) {s_z(z) * G_z1zSHASH(z1, z)}
P_z1zPilot = function(z1,z) {s_z(z) * G_z1zPilot(z1, z)}

## larval production function, using Bruno et al. model 
b_z <- function(z,b=1) {
      x=exp(z); 
      xmin=fan.area(20); 
      ifelse(x>xmin,b*x,0)
}

###########################################################
# Use the 'floor and ceiling' method to avoid eviction 
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

#####################################################
# Make the P matrix for the SHASH model 
#####################################################
m=500; L=0; U=10; L1 = 1.35; U1 = 7.9; 
out = mk_P_ceiling(m, L, U, L1, U1,P_z1zSHASH); 
Pmat = out$P; h=out$h; meshpts=out$meshpts; 
truncMesh = pmax(pmin(meshpts,U1),L1); 

######### check for eviction: are column sums < s(z)? 
plot(meshpts,apply(Pmat,2,sum)); 

points(meshpts,s_z(truncMesh),type="l",col="blue",lwd=2); 
abline(v=c(L1,U1),col="blue",lty=2); ## all is good 

#####################################################
# Make the P matrix for the pilot
#####################################################
out = mk_P_ceiling(m, L, U, L1, U1,P_z1zPilot); 
PmatPilot = out$P; 

######### check for eviction: are column sums < s(z)? 
plot(meshpts,apply(PmatPilot,2,sum)); 

points(meshpts,s_z(truncMesh),type="l",col="blue",lwd=2); 
abline(v=c(L1,U1),col="blue",lty=2); ## all is good 

######## recruitment vector 
bvec = h*recruitFun(meshpts); bvec=bvec/sum(bvec);
plot(meshpts,bvec); 

#####################################################
# Finally, compare
#####################################################
N = solve(diag(m)-Pmat);
NPilot = solve(diag(m)-PmatPilot); 

# compare remaining lifespan as a function of current size 
matplot(meshpts,cbind(apply(N,2,sum),apply(NPilot,2,sum)), xlim=c(L1,U1), type="l",lty=c(1,2)); 
sum(bvec*apply(N,2,sum)); sum(bvec*apply(NPilot,2,sum)); 

ss = N%*%bvec; 
ssPilot = NPilot%*%bvec; 
sum(ss)/sum(ssPilot); 

zq = quantile(XH$logarea.t0,c(0.05,0.50,0.95));
j5 = which.min(abs(meshpts- 1.856298))  # median size of a new recruit
j50 = which.min(abs(meshpts-zq[2]));
j95 = which.min(abs(meshpts-zq[3]));

pdf("../manuscript/figures/CoralKernelCompare_v2.pdf",height = 6, width = 4,useDingbats = F)
par(mfrow=c(2,1),yaxs="i",xaxs="i",bty="l",mar=c(4,4,2,1),mgp=c(2.2,1,0)); 
matplot(meshpts,cbind(Pmat[,j5],PmatPilot[,j5],Pmat[,j50],PmatPilot[,j50],Pmat[,j95],PmatPilot[,j95]),type="l",
col=c("black","red"),lty=c(1,2),xlab="Initial size (log area)",ylab="Subsequent size distribution",lwd=2); 
abline(h=0,col="black",lwd=2); 
add_panel_label("a") 
legend("topleft",bty="n",legend=c("SHASH model","Gaussian pilot"),col=c("black","red"),lty=c(1,2),lwd=2); 

matplot(meshpts,cbind(ss,ssPilot),type="l",lty=c(1,2),col=c("black","red"),xlab="Size (log area)",ylab="Frequency",lwd=2); 
add_panel_label("b"); 
dev.off()

## why does the SHASH model (which includes negative skew at small sizes)
## predict that seedlings grow a bit more than the Gaussian?
## I think it is because modeling skew shifted the fitted mean -- visualize this
plot(XH$logarea.t0,XH$logarea.t1,pch=1,col=alpha("black",0.25),
     xlab="size t",ylab="size t1")
points(XH$logarea.t0,predict(fitGAU)[,1],col=alpha("red",0.25),pch=16,cex=.5)
points(XH$logarea.t0,predict(fitSHASH)[,1],col=alpha("blue",0.25),pch=16,cex=.5)
## yes that's right
















