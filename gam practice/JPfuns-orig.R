###############################################################
# Functions for Jones-Pewsey (JP) distribution and relatives.  
#
# JP is a two-parameter family with skew and kurtosis parameters.  
# Our notation follows the original paper (Jones and Pewsey 2009).
# epsilon is the skew parameter (negative/positive values give
# negative and positive skew). delta > 0 is the tail-weight parameter.
# Values < 1 give fatter than Gaussian tails, values > 1 give thinner. 
# 
# JP is equivalent to SHASHo in gamlss with mu=0 and sigma=1; 
# though neither one has zero mean an unit variance. Because  
# gamlss is poorly documented, rather than relying on it we 
# self-contained code for the density dJP, quantile function qJP, and 
# random number generation rJP. 
#
# SJP is a centered and scaled version of JP, which has
# zero mean and unit variance for all values of epsilon and delta. 
#
# RSJP is the "reparameterised" SJP distribution. The parameters
# for RSJP are lambda = exp(-delta) and tau = epsilon/delta. 
# This reparameterization reduces the undesirable feature of
# JP and SJP that changes in the tail-weight parameter also
# have a large effect on the skewness, and results in more
# reliable parameter estimation.  
 
##############################################################
require(gamlss.dist); require(maxLik);  

#################################################
## Functions for the two-parameter JP distribution 
## See Jones & Pewsey 2009, p. 764 
## epsilon is real-valued, delta > 0
#################################################

## Utility function for mean and variance 
Pq = function(q) {
    fac=exp(0.25)/sqrt(8*pi); 
    t1 = besselK(x=1/4, nu = 0.5*(q+1)); 
    t2 = besselK(x=1/4, nu = 0.5*(q-1)); 
    return(fac*(t1+t2))
}


#### probability density function
dJP1 = function(x,epsilon,delta) {
        S = sinh(delta*asinh(x)- epsilon) 
        C = sqrt(1+S^2); 
        fac = 1/sqrt(2*pi*(1+x^2))
        f = fac*delta*C*exp(-S^2/2)
        return(f)
}        

#### quantile function 
qJP1 = function(p, epsilon=0, delta=1) {
    return(sinh((1/delta) * asinh(qnorm(p)) + (epsilon/delta)) )
}

#### cumulative distribution function 
pJP1 = function (q, epsilon = 0, delta = 1) {
    r <- sinh(delta * asinh(q) - epsilon)
    return(pnorm(r))
}


### random number generation 
rJP1 = function(n, epsilon=0, delta=1){
    U = runif(n); 
    return (qJP1(U,epsilon,delta))
}


### mean and variance 
JPmean = function(epsilon,delta) {
        sinh(epsilon/delta)*Pq(1/delta)
}

JPvar = function(epsilon,delta) {
        EX2 = 0.5*(cosh(2*epsilon/delta)*Pq(2/delta) -1)
        return( EX2 - JPmean(epsilon,delta)^2 ) 
}        


TESTING=FALSE; 
if(TESTING) { #--------------------------------------------------------
JPmean(0,1); JPvar(0,1);  # Gives N(0,1), so values should be 0, 1 

x=seq(-5,5,length=100); 
plot(x,dJP(x,1,2),type="l",lwd=2); 
points(x,dSHASHo(x,0,1,1,2),type="p",lty=2,col="red"); # should overplot     
    
qJP(0.26, -1, 2); qSHASHo(0.26,0,1,-1,2); 
} -----------------------------------------------------------------------

######################################################
# Functions for the standardized JP distribution, i.e. 
# centered and scaled so mean=0, variance=1 
######################################################

## probability density function 
dSJP = function(x,epsilon,delta) {
    mu = JPmean(epsilon,delta)
    sigma = sqrt(JPvar(epsilon,delta))
    return( sigma*dJP(mu + sigma*x,epsilon,delta) )
}    

## random number generation 
rSJP = function(n, epsilon=0, delta=1){
    U = runif(n); 
    X = qJP(U,epsilon,delta)
    mu = JPmean(epsilon,delta)
    sigma = sqrt(JPvar(epsilon,delta))
    return ( (X-mu)/sigma ) 
    
}

#### quantile function 
qSJP = function(p, epsilon=0, delta=1) {
    q = qJP(p,epsilon,delta); 
    mu = JPmean(epsilon,delta)
    sigma = sqrt(JPvar(epsilon,delta))
    return( (q - mu)/sigma )
}

#### cumulative distribution function 
pSJP = function (q, epsilon = 0, delta = 1) {
   mu = JPmean(epsilon,delta)
   sigma = sqrt(JPvar(epsilon,delta))
   qs = mu + sigma*q
   return(pJP(qs,epsilon,delta))
}

SJP_moments=function(epsilon,delta) {
    m3 = integrate(function(x) (x^3)*dSJP(x,epsilon,delta), -Inf, Inf)$value
    m4 = integrate(function(x) (x^4)*dSJP(x,epsilon,delta), -Inf, Inf)$value
    return(list(skew = m3, excess.kurtosis = m4/3 -1))
}    
 
    

if(TESTING) {#------------------------------------------------------------
# Testing the moments 
for(j in 1:10) {
 epsilon=rnorm(1); delta=exp(rnorm(1)); 
 u = integrate(function(x) dSJP(x,epsilon,delta),-Inf,Inf)$value # should equal 1. 
 v = integrate(function(x) x*dSJP(x,epsilon,delta), -Inf, Inf)$value   # should equal 0 
 z = integrate(function(x) (x^2)*dSJP(x,epsilon,delta), -Inf, Inf)$value   # should equal 1 
cat(round(u,digits=6), round(v,digits=6), round(z,digits=6), "\n"); 
}

## Testing the RNG -- do we recover the parameters that generated the "data"? 
par(mfrow=c(2,1)); 
X = rSJP(50000,epsilon=-0.5,delta=.25);
NLL = function(p) {
    eps=p[1]; del=p[2]; 
    -sum(log(dSJP(X,eps,del)))
}

fit=optim(par=c(0,1), NLL, control=list(maxit=50000,trace=4)); 
fit$par;     

} #------------------------------------------------------------------------------


#####################################################
# Nonparametric skew and kurtosis functions. 
# These are the same for JP and SJP
#####################################################
JP_NPskewness = function(epsilon,delta,p=0.1) {
	q = qJP(c(p,0.5,1-p),epsilon,delta)
	u = (q[3]+q[1]-2*q[2])/(q[3]-q[1]);
	return(as.numeric(u)); 
	
}	

JP_NPkurtosis=function(epsilon,delta,p=0.05) {
	q = qJP(c(p,0.25,0.75,1-p),epsilon,delta)
	qN = qnorm(c(p,0.25,0.75,1-p))
	u = (q[4]-q[1])/(q[3]-q[2]);
	uN = (qN[4]-qN[1])/(qN[3]-qN[2]);
	return (as.numeric(u/uN-1)) 
}

SJP_NPskewness = JP_NPskewness; 
SJP_NPkurtosis = JP_NPkurtosis; 


############################################################
## Reparameterized SJP distribution 
## Functions for SJP distribution in (lambda, tau) parameters
############################################################

## probability density function 
dRSJP = function(x,lambda=0,tau=0) {
    delta=exp(-tau); epsilon=lambda*delta;
    return(dSJP(x,delta,epsilon)); 
}    

## random number generation 
rRSJP = function(n, lambda=0, tau=0){
    delta=exp(-tau); epsilon=lambda*delta;
    return(rSJP(x,delta,epsilon)); 
}

#### quantile function 
qRSJP = function(p, lambda=0, tau=0) {
    delta=exp(-tau); epsilon=lambda*delta;
    return(qSJP(p,delta,epsilon)); 
}

#### cumulative distribution function 
pRSJP = function (q, lambda=0, tau=0) {
    delta=exp(-tau); epsilon=lambda*delta;
    return(pSJP(q,epsilon,delta))
}

RSJP_moments=function(lambda=0, tau=0) {
    delta=exp(-tau); epsilon=lambda*delta;
    m3 = integrate(function(x) (x^3)*dSJP(x,epsilon,delta), -Inf, Inf)$value
    m4 = integrate(function(x) (x^4)*dSJP(x,epsilon,delta), -Inf, Inf)$value
    return(list(skew = m3, excess.kurtosis = m4/3 -1))
}

RSJP_NPskewness = function(lambda=0,tau=0,p=0.1) {
    delta=exp(-tau); epsilon=lambda*delta;
    u = JP_NPskewness(epsilon,delta,p); 
	return(u); 
}	

RSJP_NPkurtosis = function(lambda=0,tau=0,p=0.05) {
    delta=exp(-tau); epsilon=lambda*delta;
    u = JP_NPkurtosis(epsilon,delta,p); 
	return(u); 
}	

############################################################
## Four-parameter Reparameterized JP distribution 
## The mean and sd parameters are the actual mean and sd
## The (lambda, tau) parameters are used for skew and kurtosis
## !! These need to be checked!! 
############################################################

## probability density function 
dCRJP = function(x, mean=0, sd=1, lambda=0, tau=0) {
    delta=exp(-tau); epsilon=lambda*delta;
    return((1/sd)*dRSJP((x-mean/sd),delta,epsilon)); 
}    

## random number generation 
rCRJP = function(n, mean=0, sd=1, lambda=0, tau=0){
    delta=exp(-tau); epsilon=lambda*delta;
    return(mu +sd*rRSJP(n,delta,epsilon)); 
}

#### quantile function 
qCRJP = function(p, mean=0, sd=1, lambda=0, tau=0) {
    delta=exp(-tau); epsilon=lambda*delta;
    return(mu + sigma*qRSJP(p,delta,epsilon)); 
}

#### cumulative distribution function 
pCRJP = function (q, mean=0, sd=1, lambda=0, tau=0) {
    delta=exp(-tau); epsilon=lambda*delta;
    qs = (z - mean)/sd; 
    return(pRSJP(qs,epsilon,delta))
}

########################################################################
## Fit parameters of SJP by maximum Likelihood, using BHHH from MaxLik
##
## Maximization uses multistart with random initial parameters,
## with input parameter nstart specifying the number of random starts. 
#######################################################################

SJPMaxlik <- function(y,nstart=10,start = c(epsilon=0,log.delta = 0), sigma.start=0.2 ) {
 
  LogLik1=function(pars,response){
    val = dSJP(response,epsilon=pars[1],delta=exp(pars[2]));  
    return(log(val)); 
  }  
    
  bestPars=numeric(2); bestMax=-10^17; 
  for(jrep in 1:nstart) {
    startj = start + sigma.start*rnorm(2); 
    fit = maxLik(logLik=LogLik1,start=startj, response=y, method="BHHH",control=list(iterlim=5000,printLevel=1),
                finalHessian=FALSE);  
    fit = maxLik(logLik=LogLik1,start=fit$estimate, response=y, method="BHHH",control=list(iterlim=5000,printLevel=1),
                finalHessian=FALSE);             
    fit = maxLik(logLik=LogLik1,start=fit$estimate, response=y, method="BHHH",control=list(iterlim=5000,printLevel=1),
                finalHessian=FALSE);     
    
    if(fit$maximum==0) fit$maximum = -10^16; 	# failed fit
    if(fit$code>2) fit$maximum = -10^16; 		# failed fit
    if(fit$maximum > bestMax)	{bestPars=fit$estimate; bestMax=fit$maximum; bestFit=fit;}
    cat(jrep,fit$maximum,"\n"); 
  }	
  
  estimate = fit$estimate; estimate[2]=exp(estimate[2]); 
            
  return(list(fit=bestFit,estimate=estimate)); 
}

########################################################################
## Fit parameters of RSJP by maximum Likelihood using BHHH from MaxLik.
##
## Maximization uses multistart with random initial parameters,
## with input parameter nstart specifying the number of random starts. 
#######################################################################

RSJP_ML <- function(y,nstart=10,start = c(lambda=0,tau=0), sigma.start=0.2 ) {
 
  LogLik1=function(pars,response){
    lambda=pars[1]; tau = pars[2]
    delta=exp(-tau); epsilon=lambda*delta;
    val = dSJP(response,epsilon,delta);  
    return(log(val)); 
  }  
    
  NLL = function(pars,response) -sum(LogLik1(pars,response)); 

  bestPars=numeric(2); bestMax=-10^17; 
  for(jrep in 1:nstart) {
    startj = start + sigma.start*rnorm(2); 
    fit = optim(startj,NLL,response=y, control=list(trace=0,maxit=20000)); 
    fit = maxLik(logLik=LogLik1,start=fit$par,response=y, method="BHHH",control=list(iterlim=5000,printLevel=1),
                finalHessian=FALSE); 
    fit = optim(fit$estimate,NLL,response=y, control=list(trace=0,maxit=20000));         
    fit = maxLik(logLik=LogLik1,start=fit$par, response=y, method="BHHH",control=list(iterlim=5000,printLevel=1),
                finalHessian=FALSE);     
    
    if(fit$maximum==0) fit$maximum = -10^16; 	# failed fit
    if(fit$code>2) fit$maximum = -10^16; 		# failed fit
    if(fit$maximum > bestMax)	{bestPars=fit$estimate; bestMax=fit$maximum; bestFit=fit;}
    cat("ML fit ", jrep,fit$maximum," <--------------------->", "\n"); 
  }	
      
  return(list(fit=bestFit,estimate=fit$estimate)); 
}






