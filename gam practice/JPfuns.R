##################################################################
## Functions for Jones-Pewsey (JP) distribution and relatives.  
##
## (1) JP(epsilon, delta) 
##     This is equivalent to SHASHo in gamlss with mu=0 and sigma=1.
##     It does NOT have mean=0 or variance=1. It is a 
##     a two-parameter family with skew and kurtosis parameters.
##   
##     Our notation follows the original paper (Jones and Pewsey 2009).
##     - epsilon is the skew parameter. Negative/positive values give
##       negative/positive skewness. 
##     - delta > 0 is the tail-weight parameter. Values < 1 
##       give fatter than Gaussian tails, values > 1 give thinner. 
## 
##     Because gamlss is poorly documented, we provide here  
##     self-contained code for the density dJP, quantile 
##     function qJP, and random number generation rJP. 
##
## (2) SJP(epsilon, delta) 
##     This is a centered and scaled version of JP. 
##     It has zero mean and unit variance, for all values of 
##     epsilon and delta. 
##
## (3) RSJP(lambda, tau) 
##     This is the "reparameterised" SJP distribution. The parameters
##     for RSJP are lambda = exp(-delta) and tau = epsilon/delta. 
##     This reparameterization reduces the undesirable feature of
##     that changes in the tail-weight parameter also have a 
##     large effect on the skewness, and results in more
##     reliable parameter estimation.
##  
##  (4) JPLS(mean,sd,lambda,tau)
##      This is the four-parameter reparameterized JP distribution,  
##      where "LS" stands for location-scale. The (lambda, tau) 
##      parameters are used for skew and kurtosis, 
##            lambda = exp(-delta), tau = epsilon/delta 
##      in the notation of Jones and Pewsey (2009). 
##      The "mean" and "sd" arguments are the actual mean and std dev. 
######################################################################      
 
 
 ###############################################################
 ## The mgcv parameterization for gam.family 'shash'
 ## mu and sigma control location and scale, epsilon determines
 ## skewness (same sign as epsilon), and delta > 0 controls tail
 ## weight. Tails are heavier than Gaussian for delta>1, lighter
 ## for delta < 1. 
 ###############################################################
dshash_mgcv = function(x,mu=0,sigma=1,epsilon=0,delta=1) {  
	z = (x - mu)/(sigma*delta); 
	Sz = sinh(delta*asinh(z) - epsilon); 
	Cz = sqrt( 1 + Sz^2 )
	fac = 1/sqrt(2*pi*(1+z^2));
	return(Cz*exp(-0.5*Sz^2)*fac/sigma)
}	

COMPARING = TRUE; 
if(COMPARING) {
mu = rnorm(1); sigma=runif(1)*5; epsilon=rnorm(1); delta=exp(rnorm(1)); 
integrate(function(x) dshash_mgcv(x,mu,sigma,epsilon,delta), -Inf, Inf)$value

m1 = integrate(function(x) x*dshash_mgcv(x,mu,sigma,epsilon,delta), -Inf, Inf)$value
m2 =  integrate(function(x) (x^2)*dshash_mgcv(x,mu,sigma,epsilon,delta), -Inf, Inf)$value
m3 = integrate(function(x) ((x-m1)^3)*dshash_mgcv(x,mu,sigma,epsilon,delta), -Inf, Inf)$value
m4 = integrate(function(x) ((x-m1)^4)*dshash_mgcv(x,mu,sigma,epsilon,delta), -Inf, Inf)$value
cat(mu, sigma, m1, sqrt(m2-m1^2), m3, m4,"\n"); 

require(gamlss.dist); 
m1 = integrate(function(x) x*dSHASHo2(x,mu,sigma,epsilon,delta), -Inf, Inf)$value
m2 =  integrate(function(x) (x^2)*dSHASHo2(x,mu,sigma,epsilon,delta), -Inf, Inf)$value
m3 = integrate(function(x) ((x-m1)^3)*dSHASHo2(x,mu,sigma,epsilon,delta), -Inf, Inf)$value
m4 = integrate(function(x) ((x-m1)^4)*dSHASHo2(x,mu,sigma,epsilon,delta), -Inf, Inf)$value
cat(mu, sigma, m1, sqrt(m2-m1^2), m3, m4,"\n"); 
} 
 
##########################################################
##              JP Distribution 
## Functions for the two-parameter JP distribution 
## See Jones & Pewsey 2009, p. 764.  
## epsilon is real-valued, delta > 0
##########################################################

## Utility function for mean and variance 
Pq = function(q) {
    fac=exp(0.25)/sqrt(8*pi); 
    t1 = besselK(x=1/4, nu = 0.5*(q+1)); 
    t2 = besselK(x=1/4, nu = 0.5*(q-1)); 
    return(fac*(t1+t2))
}


#### Probability density function
dJP = function(x,epsilon,delta) {
        S = sinh(delta*asinh(x)- epsilon) 
        C = sqrt(1+S^2); 
        fac = 1/sqrt(2*pi*(1+x^2))
        f = fac*delta*C*exp(-S^2/2)
        return(f)
}        

#### Quantile function 
qJP = function(p, epsilon=0, delta=1) {
    return(sinh((1/delta) * asinh(qnorm(p)) + (epsilon/delta)) )
}

#### Cumulative distribution function 
pJP = function (q, epsilon = 0, delta = 1) {
    r <- sinh(delta * asinh(q) - epsilon)
    return(pnorm(r))
}


### Random number generation 
rJP = function(n, epsilon=0, delta=1){
    U = runif(n); 
    return (qJP(U,epsilon,delta))
}


### Mean and variance 
JPmean = function(epsilon,delta) {
        sinh(epsilon/delta)*Pq(1/delta)
}

JPvar = function(epsilon,delta) {
        EX2 = 0.5*(cosh(2*epsilon/delta)*Pq(2/delta) -1)
        return( EX2 - JPmean(epsilon,delta)^2 ) 
}        


TESTING=FALSE; 
if(TESTING) { #---- Check mean and var; compare with dSHASHo in gamlss.dist 
  require(gamlss.dist)
JPmean(0,1); JPvar(0,1);  # Gives N(0,1), so values should be 0, 1 

x=seq(-5,5,length=100); 
plot(x,dJP(x,1,2),type="l",lwd=1); 
points(x,dSHASHo(x,0,1,1,2),type="p",lty=2,col="red"); # should overplot     
    
qJP(0.26, -1, 2); qSHASHo(0.26,0,1,-1,2); 
} #--------------------------------------------------------------------

###############################################################
#              SJP Distribution 
# Functions for the standardized JP distribution.  
# JP distribution is centered and scaled so mean=0, variance=1 
###############################################################
 
## probability density function 
dSJP = function(x,epsilon=0,delta=1) {
    mu = JPmean(epsilon,delta);
    sigma = sqrt(JPvar(epsilon,delta));
    return( sigma*dJP(mu + sigma*x,epsilon,delta) );
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
    
    
TESTING=FALSE;     
if(TESTING) {
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

} 


######################################################
# Nonparametric skew and kurtosis functions. 
# These are the same for JP and SJP distributions. 
######################################################
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


###############################################################
##              RSJP Distribution  
##       Reparameterized SJP distribution 
## Functions for SJP distribution in (lambda, tau) parameters
## where lambda controls skewness and tau controls kurtosis. 
## In terms of the parameters of the original distribution, 
## tau = -log(delta) and lambda = epsilon/delta. 
##############################################################

## probability density function 
dRSJP = function(x,lambda=0,tau=0) {
    delta=exp(-tau); epsilon=lambda*delta;
    return(dSJP(x,epsilon,delta)); 
}    


## random number generation 
rRSJP = function(n, lambda=0, tau=0){
    delta=exp(-tau); epsilon=lambda*delta;
    return(rSJP(n,epsilon,delta)); 
}

#### quantile function 
qRSJP = function(p, lambda=0, tau=0) {
    delta=exp(-tau); epsilon=lambda*delta;
    return(qSJP(p,epsilon,delta)); 
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

#########################################################################
##         JPLS Distribution  (LS = "location-scale") 
## Four-parameter JP distribution, reparameterized so that 
## (1) mu and sigma are the actual mean and standard deviation.  
## (2) The (lambda, tau) parameters control skew and kurtosis
## In terms of the parameters of the original distribution, 
## tau = -log(delta) and lambda = epsilon/delta. 
## 
## The (lambda, tau) parameterization reduces the undesirable 
## feature of the (epsilon, delta) parameterization that changes 
## in the tail-weight parameter also have a large effect on the 
## skewness, and results in more reliable parameter estimation.
#########################################################################

## probability density function 
dJPLS = function(x, mean=0, sd=1, lambda=0, tau=0) {
    delta=exp(-tau); epsilon=lambda*delta;
    return((1/sd)*dSJP((x-mean)/sd,epsilon,delta)); 
}    

## random number generation 
rJPLS = function(n, mean=0, sd=1, lambda=0, tau=0){
    delta=exp(-tau); epsilon=lambda*delta;
    return(mean +sd*rSJP(n,epsilon,delta)); 
}

#### quantile function 
qJPLS = function(p, mean=0, sd=1, lambda=0, tau=0) {
    delta=exp(-tau); epsilon=lambda*delta;
    return(mean + sd*qSJP(p,epsilon,delta)); 
}

#### cumulative distribution function 
pJPLS = function (x, mean=0, sd=1, lambda=0, tau=0) {
    delta=exp(-tau); epsilon=lambda*delta;
    xs = (x - mean)/sd; 
    return(pSJP(xs,epsilon,delta))
}
pJPLS = Vectorize(pJPLS,vectorize.args="x"); 


TESTING=FALSE;     
if(TESTING) {
# Testing the moments 
for(j in 1:10) {
 mu = rnorm(1); sigma = exp(rnorm(1)); lambda=rnorm(1); tau = rnorm(1); 
 u = integrate(function(x) dJPLS(x,mu,sigma,lambda,tau),-Inf,Inf)$value # should equal 1
 v = integrate(function(x) x*dJPLS(x,mu,sigma,lambda,tau), -Inf, Inf)$value   # should equal mu 
 z = integrate(function(x) (x^2)*dJPLS(x,mu,sigma,lambda,tau), -Inf, Inf)$value  - v^2  # should equal 1 
cat(round(u,digits=6), round(mu-v,digits=6), round(z-sigma^2,digits=6), "\n"); 
}
}

if(TESTING) {
# Testing RNG and quantiles 
for(j in 1:10) {
 mu = rnorm(1); sigma = exp(rnorm(1)); lambda=rnorm(1); tau = rnorm(1); qu = runif(1); 
 X = rJPLS(10^7,mu,sigma,lambda,tau); 
 cat(round(mu-mean(X),digits=4), round(sigma-sd(X), digits=4), round(quantile(X,qu) - qJPLS(qu, mu,sigma,lambda,tau),digits=4), "\n"); 
}
}

if(TESTING) {
# Testing CDF 
for(j in 1:10) {
 mu = rnorm(1); sigma = exp(rnorm(1)); lambda=rnorm(1); tau = rnorm(1); qu = runif(1); 
 q = qJPLS(qu, mu, sigma, lambda, tau) 
 p = pJPLS(q, mu, sigma, lambda, tau) 
 cat(round(qu-p,digits=6),  "\n"); 
}
}


########################################################################
## Fit parameters (epsilon, delta) of SJP by maximum Likelihood, 
## using BHHH from MaxLik. 
## Maximization uses multistart with random initial parameters,
## with input parameter nstart specifying the number of random starts. 
#######################################################################

require(maxLik); 

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
## Fit parameters (lambda, tau) of RSJP by maximum Likelihood 
## using BHHH from MaxLik.
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









