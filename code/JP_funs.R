####################################################################
## Functions for Jones-Pewsey (JP) distribution and relatives.  
##
## (1) JP2(epsilon, delta) 
##		The Jonew-Pewsey (2009) 2-parameter distribution family, 
## 		with epsilon determining skewness and delta determining tail 
## 		weight. This is the distribution of $X_\{epsilon,delta}$ in 
##		Jones and Pewsey (2009). It is equivalent to SHASHo2 in 
##		gamlss with parameters mu=0, sigma=1 (but its mean is not
##		zero, and its variance is not 1).  
##   
##     Our notation follows the original paper (Jones and Pewsey 2009).
##     - epsilon is the skew parameter. Negative/positive values give
##       negative/positive skewness. 
##     - delta > 0 is the tail-weight parameter. Values < 1 
##       give fatter than Gaussian tails, values > 1 give thinner. 
## 
##
## (2) JPS(mu,sigma,epsilon, delta) 
##     This is a four-parameter family in which the JP2(epsilon,delta) 
##     distribution is shifted and scaled to have mean=mu, and standard 
##     deviation=sigma, for ALL values of epsilon and delta. To our 
##	   knowledge it has not been introduced elsewhere. 
##
##
## (3) JPR(mu, sigma, lambda, tau) 
##     This is a "reparameterised" JPS distribution. The skewness
##     and kurtosis parameters are lambda = exp(-delta) and 
## 	   tau = epsilon/delta. This reparameterization reduces the 
## 	   undesirable feature of JPS that changes in the tail-weight 
## 	   parameter also have a large effect on skewness, and improves 
##     the reliability of parameter estimation.
##  
########################################################################      
 
 
##########################################################
##              JP2 Distribution 
## Functions for the two-parameter JP distribution 
## See Jones & Pewsey 2009, p. 764.  
## epsilon is real-valued, delta > 0
##########################################################

#### Probability density function
dJP2 = function(x,epsilon,delta) {
        S = sinh(delta*asinh(x)- epsilon) 
        C = sqrt(1+S^2); 
        fac = 1/sqrt(2*pi*(1+x^2))
        f = fac*delta*C*exp(-S^2/2)
        return(f)
}        

#### Quantile function 
qJP2 = function(p, epsilon=0, delta=1) {
    return(sinh( (1/delta) * asinh(qnorm(p)) + (epsilon/delta)) )
}

#### Cumulative distribution function 
pJP2 = function (x, epsilon = 0, delta = 1) {
    r <- sinh(delta * asinh(x) - epsilon)
    return(pnorm(r))
}


### Random number generation 
rJP2 = function(n, epsilon=0, delta=1){
    U = runif(n); 
    return (qJP2(U,epsilon,delta))
}

##########################################
## Mean and variance of JP2 
##########################################

## Utility function, see Jones and Pewsey (2009, p. 764)
Pq = function(q) {
    fac=exp(0.25)/sqrt(8*pi); 
    t1 = besselK(x=1/4, nu = 0.5*(q+1)); 
    t2 = besselK(x=1/4, nu = 0.5*(q-1)); 
    return(fac*(t1+t2))
}

JP2mean = function(epsilon,delta) {
        sinh(epsilon/delta)*Pq(1/delta)
}

JP2var = function(epsilon,delta) {
        EX2 = 0.5*(cosh(2*epsilon/delta)*Pq(2/delta) -1)
        return( EX2 - JP2mean(epsilon,delta)^2 ) 
}        

JP2sd = function(epsilon,delta) { sqrt(JP2var(epsilon,delta)) } 


###############################################################
##              JPS Distribution 
## Functions for the JPS distribution, four parameter family
## consisting of shifted and scaled JP2 distributions having
## mean=mu, standard deviation=sigma for all epsilon, delta  
###############################################################
 
## probability density function 
dJPS = function(x, mean=0, sd=1, epsilon=0, delta=1) {
    m = JPmean(epsilon,delta);
    s = JPsd(epsilon,delta);
	mu = mean; sigma = sd; 
	fac = s/sigma;
	return( fac*dJP2(m + fac*(x-mu)) )
} 	
	
## random number generation 
rJPS = function(n, mean=0, sd=1, epsilon=0, delta=1){
	Z = rJP2(n,epsilon,delta) 
	Z01 = (Z - JP2mean(epsilon,delta))/JP2sd(epsilon,delta); 
    return(mean + sd*Z01)
} 











#### quantile function 
qJPS = function(p, mean=0, sd=1, epsilon=0, delta=1) {
    q = qJP2(p,epsilon,delta); 
    m = JP2mean(epsilon,delta);
    s = JPsd(epsilon,delta));
    return( (q - m)/sigma );
}

#### cumulative distribution function 
pSJP = function (q, epsilon = 0, delta = 1) {
   mu = JPmean(epsilon,delta)
   sigma = sqrt(JPvar(epsilon,delta))
   qs = mu + sigma*q
   return(pJP(qs,epsilon,delta))
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









