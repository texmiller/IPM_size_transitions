###############################################################
# Functions for a Standardized Jones-Pewsey distribution, 
# a two-parameter family with mean=0 and variance=1. 
#
# The "Jones Pewsey distribution" is equivalent to SHASHo
# in gamlss with mu=0 and sigma=1, which DOES NOT have
# zero mean and unit variance.  
#
# Epsilon is the skew parameter, negative/positive values
# giving negative/positive skew. 
# It corresponds to nu in gamlss. 
#
# delta>0 is the tail weight parameter. Values < 1 give fatter
# than Gaussian tails, values > 1 give thinner. 
# It corresponds to tau in gamlss. 
##############################################################
require(gamlss.dist); require(maxLik);  

#################################################
# Functions for the original JP distribution 
#################################################
## See Jones & Pewsey 2009, p. 764 
Pq = function(q) {
    fac=exp(0.25)/sqrt(8*pi); 
    t1 = besselK(x=1/4, nu = 0.5*(q+1)); 
    t2 = besselK(x=1/4, nu = 0.5*(q-1)); 
    return(fac*(t1+t2))
}

#### epsilon is real-valued, delta > 0
JPmean = function(epsilon,delta) {
        sinh(epsilon/delta)*Pq(1/delta)
}

#### epsilon is real-valued, delta > 0
JPvar = function(epsilon,delta) {
        EX2 = 0.5*(cosh(2*epsilon/delta)*Pq(2/delta) -1)
        return( EX2 - JPmean(epsilon,delta)^2 ) 
}

#### epsilon is real-valued, delta > 0
dJP = function(x,epsilon,delta) {
        S = sinh(delta*asinh(x)- epsilon) 
        C = sqrt(1+S^2); 
        fac = 1/sqrt(2*pi*(1+x^2))
        f = fac*delta*C*exp(-S^2/2)
        return(f)
}        


qJP = function(p, epsilon=0, delta=1) {
    nu = epsilon; tau=delta; mu=0; sigma=1; 
    return(  sinh((1/tau) * asinh(qnorm(p)) + (nu/tau)) )
}

rJP = function(n, epsilon=0, delta=1){
    U = runif(n); 
    return (qJP(U,epsilon,delta))
}
    
JPmean(0,1); JPvar(0,1);  # should be zero, 1 

   
x=seq(-5,5,length=100); 
plot(x,dJP(x,1,2),type="l",lwd=2); 
points(x,dSHASHo(x,0,1,1,2),type="p",lty=2,col="red"); # should overplot     
    
qJP(0.26, -1, 2); qSHASHo(0.26,0,1,-1,2); 


###################################################
# Functions for the standardized JP distribution 
###################################################
dSJP = function(x,epsilon,delta) {
    mu = JPmean(epsilon,delta)
    sigma = sqrt(JPvar(epsilon,delta))
    return( sigma*dJP(mu + sigma*x,epsilon,delta) )
}    

rSJP = function(n, epsilon=0, delta=1){
    U = runif(n); 
    X = qJP(U,epsilon,delta)
    mu = JPmean(epsilon,delta)
    sigma = sqrt(JPvar(epsilon,delta))
    return ( (X-mu)/sigma ) 
    
}

TESTING=FALSE; 
if(TESTING) {######################
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

} ###################################


###################################################
# Nonparametric skew and kurtosis functions. 
# These are the same for JP and standardized JP
###################################################

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

########################################################################
# Fit parameters of standardized JP by maximum Likelihood
# using BHHH from MaxLik package.  
#
# Maximization uses multistart with random initial parameters,
# with input parameter nstart specifying the number of random starts. 
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





