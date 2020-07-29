###########################################################################################
# Fit a gamlss distribution family to a set of values by maximum likelihood using MaxLik.
#    y is the set of values
#    DIST is the name of the distribution family (e.g., DIST = "SEP2") 
#    density is the density function, e.g. density = dSEP2
# DIST and density have to match. There has to be a better way than requirint the user to
# give both as argument, but I haven't figured one out. 
#  
# This function takes advantage of the structure of a gamlss family, to specify start
# values for parameters based on the data, and to guarantee valid parameters through the
# use of the family's link and link-inverse functions. 
#
# It does 10-fold multistart with jittered initial parameters, and does a final fit
# from the best of them. It's not entirely unreasonable to consider that it might 
# be somewhat reliable, though 20 or 100 would be more reassuring. 
#
# The returned value is a maxLik fit, with parameters on the inverse-link scale (e.g,
# typically log(sigma) rather than sigma). 
###########################################################################################

gamlssMaxlik <- function(y,DIST,density) {
	fam = as.gamlss.family(DIST); 
	mu = mean(eval(fam$mu.initial, list(y = y)))
	sigma = mean(eval(fam$sigma.initial, list(y = y, mu = mu)))
	nu=mean(eval(fam$nu.initial, list(y = y, mu = mu, sigma = sigma)))
	tau = mean(eval(fam$tau.initial, list(y = y, mu = mu, sigma = sigma, nu = nu)))
	
	LogLik=function(pars,response){
		val = density(x=response, mu=fam$mu.linkinv(pars[1]), 
								  sigma=fam$sigma.linkinv(pars[2]), 
								  nu=fam$nu.linkinv(pars[3]),
								  tau=fam$tau.linkinv(pars[4]),log=TRUE); 
		return(val); 
}

start=c(eta.mu=fam$mu.linkfun(mu),eta.sigma=fam$sigma.linkfun(sigma),eta.nu=fam$nu.linkfun(nu),eta.tau=fam$tau.linkfun(tau)); 
bestPars=numeric(4); bestMax=-10^16; 
for(jrep in 1:10) {
	startj = start*exp(0.1*rnorm(4)); 
	fit = maxLik(logLik=LogLik,start=startj, response=y, method="BHHH",control=list(iterlim=5000,printLevel=0),
		finalHessian=FALSE); 
	if(fit$maximum==0) fit$maximum = -10^16; 	# failed fit
	if(fit$code>2) fit$maximum = -10^16; 		# failed fit
	if(fit$maximum > bestMax)	{bestPars=fit$estimate; bestMax=fit$maximum}
	cat(jrep,fit$maximum,"\n"); 
}	

fit = maxLik(logLik=LogLik,start=bestPars, response=y, method="BHHH",control=list(iterlim=5000,printLevel=0),
		finalHessian=FALSE); 
			
return(fit); 
	
}

#y = rSHASHo(500,mu=1,sigma=1,nu=2,tau=2); 
#out = gamlssMaxlik(y=y, DIST="SHASHo",density=dSHASHo);
#out$estimate; 
