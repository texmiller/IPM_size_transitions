###########################################################################################
# Fit a gamlss distribution family to a set of values by maximum likelihood using MaxLik.
#    y is the set of values
#    DIST is the name of the distribution family (e.g., DIST = "JSU")
# 
# The function can accommodate 2-, 3-, and 4-parameter gamlss distribution families
# (this covers all of the continuous gamlss families except exponential)
#  
# The function takes advantage of the structure of a gamlss family, to specify start
# values for parameters based on the data, and to guarantee valid parameters through the
# use of the family's link and link-inverse functions. The multiple uses of the
# .linkfun and .lininv for distribution parameters are (sadly) necessary to make
# sure that the optimizer can send any real numbers as arguments to the likelihood
# function, and they translate to valid distribution parameters (e.g., sigma > 0). 
#
# Maximization uses 10-fold multistart with jittered initial parameters, plus a final fit
# starting from the best of them. It's not entirely unreasonable to consider that it might 
# be somewhat reliable, though 20 or 100 would be more reassuring. 
#
# RETURNED VALUE is a maxLik() fit, with parameters on the family's link-transformed 
# scale (e.g, typically log(sigma) rather than sigma). 
#
###########################################################################################

gamlssMaxlik <- function(y,DIST) {
  out <- list()
  aics <- numeric(length=length(DIST))
  for(d in 1:length(DIST)){
  fam = as.gamlss.family(DIST[d])
  n_par <-   fam$nopar
  
  ## assume at least two parameters  
  mu = mean(eval(fam$mu.initial, list(y = y)))
  sigma = mean(eval(fam$sigma.initial, list(y = y, mu = mu)))
  start=c(eta.mu=fam$mu.linkfun(mu),eta.sigma=fam$sigma.linkfun(sigma))
  ## and maybe a third
  if(n_par==3){
    nu=mean(eval(fam$nu.initial, list(y = y, mu = mu, sigma = sigma)))
    start=c(start,eta.nu=fam$nu.linkfun(nu))
  }
  ## and maybe a fourth
  if(n_par==4){
    nu=mean(eval(fam$nu.initial, list(y = y, mu = mu, sigma = sigma)))
    tau=mean(eval(fam$tau.initial, list(y = y, mu = mu, sigma = sigma, nu = nu)))
    start=c(start,eta.nu=fam$nu.linkfun(nu),eta.tau=fam$tau.linkfun(tau))
  }
  
  LogLik1=function(pars,response){
    if(n_par==2) fun_args = list(x=response, mu=fam$mu.linkinv(pars[1]), sigma=fam$sigma.linkinv(pars[2]),log=TRUE)
    if(n_par==3) fun_args = list(x=response, mu=fam$mu.linkinv(pars[1]), 
                                   sigma=fam$sigma.linkinv(pars[2]), 
                                   nu=fam$nu.linkinv(pars[3]),log=TRUE)
    if(n_par==4) fun_args = list(x=response, mu=fam$mu.linkinv(pars[1]), 
                                   sigma=fam$sigma.linkinv(pars[2]), 
                                   nu=fam$nu.linkinv(pars[3]),
                                   tau=fam$tau.linkinv(pars[4]),log=TRUE)
    val = do.call(paste("d",DIST[d],sep=""),fun_args)
    return(val); 
  }
  
  bestPars=numeric(n_par); bestMax=-10^17; 
  for(jrep in 1:20) {
    startj = start*exp(0.1*rnorm(n_par)); 
    fit = maxLik(logLik=LogLik1,start=startj, response=y, method="BHHH",control=list(iterlim=5000,printLevel=0),
                finalHessian=FALSE);  
    
    if(fit$maximum==0) fit$maximum = -10^16; 	# failed fit
    if(fit$code>2) fit$maximum = -10^16; 		# failed fit
    if(fit$maximum > bestMax)	{bestPars=fit$estimate; bestMax=fit$maximum; bestFit=fit;}
    cat(DIST[d],jrep,fit$maximum,"\n"); 
  }	
  
  #fit = maxLik(logLik=LogLik1,start=bestPars, response=y, method="BHHH",control=list(iterlim=5000,printLevel=0),
  #              finalHessian=FALSE); 
  #if(fit$maximum==0) fit$maximum = -10^16; 	# failed fit
  #if(fit$code>2) fit$maximum = -10^16; 		# failed fit                
  fit=bestFit; 
  aics[d] <- 2*n_par - 2*fit$maximum
  out[[d]] <- fit
  }

  return(list(out=out,aics=aics)); 
}

if(FALSE) {
## can I recover parameters of a normal distribution?
y = rNO(5000,mu=1,sigma=1)
out.NO = gamlssMaxlik(y=y, DIST="NO")
out.NO$estimate[1];exp(out.NO$estimate[2])
# yes. does it favor a normal over other candidates? -- picking a few favorites
out.JSU = gamlssMaxlik(y=y, DIST="JSU")
out.SHASH = gamlssMaxlik(y=y, DIST="SHASH")
out.ST1 = gamlssMaxlik(y=y, DIST="ST1")
# yes, normal is selected
out.NO$AIC;out.JSU$AIC;out.SHASH$AIC;out.ST1$AIC
}
