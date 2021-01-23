z = log(CYIM$vol_t); z1=scaledResids;

rollMomentsNP(scaledResids); 

make_pql = function(z,pars){
    u = pars[1] + pars[2]*z + pars[3]*z^2
    lambda = -1 + 2*exp(u)/(1+exp(u)) 
    u = pars[4] + pars[5]*z + pars[6]*z^2; pz = 10*exp(u)/(1+exp(u)); 
    u = pars[7] + pars[8]*z + pars[9]*z^2; qz = 100*exp(u)/(1+exp(u));
    qz = pmax(qz,2.01/pz); 
    return(list(lambda=lambda,p=pz,q=qz)); 
 }

sgt_logLik = function(pars) {
    pql = make_pql(z,pars); 
    lik = dsgt(z1, mu = 0, sigma = 1, lambda = pql$lambda, p = pql$p, q = pql$q, mean.cent = TRUE, var.adj = TRUE, log = TRUE)
    return(lik)    
}

sgt_logLik1 = function(pars) {
    val = sum(sgt_logLik(pars))
    if(!is.finite(val)) val = -10^16
    return(val)
}    

par0 = c(0,0,0,-2,0,0,-2,0,0);     
mfit = maxLik(logLik=sgt_logLik,start=par0,method="BHHH",control=list(iterlim=1000,printLevel=2)); 

#fit=optim(par=par0, fn=sgt_logLik1,method="SANN",control=list(maxit=25000,trace=4,fnscale=-1));  

#for(j in 1:5) {
#fit1 = optim(par=fit$par, fn=sgt_logLik1,control=list(maxit=5000,trace=4,REPORT=2,fnscale=-1));  
#fit = optim(par=fit1$par, fn=sgt_logLik1,control=list(maxit=5000,trace=4,fnscale=-1));  
#}

### plot the results 
pars=fit$par
zx = seq(min(z),max(z),length=100); 
pql = make_pql(zx,pars)
    
par(mfrow=c(2,2),bty="l",mar=c(4,4,1,1));
plot(zx,pql$lambda); plot(zx,pql$p); plot(zx,pql$q); 
    
px=seq(-4,4,length=200); py=matrix(NA,length(px),3); 
pql=make_pql(quantile(z,0.15),pars); 
py[,1] = dsgt(px, mu = 0, sigma = 1, lambda = pql$lambda, p = pql$p, q = pql$q,
mean.cent = TRUE, var.adj = TRUE, log = FALSE)

pql=make_pql(quantile(z,0.5),pars); 
py[,2] = dsgt(px, mu = 0, sigma = 1, lambda = pql$lambda, p = pql$p, q = pql$q,
mean.cent = TRUE, var.adj = TRUE, log = FALSE)

pql=make_pql(quantile(z,0.9),pars); 
py[,3] = dsgt(px, mu = 0, sigma = 1, lambda = pql$lambda, p = pql$p, q = pql$q,
mean.cent = TRUE, var.adj = TRUE, log = FALSE)

matplot(px,py, type="l",col=c("red","black","blue"),lty=1); 



