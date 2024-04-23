setwd("c:/repos/IPM_size_transitions/cactus"); #Steve

require(maxLik); 
load("cactusSetup.Rdata");  # the start of running a cactus_growth_modeling file

z = log(CYIM$vol_t); z1=scaledResids;
rollMomentsNP(z,residuals(CYIM_lmer_best)); 



CYIM_gam_model <- gam(list( log(vol_t1) ~ s(log(vol_t),k=5) + s(year_t,bs="re") + s(plot,bs="re"), ~s(log(vol_t),k=5) ), data=CYIM, gamma=1.4, family=gaulss(b=0))
plot(CYIM_gam_model,scale=0); 

rollMomentsNP(z,residuals(CYIM_gam_model,type="response")); 

fitted_vals = predict(CYIM_gam_model,type="response"); 
sigma_hat = 1/fitted_vals[,2]; 
dev.new(); plot(z,sigma_hat); 


CYIM_gam_model0 <- gam(list( log(vol_t1) ~ s(log(vol_t)) + s(year_t,bs="re") + s(plot,bs="re"), ~log(vol_t)), data=CYIM, family=gaulss(b=0))
CYIM_gam_model1 <- gam(list( log(vol_t1) ~ s(log(vol_t)) + s(year_t,bs="re") + s(plot,bs="re"), ~log(vol_t) + I(log(vol_t)^2) ), data=CYIM, family=gaulss(b=0))








###########################################################################
# Function to give p and q as a function of covariate z (initial size) 
#
# We restrict the range of possible values for lambda, p, q: 
# lambda is in (-1,1), p is in (0,50), q is in (0.5,100), and pq>2. 
# For lambda this is a requirement of the sgt distribution
# For p and q, more extreme values make little difference in shape but
# create numerical problems in evaluation of the Beta function in dsgt. 
###########################################################################
make_pql = function(z,pars){
    u = pars[1] + pars[2]*z + pars[3]*z^2
    lambda = -1 + 2*exp(u)/(1+exp(u)) 
    u = pars[4] + pars[5]*z + pars[6]*z^2; pz = 50*exp(u)/(1+exp(u)); 
    u = pars[7] + pars[8]*z + pars[9]*z^2; qz = 0.5 + 95.5*exp(u)/(1+exp(u));
    # ensure pq > 2 by adjusting both if it fails. 
    r = rep(1,length(z));  
    e = which(pz*qz < 2.01); r[e]=sqrt(2.01/(pz[e]*qz[e]));
    pz=pz*r; qz = qz*r; 
    return(list(lambda=lambda,p=pz,q=qz)); 
 }

## log likelihood function for maxLik 
sgt_logLik = function(pars) {
    pql = make_pql(z,pars); 
    lik = dsgt(z1, mu = 0, sigma = 1, lambda = pql$lambda, p = pql$p, q = pql$q, mean.cent = TRUE, var.adj = TRUE, log = TRUE)
    return(lik)    
}

## log likelihood function for optim 
sgt_logLik1 = function(pars) {
    val = sum(sgt_logLik(pars))
    if(!is.finite(val)) val = -(10^16)
    return(val)
}    

par0 = c(0,0,0,-3,0,0,-3,0,0); # start with lambda=0, p and q small 
mfit = maxLik(logLik=sgt_logLik,start=par0,method="BHHH",control=list(iterlim=1000,printLevel=2)); 

require(bbmle); 
pars=mfit$estimate; names(pars) = LETTERS[1:length(pars)]; 

minuslogl = function(A,B,C,D,E,F,G,H,I) {
        pars=c(A,B,C,D,E,F,G,H,I) 
        val = -sgt_logLik1(pars); 
        return(val)
}        

start = list(A=pars[1],B=pars[2],C=pars[3],D=pars[4],E=pars[1],F=pars[2],G=pars[3],H=pars[4],I=pars[5])

### plot the results 
pars=mfit$estimate
zx = seq(min(z),max(z),length=100); 
pql = make_pql(zx,pars)
    
par(mfrow=c(2,2),bty="l",mar=c(4,4,1,1));
plot(zx,pql$lambda,xlab="Initial size",ylab="Skew parameter lambda"); 
plot(zx,pql$p,xlab="Initial size",ylab="Shape parameter p"); 
plot(zx,pql$q,xlab="Initial size",ylab="Shape parameter q"); 
    
px=seq(-4,4,length=200); py=matrix(NA,length(px),3); 
pql=make_pql(quantile(z,0.1),pars); 
py[,1] = dsgt(px, mu = 0, sigma = 1, lambda = pql$lambda, p = pql$p, q = pql$q,
mean.cent = TRUE, var.adj = TRUE, log = FALSE)

pql=make_pql(quantile(z,0.5),pars); 
py[,2] = dsgt(px, mu = 0, sigma = 1, lambda = pql$lambda, p = pql$p, q = pql$q,
mean.cent = TRUE, var.adj = TRUE, log = FALSE)

pql=make_pql(quantile(z,0.9),pars); 
py[,3] = dsgt(px, mu = 0, sigma = 1, lambda = pql$lambda, p = pql$p, q = pql$q,
mean.cent = TRUE, var.adj = TRUE, log = FALSE)

matplot(px,py, type="l",col=c("red","black","blue"),lty=1,xlab="Size t", ylab="Distribution of size t+1"); 



