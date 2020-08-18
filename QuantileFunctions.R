### A qsn function that actually works, unlike the one in the sn package. 
### Vectorized in p, but not the distribution parameters 
my.qsn = function(p,xi,omega,alpha) {
    px = seq(-20,20,length=250); 
    py = psn(px,omega=omega,alpha=alpha,tau=0);
    F = splinefun(py,px,method="monoH.FC"); 
    return(F(p))
}    

qrSEP1 = function(p,mu,sigma,nu,tau) {
    px = seq(mu-8*sigma,mu+8*sigma,length=150); 
    py = pSEP1(px,mu=mu,sigma=sigma,nu=nu,tau=tau);
    F = splinefun(py,px,method="monoH.FC"); 
    return(F(p))
}    


qrSEP2 = function(p,mu,sigma,nu,tau) {
    px = seq(mu-5*sigma,mu+5*sigma,length=150); 
    py = pSEP2(px,mu=mu,sigma=sigma,nu=nu,tau=tau);
    F = splinefun(py,px,method="monoH.FC"); 
    return(F(p))
}    

qrST1 = function(p,mu,sigma,nu,tau) {
    px = seq(mu-8*sigma,mu+8*sigma,length=150); 
    py = pST1(px,mu=mu,sigma=sigma,nu=nu,tau=tau);
    F = splinefun(py,px,method="monoH.FC"); 
    return(F(p))
}    

qrST2 = function(p,mu,sigma,nu,tau) {
    px = seq(mu-8*sigma,mu+8*sigma,length=150); 
    py = pST2(px,mu=mu,sigma=sigma,nu=nu,tau=tau);
    F = splinefun(py,px,method="monoH.FC"); 
    return(F(p))
}    

qrST3 = function(p,mu,sigma,nu,tau) {
    px = seq(mu-8*sigma,mu+8*sigma,length=150); 
    py = pST3(px,mu=mu,sigma=sigma,nu=nu,tau=tau);
    F = splinefun(py,px,method="monoH.FC"); 
    return(F(p))
}    

qrST5 = function(p,mu,sigma,nu,tau) {
    px = seq(mu-8*sigma,mu+8*sigma,length=150); 
    py = pST5(px,mu=mu,sigma=sigma,nu=nu,tau=tau);
    F = splinefun(py,px,method="monoH.FC"); 
    return(F(p))
}    
