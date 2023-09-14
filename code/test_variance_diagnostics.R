########################################################
# Try out multiple_levene_test and multiple_bs_test 
########################################################

### move to the right local directory 
tom = "C:/Users/tm9/Dropbox/github/IPM_size_transitions"
steve = "c:/repos/IPM_size_transitions" 
home = ifelse(Sys.info()["user"] == "Ellner", steve, tom)
setwd(home); setwd("code"); 

source("variance_diagnostics.R"); 
	
stopCluster(c1); 
c1<- makeCluster(8); 
registerDoParallel(c1);

N = 1000; ## number of fitted values 
nreps = 200; ## number of replicate simulations 
R = 1000; ## number of randomizations for randomization test  

pbin = pspline = numeric(nreps); 
for(jrep in 1:nreps){
	cat("Rep ", jrep, "------------------------------", "\n")
	fitted_vals = sort(2*rbeta(N,5,3)); 
	sd_vals = rep(1,N); 
	sd_vals = 1 + 0.2*sin(4*pi*fitted_vals);  
	scaled_resids = rSST(N,mu=rep(0,N), sigma = sd_vals, nu = exp(-2 + 2*fitted_vals), tau = 5); 
	out = multiple_levene_test(fitted_vals, scaled_resids, 2, 6, R) 
	pbin[jrep]=out$p_value
	out = multiple_bs_test(fitted_vals, scaled_resids, 4, 8, R) 
	pspline[jrep]=out$p_value
}

stopCluster(c1); 
par(mfrow=c(2,1));
hist(pbin,main=mean(pbin<0.05)); 
hist(pspline,main=mean(pspline<0.05)); 
