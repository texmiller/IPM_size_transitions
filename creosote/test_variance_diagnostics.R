#############################################################
##  Levene Test, binning by percentiles of the data
##  u is a factor variable to indicate bin membership 
##  p-value from randomization of scaled residuals 
#############################################################

### move to the right local directory 
tom = "C:/Users/tm9/Dropbox/github/IPM_size_transitions"
steve = "c:/repos/IPM_size_transitions" 
home = ifelse(Sys.info()["user"] == "Ellner", steve, tom)
setwd(home); setwd("creosote"); 

library(car); 
library(gamlss.dist); 
require(parallel); 
require(doParallel); 

#################################################################
## Compute minimum Levene test p-value on data and R random 
## permutations, scanning across a range of bin numbers
#################################################################
multiple_levene_test = function(fitted_vals, residuals, min_bins, max_bins, R) {
	e = order(fitted_vals); residuals = residuals[e];  
	indexx = seq_along(residuals)
	p_true=rep(NA,max_bins); 
	for(nbins in min_bins:max_bins){
		u = nbins*indexx/(1 + max(indexx)); 
		u = floor(u); u = factor(u); 
		p_true[nbins] = leveneTest(residuals,u,center=mean, trim = 0.1)$"Pr(>F)"[[1]]
	}
	p_min_true = min(p_true,na.rm=TRUE); 
	bins_min_true = which.min(p_true); 

	 
	out = foreach(j=1:R,.combine = c,.packages="car") %dopar% 
    {
	ran_resids = sample(residuals);
	p_ran = rep(NA,max_bins); 
	for(nbins in min_bins:max_bins) {
	u = nbins*indexx/(1 + max(indexx)); 
	u = floor(u); u = factor(u);  
	p_ran[nbins] = leveneTest(ran_resids, u, center=mean, trim = 0.1)$"Pr(>F)"[[1]]
	}	
	min(p_ran,na.rm=TRUE); 
}
	p_value = mean(out < p_min_true); ## how many randomizations exceed the true value? 
	return(list(p_value = p_value, p_min_true = p_min_true, bins_min_true=bins_min_true, p_min_random = out))
} 
	
#################################################################
## Compute maximum b-spline r^2 value on data and R random 
## permutations, scanning across a range of basis functions
#################################################################
multiple_bs_test = function(fitted_vals, residuals, min_basis, max_basis, R) {
	e = order(fitted_vals); residuals = residuals[e];  
	indexx = seq_along(residuals)
	rsq_true=rep(NA,max_basis); 
	for(nbasis in min_basis:max_basis){
		X = bs(fitted_vals,df=nbasis,intercept=TRUE); 
		fit_true = lm(I(residuals^2) ~ X-1) 
		rsq_true[nbasis] = summary(fit_true)$r.squared; 
	}
	rsq_max_true = max(rsq_true,na.rm=TRUE); 
	df_max_true = which.max(rsq_true); 

	 
	out = foreach(j=1:R,.combine = c,.packages="splines") %dopar% 
    {
	ran_resids = sample(residuals);
	rsq_ran = rep(NA,max_basis); 
	for(nbasis in min_basis:max_basis) {
		X = bs(fitted_vals,df=nbasis,intercept=TRUE); 
		fit_ran = lm(I(ran_resids^2) ~ X-1) 
		rsq_ran[nbasis] = summary(fit_ran)$r.squared; 
	}
	max(rsq_ran,na.rm=TRUE); 
}
	p_value = mean(out > rsq_max_true); ## how extreme is the true value? 
	return(list(p_value = p_value, rsq_max_true = rsq_max_true, df_max_true=df_max_true, rsq_max_random = out))
} 

if(TESTING) {
########################################################
# Try out multiple_levene_test and multiple_bs_test 
########################################################	
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
}