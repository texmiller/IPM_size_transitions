#####################################################################
## Code for SiZeR-inspired multiple Levene and B-spline tests 
## for nonconstant variance as a function of a covariate
#####################################################################

library(car); 
require(parallel); 
require(doParallel); 
require(splines); 

###################################################################
## Compute minimum Bartlett test p-value on data and R random 
## permutations, scanning across a range of bin numbers.
###################################################################
multiple_bartlett_test = function(fitted_vals, residuals, min_bins, max_bins, R) {
	e = order(fitted_vals); 
	fitted_vals = fitted_vals[e]; residuals = residuals[e];  
	indexx = seq_along(residuals)
	p_true=rep(NA,max_bins); 
	for(nbins in min_bins:max_bins){
		u = nbins*indexx/(1 + max(indexx)); 
		u = floor(u); u = factor(u); 
		p_true[nbins] = bartlett.test(residuals,u)$p.value
	}
	p_min_true = min(p_true,na.rm=TRUE); 
	bins_min_true = which.min(p_true); 
	
	## compute variances for most significant bin number 
	u = bins_min_true*indexx/(1 + max(indexx)); 
	u = floor(u); 
	vars = medians = numeric(max(u)); 
	for(j in 1:length(vars)) {
		e = which(u==j)
		vars[j]=var(residuals[e]);
		medians[j] = median(fitted_vals[e]); 
	}	
	
	## residual randomizations   
	out = foreach(j=1:R,.combine = c) %dopar% 
    {
	ran_resids = sample(residuals);
	p_ran = rep(NA,max_bins); 
	for(nbins in min_bins:max_bins) {
	u = nbins*indexx/(1 + max(indexx)); 
	u = floor(u); u = factor(u);  
	p_ran[nbins] = bartlett.test(ran_resids, u)$p.value
	}	
	min(p_ran,na.rm=TRUE); 
}
	p_value = mean(out < p_min_true); ## how many randomizations exceed the true value? 
	return(list(p_value = p_value, p_min_true = p_min_true, p_true = p_true, 
			bins_min_true=bins_min_true, bin_variances = vars, bin_medians=medians))
} 


#################################################################
## Compute maximum b-spline r^2 value on data and R random 
## permutations, scanning across a range of basis functions
#################################################################
multiple_bs_test = function(fitted_vals, residuals, min_basis, max_basis, R) {
	e = order(fitted_vals); 
	fitted_vals = fitted_vals[e]; residuals = residuals[e]; 
	indexx = seq_along(residuals)
	rsq_true=rep(NA,max_basis); 
	for(nbasis in min_basis:max_basis){
		X = bs(fitted_vals,df=nbasis,intercept=TRUE); 
		fit_true = lm(I(residuals^2) ~ X-1) 
		rsq_true[nbasis] = summary(fit_true)$adj.r.squared; 
	}
	rsq_max_true = max(rsq_true,na.rm=TRUE); 
	df_max_true = which.max(rsq_true); 
	
	###### make the maximum-rsquare trend 
	X = bs(fitted_vals,df=df_max_true,intercept=TRUE); 
	fit_true = lm(I(residuals^2) ~ X-1)
	
	px = seq(min(fitted_vals),max(fitted_vals),length=200); 
	X = bs(px,df=df_max_true,intercept=TRUE); 
	py = X%*%coef(fit_true); 
	
	out = foreach(j=1:R,.combine = c,.packages="splines") %dopar% 
		{
		ran_resids = sample(residuals);
		rsq_ran = rep(NA,max_basis); 
		for(nbasis in min_basis:max_basis) {
			X = bs(fitted_vals,df=nbasis,intercept=TRUE); 
			fit_ran = lm(I(ran_resids^2) ~ X-1) 
			rsq_ran[nbasis] = summary(fit_ran)$adj.r.squared; 
		}
	max(rsq_ran,na.rm=TRUE); 
	}
	p_value = mean(out > rsq_max_true); ## how extreme is the true value? 
	return(list(p_value = p_value, rsq_max_true = rsq_max_true, rsq_true = rsq_true, 
	df_max_true=df_max_true, trend_x = px, trend_y = py))
} 


#################################################################
## Compute maximum b-spline r^2 value on data and R random 
## permutations, scanning across a range of basis functions
#################################################################
multiple_bs_test = function(fitted_vals, residuals, min_basis, max_basis, R) {
	e = order(fitted_vals); 
	fitted_vals = fitted_vals[e]; residuals = residuals[e]; 
	indexx = seq_along(residuals)
	AIC_true=rep(NA,min_basis); 
	for(nbasis in min_basis:min_basis){
		X = bs(fitted_vals,df=nbasis,intercept=TRUE); 
		fit_true = lm(I(residuals^2) ~ X-1) 
		AIC_true[nbasis] = AIC(fit_true)
	}
	AIC_min_true = min(AIC_true,na.rm=TRUE); 
	df_min_true = which.min(AIC_true); 
	
	###### make the minimum-AIC trend 
	X = bs(fitted_vals,df=df_min_true,intercept=TRUE); 
	fit_true = lm(I(residuals^2) ~ X-1)
	
	px = seq(min(fitted_vals),max(fitted_vals),length=200); 
	X = bs(px,df=df_min_true,intercept=TRUE); 
	py = X%*%coef(fit_true); 
	
	out = foreach(j=1:R,.combine = c,.packages="splines") %dopar% 
		{
		ran_resids = sample(residuals);
		AIC_ran = rep(NA,min_basis); 
		for(nbasis in min_basis:min_basis) {
			X = bs(fitted_vals,df=nbasis,intercept=TRUE); 
			fit_ran = lm(I(ran_resids^2) ~ X-1) 
			AIC_ran[nbasis] = AIC(fit_ran)
		}
	min(AIC_ran,na.rm=TRUE); 
	}
	p_value = mean(out < AIC_min_true); ## how extreme is the true value? 
	return(list(p_value = p_value, AIC_min_true = AIC_min_true, AIC_true = AIC_true, 
	df_min_true=df_min_true, trend_x = px, trend_y = py))
} 

###############################################################################
## Compute minimum Levene test p-value on data and R random 
## permutations, scanning across a range of bin numbers.
## This uses the leveneTest() from car package, with 5% trimmed means 
## to give some robustness against outliers. Not used in the paper 
## because it is based on absolute rather than squared residuals, and 
## (more practically) the multiple_bartlett_test seems to perform better
## in test cases (admittedly, just a few). But what the heck, we wrote it
## so we might as well give it to you. 
###############################################################################
multiple_levene_test = function(fitted_vals, residuals, min_bins, max_bins, R) {
	e = order(fitted_vals); 
	fitted_vals = fitted_vals[e]; residuals = residuals[e];   
	indexx = seq_along(residuals)
	p_true=rep(NA,max_bins); 
	for(nbins in min_bins:max_bins){
		u = nbins*indexx/(1 + max(indexx)); 
		u = floor(u); u = factor(u); 
		p_true[nbins] = leveneTest_2(residuals,u,trim=0.05)$"Pr(>F)"[[1]]
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
	p_ran[nbins] = leveneTest_2(ran_resids, u,trim=0.05)$"Pr(>F)"[[1]]
	}	
	min(p_ran,na.rm=TRUE); 
}
	p_value = mean(out < p_min_true); ## how many randomizations exceed the true value? 
	return(list(p_value = p_value, p_min_true = p_min_true, p_true = p_true, bins_min_true=bins_min_true, p_min_random = out))
} 

leveneTest_2(y, group,...){ 
    if (!is.numeric(y)) 
        stop(deparse(substitute(y)), " is not a numeric variable")
    if (!is.factor(group)) {
        warning(deparse(substitute(group)), " coerced to factor.")
        group <- as.factor(group)
    }
    valid <- complete.cases(y, group)
    meds <- tapply(y[valid], group[valid], mean, ...)
    resp <- (y - meds[group])^2
    table <- anova(lm(resp ~ group))[, c(1, 4, 5)]
    rownames(table)[2] <- " "
    dots <- deparse(substitute(...))
    attr(table, "heading") <- paste("Levene's Test for Homogeneity of Variance (center = ", 
        deparse(substitute(center)), if (!(dots == "NULL")) 
            paste(":", dots), ")", sep = "")
    table
}
