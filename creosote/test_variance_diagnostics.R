### move to the right local directory 
tom = "C:/Users/tm9/Dropbox/github/IPM_size_transitions"
steve = "c:/repos/IPM_size_transitions" 
home = ifelse(Sys.info()["user"] == "Ellner", steve, tom)
setwd(home); setwd("creosote"); 

library(lme4)
library(mgcv)
library(maxLik)
library(bbmle)
library(car); 
library(gamlss.dist); 

N = 2000; 
## sd_vals = 1 + 0.1*sin(4*pi*fitted_vals);  
## with Gaussian^1.25 residuals
# mean(pbin<0.05); [1] 0.76
# mean(pspline<0.05); [1] 0.35

## sd_vals = rep(1,N); 
## with Gaussian^1.25 residuals
## mean(pspline<0.05); [1] 0.01
## mean(pbin<0.05); [1] 0.02

nreps = 25; pbin = pspline = numeric(100); 
for(jrep in 1:nreps){
cat("Rep ", jrep, "------------------------------", "\n")

#fitted_vals = log(rexp(N,rate=1)); fitted_vals=sort(fitted_vals); 
#fitted_vals = 2*(fitted_vals-min(fitted_vals))/(max(fitted_vals)-min(fitted_vals)); 
fitted_vals = sort(2*rbeta(N,5,3)); 
# sd_vals = 1 + 0.1*sin(5*pi*fitted_vals);  
# resids = rnorm(N); resids=resids*(abs(scaled_resids)^0.25); 
# scaled_resids = (resids/sd(resids))*sd_vals; 
sd_vals = rep(1,N); 
scaled_resids = rSST(N,mu=rep(0,N), sigma = sd_vals, nu = exp(-2 + 2*fitted_vals), tau = 5); 
par(mfrow=c(2,2)); 
hist(fitted_vals); 
plot(fitted_vals,sd_vals); 
plot(fitted_vals,scaled_resids); 
plot(1:N,scaled_resids); 

#############################################################
##  Levene Test, binning by percentiles of the data
##  u is a factor variable to indicate bin membership 
##  p-value from randomization of scaled residuals 
#############################################################

## sort the fitted values and residuals 
indexx = c(1:length(fitted_vals)); p_vals = numeric(0); 
for(nbins in 3:10){
    u = nbins*indexx/(1 + max(indexx)); 
    u = floor(u); u = factor(u); unique(u); 
    p_true = leveneTest(scaled_resids,u,center=mean, trim = 0.1)$"Pr(>F)"[[1]]
    p_ran = numeric(1000)
    for(j in 1:1000) {
        p_ran[j] = leveneTest(sample(scaled_resids),u, center=mean, trim = 0.1)$"Pr(>F)"[[1]]
    }
    hist(p_ran);  abline(v=p_true);
    vars = numeric(nbins); for(k in 1:nbins) vars[k]=var(scaled_resids[u==(k-1)]); 
    cat(nbins, mean(p_ran < p_true), "    ", signif(vars,3),"\n"); 
    p_vals = c(p_vals,mean(p_ran < p_true)); 
}
cat(jrep, p.adjust(p_vals,method="hommel"),"\n"); 
pbin[jrep]=min(p.adjust(p_vals,method="hommel"))

############################################################
## Breusch-Pagan-like test using spline covariates 
## to test for missed patterns. Randomization p-value 
## based on shuffling scaled residuals.  
############################################################
library(splines); 

p_vals = numeric(0); 
for(df in 4:10) {
    X = bs(fitted_vals,df=df,intercept=TRUE); 
    matplot(fitted_vals,X,type="l",lty=1)

    fit_true = lm(I(scaled_resids^2) ~ X-1) 
    rsq_true = summary(fit_true)$r.squared; 

    rsq_ran = numeric(1000); 
    for(j in 1:1000) {
        yj = sample(scaled_resids); 
        fit_ran = lm(I(yj^2) ~ X-1) 
        rsq_ran[j] = summary(fit_ran)$r.squared; 
    }    
    hist(log(rsq_ran)); abline(v=log(rsq_true));
    cat(df, mean(rsq_ran>rsq_true), "\n"); 
    p_vals = c(p_vals,mean(rsq_ran > rsq_true)); 
}
cat(jrep,p.adjust(p_vals,method="hommel"),"\n"); 
pspline[jrep]=min(p.adjust(p_vals,method="hommel"))
} 



