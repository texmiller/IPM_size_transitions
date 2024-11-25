############################################################################################# 
## Code for BCa bootstrap confidence interval, based on code from the "coxed" R package 
## @Manual{,
##    title = {coxed: Duration-Based Quantities of Interest for the Cox Proportional Hazards Model},
##    author = {Kropko, Jonathan and Harden, {Jeffrey J.}},
##    year = {2020},
##    note = {R package version 0.3.3},
##    url = {https://CRAN.R-project.org/package=coxed},
## }
## which in turn was adapted from code in the "mediation" R package. 
## The "coxed" package was released with a GPL-2 license, which permits copying of the
## source code, as we have done here. Further use of the code, by yourself, is governed
## by the terms of the GPL-2 license
##
## This scripts corrects several important errors in the coxed and mediation routines.
##  - the bias-correction factor z0 is correctly defined, using theta_hat instead of mean(theta)
##  - the acceleration parameter a is an input parameter. The code in coxed and mediation
##		calculates a from bootstrap parameter estimates, which is incorrect. The simplest
##		approximation used jackknifed parameter estimates. 
##
## If given a valid estimate of a, this code implements the BCa interval as described in
##   Thomas J. DiCiccio and Bradley Efron (1996). Bootstrap Confidence Intervals. 
##   Statistical Science 11(3): 189–228. 
##   Bradley Efron and Robert J. Tibshirani (1994) An Introduction to the Bootstrap. 
##	 Chapman & Hall/CRC, New York
## 
## If the acceleration parameter a is estimated from jackknifed parameter estimates
## (eqn. (14.15), pg. 186 in Efron & Tibshirani) then it is exactly equivalent to
## bcanon() in the bootstrap package. 
## @Manual{,
##    title = {bootstrap: Functions for the Book "An Introduction to the Bootstrap"},
##    author = {S original and from StatLib and by Rob Tibshirani. R port by Friedrich Leisch.},
##    year = {2019},
##    note = {R package version 2019.6},
##    url = {https://CRAN.R-project.org/package=bootstrap},
## }
##
## With the default value a=0, the resulting intervals are bias-corrected by not "accelerated"
## (i.e., there is no adjustment for skewness) and therefore not second-order accurate.  
################################################################################################
bca = function (theta, theta_hat, a=0, conf.level = 0.95) {
    low <- (1 - conf.level)/2
    high <- 1 - low
    sims <- length(theta)
    z.inv <- mean(theta < as.numeric(theta_hat)); 
    z <- qnorm(z.inv)
    lower.inv <- pnorm(z + (z + qnorm(low))/(1 - a * (z + qnorm(low))))
    lower <- quantile(theta, lower.inv, names = FALSE)
    upper.inv <- pnorm(z + (z + qnorm(high))/(1 - a * (z + qnorm(high))))
    upper <- quantile(theta, upper.inv, names = FALSE)
    return(c(lower, upper))
}
