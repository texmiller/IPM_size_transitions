############################################################################################# 
## Code for BCa bootstrap confidence interval, copied from the "coxed" R package 
## @Manual{,
##    title = {coxed: Duration-Based Quantities of Interest for the Cox Proportional Hazards Model},
##    author = {Kropko, Jonathan and Harden, {Jeffrey J.}},
##    year = {2020},
##    note = {R package version 0.3.3},
##    url = {https://CRAN.R-project.org/package=coxed},
## }
## which in turn was adapted from code in the "mediation" R package. 
##
## The "coxed" package was released with a GPL-2 license, which permits copying of the
## source code, as we have done here. Further use of the code, by yourself, is governed
## by the terms of the GPL-2 license. 
##
## This code implements the simplest form of the BCa interval, as described in 
##   Thomas J. DiCiccio and Bradley Efron (1996). Bootstrap Confidence Intervals. 
##   Statistical Science 11(3): 189–228. 
##
## It is exactly equivalent to bcanon() in the bootstrap package. 
## @Manual{,
##    title = {bootstrap: Functions for the Book "An Introduction to the Bootstrap"},
##    author = {S original and from StatLib and by Rob Tibshirani. R port by Friedrich Leisch.},
##    year = {2019},
##    note = {R package version 2019.6},
##    url = {https://CRAN.R-project.org/package=bootstrap},
## }
################################################################################################
bca = function (theta, theta_hat, conf.level = 0.95) {
    low <- (1 - conf.level)/2
    high <- 1 - low
    sims <- length(theta)
    z.inv <- length(theta[theta < theta_hat])/sims
    z <- qnorm(z.inv)
    U <- (sims - 1) * (mean(theta, na.rm = TRUE) - theta)
    top <- sum(U^3)
    under <- 6 * (sum(U^2))^{3/2}
    a <- top/under
    lower.inv <- pnorm(z + (z + qnorm(low))/(1 - a * (z + qnorm(low))))
    lower <- quantile(theta, lower.inv, names = FALSE)
    upper.inv <- pnorm(z + (z + qnorm(high))/(1 - a * (z + qnorm(high))))
    upper <- quantile(theta, upper.inv, names = FALSE)
    return(c(lower, upper))
}
