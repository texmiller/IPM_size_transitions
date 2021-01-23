library(mgcv); 
n_obs<-100
n_levels<-10
x1 = runif(n_obs,-2,2)
rfx<-rnorm(n_levels,0,3)
mu = rep(rfx,(n_obs/n_levels)) + 2.3*x1 - 0.5*x1^2 + 1.9*x1^3
y = rnorm(n_obs,mu,5)
practice.dat <- data.frame(x1=x1,
                           rfx=as.factor(rep(letters[1:n_levels],(n_obs/n_levels))),
                           y=y)

testgam <- gam(y~s(x1)+s(rfx,bs="re"),data=practice.dat)
Xp <- predict.gam(testgam,type = "lpmatrix")[,2:10]  # first column is the intercept term (all 1's) 

B <- smoothCon(s(x1,k=10),data=practice.dat,absorb.cons=TRUE)[[1]]$X; 

plot(B,Xp); # identical 

 
 


