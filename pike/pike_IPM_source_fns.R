# IPM ---------------------------------------------------------------------
## GROWTH - SHASH
gxy_SHASH<-function(x,y){
  xb=pmin(pmax(x,min(pike_final$log_t0)),max(pike_final$log_t1)) 
  pred=predict(pike_shash,
               newdata = data.frame(log_t0=xb))
  return(dSHASHo2(x=y, 
                  mu=pred[,1],
                  sigma = exp(pred[,2]), 
                  nu = pred[,3], 
                  tau = exp(pred[,4])))
}

## GROWTH - Gaussian
gxy_GAU<-function(x,y){
  xb=pmin(pmax(x,min(pike_final$log_t0)),max(pike_final$log_t1)) 
  pred = predict(pike_gau,newdata = data.frame(log_t0=xb))
  return(dnorm(y,mean=pred[,1],sd=exp(pred[,2])))
}

## SURVIVAL
sx<-function(x){
  xb=pmin(pmax(x,min(pike_final$log_t0)),max(pike_final$log_t1)) 
  pred = predict(pike_surv_gam,newdata = data.frame(log_t0=xb),type="response")
  return(pred)
}

## COMBINED GROWTH_SURVIVAL
pxy <- function(x,y,dist){
  result <- sx(x)*do.call(paste0("gxy_",dist),list(x,y))
  return(result)
}

## FERTILITY
## these parameter values come from Table 2
r<-0.9 #fertilization probability
q<-0.5 #fraction female
S1<-0.000623 #survival from fertilized egg to 1yo

fx<-function(x){
  xb=pmin(pmax(x,min(pike_final$log_t0)),max(pike_final$log_t1)) 
  eggs=predict(pike_fert_gam,newdata = data.frame(log_t0=xb),type="response")
  return(eggs*r*q*S1)  
}

#SIZE DISTRIBUTION OF RECRUITS
#recruit.size<-function(y){
#  dnorm(x=y,
#        mean=pike_recruitsize_gau$coefficients[1],
#        sd=exp(pike_recruitsize_gau$coefficients[2]))
#}

recruit.size<-function(y){
  dSHASHo2(y,
           pike_recruitsize_shash$coefficients[1],
           exp(pike_recruitsize_shash$coefficients[2]),
           pike_recruitsize_shash$coefficients[3],
           exp(pike_recruitsize_shash$coefficients[4]),
  )
}

##COMBINED FERTILITY/RECRUITMENT
fxy<-function(x,y){
  return(fx(x)*recruit.size(y))
}

##PUT IT ALL TOGETHER
## defaults here come from experimentation in the basement
ApproxMatrix <- function(lower,upper,ext.lower=0,ext.upper=0.1,mat.size=800,dist){
  # Matrix size and size extensions (upper and lower integration limits)
  n <- mat.size
  L <- lower - ext.lower
  U <- upper + ext.upper
  # Bin size for n bins
  h <- (U - L)/n
  # Lower boundaries of bins 
  b <- L + c(0:n)*h
  # Bin midpoints
  y <- 0.5*(b[1:n] + b[2:(n + 1)])
  
  # Growth/Survival matrix
  Pmat <- t(outer(y, y, pxy, dist=dist)) * h 
  # Fertility/Recruitment matrix
  Fmat <- t(outer(y, y, fxy)) * h 
  # Put it all together
  IPMmat <- Pmat + Fmat
  
  return(list(IPMmat = IPMmat, Fmat = Fmat, Pmat = Pmat, meshpts = y))
}
