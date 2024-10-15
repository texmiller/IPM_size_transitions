# function for converting cactus size measurements to volume
volume <- function(h, w, p){
  (1/3)*pi*h*(((w + p)/2)/2)^2
}
invlogit<-function(x){exp(x)/(1+exp(x))}


## IPM source functions
## SURVIVAL
sx<-function(x,exclude,year=2004,plot=1){
  xb=pmin(pmax(x,min.size),max.size)
  pred=predict(surv_mod,newdata = data.frame(logvol_t=xb,plot=plot,year_t=year),
               exclude=paste0(exclude))
  return(invlogit(pred))
}
## GROWTH - SHASH
gxy_SHASH<-function(x,y,exclude,year=2004,plot=1){
  xb=pmin(pmax(x,min.size),max.size) #Transforms all values below/above limits in min/max size
  pred=predict(grow_SHASH,
               newdata = data.frame(logvol_t=xb,plot=plot,year_t=year),
               exclude=paste0(exclude))
  return(dSHASHo2(x=y, 
                  mu=pred[,1],
                  sigma = exp(pred[,2]), 
                  nu = pred[,3], 
                  tau = exp(pred[,4])))
}
## GROWTH - Gaussian
gxy_GAU<-function(x,y,exclude,year=2004,plot=1){
  xb=pmin(pmax(x,min.size),max.size) #Transforms all values below/above limits in min/max size
  pred = predict(grow_GAU,
                 newdata = data.frame(logvol_t=xb,plot=plot,year_t=year),
                 exclude=paste0(exclude))
  return(dnorm(y,mean=pred[,1],sd=exp(pred[,2])))
}
## COMBINED GROWTH_SURVIVAL
pxy <- function(x,y,dist,exclude,year=2004,plot=1){
  result <- sx(x,exclude,year,plot)*do.call(paste0("gxy_",dist),list(x,y,exclude,year,plot))
  return(result)
}
#PR FLOWERING
flow.x <- function(x,exclude,year=2004,plot=1){
  xb=pmin(pmax(x,min.size),max.size)
  pred=predict(flow_mod,
               newdata = data.frame(logvol_t=xb,plot=plot,year_t=year),
               exclude=paste0(exclude))
  return(invlogit(pred))
}
##FLOWERBUD PRODUCTION BY FLOWERING PLANTS
fert.x <- function(x,exclude,year=2004,plot=1){
  xb=pmin(pmax(x,min.size),max.size)
  pred=predict(fert_mod,
               newdata = data.frame(logvol_t=xb,plot=plot,year_t=year),
               exclude=paste0(exclude))
  return(exp(pred))
}
##SEED BANK CONTRIBUTION OF X-SIZED PLANTS
fx<-function(x,exclude,year=2004,plot=1){
  return(flow.x(x,exclude,year,plot)*fert.x(x,exclude,year,plot)*seeds_per_fruit$seeds_per_fruit*seed_survival$seed_survival)  
}
#SIZE DISTRIBUTION OF RECRUITS
recruit.size<-function(y){
  dnorm(x=y,mean=seedling_size$mean_size,sd=seedling_size$sd_size)
}

#PUT IT ALL TOGETHER
bigmatrix<-function(lower.extension = 0,upper.extension = 0,
                    mat.size,dist,exclude,year=2004,plot=1){
  n<-mat.size
  L<-min.size + lower.extension
  U<-max.size + upper.extension
  #these are the upper and lower integration limits
  h<-(U-L)/n                   #Bin size
  b<-L+c(0:n)*h;               #Lower boundaries of bins 
  y<-0.5*(b[1:n]+b[2:(n+1)]);  #Bins' midpoints
  # Fertility matrix
  Fmat<-matrix(0,(n+2),(n+2))
  # 1-yo banked seeds go in top row
  Fmat[1,3:(n+2)]<-fx(y,exclude,year,plot)
  # Growth/survival transition matrix
  Tmat<-matrix(0,(n+2),(n+2))
  # Graduation to 2-yo seed bank = pr(not germinating as 1-yo)
  Tmat[2,1]<-(1-germination$germ1)
  # Graduation from 1-yo bank to cts size = germination * size distn * pre-census survival
  Tmat[3:(n+2),1]<- germination$germ1 * precensus_survival$precensus_survival * recruit.size(y) * h   
  # Graduation from 2-yo bank to cts size = germination * size distn * pre-census survival
  Tmat[3:(n+2),2]<- germination$germ2 * precensus_survival$precensus_survival * recruit.size(y) * h  
  # Growth/survival transitions among cts sizes
  Tmat[3:(n+2),3:(n+2)]<-t(outer(y,y,pxy,dist=dist,exclude=exclude,year=year,plot=plot)) * h 
  # Put it all together
  IPMmat<-Fmat+Tmat     #Full Kernel is simply a summation ot fertility
  #and transition matrix
  return(list(IPMmat=IPMmat,Fmat=Fmat,Tmat=Tmat,meshpts=y))
}