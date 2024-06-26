invlogit<-function(x){exp(x)/(1+exp(x))}

# VITAL RATE FUNCTIONS ----------------------------------------------------

## GROWTH - SHASH
gxy_SHASH<-function(x,y,params,plot,year){
  xb=pmin(pmax(x,params$min.size),params$max.size) #Transforms all values below/above limits in min/max size
  pred=predict(CYIM_gam_shash,
               newdata = data.frame(logvol_t=xb,plot=plot,year_t=year))
  return(dSHASHo2(x=y, 
              mu=pred[,1],
              sigma = exp(pred[,2]), 
              nu = pred[,3], 
              tau = exp(pred[,4])))
}

## GROWTH - Gaussian
gxy_GAU<-function(x,y,params){
  xb=pmin(pmax(x,params$min.size),params$max.size) #Transforms all values below/above limits in min/max size
  pred = predict(CYIM_grow_m1,
               newdata = data.frame(logvol_t=xb,plot="1",year_t="2004"),
               exclude=c("s(plot)","s(year_t)"))
  return(dnorm(y,mean=pred[,1],sd=exp(pred[,2])))
}

## SURVIVAL
sx<-function(x,params){
  xb=pmin(pmax(x,params$min.size),params$max.size)
  p.surv<-params$surv.mu + params$surv.bsize*xb
  return(invlogit(p.surv))
}

## COMBINED GROWTH_SURVIVAL
pxy <- function(x,y,params,dist){
  result <- sx(x,params)*do.call(paste0("gxy_",dist),list(x,y,params))
  return(result)
}

#PRODUCTION OF 1-YO SEEDS IN THE SEED BANK FROM X-SIZED MOMS
flow.x <- function(x,params){
  xb=pmin(pmax(x,params$min.size),params$max.size)
  p.flow<-params$flow.mu + params$flow.bsize*xb
  return(invlogit(p.flow))
}

fert.x <- function(x,params){
  xb=pmin(pmax(x,params$min.size),params$max.size)
  nfruits<-params$fert.mu + params$fert.bsize*xb
  return(exp(nfruits))
}

fx<-function(x,params){
  return(flow.x(x,params)*fert.x(x,params)*params$mu_spf*params$seedsurv)  
}

#SIZE DISTRIBUTION OF RECRUITS
recruit.size<-function(y,params){
  dnorm(x=y,mean=params$mu_sdlgsize,sd=params$sigma_sdlgsize)
}

#PUT IT ALL TOGETHER
bigmatrix<-function(params,
                    lower.extension = 0, ## I'll need to extend lower and upper beyond true size limits
                    upper.extension = 0,
                    mat.size, ## matrix dimensions
                    dist){
  
  n<-mat.size
  L<-params$min.size + lower.extension
  U<-params$max.size + upper.extension
  #these are the upper and lower integration limits
  h<-(U-L)/n                   #Bin size
  b<-L+c(0:n)*h;               #Lower boundaries of bins 
  y<-0.5*(b[1:n]+b[2:(n+1)]);  #Bins' midpoints
  
  # Fertility matrix
  Fmat<-matrix(0,(n+2),(n+2))
  
  # 1-yo banked seeds go in top row
  Fmat[1,3:(n+2)]<-fx(y,params)
  
  # Growth/survival transition matrix
  Tmat<-matrix(0,(n+2),(n+2))
  
  # Graduation to 2-yo seed bank = pr(not germinating as 1-yo)
  Tmat[2,1]<-(1-params$germ1)
  
  # Graduation from 1-yo bank to cts size = germination * size distn * pre-census survival
  Tmat[3:(n+2),1]<- params$germ1 * params$precenus_surv * recruit.size(y,params) * h   
  
  # Graduation from 2-yo bank to cts size = germination * size distn * pre-census survival
  Tmat[3:(n+2),2]<- params$germ2 * params$precenus_surv * recruit.size(y,params) * h  
  
  # Growth/survival transitions among cts sizes
  Tmat[3:(n+2),3:(n+2)]<-t(outer(y,y,pxy,params=params,dist=dist)) * h 
  
  # Put it all together
  IPMmat<-Fmat+Tmat     #Full Kernel is simply a summation ot fertility
  #and transition matrix
  return(list(IPMmat=IPMmat,Fmat=Fmat,Tmat=Tmat,meshpts=y))
}

# lambdaS function##########################################################
lambdaSim<-function(mat_list, ## a list of transition matrices, each corresponding to a study year
                    max_yrs=1000 ## how many years the simulation runs (arbitrarily large)
){
  ## grab the dimension of the projection matrix
  matdim<-dim(mat_list[[1]])[1]
  ## grab the number of study years / matrices we have available
  n_years <- length(mat_list)
  ## vector that will hold year-by-year growth rates
  rtracker <- rep(0,max_yrs)
  ## initial vector of population structure -- note that this sums to one, which will be convenient
  n0 <- rep(1/matdim,matdim)
  for(t in 1:max_yrs){ #Start loop
    ## for each year, randomly sample one of the matrices
    A_t <- mat_list[[sample.int(n=n_years,size=1)]]
    ## project the population one step forward
    n0 <- A_t %*% n0
    ## total population size after one year of growth
    N  <- sum(n0)
    ## calculate r as log(N_t+1 / N_t), note that here N_t=1
    rtracker[t]<-log(N)
    ## rescale population vector to sum to one, so the same trick works again next time step
    n0 <-n0/N
  }
  #discard first 10% of time series
  burnin    <- round(max_yrs*0.1)
  #Finish and return
  log_lambdaS <- mean(rtracker[-c(1:burnin)])
  lambdaS<-exp(log_lambdaS)
  return(list(log_lambdaS=log_lambdaS,lambdaS=lambdaS))
}