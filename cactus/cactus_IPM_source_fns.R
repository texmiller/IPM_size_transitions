invlogit<-function(x){exp(x)/(1+exp(x))}

# VITAL RATE FUNCTIONS ----------------------------------------------------
## GROWTH - skewed t
gxy_ST<-function(x,y,params){
  xb=pmin(pmax(x,params$min.size),params$max.size) #Transforms all values below/above limits in min/max size
  grow_mu <- params$grow.mu + params$grow.bsize * xb
  return(dST3(x=y, 
              mu=mu_size,
              sigma = exp(params$sigma_b0 + params$sigma_b1*mu_size + params$sigma_b2*mu_size^2), 
              nu = exp(params$nu_b0 + params$nu_b1*mu_size), 
              tau = exp(params$tau_b0 + params$tau_b1*mu_size + params$tau_b2*mu_size^2)))
  
}

## GROWTH - Gaussian
gxy_norm<-function(x,y,params){
  xb=pmin(pmax(x,params$min.size),params$max.size) #Transforms all values below/above limits in min/max size
  grow_mu <- params$grow.mu.norm + params$grow.bsize.norm * xb
  grow_sd <- exp(params$grow.sd.b0 + params$grow.sd.b1*grow_mu)
  return(dnorm(y,mean=grow_mu,sd=growth_sd))
    
}

## SURVIVAL
sx<-function(x,params){
  xb=pmin(pmax(x,params$min.size),params$max.size)
  p.surv<-params$surv.mu + params$surv.bsize*xb
  return(invlogit(p.surv))
}

## COMBINED GROWTH_SURVIVAL
pxy <- function(x,y,params,dist){
  result <- ifelse(dist=="ST",sx(x,params)*gxy_ST(x,y,params),sx(x,params)*gxy_norm(x,y,params))
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