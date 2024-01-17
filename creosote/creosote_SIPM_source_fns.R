



# SIPM --------------------------------------------------------------------


boot.on = FALSE
boot.num <- 1
boot.switch <- TRUE
ran.seeds <- sample.int(100000,size=boot.num)
# "01_SeedVelocities"
# Calculate terminal velocities of seeds and fit to distributions
source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/01_SeedVelocities.R")
# "02_WindSpeeds"
# Load in wind speeds and fit to distributions
source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/02_WindSpeeds.R")
# "03_Dispersal"
# Construct dispersal kernel functions for seeds
source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/03_Dispersal.R")
# "06_BootRes"
# Run resampling subroutine for wind speeds, terminal velocities, and demography
source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/06_BootRes.R")


## Function to compute the WALD mgf
WALDmgf <- function(s,nu,lambda) {
  t1 <- (lambda/nu) 
  t2 <- 2*(nu^2)*s/lambda 
  mgf <- exp(t1*(1-sqrt(1-t2)))
  return(mgf)
}    

## Function to compute the marginalize WALD mgf
margWALDmgf <- function(s,nu,lambda) {
  (1/pi)*integrate(function(q) WALDmgf(s*cos(q),nu,lambda),0,pi)$value
}
WALD_par <- function(h=0.15,elas="none",sens="none",mat){
  
  # Fit equation to convert volume to height for dispersal kernel use
  LATR_full %>%
    dplyr::select(max.ht_t, volume_t) %>% 
    drop_na(max.ht_t, volume_t) %>% 
    rename("h" = max.ht_t, "v" = volume_t) %>% 
    arrange(v) %>% 
    nlsLM(h ~ A*v^(1/3),
          start = list(A = 0), data = .) %>% 
    coef() %>% 
    as.numeric() -> A
  
  # Function converting volume to height (embedded here bc LATR_full will change with bootstrap iterations)
  # returns height in centimeters
  vol.to.height <- function(v){A*v^(1/3)}
  
  # size vector (log(volume))
  zvals <- mat$meshpts
  ## eviction problem!! if the zval is below true min or above true max, set to true min and true max
  zvals[zvals<LATR_size_bounds$min_size]=LATR_size_bounds$min_size
  zvals[zvals>=LATR_size_bounds$max_size]=LATR_size_bounds$max_size
  
  # Vector of heights across which dispersal kernel will be evaluated
  heights <- sapply(exp(zvals), vol.to.height)/100
  WALD.par <- vector("list",length(heights))
  WALD.par[heights>=h] <- lapply(heights[heights>=h],WALD.b.tom,elas=elas,sens=sens)
  
  return(list(heights=heights,WALD.par=WALD.par))
}
WALD_samples<-function(N,h=0.15,elas="none",sens="none",seed=NULL,params){
  r=matrix(0,nrow=N,ncol=length(params$heights))
  r[,params$heights>h]=sapply(params$heights[params$heights>h],WALD.f.e.h.tom,n=N,elas=elas,sens=sens,seed=seed)
  alpha <- matrix(runif(N*length(params$heights),0,2*pi),dim(r))
  X=r*cos(alpha)
  return(X)
}
## Function to compute wave speeds c(s)
cs <- function(s,h=0.15,emp=F,params,D.samples,mat) {
  # survival-growth matrix
  P <- mat$Pmat 
  # fertility matrix
  Fs <- mat$Fmat 
  for(j in 1:length(params$heights)){
    if(params$heights[j]<h){next}
    if(emp==F){Fs[,j] <- Fs[,j]*margWALDmgf(s,nu=params$WALD.par[[j]]$nu,lambda=params$WALD.par[[j]]$lambda)}
    if(emp==T){Fs[,j] <- Fs[,j]*empiricalWALDmgf(s,D.samples[,j])}
  }
  Hs <- P+Fs 
  L1 = abs(eigen(Hs)$values[1]); 
  return((1/s)*log(L1))  #
}
