

# IPM functions -----------------------------------------------------------

# Growth from size x to y at density d, using best GAUSSIAN
growth.GAU <- function(x, y, d){
  xb = pmin(pmax(x, LATR_size_bounds$min_size), LATR_size_bounds$max_size)
  mu = predict(LATR_GAU_best,
               newdata = data.frame(dens_scaled = d, log_volume_t = xb, unique.transect = "1.FPS"),
               re.form = NA)
  sigma = exp(GAU_sd_coef$estimate[1]+GAU_sd_coef$estimate[2]*mu)
  return(dnorm(y,mu,sigma))
}

# Growth from size x to y at density d, using JSU
growth.JSU <- function(x, y, d){
  xb = pmin(pmax(x, LATR_size_bounds$min_size), LATR_size_bounds$max_size)
  mu = predict(LATR_GAU_best,
               newdata = data.frame(dens_scaled = d, log_volume_t = xb, unique.transect = "1.FPS"),
               re.form = NA)
  return(dJSU(y,
              mu=mu,
              sigma=exp(GAU_sd_coef$estimate[1]+GAU_sd_coef$estimate[2]*mu),
              nu=JSUout$estimate[3]+JSUout$estimate[4]*mu,
              tau=exp(JSUout$estimate[5]+JSUout$estimate[6]*mu)))
}

# Survival of size x at density d using best GAM
# For nnaturally occuring plants (transplant = FALSE)
survival <- function(x, d){
  xb = pmin(pmax(x, LATR_size_bounds$min_size), LATR_size_bounds$max_size)
  pred <- predict.gam(LATR_surv_best,
                      newdata = data.frame(dens_scaled = d, log_volume_t = xb, transplant = FALSE,
                                           unique.transect = "1.FPS"),
                      type = "response",
                      exclude = "s(unique.transect)")
  return(pred)
}

# Combined growth and survival at density d
growsurv <- function(x, y, d, dist){
  survival(x, d) * do.call(paste0("growth.",dist),list(x, y, d))
}

# Flowering at size x and density d using best GAM
flower <- function(x, d){
  xb = pmin(pmax(x, LATR_size_bounds$min_size), LATR_size_bounds$max_size)
  pred <- predict.gam(LATR_flower_best,
                      newdata = data.frame(dens_scaled = d, log_volume_t = xb, unique.transect = "1.FPS"),
                      type = "response",
                      exclude = "s(unique.transect)")
  return(pred)
}

# Seed production (fruits * seeds/fruit) at size x and density d using best GAM
seeds <- function(x, d, seeds.per.fruit = 5){
  xb = pmin(pmax(x, LATR_size_bounds$min_size), LATR_size_bounds$max_size)
  pred <- predict.gam(LATR_fruits_best,
                      newdata = data.frame(dens_scaled = d, log_volume_t = xb, unique.transect = "1.FPS"),
                      type = "response",
                      exclude = "s(unique.transect)")
  return(pred*seeds.per.fruit)
}

# Seed-to-Seedling recruitment probability at density d
recruitment <- function(d){
  pred <- predict.gam(LATR_recruit_best,
                      newdata = data.frame(dens_scaled = d, unique.transect = "1.FPS"),
                      type = "response",
                      exclude = "s(unique.transect)")
  return(pred[1])
}

# Recruit size distribution at size y
recruitsize <- function(y,d){
  lpmat <- predict.gam(LATR_recruitsize_best,
                       newdata = data.frame(dens_scaled = d, unique.transect = "1.FPS"),
                       type = "lpmatrix",
                       exclude = "s(unique.transect)")
  recruitsize_mu <- lpmat[, 1:(recruit_size_sd_index-1)] %*% coef(LATR_recruitsize_best)[1:(recruit_size_sd_index-1)]
  recruitsize_sigma <- exp(lpmat[, recruit_size_sd_index:recruit_size_coef_length] %*% coef(LATR_recruitsize_best)[recruit_size_sd_index:recruit_size_coef_length])
  return(dnorm(x = y, mean = recruitsize_mu, sd = recruitsize_sigma))
}

# Combined flowering, fertility, and recruitment
fertrecruit <- function(x, y, d){
  flower(x,d)*seeds(x,d)*recruitment(d)*recruitsize(y,d)
}

# Put it all together; projection matrix is a function of weighted density (dens)
# We need a large lower extension because growth variance (gaussian) is greater for smaller plants
ApproxMatrix <- function(dens,ext.lower=lower.extension,ext.upper=upper.extension,
                         min.size=LATR_size_bounds$min_size,max.size=LATR_size_bounds$max_size,
                         mat.size=matdim,dist){
  
  # Matrix size and size extensions (upper and lower integration limits)
  n <- mat.size
  L <- min.size + ext.lower
  U <- max.size + ext.upper
  
  # Bin size for n bins
  h <- (U - L)/n
  
  # Lower boundaries of bins 
  b <- L + c(0:n)*h
  
  # Bin midpoints
  y <- 0.5*(b[1:n] + b[2:(n + 1)])
  
  # Growth/Survival matrix
  Pmat <- t(outer(y, y, growsurv, d = dens, dist=dist)) * h 
  
  # Fertility/Recruiment matrix
  Fmat <- t(outer(y, y, fertrecruit, d = dens)) * h 
  
  # Put it all together
  IPMmat <- Pmat + Fmat
  
  #and transition matrix
  return(list(IPMmat = IPMmat, Fmat = Fmat, Pmat = Pmat, meshpts = y))
}

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
