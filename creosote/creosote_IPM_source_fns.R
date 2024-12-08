# IPM functions -----------------------------------------------------------
# using fitted model objects

# Growth from size x to y at density d, using best GAUSSIAN
growth.GAU <- function(x, y, d){
  xb = pmin(pmax(x, LATR_size_bounds$min_size), LATR_size_bounds$max_size)
  ## need expected value before I can use predict.gam because fitted_val is an input
  pred_mat = predict.gam(LATR_GAU_best,
                         newdata = data.frame(dens_scaled = d, log_volume_t = xb, fitted_vals = 0, unique.transect = "1.FPS"),
                         type = "lpmatrix",
                         exclude = "s(unique.transect)")
  mu = pred_mat[,1:19]%*%coef(LATR_GAU_best)[1:19]
  pred = predict.gam(LATR_GAU_best,
                     newdata = data.frame(dens_scaled = d, log_volume_t = xb, fitted_vals = mu, unique.transect = "1.FPS"),
                     type = "response",
                     exclude = "s(unique.transect)")
  return(dnorm(y,mean=pred[,1],1/pred[,2]))
}

# Growth from size x to y at density d, using JSU
growth.JSU <- function(x, y, d){
  xb = pmin(pmax(x, LATR_size_bounds$min_size), LATR_size_bounds$max_size)
  pred_mat = predict.gam(LATR_GAU_best,
                         newdata = data.frame(dens_scaled = d, log_volume_t = xb, fitted_vals = 0, unique.transect = "1.FPS"),
                         type = "lpmatrix",
                         exclude = "s(unique.transect)")
  mu = pred_mat[,1:19]%*%coef(LATR_GAU_best)[1:19]
  pred = predict.gam(LATR_GAU_best,
                     newdata = data.frame(dens_scaled = d, log_volume_t = xb, fitted_vals = mu, unique.transect = "1.FPS"),
                     type = "response",
                     exclude = "s(unique.transect)")
  return(dJSU(y,
              mu=pred[,1],
              sigma=1/pred[,2],
              nu=JSUout$estimate[1]+JSUout$estimate[2]*pred[,1],
              tau=exp(JSUout$estimate[3]+JSUout$estimate[4]*pred[,1])))
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
# We need a large lower extension because growth variance is greater for smaller plants
ApproxMatrix <- function(dens,ext.lower,ext.upper,
                         min.size=LATR_size_bounds$min_size,max.size=LATR_size_bounds$max_size,
                         mat.size,dist){
  
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

