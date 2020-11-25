## Source functions for creosote IPM, using gams fit in creosote_growth_modeling_gammit.R
## These functions are size (x) and density (d), with no other spatial or temporal random effects
## d is the weighted density of the window (sum of log volume)

## Growth -- Gaussian using best gam
growth_fn_norm <- function(x,y,d){
  lpmat <- predict.gam(LATR_gam_model,
                      newdata = data.frame(
                        weighted.dens = d,
                        log_volume_t = x,
                        unique.transect="foo"),
                      type="lpmatrix",
                      exclude = "s(unique.transect)")
  ## linear predictor for mean and log sigma -- need to update so these indices are not hard-coded but for now they work
  grow_mu <- lpmat[,1:19]%*%coef(LATR_gam_model)[1:19]
  grow_sigma <- exp(lpmat[,32:50]%*%coef(LATR_gam_model)[32:50])
  return(dnorm(y,mean=grow_mu,sd=grow_sigma))
}

## Survival-- prediction from naturally occuring plants
survival_fn <- function(x,d){
  pred <- predict.gam(LATR_surv_best,
              newdata = data.frame(
                weighted.dens = d,
                log_volume_t = x,
                transplant=F,
                unique.transect="foo"),
              exclude = "s(unique.transect)")
  return(invlogit(pred))
}

## Flowering
flower_fn <- function(x,d){
  pred <- predict.gam(LATR_flower_best,
              newdata = data.frame(
                weighted.dens = d,
                log_volume_t = x,
                unique.transect="foo"),
              exclude = "s(unique.transect)")
  return(invlogit(pred))
}

## Fruits
fruits_fn <- function(x,d){
  pred <- predict.gam(LATR_fruits_best,
                      newdata = data.frame(
                        weighted.dens = d,
                        log_volume_t = x,
                        unique.transect="foo"),
                      exclude = "s(unique.transect)")
  return(exp(pred))
}
