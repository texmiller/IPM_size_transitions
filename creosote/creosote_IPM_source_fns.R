## Source functions for creosote IPM, using gams fit in creosote_growth_modeling_gammit.R
## These functions are size (x) and density (d), with no other spatial or temporal random effects
## d is the weighted density of the window (sum of log volume)

## Growth
#growth_fn <- 

## Survival-- prediction from naturally occuring plants
survival_fn <- function(x,d){
  pred <- predict.gam(LATR_surv_best,
              newdata = data.frame(
                d.stand = d,
                log_volume_t = x,
                transplant=F,
                unique.transect="foo"),
              exclude = "s(unique.transect)")
  return(invlogit(pred))
}
survival_fn(x=5,d=0)

## Flowering
flower_fn <- function(x,d){
  pred <- predict.gam(LATR_flower_best,
              newdata = data.frame(
                d.stand = d,
                log_volume_t = x,
                unique.transect="foo"),
              exclude = "s(unique.transect)")
  return(invlogit(pred))
}

flower_fn(x=5,d=0)
