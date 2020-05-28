##############################################################################
# Growth kernel modeling for Cylindropuntia imbricata - modeled off of Steve's PSSPgrowthModeling-ML.R script. 
#
# Data collection described in Miller et al. 2009, Ohm and Miller 2014, Czachura and Miller 2020, and elsewhere  
# This data set includes 2018 data, which has not been previously published. 
#
# Last update: 38 May, 2020
#
##############################################################################

rm(list=ls(all=TRUE));
setwd("./cactus"); 

require(car); require(lme4); require(zoo); require(moments); require(mgcv); 
require(gamlss); require(gamlss.tr); require(AICcmodavg); 
require(lmerTest); require(tidyverse); require(maxLik); 

# Steve's diagnostics functions
source("../Diagnostics.R") 
# function for converting cactus size measurements to volume
volume <- function(h, w, p){
  (1/3)*pi*h*(((w + p)/2)/2)^2
}


##############################################################
# 1. Read in the data and pare down data set to just size transitions
##############################################################

cactus<-read_csv("cholla_demography_20042018_EDI.csv") %>% 
  ## filter out transplants and drop seed addition plots (which start with letter H)
  ## also drop Year_t==2018 because I don't have 2019 size data (not entered yet). 2018 data still included in 2017-2018 transition
    filter(Transplant == 0,
         str_sub(Plot,1,1)!="H",
         Year_t!=2018) %>% 
  ## convert height, max width, perp width to volume of cone, take natural log
  mutate(vol_t = volume(Height_t,Width_t,Perp_t),
         vol_t1 = volume(Height_t1,Width_t1,Perp_t1),
         plot = as.factor(Plot),
         year_t = as.factor(Year_t)) %>%
  select(year_t,plot,vol_t,vol_t1) %>% 
  filter() %>% 
  ## sort by initial size
  arrange(vol_t) %>% 
  ## drop rows with NAs
  drop_na()

## how much coverage do I have across years and plots?
table(cactus$plot,cactus$year_t) ## first four years has fewer plants from fewer plots


