## Purpose: prep growth data for analysis
## hand-picking four data sets that appear easy to work with
## and represent various taxonomic groups
library(tidyverse)

## data set 1: gorgonian coral from Peterson et al. (data originator: Linares)
## three sites, each with 3-5 years (1999-2004)
## size is max height, units not given in Peterson et al. (check Linares papers)
coral <- read.csv("coral/Gorgonian_raw_data.csv")
colnames(coral)[6:8] <- c("Size"="t0", "Sizenext"="t1", "Survnext"= "survival") # rename column header
## drop rows without t0 and t1 size data
coral %>% filter(!is.na(t0) & !is.na(t1))-> coral_final
## plot all sites and years
plot(coral$t0,coral$t1) 
plot(log(coral$t0),log(coral$t1)) 

## data set 2:  my cactus data
## spatial (plots) and temporal (2004-2018) structure
## three size measurements, combined into a measure of volume (cm^3)
## this is the data file submitted to EDI
volume <- function(h, w, p){
  (1/3)*pi*h*(((w + p)/2)/2)^2
}
cactus <- read.csv("cactus/cholla_demography_20042018_EDI.csv")
cactus$t0 <- log(volume(cactus$Height_t,cactus$Width_t,cactus$Perp_t))
cactus$t1 <- log(volume(cactus$Height_t1,cactus$Width_t1,cactus$Perp_t1))
## drop rows without t0 and t1 size data
cactus %>% filter(!is.na(t0) & !is.na(t1))-> cactus_final
plot(cactus$t0,cactus$t1) 

## data set 3: pike (fish) data from Stubberud et al. 2019
## source paper has two sexes but this data set is female-only
## I am a little confused about this data set because the metadata
## describe destructive sampling but the actual data suggest recaptures
## original paper describes stricyly non-negative size changes
## they modeled with lognormal distn and size-dependent variance,
## truncated to prevent shrinkage
## included random year effects
pike <- read.csv("pike/data/PikeGrowthData1944_1995.csv") %>% 
  arrange(Year,Ind)
## this will need some data manipulation, starting with wide format
pike_wide <- pike %>% select(-RowID) %>% 
  pivot_wider(names_from = Year,values_from = Length)
## make sure that the column names are ordered years with no missing year
as.numeric(names(pike_wide)[-1]) == 1944:1995 #good
##now stack transition years
pike_trans_year <- pike_wide[,(1:3)]
pike_trans_year$year <- as.numeric(names(pike_wide)[2])
names(pike_trans_year)[2:3]<-c("t0","t1")
for(i in 3:(ncol(pike_wide)-1)){
  hold <- pike_wide[,c(1,i,(i+1))]
  hold$year <- as.numeric(names(pike_wide)[i])
  names(hold)[2:3]<-c("t0","t1")
  pike_trans_year <- bind_rows(pike_trans_year,hold)
}
## drop rows without t0 and t1 size data
pike_trans_year %>% filter(!is.na(t0) & !is.na(t1))-> pike_final
plot(pike_final$t0,pike_final$t1)
abline(0,1)
## add temperature data 
temps <- read.csv("pike/data/LakeTemp1944_2002.csv")
pike_final$temperature <- rep(NA,nrow(pike_final)); 
years = unique(pike_final$year); nyears=length(years); 
for(j in 1:nyears) {
	theYear = years[j];
	tj = which(temps$Year==theYear);
	tempj = temps$Temperature[tj]
	pj = which(pike_final$year==theYear);
	pike_final$temperature[pj]=tempj;
}	
	
## data set 4: voles from van Benthem et al.
## these data are from males only
## original study included fixed effects of population phase and reproductive stage
## and random effects of month and individual
## note that growth here is month to month
voles <- read.csv("voles/Males.csv")
voles$year <- as.numeric(substr(voles$date,start=1,stop=4))
voles$month <- as.numeric(substr(voles$date,start=6,stop=7))
#drop individuals with only one capture
voles_recap <- voles %>% 
  arrange(id,year,month) %>% 
  group_by(id) %>% 
  summarise(n_cap = n()) %>% 
  filter(n_cap > 1)
## convert to get size changes between consecutive months
voles_wide <- voles %>% 
  select(-X,-date,-session,-repro) %>% 
  filter(id %in% voles_recap$id) %>% 
  arrange(year,month) %>% 
  mutate(year_month = interaction(year,month)) %>% 
  select(-year,-month) %>% 
  pivot_wider(names_from = year_month,values_from = mass)
## reformat for growth over adjecent months
voles_trans_month <- voles_wide[,(1:5)]
voles_trans_month$year <- as.numeric(substr(names(voles_wide)[4],1,4))
voles_trans_month$month <- as.numeric(substr(names(voles_wide)[4],6,7))
names(voles_trans_month)[4:5] <- c("t0","t1")
## now loop over columns and stack up
for(i in 5:(ncol(voles_wide)-1)){
  hold <- voles_wide[,c(1:3,i,(i+1))]
  hold_year <- as.numeric(substr(names(hold)[4],1,4))
  hold_month1 <- as.numeric(substr(names(hold)[4],6,7))
  hold_month2 <- as.numeric(substr(names(hold)[5],6,7))
  ## check that columns i and i+1 are adjacent months
  if(hold_month2 == hold_month1+1 |
     {hold_month2==1 & hold_month1==12}){
    hold$year <- hold_year
    hold$month <- hold_month1
    names(hold)[4:5]<-c("t0","t1")
    voles_trans_month <- bind_rows(voles_trans_month,hold)
  }
}
## drop rows without t0 and t1 size data
voles_trans_month %>% filter(!is.na(t0) & !is.na(t1))-> voles_final
plot(voles_final$t0,voles_final$t1)
abline(0,1)

## save list of four data sets as rds file
growth_dat <- list(coral=coral_final,cactus=cactus_final,pike=pike_final,voles=voles_final)
write_rds(growth_dat,"growth_dat.rds")
