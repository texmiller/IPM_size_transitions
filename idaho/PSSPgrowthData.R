##############################################################################
# Create .csv file with data on Pseudoroegneria spicata growth, USSES Idaho.
# Historical and modern data from by Adler et al. removals experiments paper
# Modern data are just the Control and shrub-removal treaments.  
# Seedlings need to be modeled separately, and are excluded here. 
#
# Original: PBA March 2016
# Modified for IPM_size_transitions by SPE, April 2020 
#
# NOTE: this is PSSP-specific; some shields in the original code against 
# bringing in removals plots data for other species have been removed. 
##############################################################################

rm(list=ls(all=TRUE));
setwd("c:/repos/IPM_size_transitions/idaho"); 

source("Utilities.R");

##############################################################
#  1. Import data and calculate W's  - ARTR
##############################################################
sppList <- c("ARTR","HECO","POSE","PSSP","allcov","allpts")
doSpp <- "PSSP"   # !!!!!!!!!!!!!  Don't change this

dataDir1 <- "c:/repos/IPM_size_transitions/idaho/PSSP/legacy_data/"
dataDir2 <- "c:/repos/IPM_size_transitions/idaho/PSSP/modern_data/"
nonCompLength.s=5 #Number of columns in SppData that are not measures of competitors 

# set up distance weights------------------------------------------------
dists <- read.csv("PSSP/modern_data/IdahoModDistanceWeights_noExptl.csv");
dists$allcov <- rowMeans(dists[,1:4])  # for "other" polygons use average of big 4
dists$allpts <- dists$POSE  # set forb dist wts = smallest grass (POSE)

#########################################
# import old data------------------------
#########################################
D1 <- fetchGdat(doSpp=doSpp,speciesList=sppList,datadir=dataDir1,distWts=dists)
D1$Treatment <- "Control"

#########################################
# import modern data---------------------
#########################################
D2 <- fetchGdat(doSpp=doSpp,speciesList=sppList,datadir=dataDir2,distWts=dists)

# merge in treatment identity for each quadrat 
tmp <- read.csv(paste(dataDir2,"quad_info.csv",sep=""))
tmp <- tmp[,c("quad","Treatment")]
D2 <- merge(D2,tmp, all.x=T)

# limit to control and removals plots 
keep <- which(is.element(D2$Treatment,c("Control","No_shrub")))
D2 <- D2[keep,]; 

#######################################################################
# 2. Merge, clean up, and save 
#######################################################################
allD = rbind(D1,D2); 
allD$year[allD$year<2000] <- allD$year[allD$year<2000] + 1900

## eliminate likely recording errors: large plants that became tiny in one year. 
#tmp <- which((allD$area.t0>100)&(allD$area.t1 < 0.26)); 
#allD <- allD[-tmp,]; 

## get rid of seedlings 
allD <- trimQuadrats(allD)$data;
allD <- subset(allD,age>1); 

plot(logarea.t1~logarea.t0,data=allD); 

allD$Treatment[allD$Treatment=="Control" & allD$year>2000] <- "ControlModern"
allD$year <- as.factor(allD$year)

cols <- c("area.t0","area.t1","logarea.t0","logarea.t1","year","Group","age","W.ARTR","W.HECO","W.POSE","W.PSSP","W.allcov","W.allpts","Treatment")     
allD <- allD[,cols]; 
e = order(allD$area.t0); allD <- allD[e,]; 

write.csv(allD,file="PSSP_growth_data.csv"); 
