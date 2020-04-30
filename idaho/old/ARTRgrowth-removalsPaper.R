# PBA March 2016

# call from removal_analysis_wrapper.r

#########################################
#  1. Import data and calculate W's
#########################################

doSpp <- "ARTR"
sppList <- c("ARTR","HECO","POSE","PSSP","allcov","allpts")
dataDir1 <- paste(root,"/ExperimentTests/data/idaho",sep="")
dataDir2 <- paste(root,"/ExperimentTests/data/idaho_modern",sep="")
nonCompLength.s=5 #Number of columns in SppData that are not measures of competitors 

# set up distance weights------------------------------------------------
#dists <- read.csv(paste(dataDir2,"/speciesdata/IdahoModDistanceWeights_noExptl.csv",sep=""));
dists$allcov <- rowMeans(dists[,1:4])  # for "other" polygons use average of big 4
dists$allpts <- dists$POSE  # set forb dist wts = smallest grass (POSE)

# import old data--------------------------------------------------------

source("fetchGrowthData.r")

D1 <- fetchGdat(doSpp=doSpp,speciesList=sppList,datadir=dataDir1,distWts=dists)
D1$Treatment <- "Control"

# import modern data--------------------------------------------------------

D2 <- fetchGdat(doSpp=doSpp,speciesList=sppList,datadir=dataDir2,distWts=dists)

# merge in treatment data
tmp <- read.csv(paste(dataDir2,"/quad_info.csv",sep=""))
tmp <- tmp[,c("quad","Treatment")]
D2 <- merge(D2,tmp, all.x=T)

# account for removal in baseline years
if(doSpp!="ARTR"){
  ii <- which(D2$year>=2011 & D2$Treatment=="No_shrub")
  D2$W.ARTR[ii] <- 0
}else{
  ii <- which(D2$year>=2011 & D2$Treatment=="No_grass")
   D2$W.HECO[ii] <- 0 ; D2$W.POSE[ii] <- 0 ; D2$W.PSSP[ii] <- 0
}

# combine old and modern
allD <- rbind(D1,D2)
rm(D1,D2,tmp)

# clean up dataset ----------------------------------------------
allD$year[allD$year<2000] <- allD$year[allD$year<2000] + 1900

if(doSpp=="ARTR"){
  keep <- which(is.element(allD$Treatment,c("Control","No_grass")))
}else{
  keep <- which(is.element(allD$Treatment,c("Control","No_shrub")))
}
allD <- allD[keep,]

# remove outliers (large plants that obviously do not turn into tiny plants)
tmp=which(allD$quad=="Q23" & allD$year==1945 & allD$trackID==67)
tmp=c(tmp,which(allD$quad=="Q12" & allD$year==1955 & allD$trackID==25))
tmp=c(tmp,which(allD$quad=="Q26" & allD$year==1945 & allD$trackID==73))
allD=allD[-tmp,]

#########################################
#  2. Fit models
#########################################

# set up indicator variables
allD$Treatment2 <- allD$Treatment
allD$Treatment2[allD$year>2000] <- "Modern"
allD$Treatment3 <- allD$Treatment
allD$Treatment3[allD$Treatment=="Control" & allD$year>2000] <- "ControlModern"

allD$year <- as.factor(allD$year)

# use lmer
# #library(lme4)
# # no treatment effect
# m0 <- lmer(logarea.t1~logarea.t0+W.ARTR + W.HECO + W.POSE + W.PSSP+  W.allcov + W.allpts +
#              (1|Group)+(logarea.t0|year),data=allD) 
# # treatment effect
# m1 <- lmer(logarea.t1~logarea.t0+Treatment+W.ARTR + W.HECO + W.POSE + W.PSSP+  W.allcov + W.allpts +
#              (1|Group)+(logarea.t0|year),data=allD) 
# 

# use INLA
# Set up ID variables for INLA random effects
allD$GroupID <- as.numeric(allD$Group)
allD$yearID <- 100+as.numeric(allD$year) # for random year offset on intercept

# Basic treatment effect (intercept) model
m1 <- inla(logarea.t1 ~ logarea.t0 + Treatment + W.ARTR + W.HECO + W.POSE + W.PSSP + W.allcov + W.allpts +
  f(yearID, model="iid", prior="normal",param=c(0,0.001))+
  f(GroupID, model="iid", prior="normal",param=c(0,0.001))+
  f(year, logarea.t0, model="iid", prior="normal",param=c(0,0.001)), data=allD,
  family=c("gaussian"), verbose=FALSE,
  control.predictor = list(link = 1),control.compute=list(dic=T,mlik=T),
  control.inla = list(h = 1e-10),Ntrials=rep(1,nrow(allD)))

# # compare to other models using lmer
# m0.lmer <- lmer(logarea.t1~logarea.t0 + W.ARTR + W.HECO + W.POSE + W.PSSP+ W.allcov + W.allpts+
#              (1|Group)+(logarea.t0|year),data=allD) 
# m1.lmer <- lmer(logarea.t1~logarea.t0+ Treatment + W.HECO + W.POSE + W.PSSP+ W.ARTR + W.allcov + W.allpts+
#              (1|Group)+(logarea.t0|year),data=allD) 
# m2.lmer <- lmer(logarea.t1~logarea.t0+ Treatment + W.HECO + W.POSE + W.PSSP+ W.ARTR + W.allcov + W.allpts+
#            W.ARTR:Treatment+  
#              (1|Group)+(logarea.t0|year),data=allD) 
# print(c(AIC(m0.lmer),AIC(m1.lmer),AIC(m2.lmer))) # simplest model (no treatment effects) is best
# 
# # fit better model
# m.best <- inla(logarea.t1 ~ logarea.t0 + W.ARTR + W.HECO + W.POSE + W.PSSP + W.allcov + W.allpts +
#   f(yearID, model="iid", prior="normal",param=c(0,0.001))+
#   f(GroupID, model="iid", prior="normal",param=c(0,0.001))+
#   f(year, logarea.t0, model="iid", prior="normal",param=c(0,0.001)), data=allD,
#   family=c("gaussian"), verbose=FALSE,
#   control.predictor = list(link = 1),control.compute=list(dic=T,mlik=T),
#   control.inla = list(h = 1e-10),Ntrials=rep(1,nrow(allD)))
