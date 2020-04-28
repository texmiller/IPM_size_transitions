# From repos/ExperimentTests/removals/growth 
# Original by PBA, based on earlier scripts by SPE and BT 
# Modified by SPE for IPM_size_transitions, April 2020 

fetchGdat <- function(doSpp,speciesList,datadir,distWts){

  growDfile=paste(datadir,"growDnoNA.csv",sep="")
  growD=read.csv(file=growDfile)
  D1=growD[growD$allEdge==0,];
  D1$year <- D1$year
  D1$logarea.t0=log(D1$area.t0)
  D1$logarea.t1=log(D1$area.t1)
  D1$quad=as.character(D1$quad)
  
  # import neighbor data
  ringD <- read.csv(paste(datadir,doSpp,"_nbhood_rings.csv",sep=""))
  tmpD <- read.csv(paste(datadir,doSpp,"_nbhood_rings_allothers.csv",sep=""))
  ringD<-merge(ringD,tmpD)
  ringD$year<-ringD$year
  
  # merge D with ringD (D contains fewer rows)
  D1<-merge(D1,ringD,by.x=c("quad","year","trackID"),by.y=c("quad","year","genetID"))
  D1=D1[order(D1$X),]
  rm(ringD,growD)
  row.names(D1) <- NULL  
  
  # calculate W's (MAKE SURE NOT TO REORDER D!)
  W <- matrix(NA,NROW(D1),length(speciesList))
  colnames(W) <- paste("W.",speciesList,sep="")

  # do big 4
  for(iSpp in 1:4){
    neighborCols=which(substr(names(D1),1,4)==speciesList[iSpp]) # pull out annulus data for the focal species 
    dist_wts <- distWts[,paste0(speciesList[iSpp])]
    C <- data.matrix(D1[,neighborCols]) #matrix of conspecific areas in the annuli 
    W[,iSpp] <- C%*%dist_wts 
  }
  
  # do allcov and allpts
  for(iSpp in 5:6){
    neighborCols=which(substr(names(D1),1,6)==speciesList[iSpp]) # pull out annulus data for the focal species 
    dist_wts <- distWts[,paste0(speciesList[iSpp])]
    C <- data.matrix(D1[,neighborCols]) #matrix of conspecific areas in the annuli 
    W[,iSpp] <- C%*%dist_wts 
  }
  
  #format
  D1 <- D1[,c("quad","year","trackID","age","distEdgeMin","allEdge","QuadName","Grazing","Group","area.t0","area.t1","logarea.t0","logarea.t1","species")]
  D1 <- cbind(D1,W)
  
  return(D1)

}

##########################################################################
# This function removes individuals who would be flagged as seedlings
# based on the rule: (age==1)&(size <= log(.25)) ==> seedling 
# but who might not be seedlings because their age is > 1
#
# Individuals larger than 0.25 are assumed to be older (true age > 1) 
# regardless of what age is recorded for them. 
##########################################################################

trimQuadrats <- function(data,skip.years=TRUE) {  
  quadrats=unique(data$quad); doubtful=0; 
  for(j in 1:length(quadrats)) {
    qj = quadrats[j]; years = unique(data$year[data$quad == qj]); 
    year1 = min(years)
    novice = which((data$year==year1)&(data$age==1)&(data$quad==qj)&(data$logarea.t0<log(0.25001)));  
    if(length(novice)>0) {
        data=data[-novice,];  #cat(qj,year1, length(novice), "\n"); 
        doubtful=doubtful+length(novice)
    } 

   
    if(skip.years){ # trim dubious ages after a gap in the data for a quadrat
    skipyr = which(diff(years)>1); 
    if(length(skipyr)>0) {
      skipyr=years[1+skipyr]; 
      for(k in 1:length(skipyr)){
          novice = which((data$year==skipyr[k])&(data$age==1)&(data$seedling==0)&(data$quad==qj)&(data$logarea.t0<log(0.25001)))  
          if(length(novice)>0) {
              data=data[-novice,];  #cat(qj,skipyr[k], length(novice), "\n"); 
              doubtful=doubtful+length(novice);
          }
      } 
    }
  }
  }
  likely.seedlings=which((data$age==1)&(data$logarea.t0<=log(0.25001)))    
  return(list(data=data,doubtful=doubtful,seedlings=length(likely.seedlings))); 
}



