###### This script fits discretely estimated matrix models and IPMs based on continuous
###### vital rate functions for a range of bin numbers and discretization methods
###### using the data for alpine bistort (Polygonum viviparum). 
###### This version of the script fits with even class widths, the recommended approach 
###### fitting models. In the main text, results of the script shown here are included
###### in Figures 10, 11, and 12. 

setwd("c:/repos/ipm_size_transitions/Doak/SM"); 
graphics.off()
rm(list=ls())
library(dplyr)
library(popbio)
library(MuMIn)
library(binr)
library(matrixStats)

#Remember to set directory to access input files

### Specify what range of classes/bins to evaluate -----------------------
bin.num <- c(3,4,5,6, 8, 10,seq(from = 15, to = 100, by = 10))


### Load, manipulate, and sort input data --------------

# Define mean and SD of bulbling sizes: 
allsdszs= as.numeric(read.csv('bulblingsizes.csv', header=FALSE)[,1])
mnsldgsz=mean(allsdszs) # mean bulbling size
sdsldgsz = sqrt(var(allsdszs)) # sd bulbling size

# bulblings seen in next year per bulbil produced
bulblings_per_bulbil = 0.00676 # this is mean for all years and populations

# read in full data set and clean: 
# alldata has columns of: szs0, szs1, bulbils produced in year 0
alldata=as.matrix(read.csv('BNdataforIPMsimple.csv', header=TRUE))
alldata=as.data.frame(alldata)
alldata=alldata[which(is.na(alldata[,1])==FALSE),] # remove plants missing or dead in year 0
alldata$survival = 1
alldata$survival[which(alldata$szs1==0)] = 0 # plants with 0 size in year 1 are dead
alldata$szs1[which(alldata$szs1==0)] = NA # set the size of dead plants to NA
alldata=alldata[which(alldata$szs0 !=0),]

# define minimum and maximum sizes for modeling and enforce limits on size:
minsize <- min(c(alldata[,1],alldata[,2]), na.rm=TRUE) -0.1
maxsize=60 # this is set manually to not have a large size classes with few to no plants
alldata[which(alldata[,1]>maxsize),1]=maxsize
alldata[which(alldata[,2]>maxsize),2]=maxsize 

# set the survival of bulblings as survial for all plants <= average bulbling size:
alltinys=alldata[alldata$szs0<=mnsldgsz,]
surv_bulbings=sum(alltinys$survival)/length(alltinys$survival)

# redefine column names for compatability: 
colnames(alldata)=c('t0','t1','bulbs0','survival')

# define a variable for reproducing or not, and make separate data frame for reproductive analysis:
alldata$reproyesno = alldata$bulbs0 # number of bulbils produced
alldata$reproyesno[alldata$reproyesno >0]=1 # so, turn into yes no variable
sizeforrepro=alldata
sizeforrepro=sizeforrepro[sizeforrepro$reproyesno==1,]

#generate empirical density function for median size estimation
#this is a set of smoothed values that can then be used with weightedMedian
#in the matrixStats package to get an estimated median size for each class.
pdfsz=density(alldata$t0, n=1024, cut=0)
pdfsz2=cbind(pdfsz$x,pdfsz$y)


#### lifespan calculation function################
# the age at which a new recruit has <1% probability of still being alive
lifespan <- function(nx){
  nclasses=dim(nx)[1]
  vec=c(100,rep(0,(nclasses-1)))  
  nx[1,]=0
  jj=1
  while (sum(vec)>1){
    vec=nx%*%vec
    jj=jj+1
  }
  return(jj)
}
#############################################


### BEGIN THE ESTIMATION FOR EACH MODEL TYPE: MATRIX = DVR IN TEXT, WHILE IPM MEAN AND IPM MEDIAN = CVR-MIDPOINT AND CVR-MEDIAN IN TEXT. 

### Matrix Models ------------------------


## make storage structures
lambdas_matrix <- rep(NA, length(bin.num))
dampratio_matrix = rep(NA, length(bin.num))
lifespan_matrix=rep(NA, length(bin.num))

mincounts=NULL #tracker for minimum counts of individuals within class

for(i in 1:length(bin.num)){ # loop over different numbers of bins
  
 # making classes: 
ss=as.numeric(alldata$t0)
vec.bin = c(minsize, minsize+1:bin.num[i]*(maxsize-minsize)*(1/bin.num[i])) # evenly distributed size classes  
nums=hist(ss,breaks=vec.bin, plot=FALSE)$counts 
mincounts=c(mincounts,min(nums)) # minimum number of individuals in a class
  
if (min(nums)>2) {  
  
  ## Initialize storage
  n.bin <- length(vec.bin)-1                  
  n <- rep(NA, n.bin)                         # count of individuals per bin
  medians <- rep(NA, n.bin)                   # median size per bin for F
  surv <- rep(NA, n.bin)                      # survivorship for each class
  grow <- matrix(NA, n.bin, n.bin)            # store growth probabilites for each class
  reproduction <- rep(NA, n.bin)              # store reproduction for each class
  
  totnums = 0 # this is monitoring for errors
  # bin survival, growth
  for(j in 1:(length(vec.bin)-1)){
    # set limits for subset according to bin breaks
    bounds <- c(vec.bin[j], vec.bin[j+1])
    # subset data according to bounds
    subset <- alldata[alldata$t0 > bounds[1] & alldata$t0 <= bounds[2],]
    # store number and median size of individuals in this bin for future reference
    n[j] <- length(subset$t0)
    medians[j] <- median(subset$t0)
    # calculate survivorship for this class
    surv[j] <- sum(subset$survival) / length(subset$t0)
    # store histo as object, to access counts per bin
    histo <- hist(subset$t1, breaks = vec.bin, plot = FALSE)
    # returns the number of individuals of a certain size class
    # to give the transition frequencies to each size class in year 1
    grow[,j] <- histo$counts/length(subset$t0[subset$survival==1]) 
    reproduction[j] <- mean(subset$bulbs0) # bulbils produced per plant in start yr
    
    totnums = totnums + sum(histo$counts)
  }
 
   # make a vector of the probability of seedling sizes: 
  sdlggrow=hist(allsdszs , breaks = vec.bin, plot = FALSE)$counts/length(allsdszs)
  
  M1 <- matrix(NA, n.bin, n.bin)   # initialize survival/growth matrix
  M <- matrix(0, (n.bin+1), (n.bin+1))# initialize projection matrix

    # populate projection matrix
  for(j in 1:length(surv)){
    M1[,j] <- surv[j] * grow[,j]
    M1[,j] <- (surv[j] * grow[,j])
  }
    
    # add lines for the creation of bulblings and their transition to first size class
    M[2:(n.bin+1), 2:(n.bin+1)] = M1
    M[2:(n.bin+1),1] = surv_bulbings*sdlggrow # survival and growth of bulblings
    M[1,2:(n.bin+1)] = reproduction*bulblings_per_bulbil  # creation of bulblings
  
  
  lambdas_matrix[i] <- lambda(M) # calls popbio functions 'lambda' and damping.ratio
  dampratio_matrix[i]=damping.ratio(M)
  lifespan_matrix[i] =lifespan(M)
  
  print(c(i,totnums)) 
} else {
  lambdas_matrix[i]=NA
  dampratio_matrix[i] = NA
  lifespan_matrix[i] = NA
}
}

######################################
### IPM using midpoint size-----------------------------------------------------
# NOTE THAT ALL THE STEPS OF ANALYSIS ARE REPEATED FOR THE MIDPOINT AND MEDIAN APPROACHES, SO THAT THEY ARE INPEPENDENT OF ONE ANOTHER

# adding in columns for squared size effects for fitting vital rate models
alldata$t0sq <- alldata$t0^2           
sizeforrepro$t0sq = sizeforrepro$t0^2

# fitting functions for different vital rates (survival, growth mean, growth variance, reproduction (yes, no) and amount of reproduction if reproducing)

# prob (survival): compare linear function of size, quadratic function of size
sur_models <- list(glm(survival~ t0, family= "binomial",data = alldata),
                   glm(survival~ t0 + t0sq , family= "binomial",data =alldata))
# gives the value of the lowest AICc
min_AIC <- min(AICc(sur_models[[1]], sur_models[[2]])$AICc)
# gives you the info for the best-fit model
bestsur <- sur_models[[which(min_AIC == AICc(sur_models[[1]], sur_models[[2]])$AICc)]] 

# growth mean: compare linear function of size, quadratic function of size, power function (A+B*(size^C))
growth_models <- list(nls(t1~ a+ b*t0, data = alldata, start=list(a= 1, b=1)),
                      nls(t1~ a+ b*t0 + c*t0sq, data = alldata, start=list(a= 1, b=1,c=1)),
                      nls(t1 ~ a + b*(t0^c), data= alldata, start=list(a= 1, b=1,c=1)))
# gives the value of the lowest AICc
min_AIC <- min(AICc(growth_models[[1]], growth_models[[2]],growth_models[[3]])$AICc)
# gives you the info for the best-fit model
bestgrowth <- growth_models[[which(min_AIC == AICc(growth_models[[1]], growth_models[[2]], growth_models[[3]])$AICc)]]

# getting residuals of growth to estimate growth variance
alldata$growth_residuals <- NA
alldata$growth_residuals[which(!is.na(alldata$t1) & !is.na(alldata$t0))] <- summary(bestgrowth)$residuals^2

# yes/no reproduction: compare linear function of size, quadratic function of size
reproyesno_models <- list(glm(reproyesno~ t0, family= "binomial",data = alldata),
                   glm(reproyesno~ t0 + t0sq , family= "binomial",data =alldata))
min_AIC <- min(AICc(reproyesno_models[[1]], reproyesno_models[[2]])$AICc)
bestreproyesno <- reproyesno_models[[which(min_AIC == AICc(reproyesno_models[[1]], reproyesno_models[[2]])$AICc)]] 

# number of bulbils if reproducing: compare linear function of size, quadratic function of size 
rep_models <- list(lm(bulbs0~ t0, data = sizeforrepro),
                    lm(bulbs0~ t0+ t0sq, data = sizeforrepro)) 
min_AIC <- min(AICc(rep_models[[1]], rep_models[[2]])$AICc)
bestrep <- rep_models[[which(min_AIC == AICc(rep_models[[1]], rep_models[[2]])$AICc)]] 

# variance in growth: uses best-fit model for growth mean: 
# compare linear function of size, quadratic function of size
# constrain to 0 intercept to prevent predictions of negative variance
var_models <- list(glm(growth_residuals~ t0-1, data= alldata), glm(growth_residuals~ t0 + t0sq-1, data= alldata))
min_AIC <- min(AICc(var_models[[1]], var_models[[2]])$AICc) 
bestvar <- var_models[[which(min_AIC == AICc(var_models[[1]], var_models[[2]])$AICc)]] 

# make models with varying bin numbers

lambdas_ipm_mean <- vector("numeric", length= length(bin.num))
dampratio_mean = vector("numeric", length= length(bin.num))
lifespan_mean = vector("numeric", length= length(bin.num))

truebinsizes= matrix(0,length(bin.num),1)

for (i in 1:length(bin.num)){

  vec.bin = c(minsize, minsize+1:bin.num[i]*(maxsize-minsize)*(1/bin.num[i])) # evenly distributed size classes
  binmids <- 0.5*(vec.bin[2:length(vec.bin)] + vec.bin[1:(length(vec.bin)-1)]) # mesh points are the midpoints of each class
  n.bin = length(binmids)
  truebinsizes[i] = n.bin  
  
  # constructing the matrices
  indata <- as.data.frame(cbind(binmids, binmids^2))
  names(indata) <- c("t0", "t0sq")
  
  
  # vital rate predictions from best models:
  sur_vals <- predict(bestsur,indata, type='response')
  reproyesnovals=predict(bestreproyesno, indata, type='response')
  rep_vals <- predict(bestrep, indata, type='response')
  rep_vals[rep_vals<0] <- 0
  growth_vals <- predict(bestgrowth,indata, type='response') # growth mean
  var_vals <- predict(bestvar,indata, type='response') # growth variance
  gmx <- matrix(data=NA, nrow=n.bin, ncol=n.bin)
   
  # growth probabilities using cdf function
  for (ss in 1:(n.bin)) {
    growcdf <- pnorm(vec.bin,growth_vals[ss],sqrt(var_vals[ss])) # get the cdf at bin edges
    grows <- growcdf[2:length(vec.bin)]-growcdf[1:(length(vec.bin)-1)] # take the difference between bin edges
    if(sum(grows)>0){grows <- grows/sum(grows) # avoid eviction for re-normalizing to sum to 1
                     gmx[,ss] <- grows} else {gmx[,ss] <- NA} 
    # this if statement breaks the code (puts NA's into the matrix) if the sum of the PDF is zero (which happens if all the probability is outside of the size bounds)
  } # end ss loop
  
  # make the surv*growth matrix
  survgmx <- gmx*t(matrix( rep(sur_vals,(n.bin)),(n.bin))) # survs* growth
  mx1 <- survgmx # growth and survival matrix, without the reproduction
  
  #vector of reproductive rates (probability of reproducing * number of bulbils if reproducing)
  reprow <- rep_vals*reproyesnovals  
  
  # making estimated seedling sizes using the cdf function
  sdlgszcdf=pnorm(vec.bin,mnsldgsz,sdsldgsz) # get the cdf at bin edges
  sdlgszprobs=sdlgszcdf[2:length(vec.bin)]-sdlgszcdf[1:(length(vec.bin)-1)] # take the difference in cdfs between bin edges
  sdlgszprobs=sdlgszprobs/sum(sdlgszprobs) # avoid eviction by re-normalizing to sum to 1
  
  mx <- matrix(0, (n.bin+1), (n.bin+1))
  mx[2:(n.bin+1), 2:(n.bin+1)] = mx1
  mx[2:(n.bin+1),1] = surv_bulbings*sdlgszprobs # survival and growth of bulblings
  mx[1,2:(n.bin+1)] = reprow*bulblings_per_bulbil  # creation of new bulblings
 
  lambdas_ipm_mean[i] <- Re(eigen(mx)$values[1])
  dampratio_mean[i]= damping.ratio(mx)
  lifespan_mean[i] =lifespan(mx)
 print(i)
}


#######################################
######################################
### IPM with median sizes-----------------------------------------------------

# adding in columns for squared size effects
alldata$t0sq <- alldata$t0^2           
sizeforrepro$t0sq = sizeforrepro$t0^2

# fitting functions for different vital rates (survival, growth, growth variance, reproduction (yes, no) and amount of reproduction if reproducing)

# prob (survival): compare linear function of size, quadratic function of size
sur_models <- list(glm(survival~ t0, family= "binomial",data = alldata),
                   glm(survival~ t0 + t0sq , family= "binomial",data =alldata))
# gives the value of the lowest AICc
min_AIC <- min(AICc(sur_models[[1]], sur_models[[2]])$AICc)
# gives you the info for the best-fit model
bestsur <- sur_models[[which(min_AIC == AICc(sur_models[[1]], sur_models[[2]])$AICc)]] 

# growth mean: compare linear function of size, quadratic function of size, power function (A+B*(size^C))
growth_models <- list(nls(t1~ a+ b*t0, data = alldata, start=list(a= 1, b=1)),
                      nls(t1~ a+ b*t0 + c*t0sq, data = alldata, start=list(a= 1, b=1,c=1)),
                      nls(t1 ~ a + b*(t0^c), data= alldata, start=list(a= 1, b=1,c=1)))
# gives the value of the lowest AICc
min_AIC <- min(AICc(growth_models[[1]], growth_models[[2]],growth_models[[3]])$AICc)
# gives you the info for the best-fit model
bestgrowth <- growth_models[[which(min_AIC == AICc(growth_models[[1]], growth_models[[2]], growth_models[[3]])$AICc)]]

# getting residuals of growth to estimate growth variance
alldata$growth_residuals <- NA
alldata$growth_residuals[which(!is.na(alldata$t1) & !is.na(alldata$t0))] <- summary(bestgrowth)$residuals^2

# yes/no reproduction: compare linear function of size, quadratic function of size
reproyesno_models <- list(glm(reproyesno~ t0, family= "binomial",data = alldata),
                          glm(reproyesno~ t0 + t0sq , family= "binomial",data =alldata))
min_AIC <- min(AICc(reproyesno_models[[1]], reproyesno_models[[2]])$AICc)
bestreproyesno <- reproyesno_models[[which(min_AIC == AICc(reproyesno_models[[1]], reproyesno_models[[2]])$AICc)]] 

# number of bulbils if reproducing: compare linear function of size, quadratic function of size, 
rep_models <- list(lm(bulbs0~ t0, data = sizeforrepro),
                   lm(bulbs0~ t0+ t0sq, data = sizeforrepro)) 
min_AIC <- min(AICc(rep_models[[1]], rep_models[[2]])$AICc)
bestrep <- rep_models[[which(min_AIC == AICc(rep_models[[1]], rep_models[[2]])$AICc)]] 

# variance in growth: uses best-fit model for growth: compare linear function of size, quadratic function of size
# constrain to 0 intercept to prevent prediction of negative variance
var_models <- list(glm(growth_residuals~ t0-1, data= alldata), glm(growth_residuals~ t0 + t0sq-1, data= alldata))
min_AIC <- min(AICc(var_models[[1]], var_models[[2]])$AICc) 
bestvar <- var_models[[which(min_AIC == AICc(var_models[[1]], var_models[[2]])$AICc)]] 

# make models with varying bin numbers

lambdas_ipm_median <- vector("numeric", length= length(bin.num))
dampratio_median <- vector("numeric", length= length(bin.num))
lifespan_median = vector("numeric", length= length(bin.num))

truebinsizes= matrix(0,length(bin.num),1)

for (i in 1:length(bin.num)){

 vec.bin = c(minsize, minsize+1:bin.num[i]*(maxsize-minsize)*(1/bin.num[i])) # evenly distributed size classes
  
  nums=hist(alldata$t0,breaks=vec.bin, plot=FALSE)$counts 
  
# getting estimated median size in each bin:
# these will be used as the mesh points to evaluate vital rate functions   
  binmids =   rep(NA, length(vec.bin)-1)
  for(jj in 1:(length(vec.bin)-1)){
    # set limits for subset according to bin breaks
    bounds <- c(vec.bin[jj], vec.bin[jj+1])
    # subset data according to bounds
     subsetszs <- pdfsz2[which(pdfsz2[,1] >= bounds[1] & pdfsz2[,1] < bounds[2]),]
     # get weighted median estimates:
    binmids[jj] <- weightedMedian(subsetszs[,1],subsetszs[,2]) 
  }
  
  # if there is failure in the median estimation, set the size for these bins to the midpoint: 
  meanbinmids <- 0.5*(vec.bin[2:length(vec.bin)] + vec.bin[1:(length(vec.bin)-1)])
  binmids[is.na(binmids)]=meanbinmids[is.na(binmids)]
 
  n.bin = length(binmids)
  truebinsizes[i] = n.bin  
  
  # constructing the matrices
  indata <- as.data.frame(cbind(binmids, binmids^2))
  names(indata) <- c("t0", "t0sq")
  
  # vital rate predictions from best models:
  sur_vals <- predict(bestsur,indata, type='response')
  reproyesnovals=predict(bestreproyesno, indata, type='response')
  rep_vals <- predict(bestrep, indata, type='response')
  rep_vals[rep_vals<0] <- 0
  growth_vals <- predict(bestgrowth,indata, type='response') # growth mean
  var_vals <- predict(bestvar,indata, type='response') # growth variance
  gmx <- matrix(data=NA, nrow=n.bin, ncol=n.bin)
  
  # growth probabilities using cdf function
  for (ss in 1:(n.bin)) {
    growcdf <- pnorm(vec.bin,growth_vals[ss],sqrt(var_vals[ss])) # get the cdf at bin edges
    grows <- growcdf[2:length(vec.bin)]-growcdf[1:(length(vec.bin)-1)] # take the difference in cdf values between bin edges
    if(sum(grows)>0){grows <- grows/sum(grows) # avoid eviction by re-normalizing to sum to 1
    gmx[,ss] <- grows} else {gmx[,ss] <- NA} 
    # this if statement breaks the code (puts NA's into the matrix) if the sum of the PDF is zero (which happens if all the probability is outside of the size bounds)
  } # end ss loop
  
  # make the surv*growth matrix
  survgmx <- gmx*t(matrix( rep(sur_vals,(n.bin)),(n.bin))) # survs* growth
  mx1 <- survgmx # growth and survival matrix, without the repro
  
  #vector of reproductive rates (probability of reproduction * number of bulbils if reproducing)
  reprow <- rep_vals*reproyesnovals 
  
  # making estimated seedling sizes using the cdf function
  sdlgszcdf=pnorm(vec.bin,mnsldgsz,sdsldgsz) # get the cdf at bin edges
  sdlgszprobs=sdlgszcdf[2:length(vec.bin)]-sdlgszcdf[1:(length(vec.bin)-1)] # take the difference in the cdf between bin edges
  sdlgszprobs=sdlgszprobs/sum(sdlgszprobs) # avoid eviction by re-normalizing to sum to 1
  
  mx <- matrix(0, (n.bin+1), (n.bin+1))
  mx[2:(n.bin+1), 2:(n.bin+1)] = mx1
  mx[2:(n.bin+1),1] = surv_bulbings*sdlgszprobs # bulbling survival and growth
  mx[1,2:(n.bin+1)] = reprow*bulblings_per_bulbil  # bulbling creation
  
  
    lambdas_ipm_median[i] <- Re(eigen(mx)$values[1])
    dampratio_median[i]= damping.ratio(mx)
    lifespan_median[i] =lifespan(mx)
  print(i)
}

bistort.output.even <- as.data.frame(cbind(bin.num,lambdas_matrix, lambdas_ipm_mean,lambdas_ipm_median, dampratio_matrix,dampratio_mean,dampratio_median,lifespan_matrix,lifespan_mean,lifespan_median))
names(bistort.output.even) <- c("bin.num", 'lam.mx','lam,mn','lam.med','damp.mx','damp.mn','damp.med','life.mx','life.mn','life.med')

#Plotting########################################
# plotting lambda
dev.new(width=4, height = 8,noRStudioGD = TRUE)
par(mfrow=c(3,1))#,pty= "s")
plot(truebinsizes,lambdas_matrix, type = 'b', xlab = "# of bins in model", main="Bistorts with even class breaks", col='blue', ylim=c(0.85,1.05))
points(truebinsizes,lambdas_ipm_median, type = 'b', xlab = "# of bins in model", col ='red', lty=3, pch=8)
points(truebinsizes,lambdas_ipm_mean, type = 'b', xlab = "# of bins in model", col ='red', lty=1, pch=1)
linenames = expression('DVR','CVR-median','CVR-midpoint')
legend("bottomright",linenames, lty=c(1,3,1), lwd=2, col=c('blue','red','red'),pch=c(1,8,1), cex=1.2)

# plotting lifespan
plot(truebinsizes,lifespan_matrix, type = 'b', xlab = "# of bins in model",  ylab = "lifespan ",col='blue', ylim=c(1,200))
points(truebinsizes,lifespan_median, type = 'b', xlab = "# of bins in model", col ='red', lty=3, pch=8)
linenames = expression(matrix,ipm)
points(truebinsizes,lifespan_mean, type = 'b', xlab = "# of bins in model", col ='red', lty=1, pch=1)

# plotting dampingratio
plot(truebinsizes,dampratio_matrix, type = 'b', xlab = "# of bins in model", ylab = "damping ratio",col='blue', ylim=c(1,3))
points(truebinsizes,dampratio_median, type = 'b', xlab = "# of bins in model", col ='red', lty=3, pch=8)
linenames = expression(matrix,ipm)
points(truebinsizes,dampratio_mean, type = 'b', xlab = "# of bins in model", col ='red', lty=1, pch=1)




