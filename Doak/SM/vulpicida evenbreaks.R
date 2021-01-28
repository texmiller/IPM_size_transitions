###### This script fits discretely estimated matrix models and IPMs based on continuous
###### vital rate functions for a range of bin numbers and discretization methods
###### using the data for Vulpicida pinastri, and epiphytic lichen. 
###### This version of the script fits with even class widths, the recommended approach 
###### fitting models. In the main text, results of the script shown here are included
###### in Figures 10, 11, and 12. 

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

dat <- read.csv("lichen_demo.csv") # main data file
recruitsurv = 0.4 #set to givelambda ~= 1
recruitszlims=c(0.124,0.50) # assign new thalli to this size range, using uniform distr.

dat.size <- as.data.frame(dat[, 5:6])
colnames(dat.size) <- c("sqrtszt"="t0", "sqrtszt1"="t1") # rename column header
size <- arrange(dat.size, t0) 
size$repro <- size$t0 * 0.047*2*(pi^0.5) # reproductive rate is a fn of size: See Shriver et al. paper for details. 
size$survival <- 1
size$survival[which(size$t1==0)] <- 0
size$repro[which(size$t1==0)] = 0
size$t1[which(size$t1==0)] <- NA
minsize <- 0 
maxsize=7 # this is set manually to not have a large category at end that has no plants.
size$t0[which(size$t0>maxsize)]=maxsize-0.1
size$t1[which(size$t1>maxsize)]=maxsize-0.1

#generate empirical density function for median size estimation
#this is a set of smoothed values that can then be used with weightedMedian
#in the matrixStats package to get an estimated median size for each class.
pdfsz=density(size$t0, n=1024, cut=0)
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
}#############################################


### BEGIN THE ESTIMATION FOR EACH MODEL TYPE: MATRIX = DVR IN TEXT, WHILE IPM MEAN AND IPM MEDIAN = CVR-MIDPOINT AND CVR-MEDIAN IN TEXT. 

### Matrix Models ------------------------

## make storage structures
lambdas_matrix <- rep(NA, length(bin.num))
dampratio_matrix = rep(NA, length(bin.num))
lifespan_matrix=rep(NA, length(bin.num))

mincounts=NULL

for(i in 1:length(bin.num)){

  ss=as.numeric(size$t0)
  vec.bin = c(minsize, minsize+1:bin.num[i]*(maxsize-minsize)*(1/bin.num[i]))  
  nums=hist(ss,breaks=vec.bin, plot=FALSE)$counts 
  mincounts=c(mincounts,min(nums))
  
  if (min(nums)>2) {  
  ## Initialize storage
  n.bin <- length(vec.bin)-1                  
  n <- rep(NA, n.bin)                         # count of indivs per bin
  medians <- rep(NA, n.bin)                    # median size per bin for F
  surv <- rep(NA, n.bin)                      # survivorship for each class
  grow <- matrix(NA, n.bin, n.bin)            # store growth probabilities for each class
  reproduction <- rep(NA, n.bin)
  
  totnums = 0 # this is monitoring for errors
  # bin data to estimate survival, growth
  for(j in 1:(length(vec.bin)-1)){
    bounds <- c(vec.bin[j], vec.bin[j+1])
   subset <- size[size$t0 > bounds[1] & size$t0 <= bounds[2],]
     n[j] <- length(subset$t0)
    medians[j] <- median(subset$t0)
    # calculate survivorship for this class
    surv[j] <- sum(subset$survival) / length(subset$t0)
    histo <- hist(subset$t1, breaks = vec.bin, plot = FALSE)
    # returns the number of individuals of a certain size class
    grow[,j] <- histo$counts/length(subset$t0[subset$survival==1]) 
    reproduction[j] <- mean(subset$repro)
    
    totnums = totnums + sum(histo$counts)
  }
  
  M1 <- matrix(NA, n.bin, n.bin)  # initiate projection matrix
  M <- matrix(0, (n.bin+1), (n.bin+1)) # populate projection matrix

  for(j in 1:length(surv)){
    M1[,j] <- surv[j] * grow[,j]
    M1[,j] <- (surv[j] * grow[,j])
  }
    
    # add lines for the creation of new thalli and their transition to first size class
    M[2:(n.bin+1), 2:(n.bin+1)] = M1
    
    sdlgszcdf=punif(vec.bin,min=recruitszlims[1],max=recruitszlims[2])
    sdlgszprobs=sdlgszcdf[2:length(vec.bin)]-sdlgszcdf[1:(length(vec.bin)-1)]
    sdlgszprobs=sdlgszprobs/sum(sdlgszprobs)
    M[2:(n.bin+1),1] = recruitsurv*sdlgszprobs
    M[1,2:(n.bin+1)] = reproduction*surv 
  
  lambdas_matrix[i] <- lambda(M) 
  dampratio_matrix[i]=damping.ratio(M)
  lifespan_matrix[i] =lifespan(M)
   } else {
    lambdas_matrix[i]=NA
    dampratio_matrix[i]=NA
    lifespan_matrix[i] =NA
   }
   }


######################################
### IPM using midpoint size-----------------------------------------------------
# NOTE THAT ALL THE STEPS OF ANALYSIS ARE REPEATED FOR THE MIDPOINT AND MEDIAN APPROACHES, SO THAT THEY ARE INPEPENDENT OF ONE ANOTHER

# adding in column for squared size effects for fitting vital rate models
size$t0sq <- size$t0^2             
# fitting functions for different vital rates (survival, growth, reproduction)

# prob (survival): 
sur_models <- list(glm(survival~ t0, family= "binomial",data = size),
                   glm(survival~ t0 + t0sq , family= "binomial",data =size))
# gives the value of the lowest AICc
min_AIC <- min(AICc(sur_models[[1]], sur_models[[2]])$AICc)
# gives info for the best-fit model
bestsur <- sur_models[[which(min_AIC == AICc(sur_models[[1]], sur_models[[2]])$AICc)]] 

# growth
growth_models <- list(nls(t1~ a+ b*t0, data = size, start=list(a= 1, b=1)),
                      nls(t1~ a+ b*t0 + c*t0sq, data = size, start=list(a= 1, b=1,c=1)),
                      nls(t1 ~ a + b*(t0^c), data= size, start=list(a= 1, b=1,c=1)))
min_AIC <- min(AICc(growth_models[[1]], growth_models[[2]],growth_models[[3]])$AICc)
bestgrowth <- growth_models[[which(min_AIC == AICc(growth_models[[1]], growth_models[[2]], growth_models[[3]])$AICc)]]

# obtain residuals of growth
size$growth_residuals <- NA
size$growth_residuals[which(!is.na(size$t1) & !is.na(size$t0))] <- summary(bestgrowth)$residuals^2

# reproduction: linear function of size, quadratic function of size, power function (A+B*(size^C))
rep_models <- list(lm(repro~ t0, data = size),
                    lm(repro~ t0+ t0sq, data = size)) 
min_AIC <- min(AICc(rep_models[[1]], rep_models[[2]]))
bestrep <- rep_models[[which(min_AIC == AICc(rep_models[[1]], rep_models[[2]])$AICc)]] 

# variance in growth: 
var_models <- list(glm(growth_residuals~ t0-1, data= size), glm(growth_residuals~ t0 + t0sq-1, data= size))
min_AIC <- min(AICc(var_models[[1]], var_models[[2]])$AICc) 
bestvar <- var_models[[which(min_AIC == AICc(var_models[[1]], var_models[[2]])$AICc)]]

# make models with varying bin numbers

lambdas_ipm_mean <- vector("numeric", length= length(bin.num))
dampratio_mean = vector("numeric", length= length(bin.num))
lifespan_mean = vector("numeric", length= length(bin.num))

truebinsizes= matrix(0,length(bin.num),1)

for (i in 1:length(bin.num)){

  
  ss=as.numeric(size$t0)
  vec.bin = c(minsize, minsize+1:bin.num[i]*(maxsize-minsize)*(1/bin.num[i])) 
  nums=hist(ss,breaks=vec.bin, plot=FALSE)$counts 
  binmids <- 0.5*(vec.bin[2:length(vec.bin)] + vec.bin[1:(length(vec.bin)-1)])
  n.bin = length(binmids)
  truebinsizes[i] = n.bin  
  
  # constructing matrix models
  indata <- as.data.frame(cbind(binmids, binmids^2))
  names(indata) <- c("t0", "t0sq")
  
  # vital rate predictions from best models:
  sur_vals <- predict(bestsur,indata, type='response')
  rep_vals <- predict(bestrep, indata, type='response')
  rep_vals[rep_vals<0] <- 0
  growth_vals <- predict(bestgrowth,indata, type='response')
  var_vals <- predict(bestvar,indata, type='response')
  gmx <- matrix(data=NA, nrow=n.bin, ncol=n.bin)
   
  # growth probabilities using cdf function
  for (ss in 1:(n.bin)) {
    growcdf <- pnorm(vec.bin,growth_vals[ss],sqrt(var_vals[ss]))
    grows <- growcdf[2:length(vec.bin)]-growcdf[1:(length(vec.bin)-1)]
    if(sum(grows)>0){grows <- grows/sum(grows)
                     gmx[,ss] <- grows} else {gmx[,ss] <- NA} 
  } # end ss loop
  
  # make the surv*growth mx
  survgmx <- gmx*t(matrix( rep(sur_vals,(n.bin)),(n.bin))) # survs* growth
  reprow <- rep_vals*sur_vals
 
  mx1 <- survgmx # growth and survival, without repro
  
  mx <- matrix(0, (n.bin+1), (n.bin+1))
   mx[2:(n.bin+1), 2:(n.bin+1)] = mx1
   
  # this section assigns the new recruits to the range of estimated sizes, uniform with a min and max
     sdlgszcdf=punif(vec.bin,min=recruitszlims[1],max=recruitszlims[2])
    sdlgszprobs=sdlgszcdf[2:length(vec.bin)]-sdlgszcdf[1:(length(vec.bin)-1)]
    sdlgszprobs=sdlgszprobs/sum(sdlgszprobs)
   
  mx[2:(n.bin+1),1] = recruitsurv*sdlgszprobs
  mx[1,2:(n.bin+1)] = reprow 
  
  lambdas_ipm_mean[i] <- Re(eigen(mx)$values[1])
  dampratio_mean[i]= damping.ratio(mx)
  lifespan_mean[i] =lifespan(mx)
  print(i)
}



#############################################
################################################
### IPM with median sizes t

# adding in columns for squared size effects
size$t0sq <- size$t0^2 

# fitting functions for different vital rates (survival, growth, reproduction)

# prob (survival): 
sur_models <- list(glm(survival~ t0, family= "binomial",data = size),
                   glm(survival~ t0 + t0sq , family= "binomial",data =size))
min_AIC <- min(AICc(sur_models[[1]], sur_models[[2]])$AICc)
bestsur <- sur_models[[which(min_AIC == AICc(sur_models[[1]], sur_models[[2]])$AICc)]] 

# growth: 
growth_models <- list(nls(t1~ a+ b*t0, data = size, start=list(a= 1, b=1)),
                      nls(t1~ a+ b*t0 + c*t0sq, data = size, start=list(a= 1, b=1,c=1)),
                      nls(t1 ~ a + b*(t0^c), data= size, start=list(a= 1, b=1,c=1)))
min_AIC <- min(AICc(growth_models[[1]], growth_models[[2]],growth_models[[3]])$AICc)
bestgrowth <- growth_models[[which(min_AIC == AICc(growth_models[[1]], growth_models[[2]], growth_models[[3]])$AICc)]]

# get residuals of growth:
size$growth_residuals <- NA
size$growth_residuals[which(!is.na(size$t1) & !is.na(size$t0))] <- summary(bestgrowth)$residuals^2

# reproduction: 
rep_models <- list(lm(repro~ t0, data = size),
                   lm(repro~ t0+ t0sq, data = size)) 
min_AIC <- min(AICc(rep_models[[1]], rep_models[[2]]))
bestrep <- rep_models[[which(min_AIC == AICc(rep_models[[1]], rep_models[[2]])$AICc)]] 

# variance in growth: 
var_models <- list(glm(growth_residuals~ t0-1, data= size), glm(growth_residuals~ t0 + t0sq-1, data= size))
min_AIC <- min(AICc(var_models[[1]], var_models[[2]])$AICc) 
bestvar <- var_models[[which(min_AIC == AICc(var_models[[1]], var_models[[2]])$AICc)]] 

# make models with varying bin numbers

lambdas_ipm_median <- vector("numeric", length= length(bin.num))
dampratio_median <- vector("numeric", length= length(bin.num))
lifespan_median = vector("numeric", length= length(bin.num))

truebinsizes= matrix(0,length(bin.num),1)

for (i in 1:length(bin.num)){
  
  ss=as.numeric(size$t0)
  vec.bin = c(minsize, minsize+1:bin.num[i]*(maxsize-minsize)*(1/bin.num[i])) 
  nums=hist(ss,breaks=vec.bin, plot=FALSE)$counts 
  binmids =   rep(NA, length(vec.bin)-1)
  for(jj in 1:(length(vec.bin)-1)){
    # set limits for subset according to bin breaks
    bounds <- c(vec.bin[jj], vec.bin[jj+1])
    # subset data according to bounds
    subsetszs <- pdfsz2[which(pdfsz2[,1] >= bounds[1] & pdfsz2[,1] < bounds[2]),]
    if (length(which(pdfsz2[,1] >= bounds[1] & pdfsz2[,1] < bounds[2]))>1) { binmids[jj] <- weightedMedian(subsetszs[,1],subsetszs[,2])} else {binmids[jj]=subsetszs[1]}
    
  }
  meanbinmids <- 0.5*(vec.bin[2:length(vec.bin)] + vec.bin[1:(length(vec.bin)-1)])
  binmids[is.na(binmids)]=meanbinmids[is.na(binmids)]

  n.bin = length(binmids)
  truebinsizes[i] = n.bin  
  
  # constructing the matrices
  indata <- as.data.frame(cbind(binmids, binmids^2))
  names(indata) <- c("t0", "t0sq")
  
  # vital rate predictions from best models:
  sur_vals <- predict(bestsur,indata, type='response')
  rep_vals <- predict(bestrep, indata, type='response')
  rep_vals[rep_vals<0] <- 0
  growth_vals <- predict(bestgrowth,indata, type='response')
  var_vals <- predict(bestvar,indata, type='response')
  gmx <- matrix(data=NA, nrow=n.bin, ncol=n.bin)
  
  # growth probs using cdf fn
  for (ss in 1:(n.bin)) {
    growcdf <- pnorm(vec.bin,growth_vals[ss],sqrt(var_vals[ss]))
    grows <- growcdf[2:length(vec.bin)]-growcdf[1:(length(vec.bin)-1)]
    if(sum(grows)>0){grows <- grows/sum(grows)
    gmx[,ss] <- grows} else {gmx[,ss] <- NA} 
      } # end ss loop
  
  # make the surv*growth mx
  survgmx <- gmx*t(matrix( rep(sur_vals,(n.bin)),(n.bin))) # survs* growth
  reprow <- rep_vals*sur_vals
  
  mx1 <- survgmx # growth and survival, without repro
  
  mx <- matrix(0, (n.bin+1), (n.bin+1))
  mx[2:(n.bin+1), 2:(n.bin+1)] = mx1
  
  # this section assigns the new recruits to the range of estimated sizes, uniform with a min and max
  sdlgszcdf=punif(vec.bin,min=recruitszlims[1],max=recruitszlims[2])
  sdlgszprobs=sdlgszcdf[2:length(vec.bin)]-sdlgszcdf[1:(length(vec.bin)-1)]
  sdlgszprobs=sdlgszprobs/sum(sdlgszprobs)
  
  mx[2:(n.bin+1),1] = recruitsurv*sdlgszprobs
  mx[1,2:(n.bin+1)] = reprow #  or recruit to it.

  lambdas_ipm_median[i] <- Re(eigen(mx)$values[1])
  dampratio_median[i]= damping.ratio(mx)
  lifespan_median[i] =lifespan(mx)
  print(i)
} # end loop 


vulpicida.output.even <- as.data.frame(cbind(bin.num,lambdas_matrix, lambdas_ipm_mean,lambdas_ipm_median, dampratio_matrix,dampratio_mean,dampratio_median,lifespan_matrix,lifespan_mean,lifespan_median))
names(vulpicida.output.even) <- c("bin.num", 'lam.mx','lam.mn','lam.med','damp.mx','damp.mn','damp.med','life.mx','life.mn','life.med')


#Plotting########################################
# plotting lambda
dev.new(width=4, height = 8,noRStudioGD = TRUE)
par(mfrow=c(3,1))#,pty= "s")
plot(truebinsizes,lambdas_matrix, type = 'b', xlab = "# of bins in model", main="Vulpicida with even class breaks", col='blue', ylim=c(0.85,1.1))
points(truebinsizes,lambdas_ipm_median, type = 'b', xlab = "# of bins in model", col ='red', lty=3, pch=8)
points(truebinsizes,lambdas_ipm_mean, type = 'b', xlab = "# of bins in model", col ='red', lty=1, pch=1)
linenames = expression('DVR','CVR-median','CVR-midpoint')
legend("bottomright",linenames, lty=c(1,3,1), lwd=2, col=c('blue','red','red'),pch=c(1,8,1), cex=1.2)

# plotting lifespan
plot(truebinsizes,lifespan_matrix, type = 'b', xlab = "# of bins in model",  ylab = "lifespan ",col='blue', ylim=c(1,100))
points(truebinsizes,lifespan_median, type = 'b', xlab = "# of bins in model", col ='red', lty=3, pch=8)
linenames = expression(matrix,ipm)
points(truebinsizes,lifespan_mean, type = 'b', xlab = "# of bins in model", col ='red', lty=1, pch=1)

# plotting dampingratio
plot(truebinsizes,dampratio_matrix, type = 'b', xlab = "# of bins in model", ylab = "damping ratio",col='blue', ylim=c(1,2))
points(truebinsizes,dampratio_median, type = 'b', xlab = "# of bins in model", col ='red', lty=3, pch=8)
linenames = expression(matrix,ipm)
points(truebinsizes,dampratio_mean, type = 'b', xlab = "# of bins in model", col ='red', lty=1, pch=1)




