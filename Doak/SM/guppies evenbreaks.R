###### This script fits discretely estimated matrix models and IPMs based on continuous
###### vital rate functions for a range of bin numbers and discretization methods
###### using the data for Trinidadian guppies. 
###### This version of the script fits with even class widths.
###### In the main text, results of the script shown here are included
###### in Figures 10, 11, and 12. 

graphics.off()
rm(list=ls())
library(dplyr)
library(popbio)
library(MuMIn)
library(binr)
library(matrixStats)
setwd("C:/Users/Dan Doak/Desktop/CurrentProjects/IPMmodels/code for guppy")

### Specify what range of classes/bins to evaluate -----------------------
bin.num <- c(3,4,5,6, 8, 10, 15, 20, 25, 30)

### Load, manipulate, and sort data --------------

#offspring size data:
allsdszs= as.numeric(read.csv('guppy offspring masses.csv', header=TRUE)[,1])

#main data file:
alldata=read.csv('guppiescleanedC.csv', header=TRUE)
#remove data that do not correspond to the intervals for whcih reproductive data exist
alldata=alldata[which(alldata$Capture.event>4 & alldata$Capture.event<15),]
alldata=alldata[which(alldata$Sex.Stage=='F'),]
# purge delta masses > 0.20, which are unrealistically high
deltamass=abs(alldata$mass1-alldata$mass0)
deltamass[is.na(deltamass)]=0
alldata=alldata[which(deltamass < 0.2),]

# check and clean out NA rows
alldata=as.data.frame(alldata)
alldata$mass0[which(alldata$mass0 == 0)] = NA
alldata$mass1[which(alldata$mass1 == 0)] = NA
alldata=alldata[which(is.na(alldata$mass0)==FALSE),]

# set min and max sizes
minsize <- min(c(allsdszs,alldata$mass0,alldata$mass1), na.rm=TRUE) - 0.001
maxsize <- max(c(allsdszs,alldata$mass0,alldata$mass1), na.rm=TRUE)+0.001

# find the mean predicted survival of new born fish to use in all models:
survmod=glm(surv~mass0 +I(mass0^2), family= "binomial",data=alldata)
invals=as.data.frame(allsdszs)#  as.data.frame(mean(allsdszs))
colnames(invals)='mass0'
survbabys= predict(survmod,invals, type='response')  
survbaby=mean(survbabys) 

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
    #print(sum(vec))
  }
  return(jj)
}#############################################

# the variable 'size' stores the main demographic data for use below 
size=alldata[,c(4,6,7,8)] # keep just columns of size, survival, and offspring
colnames(size)=c('t0','survival','t1','repro')
size$reproyesno = size$repro
size$reproyesno[size$reproyesno >0]=1 # so, turn into yes no variable
sizeforrepro=size
sizeforrepro=sizeforrepro[sizeforrepro$reproyesno==1,]

##### size density estimation for median size estimation
pdfsz=density(size$t0, n=1024, cut=0)
pdfsz2=cbind(pdfsz$x,pdfsz$y)
# this is a set of smoothed values that can be used with weightedMedian in the matrixStats package to generate a median size for each class.



#############################################

### BEGIN THE ESTIMATION FOR EACH MODEL TYPE: MATRIX = DVR IN TEXT, WHILE IPM MEAN AND IPM MEDIAN = VRF-MIDPOINT AND VRF-MEDIAN IN TEXT. 

### Matrix Models ------------------------

## make storage structures
lambdas_matrix <- rep(NA, length(bin.num))
dampratio_matrix = rep(NA, length(bin.num))
lifespan_matrix=rep(NA, length(bin.num))

mincounts=NULL

for(i in 1:length(bin.num)){# loop over different numbers of bins
  
# making classes: 
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
  # bin, survival, growth
  for(j in 1:(length(vec.bin)-1)){
    # set limits for subset according to bin breaks
    bounds <- c(vec.bin[j], vec.bin[j+1])
    # subset data according to bounds
    subset <- size[size$t0 > bounds[1] & size$t0 <= bounds[2],]
    # store number of individuals in this bin for future reference
    n[j] <- length(subset$t0)
    medians[j] <- median(subset$t0)
    # calculate survivorship for this class
    surv[j] <- sum(subset$survival) / length(subset$t0)
    # store histo as object, to access counts per bin
    histo <- hist(subset$t1, breaks = vec.bin, plot = FALSE)
    # returns the number of individuals of a certain size class
    grow[,j] <- histo$counts/sum(histo$counts) 
    reproduction[j] <- mean(subset$repro) # babies produced per mom and seen at two months
    
    totnums = totnums + sum(histo$counts)
  }
  
  # make a vector of the prob of recruit sizes: 
  sdlggrow=hist(allsdszs , breaks = vec.bin, plot = FALSE)$counts/length(allsdszs)
  
  M1 <- matrix(NA, n.bin, n.bin)   # initialize survival/growth matrix
  M <- matrix(0, (n.bin+1), (n.bin+1)) # initialize projection matrix
 
   # populate projection matrix
  for(j in 1:length(surv)){
    M1[,j] <- surv[j] * grow[,j]
    M1[,j] <- (surv[j] * grow[,j])
  }
    
    # add lines for the creation of recruits and their transition to first size class
    M[2:(n.bin+1), 2:(n.bin+1)] = M1
    M[2:(n.bin+1),1] = survbaby*sdlggrow 
    M[1,2:(n.bin+1)] = reproduction/(survbaby*2) # this accounts for the survival in the first two years and also the discounting of male offspring
  
  lambdas_matrix[i] <- lambda(M) # calls popbio f'n 'lambda' and 'damping.ratio'
  dampratio_matrix[i]=damping.ratio(M) 
  lifespan_matrix[i] =lifespan(M)
} else { #if the matrix could not be made due to low sample sizes per class
  lambdas_matrix[i]=NA
  dampratio_matrix[i] = NA
  lifespan_matrix[i] = NA
}
}


######################################
### IPM using midpoint size-----------------------------------------------------
# NOTE THAT ALL THE STEPS OF ANALYSIS ARE REPEATED FOR THE MIDPOINT AND MEDIAN APPROACHES, SO THAT THEY ARE INPEPENDENT OF ONE ANOTHER

# adding in columns for squared size effects for fitting vital rate models
size$t0sq <- size$t0^2           
sizeforrepro$t0sq = sizeforrepro$t0^2

# fitting functions for different vital rates (survival, growth, reproduction)

# prob (survival): linear function of size, quadratic function of size
sur_models <- list(glm(survival~ t0, family= "binomial",data = size),
                   glm(survival~ t0 + t0sq , family= "binomial",data =size))
# gives the value of the lowest AICc
min_AIC <- min(AICc(sur_models[[1]], sur_models[[2]])$AICc)
# gives the info for the best-fit model
bestsur <- sur_models[[which(min_AIC == AICc(sur_models[[1]], sur_models[[2]])$AICc)]] 

# growth: linear function of size, quadratic function of size, power function (A+B*(size^C))
growth_models <- list(nls(t1~ a+ b*t0, data = size, start=list(a= 1, b=1)),
                      nls(t1~ a+ b*t0 + c*t0sq, data = size, start=list(a= 1, b=1,c=1)),
                      nls(t1 ~ a + b*(t0^c), data= size, start=list(a= 1, b=1,c=1)))
# gives the value of the lowest AICc
min_AIC <- min(AICc(growth_models[[1]], growth_models[[2]],growth_models[[3]])$AICc)
# gives info for the best-fit model
bestgrowth <- growth_models[[which(min_AIC == AICc(growth_models[[1]], growth_models[[2]], growth_models[[3]])$AICc)]]

# getting residuals of growth to estimate growth variance
size$growth_residuals <- NA
size$growth_residuals[which(!is.na(size$t1) & !is.na(size$t0))] <- summary(bestgrowth)$residuals^2

# yes/no reproduction: compare linear function of size, quadratic function of size
reproyesno_models <- list(glm(reproyesno~ t0, family= "binomial",data = size),
                   glm(reproyesno~ t0 + t0sq , family= "binomial",data =size))
# gives the value of the lowest AICc
min_AIC <- min(AICc(reproyesno_models[[1]], reproyesno_models[[2]])$AICc)
# gives info for the best-fit model
bestreproyesno <- reproyesno_models[[which(min_AIC == AICc(reproyesno_models[[1]], reproyesno_models[[2]])$AICc)]] 

# reproduction if reproducing: linear function of size, quadratic function of size 
rep_models <- list(lm(repro~ t0, data = sizeforrepro),
                    lm(repro~ t0+ t0sq, data = sizeforrepro))
# gives the value of the lowest AICc
min_AIC <- min(AICc(rep_models[[1]], rep_models[[2]])$AICc)

bestrep <- rep_models[[which(min_AIC == AICc(rep_models[[1]], rep_models[[2]])$AICc)]] 

# variance in growth: uses best-fit model for growth mean: 
# compare linear function of size, quadratic function of size
# constrain to 0 intercept to prevent predictions of negative variance
var_models <- list(glm(growth_residuals~ t0-1, data= size), glm(growth_residuals~ t0 + t0sq-1, data= size))
# gives the value of the lowest AICc
min_AIC <- min(AICc(var_models[[1]], var_models[[2]])$AICc) 
# gives info for the best-fit model
bestvar <- var_models[[which(min_AIC == AICc(var_models[[1]], var_models[[2]])$AICc)]] 

# make models with varying bin numbers
lambdas_ipm_mean <- vector("numeric", length= length(bin.num))
dampratio_mean = vector("numeric", length= length(bin.num))
lifespan_mean = vector("numeric", length= length(bin.num))

truebinsizes= matrix(0,length(bin.num),1)

for (i in 1:length(bin.num)){

  vec.bin = c(minsize, minsize+1:bin.num[i]*(maxsize-minsize)*(1/bin.num[i])) 
  binmids <- 0.5*(vec.bin[2:length(vec.bin)] + vec.bin[1:(length(vec.bin)-1)])
  n.bin = length(binmids)
  truebinsizes[i] = n.bin  
  
  # constructing matrix models
  indata <- as.data.frame(cbind(binmids, binmids^2))
  names(indata) <- c("t0", "t0sq")
  
  # vital rate predictions from best models:
  sur_vals <- predict(bestsur,indata, type='response')
  reproyesnovals=predict(bestreproyesno, indata, type='response')
  rep_vals <- predict(bestrep, indata, type='response')
  rep_vals[rep_vals<0] <- 0
  growth_vals <- predict(bestgrowth,indata, type='response')
  var_vals <- predict(bestvar,indata, type='response')
  var_vals[which(var_vals<0)]=10^(-10)
  gmx <- matrix(data=NA, nrow=n.bin, ncol=n.bin)
   
  #obtain estimates of probability of new fish entering each size class: 
  sdecdf= ecdf(allsdszs)
  sdlgszcdf=sdecdf(vec.bin)
  sdlgszprobs=sdlgszcdf[2:length(vec.bin)]-sdlgszcdf[1:(length(vec.bin)-1)]
  sdlgszprobs=sdlgszprobs/sum(sdlgszprobs)
  
  # growth probs using cdf fn
  for (ss in 1:(n.bin)) {
    growcdf <- pnorm(vec.bin,growth_vals[ss],sqrt(var_vals[ss]))
    grows <- growcdf[2:length(vec.bin)]-growcdf[1:(length(vec.bin)-1)]
    if(sum(grows)>0){grows <- grows/sum(grows)
                     gmx[,ss] <- grows} else {gmx[,ss] <- NA} 
    # this if statement breaks the code (puts NA's into the matrix) if the sum of the PDF is zero (which happens if all the probability is outside of the size bounds)
  } # end ss loop
  
  # make the surv*growth mx
  survgmx <- gmx*t(matrix( rep(sur_vals,(n.bin)),(n.bin))) # survs* growth
  reprow <- rep_vals*reproyesnovals/(2*survbaby) #production of one month olds, discounted for males 
  
  mx1 <- survgmx # growth and survival, without the repro
  
  mx <- matrix(0, (n.bin+1), (n.bin+1))
   mx[2:(n.bin+1), 2:(n.bin+1)] = mx1
   mx[2:(n.bin+1),1] = sdlgszprobs*survbaby  
   mx[1,2:(n.bin+1)] = reprow
   
 
  lambdas_ipm_mean[i] <- Re(eigen(mx)$values[1])
  dampratio_mean[i]= damping.ratio(mx)
  lifespan_mean[i] =lifespan(mx)
 print(i)
}


#######################################
######################################

### IPM with median sizes-----------------------------------------------------

# adding in columns for squared size effects for fitting vital rate models
size$t0sq <- size$t0^2           
sizeforrepro$t0sq = sizeforrepro$t0^2

# fitting functions for different vital rates (survival, growth, reproduction)

# prob (survival): linear function of size, quadratic function of size
sur_models <- list(glm(survival~ t0, family= "binomial",data = size),
                   glm(survival~ t0 + t0sq , family= "binomial",data =size))
# gives the value of the lowest AICc
min_AIC <- min(AICc(sur_models[[1]], sur_models[[2]])$AICc)
# gives the info for the best-fit model
bestsur <- sur_models[[which(min_AIC == AICc(sur_models[[1]], sur_models[[2]])$AICc)]] 

# growth: linear function of size, quadratic function of size, power function (A+B*(size^C))
growth_models <- list(nls(t1~ a+ b*t0, data = size, start=list(a= 1, b=1)),
                      nls(t1~ a+ b*t0 + c*t0sq, data = size, start=list(a= 1, b=1,c=1)),
                      nls(t1 ~ a + b*(t0^c), data= size, start=list(a= 1, b=1,c=1)))
# gives the value of the lowest AICc
min_AIC <- min(AICc(growth_models[[1]], growth_models[[2]],growth_models[[3]])$AICc)
# gives info for the best-fit model
bestgrowth <- growth_models[[which(min_AIC == AICc(growth_models[[1]], growth_models[[2]], growth_models[[3]])$AICc)]]

# getting residuals of growth to estimate growth variance
size$growth_residuals <- NA
size$growth_residuals[which(!is.na(size$t1) & !is.na(size$t0))] <- summary(bestgrowth)$residuals^2

# yes/no reproduction: compare linear function of size, quadratic function of size
reproyesno_models <- list(glm(reproyesno~ t0, family= "binomial",data = size),
                          glm(reproyesno~ t0 + t0sq , family= "binomial",data =size))
# gives the value of the lowest AICc
min_AIC <- min(AICc(reproyesno_models[[1]], reproyesno_models[[2]])$AICc)
# gives info for the best-fit model
bestreproyesno <- reproyesno_models[[which(min_AIC == AICc(reproyesno_models[[1]], reproyesno_models[[2]])$AICc)]] 

# reproduction if reproducing: linear function of size, quadratic function of size 
rep_models <- list(lm(repro~ t0, data = sizeforrepro),
                   lm(repro~ t0+ t0sq, data = sizeforrepro))
# gives the value of the lowest AICc
min_AIC <- min(AICc(rep_models[[1]], rep_models[[2]])$AICc)

bestrep <- rep_models[[which(min_AIC == AICc(rep_models[[1]], rep_models[[2]])$AICc)]] 

# variance in growth: uses best-fit model for growth mean: 
# compare linear function of size, quadratic function of size
# constrain to 0 intercept to prevent predictions of negative variance
var_models <- list(glm(growth_residuals~ t0-1, data= size), glm(growth_residuals~ t0 + t0sq-1, data= size))
# gives the value of the lowest AICc
min_AIC <- min(AICc(var_models[[1]], var_models[[2]])$AICc) 
# gives info for the best-fit model
bestvar <- var_models[[which(min_AIC == AICc(var_models[[1]], var_models[[2]])$AICc)]] 

# make models with varying bin numbers
lambdas_ipm_median <- vector("numeric", length= length(bin.num))
dampratio_median <- vector("numeric", length= length(bin.num))
lifespan_median = vector("numeric", length= length(bin.num))

truebinsizes= matrix(0,length(bin.num),1)

for (i in 1:length(bin.num)){

 vec.bin = c(minsize, minsize+1:bin.num[i]*(maxsize-minsize)*(1/bin.num[i])) 
 nums=hist(size$t0,breaks=vec.bin, plot=FALSE)$counts 
 binmids =   rep(NA, length(vec.bin)-1)
 for(jj in 1:(length(vec.bin)-1)){
    # set limits for subset according to bin breaks
    bounds <- c(vec.bin[jj], vec.bin[jj+1])
    # subset data according to bounds
     subsetszs <- pdfsz2[which(pdfsz2[,1] >= bounds[1] & pdfsz2[,1] < bounds[2]),]
     binmids[jj] <- weightedMedian(subsetszs[,1],subsetszs[,2])
 }
 
 # if there is failure in the median estimation, set the size for these bins to the midpoint:
  meanbinmids <- 0.5*(vec.bin[2:length(vec.bin)] + vec.bin[1:(length(vec.bin)-1)])
  binmids[is.na(binmids)]=meanbinmids[is.na(binmids)]
 
  n.bin = length(binmids)
  truebinsizes[i] = n.bin  
  
  # constructing matrix models
  indata <- as.data.frame(cbind(binmids, binmids^2))
  names(indata) <- c("t0", "t0sq")
  
  # vital rate predictions from best models:
  sur_vals <- predict(bestsur,indata, type='response')
  reproyesnovals=predict(bestreproyesno, indata, type='response')
  rep_vals <- predict(bestrep, indata, type='response')
  rep_vals[rep_vals<0] <- 0
  growth_vals <- predict(bestgrowth,indata, type='response')
  var_vals <- predict(bestvar,indata, type='response')
  var_vals[which(var_vals<0)]=10^(-10)
  gmx <- matrix(data=NA, nrow=n.bin, ncol=n.bin)
  
  #obtain estimates of probability of new fish entering each size class: 
  sdecdf= ecdf(allsdszs)
  sdlgszcdf=sdecdf(vec.bin)
  sdlgszprobs=sdlgszcdf[2:length(vec.bin)]-sdlgszcdf[1:(length(vec.bin)-1)]
  sdlgszprobs=sdlgszprobs/sum(sdlgszprobs)
  
  # growth probs using cdf fn
  for (ss in 1:(n.bin)) {
    growcdf <- pnorm(vec.bin,growth_vals[ss],sqrt(var_vals[ss]))
    grows <- growcdf[2:length(vec.bin)]-growcdf[1:(length(vec.bin)-1)]
    if(sum(grows)>0){grows <- grows/sum(grows)
    gmx[,ss] <- grows} else {gmx[,ss] <- NA} 
    # this if statement breaks the code (puts NA's into the matrix) if the sum of the PDF is zero (which happens if all the probability is outside of the size bounds)
  } # end ss loop
  
  # make the surv*growth mx
  survgmx <- gmx*t(matrix( rep(sur_vals,(n.bin)),(n.bin))) # survs* growth
  reprow <- rep_vals*reproyesnovals/(2*survbaby) #production of one month olds, discounted for males 
  
  mx1 <- survgmx # growth and survival, without the repro
  
  mx <- matrix(0, (n.bin+1), (n.bin+1))
  mx[2:(n.bin+1), 2:(n.bin+1)] = mx1
  mx[2:(n.bin+1),1] = sdlgszprobs*survbaby  
  mx[1,2:(n.bin+1)] = reprow
  
  lambdas_ipm_median[i] <- Re(eigen(mx)$values[1])
  dampratio_median[i]= damping.ratio(mx)
  lifespan_median[i] =lifespan(mx)
  print(i)
}
mxmedian=mx

guppy.output.even <- as.data.frame(cbind(bin.num,lambdas_matrix, lambdas_ipm_mean,lambdas_ipm_median, dampratio_matrix,dampratio_mean,dampratio_median,lifespan_matrix,lifespan_mean,lifespan_median))
names(guppy.output.even) <- c("bin.num", 'lam.mx','lam,mn','lam.med','damp.mx','damp.mn','damp.med','life.mx','life.mn','life.med')


#Plotting########################################
# plotting lambda
dev.new(width=4, height = 8,noRStudioGD = TRUE)
par(mfrow=c(3,1))#,pty= "s")
plot(truebinsizes,lambdas_matrix, type = 'b', xlab = "# of bins in model", main="Guppies with even class breaks", col='blue', ylim=c(0.8,1.0))
points(truebinsizes,lambdas_ipm_median, type = 'b', xlab = "# of bins in model", col ='red', lty=3, pch=8)
points(truebinsizes,lambdas_ipm_mean, type = 'b', xlab = "# of bins in model", col ='red', lty=1, pch=1)
linenames = expression('DVR','CVR-median','CVR-midpoint')
legend("bottomright",linenames, lty=c(1,3,1), lwd=2, col=c('blue','red','red'),pch=c(1,8,1), cex=1.2)

# plotting lifespan
plot(truebinsizes,lifespan_matrix, type = 'b', xlab = "# of bins in model",  ylab = "lifespan ",col='blue', ylim=c(1,40))
points(truebinsizes,lifespan_median, type = 'b', xlab = "# of bins in model", col ='red', lty=3, pch=8)
linenames = expression(matrix,ipm)
points(truebinsizes,lifespan_mean, type = 'b', xlab = "# of bins in model", col ='red', lty=1, pch=1)

# plotting dampingratio
plot(truebinsizes,dampratio_matrix, type = 'b', xlab = "# of bins in model", ylab = "damping ratio",col='blue', ylim=c(1,3))
points(truebinsizes,dampratio_median, type = 'b', xlab = "# of bins in model", col ='red', lty=3, pch=8)
linenames = expression(matrix,ipm)
points(truebinsizes,dampratio_mean, type = 'b', xlab = "# of bins in model", col ='red', lty=1, pch=1)
