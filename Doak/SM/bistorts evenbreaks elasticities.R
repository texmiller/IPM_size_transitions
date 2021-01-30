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
bin.num <- seq(from = 20, to = 100, by=20); bin.num=round(bin.num); 
nSizes = length(bin.num); 
kernels = list(nSizes); 
all_meshpts = list(nSizes); 


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
 
  kernels[[i]]= mx; all_meshpts[[i]] = binmids; 
  Pi = mx; Pi[1,]=0; N=solve(diag(ncol(Pi))-Pi);
 cat(i,sum(N[,1]),"\n")
}


allElas = matrix(NA,400,nSizes)
K = kernels[[1]]; eK = eigen(K); 
lambda=abs(eK$values[1]); 
w = eK$vectors[,1]; w = abs(w)/sum(abs(w)); 
v = eigen(t(K))$vectors[,1]; v=abs(v); 
sens = outer(v,w)/sum(v*w); 
elas = (K/lambda)*sens; 
allElas[,1]=matrix(elas,400,1); 


for(js in 1:nSizes) {
    K = kernels[[js]]; eK = eigen(K); 
    lambda=abs(eK$values[1]); 
    w = eK$vectors[,1]; w = abs(w)/sum(abs(w)); 
    v = eigen(t(K))$vectors[,1]; v=abs(v); 
    sens = outer(v,w)/sum(v*w); 
    elas = (K/lambda)*sens; 
    binElas = matrix(NA,20,20)
    for(i in 1:20){
    for(j in 1:20){
            i.end = i*js; i.start = (i-1)*js + 1; 
            j.end = j*js; j.start = (j-1)*js + 1; 
            ebin = elas[i.start:i.end,j.start:j.end]
            binElas[i,j]=sum(ebin)
    } }       
    allElas[,js]=matrix(binElas,400,1); 
}

require(viridis); colors=magma(nSizes); 
K = kernels[[1]]; eK = eigen(K); 
lambda=abs(eK$values[1]); w = eK$vectors[,1]; w = abs(w)/sum(abs(w)); 
v = eigen(t(K))$vectors[,1]; v=abs(v); 
sens = outer(v,w)/sum(v*w); 
elas = (K/lambda)*sens; 
meshpts = all_meshpts[[1]]; 
h=diff(meshpts)[1]; 
intElas = apply(elas,2,sum);  
plot(c(meshpts[1]-h,meshpts), intElas,type="l",xlab="Size", ylab="Integrated elasticity",col=colors[1],ylim=c(0,1.2*max(intElas)));  
z0 = meshpts[1]-h; 

for(j in 2:nSizes) {
    K = kernels[[j]]; eK = eigen(K); 
    lambda=abs(eK$values[1]); w = eK$vectors[,1]; w = abs(w)/sum(abs(w)); 
    v = eigen(t(K))$vectors[,1]; v=abs(v); 
    sens = outer(v,w)/sum(v*w); 
    elas = (K/lambda)*sens; 
    meshpts = all_meshpts[[j]]; 
    h=diff(meshpts)[1]; 
    intElas = apply(elas,2,sum); intElas[-1]=intElas[-1]*length(meshpts)/length(all_meshpts[[1]]); 
    points(c(z0,meshpts),intElas,type="l",xlab="Size", ylab="Integrated elasticity",col=colors[j]);  
}



