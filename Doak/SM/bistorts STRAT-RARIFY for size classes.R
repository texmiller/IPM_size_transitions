###### This script fits discretely estimated matrix models and median-size based IPMs based on continuous
###### vital rate functions for rarefied data and using a range of size class numbers, 
###### to explore issues with low sample sizes per class, 
###### using data for alpine bistort (Polygonum viviparum). 
###### This version of the script fits with even class widths, the recommended approach. 
###### Results of the script shown here are included
###### in Figures 15 and S21.

graphics.off()
rm(list=ls())
library(dplyr)
library(popbio)
library(MuMIn)
library(binr)
library(beanplot)
library(matrixStats)
library(splitstackshape)

#Remember to set directory to access input files

### Specify what range of classes/bins to evaluate -----------------------
allbinnums=rev(c(5,10,15,20,25,30,35,40,45,50))

reps = 200 # number of reps at each ratification number

### Load, manipulate, and sort data --------------

allsdszs= as.numeric(read.csv('bulblingsizes.csv', header=FALSE)[,1])
mnsldgsz=mean(allsdszs)
sdsldgsz = sqrt(var(allsdszs))

alldata=as.matrix(read.csv('BNdataforIPMsimple.csv', header=TRUE))
alldata=as.data.frame(alldata)
alldata=alldata[which(is.na(alldata[,1])==FALSE),]

alldata$survival = 1
alldata$survival[which(alldata$szs1==0)] = 0
alldata$szs1[which(alldata$szs1==0)] = NA
alldata=alldata[which(alldata$szs0 !=0),]

minsize <- min(c(alldata[,1],alldata[,2]), na.rm=TRUE)

bulblings_per_bulbil = 0.00676 

#get total sample size for stratified rarefaction:
totsample=length(alldata[,1])
sampsizes=round(totsample/c(16))
sampfraction = sampsizes/totsample  
sampfraction=sort(sampfraction, decreasing=FALSE)

#### find the class boundaries and make a factor to define size classes for stratified bootstrap:
szs=sort(alldata$szs0)
szmin=min(szs)-0.1
szmax=max(szs)+0.1
sz0595=quantile(szs,c(0.05,0.95))
sbrks=  seq(from=sz0595[1],to=sz0595[2],length.out=19) # these are the breaks for making 20 class models

alldata$szsfactor=cut(alldata$szs0,breaks=c(szmin,sbrks,szmax), labels = FALSE)# 
sbrks=c(szmin,sbrks,szmax)

bulblingsz=sqrt(4.388*4)
alltinys=alldata[alldata$szs0<=bulblingsz,]
surv_bulbings=sum(alltinys$survival)/length(alltinys$survival)

minsize <- min(c(alldata[,1],alldata[,2]), na.rm=TRUE)
maxsize <- max(c(alldata[,1],alldata[,2]), na.rm=TRUE) +0.1
maxsize=60 # this is set manually to not have a large category at end that has no plants.
alldata[which(alldata[,1]>maxsize),1]=maxsize
alldata[which(alldata[,2]>maxsize),2]=maxsize

####################################################
lifespan <- function(nx){
  nclasses=dim(nx)[1]
  vec=c(100,rep(0,(nclasses-1)))
  nx[1,]=0
  jj=1
  while (sum(vec)>1 & jj<1001){
    vec=nx%*%vec
    jj=jj+1
    #print(sum(vec))
  }
  return(jj)
}
#############################################

# Make storage structures: 
all_lams_ipm = NULL 
all_lams_mx = NULL
alldat = alldata # to store the full set of data
alloutdata = NULL

#loop across the number of size classes used:
for (bin.num in allbinnums) {

  ## Loops which create a stratified subset of the data with a given sample size: 
  ## Note that as run for paper, there is only one sample size used
for (sams in sampfraction) {  
  for (rr in 1:reps) {
    alldata=stratified(alldat, "szsfactor",sams,replace=TRUE)
    print(c(totsample*sams,rr, length(alldata$szs0)))
    
size=alldata
colnames(size)=c('t0','t1','bulbs0','survival', 'szsfactor')

size$reproyesno = size$bulbs0
size$reproyesno[size$reproyesno >0]=1 
sizeforrepro=size
sizeforrepro=sizeforrepro[sizeforrepro$reproyesno==1,]

#####new size density estimation for median size estimation
pdfsz=density(size$t0, n=1024, cut=0)
pdfsz2=cbind(pdfsz$x,pdfsz$y)


### Matrix Model ------------------------

## make storage structures
lambdas_matrix <- rep(NA, length(bin.num))
dampratio_matrix = rep(NA, length(bin.num))
lifespan_matrix=rep(NA, length(bin.num))

sz0595now=quantile(alldata$szs0,c(1/bin.num,(1-1/bin.num))) # make these quantiles match the bin numbers
sbrksnow=  seq(from=sz0595now[1],to=sz0595now[2],length.out=(bin.num-1)) 

  vec.bin = c(szmin,sbrksnow,szmax) 
  ss=as.numeric(size$t0)
  nums=hist(ss,breaks=vec.bin, plot=FALSE)$counts 


  ## Initialize storage
  n.bin <- length(vec.bin)-1                  # this is a workaround
  n <- rep(NA, n.bin)                         # count of indvs per bin
  medians <- rep(NA, n.bin)                    # median size per bin for F
  surv <- rep(NA, n.bin)                      # survivorship for each class
  grow <- matrix(NA, n.bin, n.bin)            # store growth probabilities for each class
  reproduction <- rep(NA, n.bin)
  
  totnums = 0 # this is monitoring for errors
  # bin, survival, growth
  for(j in 1:(length(vec.bin)-1)){
    bounds <- c(vec.bin[j], vec.bin[j+1])
    subset <- size[size$t0 > bounds[1] & size$t0 <= bounds[2],]
    n[j] <- length(subset$t0)
    medians[j] <- median(subset$t0)
    surv[j] <- sum(subset$survival) / length(subset$t0)
    histo <- hist(subset$t1, breaks = vec.bin, plot = FALSE)
    grow[,j] <- histo$counts/length(subset$t0[subset$survival==1]) 
    reproduction[j] <- mean(subset$bulbs0) # bulbils produced per plant in start yr
    totnums = totnums + sum(histo$counts)
  }
  
  # make a vector of the prob of seedling sizes: 
  sdlggrow=hist(allsdszs , breaks = vec.bin, plot = FALSE)$counts/length(allsdszs)
  
  M1 <- matrix(NA, n.bin, n.bin)   # initiate projection matrix
  M <- matrix(0, (n.bin+1), (n.bin+1))

  for(j in 1:length(surv)){
    M1[,j] <- surv[j] * grow[,j]
    M1[,j] <- (surv[j] * grow[,j])
  }
    
    M[2:(n.bin+1), 2:(n.bin+1)] = M1
    M[2:(n.bin+1),1] = surv_bulbings*sdlggrow 
    M[1,2:(n.bin+1)] = reproduction*bulblings_per_bulbil  
  
  if (length(which(is.na(M)==TRUE))==0) {
  lambdas_matrix <- lambda(M) # calls popbio f'n 'lambda'
  lifespan_matrix <- lifespan(M) # calls popbio f'n 'lambda'
  dampratio_matrix=damping.ratio(M) } else {
    lambdas_matrix <- NA 
    lifespan_matrix <- NA 
    dampratio_matrix=NA } 
    
  alloutdata = rbind(alloutdata,cbind(sams,rr,bin.num,0,  lambdas_matrix,lifespan_matrix,dampratio_matrix,  dim(M)[1],min(nums), median(nums), length(nums[nums<6]),(length(which(is.na(M)==TRUE))==0) )  )
  

  #print(c(i,totnums))
#}


#IPM with 20 sizes-----------------------------------------------------

  
size$t0sq <- size$t0^2           
sizeforrepro$t0sq = sizeforrepro$t0^2
# fitting functions for different vital rates (survival, growth, reproduction)

# prob (survival): 
sur_models <- list(glm(survival~ t0, family= "binomial",data = size),
                   glm(survival~ t0 + t0sq , family= "binomial",data =size))
min_AIC <- min(AICc(sur_models[[1]], sur_models[[2]])$AICc)
bestsur <- sur_models[[which(min_AIC == AICc(sur_models[[1]], sur_models[[2]])$AICc)]] 

# growth: 
growth_models <- list(nls(t1~ a+ b*t0, data = size, start=list(a= 1, b=1), control=nls.control( maxiter = 200,warnOnly = TRUE)),
                      nls(t1~ a+ b*t0 + c*t0sq, data = size, start=list(a= 1, b=1,c=1),control=nls.control( maxiter = 200,warnOnly = TRUE)),
                      nls(t1 ~ a + b*(t0^c), data= size, start=list(a= 1, b=1,c=1),control=nls.control( maxiter = 200,warnOnly = TRUE)))
min_AIC <- min(AICc(growth_models[[1]], growth_models[[2]],growth_models[[3]])$AICc)
bestgrowth <- growth_models[[which(min_AIC == AICc(growth_models[[1]], growth_models[[2]], growth_models[[3]])$AICc)]]

# getting residuals of growth
size$growth_residuals <- NA
size$growth_residuals[which(!is.na(size$t1) & !is.na(size$t0))] <- summary(bestgrowth)$residuals^2


reproyesno_models <- list(glm(reproyesno~ t0, family= "binomial",data = size),
                          glm(reproyesno~ t0 + t0sq , family= "binomial",data =size))
min_AIC <- min(AICc(reproyesno_models[[1]], reproyesno_models[[2]])$AICc)
bestreproyesno <- reproyesno_models[[which(min_AIC == AICc(reproyesno_models[[1]], reproyesno_models[[2]])$AICc)]] 

# reproduction if reproducing: 
rep_models <- list(lm(bulbs0~ t0, data = sizeforrepro),
                   lm(bulbs0~ t0+ t0sq, data = sizeforrepro)) 
min_AIC <- min(AICc(rep_models[[1]], rep_models[[2]])$AICc)
bestrep <- rep_models[[which(min_AIC == AICc(rep_models[[1]], rep_models[[2]])$AICc)]] 

# variance in growth: 
var_models <- list(glm(growth_residuals~ t0-1, data= size), glm(growth_residuals~ t0 + t0sq-1, data= size))
min_AIC <- min(AICc(var_models[[1]], var_models[[2]])$AICc) 
bestvar <- var_models[[which(min_AIC == AICc(var_models[[1]], var_models[[2]])$AICc)]] 

sz0595now=quantile(alldata$szs0,c(1/bin.num,(1-1/bin.num))) # make these quantiles match the bin numbers
sbrksnow=  seq(from=sz0595now[1],to=sz0595now[2],length.out=(bin.num-1)) 

vec.bin = c(szmin,sbrksnow,szmax) 

  binmids =   rep(NA, length(vec.bin)-1)
  for(jj in 1:(length(vec.bin)-1)){
    bounds <- c(vec.bin[jj], vec.bin[jj+1])
    subsetszs <- pdfsz2[which(pdfsz2[,1] >= bounds[1] & pdfsz2[,1] < bounds[2]),]
    if (length(which(pdfsz2[,1] >= bounds[1] & pdfsz2[,1] < bounds[2]))>1) { binmids[jj] <- weightedMedian(subsetszs[,1],subsetszs[,2])} else {binmids[jj]=subsetszs[1]}
    
  }
  
  meanbinmids <- 0.5*(vec.bin[2:length(vec.bin)] + vec.bin[1:(length(vec.bin)-1)])
  binmids[is.na(binmids)]=meanbinmids[is.na(binmids)]
  n.bin = length(binmids)
  truebinsizes= n.bin  
  
  # constructing matrices
  indata <- as.data.frame(cbind(binmids, binmids^2))
  names(indata) <- c("t0", "t0sq")
  sur_vals <- predict(bestsur,indata, type='response')
  reproyesnovals=predict(bestreproyesno, indata, type='response')
  rep_vals <- predict(bestrep, indata, type='response')
  rep_vals[rep_vals<0] <- 0
  growth_vals <- predict(bestgrowth,indata, type='response')
  growth_vals[which(growth_vals<=0)]= szmin
  var_vals <- predict(bestvar,indata, type='response')
  var_vals[which(var_vals<0)] = 0.001
  gmx <- matrix(data=NA, nrow=n.bin, ncol=n.bin)
  
  sdlgszcdf=pnorm(vec.bin,mnsldgsz,sdsldgsz)
  sdlgszprobs=sdlgszcdf[2:length(vec.bin)]-sdlgszcdf[1:(length(vec.bin)-1)]
  sdlgszprobs=sdlgszprobs/sum(sdlgszprobs)
  
  # growth probs 
  for (ss in 1:(n.bin)) {
    growcdf <- pnorm(vec.bin,growth_vals[ss],sqrt(var_vals[ss]))
    grows <- growcdf[2:length(vec.bin)]-growcdf[1:(length(vec.bin)-1)]
    if(sum(grows)>0){grows <- grows/sum(grows)
    gmx[,ss] <- grows} else {gmx[,ss] <- NA} 
  } # end ss loop
  
  survgmx <- gmx*t(matrix( rep(sur_vals,(n.bin)),(n.bin))) # survs* growth
  reprow <- rep_vals*reproyesnovals 
  mx1 <- survgmx # growth and survival
  mx <- matrix(0, (n.bin+1), (n.bin+1))
  mx[2:(n.bin+1), 2:(n.bin+1)] = mx1
  mx[2:(n.bin+1),1] = surv_bulbings*sdlgszprobs
  mx[1,2:(n.bin+1)] = reprow*bulblings_per_bulbil  
  
  # this says that if the matrix model can't be fit, do not store the data from the ipm fit either
  if (length(which(is.na(M)==TRUE))==0) {
    lambdas_ipm <- Re(eigen(mx)$values[1])
    dampratio_median= damping.ratio(mx)
    lifespan_median =lifespan(mx) } else {
      lambdas_ipm <- NA # calls popbio f'n 'lambda'
      lifespan_median <- NA # calls popbio f'n 'lambda'
      dampratio_median=NA } 

  alloutdata = rbind(alloutdata,cbind(sams,rr,bin.num,1,  lambdas_ipm,lifespan_median,dampratio_median,dim(mx)[1],min(nums), median(nums), length(nums[nums<6]),(length(which(is.na(M)==TRUE))==0) ))


  }} # end loops for sample size and rep# 
  
} # end loop for different class nums


colnames(alloutdata)=c('num','rep','classnum','model','lam','lifespan','damping', 'mxdim', 'minnum', 'mediannum', 'n.lessthan6','n.zeros')

#################################################
#Plotting########################################


dev.new(width=5, height = 8,noRStudioGD = TRUE)
par(mfrow=c(3,1))

alloutdata2020=as.data.frame(alloutdata)
samp=round(alloutdata2020[,1]*totsample)
lams = alloutdata2020$lam
lifespan=alloutdata2020$lifespan
dampingratios=alloutdata2020$damping
mxipm=alloutdata2020$model
classes = alloutdata2020$classnum

  beanplot(lams ~ mxipm*classes, ll = 0.04, main = c('Bistort'), ylim=c(0.95,1.05), xlab='class number', ylab = "lambda", side = "both",border = NA, col = list("black", c("grey", "white")),what = c(FALSE, TRUE, TRUE, FALSE))
  legend("topright", fill = c("black", "grey"), legend = c("Discrete", "Continuous"))

  beanplot(dampingratios ~ mxipm*classes, ll = 0.04, main = c('Bistort'), xlab='class number', ylab = "damping ratio", side = "both",border = NA, col = list("black", c("grey", "white")),what = c(FALSE, TRUE, TRUE, FALSE))
  legend("topright", fill = c("black", "grey"))
  
 beanplot(lifespan ~ mxipm*classes, ll = 0.04, xlab='class number', ylab = "lifespan", side = "both",border = NA, col = list("black", c("grey", "white")),what = c(FALSE, TRUE, TRUE, FALSE))
  
 
