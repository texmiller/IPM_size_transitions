###### This script fits discretely estimated matrix models and IPMs based on continuous
###### vital rate functions for a range of bin numbers and discretization methods
###### using the data for Borderea chouardii. 
###### This version of the script fits with even class widths, the recommended approach 
###### fitting models. In the main text, results of the script shown here are included
###### in Figures 10, 11, and 12. 

#Borderea notes:for Borderea, the first class is seeds, then the first column is seeds becoming small plants (seedlings, and can be done using the size distr of seedling plants from the data

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

### Load, manipulate, and sort data --------------
dat <- read.csv("allBCdata.csv")
dat$stat0[dat$stat0=="Seedl"] = 'seedl'
dat=dat[which(is.na(dat$sz0)==FALSE),]

seedrates = c(0.500,0.030) # the entries in first column of mx, from Garcia: these are seeds-> seeds, and surv of seeds to be new plants the next yr
seedsperfrt=4.5 # rough average from Garcia data

dat.size <- as.data.frame(dat[, c(6,10,11)])
colnames(dat.size) <- c("sz0", "sz1", "repro") 
size <- arrange(dat.size, sz0) # arrange initial column in ascending orders

size$survival <- 1
size$survival[which(is.na(size$sz1)==TRUE)] <- 0

size$reproyesno= size$repro # a variable to have yes/no reproduction
size$reproyesno[size$reproyesno>1]=1
sizeforrepro=size[which(size$repro>0),]

# get mean and variance of sdling sizes:
mnsldgsz=mean(dat$sz0[which(dat$stat0=="seedl")], na.rm=TRUE)
sdsldgsz=  sd(dat$sz0[which(dat$stat0=="seedl")])
allsdszs=dat$sz0[which(dat$stat0=="seedl")]

minsize=2 # set manually to get rid of an outlier; important for even size classes
maxsize <- max(c(size$sz0,size$sz1), na.rm=TRUE)+0.1
size[which(size$sz0 < minsize),1] = minsize
size[which(size$sz1 < minsize),2] = minsize

size$t0=size$sz0 # converting names for consistency to other species code
size$t1=size$sz1
if(length(size$t0[which(is.na(size$t0))])>0){ 
  warning("remove NA's in sizes")
}

size$reproyesno= size$repro # a variable to have yes/no reproduction
size$reproyesno[size$reproyesno>1]=1
sizeforrepro=size[which(size$repro>0),]

#generate empirical density function for median size estimation
pdfsz=density(size$t0, n=1024, cut=0)
pdfsz2=cbind(pdfsz$x,pdfsz$y)

####################################################
lifespan <- function(nx){
  nclasses=dim(nx)[1]
  vec=c(100,rep(0,(nclasses-1)))
  # nx[1,]=0
  # nx[,1]=nx[,1]/sum(nx[,1]) # this is only for borderea, to get to seedlings
  nx[1,2:length(nx[1,])]=0
  jj=1
  while (sum(vec)>1){
    vec=nx%*%vec
    jj=jj+1
    #print(sum(vec))
  }
  return(jj)
}#############################################



### BEGIN THE ESTIMATION FOR EACH MODEL TYPE: MATRIX = DVR IN TEXT, WHILE IPM MEAN AND IPM MEDIAN = CVR-MIDPOINT AND CVR-MEDIAN IN TEXT. 

### Matrix Models ------------------------


lambdas_matrix <- rep(NA, length(bin.num))
dampratio_matrix = rep(NA, length(bin.num))
lifespan_matrix=rep(NA, length(bin.num))
classNdataout = NULL # this stores the new results on sample size per class

mincounts=NULL

for(i in 1:length(bin.num)){

  ss=as.numeric(size$t0)
  vec.bin = c(minsize, minsize+1:bin.num[i]*(maxsize-minsize)*(1/bin.num[i]))  
  nums=hist(ss,breaks=vec.bin, plot=FALSE)$counts 
  mincounts=c(mincounts,min(nums))
  
  if (min(nums)>2) { 
    
  ## Initialize storage
  n.bin <- length(vec.bin)-1                 
  n <- rep(NA, n.bin)                         # count of indvs per bin
  medians <- rep(NA, n.bin)                    # median size per bin for F
  surv <- rep(NA, n.bin)                      # survivorship for each class
  grow <- matrix(NA, n.bin, n.bin)            # store growth probabilities for each class
  reproduction <- rep(NA, n.bin)
  
  totnums = 0 # this is monitoring for errors
  # bin survival, growth
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
    histo <- hist(subset$t1, breaks = vec.bin, plot = FALSE)
    grow[,j] <- histo$counts/length(subset$t0[subset$survival==1]) 
    reproduction[j] <- mean(subset$repro) # note that repro = frt # for Borderea
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
    
    # add lines for the creation of seeds and their transition to first size class
    M[2:(n.bin+1), 2:(n.bin+1)] = M1
    M[,1] = c(seedrates[1], seedrates[2]* sdlggrow)
    M[1,2:(n.bin+1)] = surv*reproduction*seedsperfrt #  
  
  
  lambdas_matrix[i] <- lambda(M) 
  lifespan_matrix[i] <- lifespan(M) 
  dampratio_matrix[i]=damping.ratio(M)
  classNdataout= rbind(classNdataout,c(bin.num[i],min(nums), median(nums), length(nums[nums<6])))
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
size$t0sq <- size$t0^2           
sizeforrepro$t0sq = sizeforrepro$t0^2
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

# getting residuals of growth
size$growth_residuals <- NA
size$growth_residuals[which(!is.na(size$t1) & !is.na(size$t0))] <- summary(bestgrowth)$residuals^2


reproyesno_models <- list(glm(reproyesno~ t0, family= "binomial",data = size),
                   glm(reproyesno~ t0 + t0sq , family= "binomial",data =size))
min_AIC <- min(AICc(reproyesno_models[[1]], reproyesno_models[[2]])$AICc)
bestreproyesno <- reproyesno_models[[which(min_AIC == AICc(reproyesno_models[[1]], reproyesno_models[[2]])$AICc)]] 

# reproduction if reproducing:  
rep_models <- list(lm(repro~ t0, data = sizeforrepro),
                    lm(repro~ t0+ t0sq, data = sizeforrepro)) 
min_AIC <- min(AICc(rep_models[[1]], rep_models[[2]])$AICc)
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

  vec.bin = c(minsize, minsize+1:bin.num[i]*(maxsize-minsize)*(1/bin.num[i])) 
 binmids <- 0.5*(vec.bin[2:length(vec.bin)] + vec.bin[1:(length(vec.bin)-1)])
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
  growth_vals <- predict(bestgrowth,indata, type='response')
  var_vals <- predict(bestvar,indata, type='response')
  var_vals[var_vals<0]=0.0001
  gmx <- matrix(data=NA, nrow=n.bin, ncol=n.bin)
  sdlgszcdf=pnorm(vec.bin,mnsldgsz,sdsldgsz)
  sdlgszprobs=sdlgszcdf[2:length(vec.bin)]-sdlgszcdf[1:(length(vec.bin)-1)]
  sdlgszprobs=sdlgszprobs/sum(sdlgszprobs)
  
  # growth probs using cdf fn
  for (ss in 1:(n.bin)) {
    growcdf <- pnorm(vec.bin,growth_vals[ss],sqrt(var_vals[ss]))
    grows <- growcdf[2:length(vec.bin)]-growcdf[1:(length(vec.bin)-1)]
    if(sum(grows)>0){grows <- grows/sum(grows)
                     gmx[,ss] <- grows} else {gmx[,ss] <- NA} 
  } # end ss loop
  
  # make the surv*growth mx
  survgmx <- gmx*t(matrix( rep(sur_vals,(n.bin)),(n.bin))) # survs* growth
  reprow <- rep_vals*reproyesnovals
 
  mx1 <- survgmx # growth and survival, without repro
  
  mx <- matrix(0, (n.bin+1), (n.bin+1))
   mx[2:(n.bin+1), 2:(n.bin+1)] = mx1
   mx[,1] = c(seedrates[1], seedrates[2]* sdlgszprobs)
   
 mx[1,2:(n.bin+1)] = reprow*seedsperfrt 
 
  lambdas_ipm_mean[i] <- Re(eigen(mx)$values[1])
  dampratio_mean[i]= damping.ratio(mx)
  lifespan_mean[i] =lifespan(mx)
  print(i)
}


#######################################
######################################
### IPM with median sizes-----------------------------------------------------

# adding in columns for squared size effects
size$t0sq <- size$t0^2           
sizeforrepro$t0sq = sizeforrepro$t0^2
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

# getting residuals of growth
size$growth_residuals <- NA
size$growth_residuals[which(!is.na(size$t1) & !is.na(size$t0))] <- summary(bestgrowth)$residuals^2

reproyesno_models <- list(glm(reproyesno~ t0, family= "binomial",data = size),
                          glm(reproyesno~ t0 + t0sq , family= "binomial",data =size))
min_AIC <- min(AICc(reproyesno_models[[1]], reproyesno_models[[2]])$AICc)
bestreproyesno <- reproyesno_models[[which(min_AIC == AICc(reproyesno_models[[1]], reproyesno_models[[2]])$AICc)]] 

# reproduction if reproducing: 
rep_models <- list(lm(repro~ t0, data = sizeforrepro),
                   lm(repro~ t0+ t0sq, data = sizeforrepro)) 
min_AIC <- min(AICc(rep_models[[1]], rep_models[[2]])$AICc)
bestrep <- rep_models[[which(min_AIC == AICc(rep_models[[1]], rep_models[[2]])$AICc)]]
# variance in growth: 
var_models <- list(glm(growth_residuals~ t0-1, data= size), glm(growth_residuals~ t0 + t0sq-1, data= size))
min_AIC <- min(AICc(var_models[[1]], var_models[[2]])$AICc) 
bestvar <- var_models[[which(min_AIC == AICc(var_models[[1]], var_models[[2]])$AICc)]] 

# vary bin sizes
lambdas_ipm_median <- vector("numeric", length= length(bin.num))
dampratio_median <- vector("numeric", length= length(bin.num))
lifespan_median = vector("numeric", length= length(bin.num))
truebinsizes= matrix(0,length(bin.num),1)

for (i in 1:length(bin.num)){

  vec.bin = c(minsize, minsize+1:bin.num[i]*(maxsize-minsize)*(1/bin.num[i])) 
  nums=hist(size$t0,breaks=vec.bin, plot=FALSE)$counts 
  binmids =   rep(NA, length(vec.bin)-1)
  
  for(jj in 1:(length(vec.bin)-1)){
    bounds <- c(vec.bin[jj], vec.bin[jj+1])
    subsetszs <- pdfsz2[which(pdfsz2[,1] >= bounds[1] & pdfsz2[,1] < bounds[2]),]
   if (length(which(pdfsz2[,1] >= bounds[1] & pdfsz2[,1] < bounds[2]))>1) { binmids[jj] <- weightedMedian(subsetszs[,1],subsetszs[,2])} else {binmids[jj]=subsetszs[1]}
  }
  
  meanbinmids <- 0.5*(vec.bin[2:length(vec.bin)] + vec.bin[1:(length(vec.bin)-1)])
  binmids[is.na(binmids)]=meanbinmids[is.na(binmids)]
  n.bin = length(binmids)
  truebinsizes[i] = n.bin  
  
  # constructing matrices
  
  indata <- as.data.frame(cbind(binmids, binmids^2))
  names(indata) <- c("t0", "t0sq")
  
  sur_vals <- predict(bestsur,indata, type='response')
  reproyesnovals=predict(bestreproyesno, indata, type='response')
  rep_vals <- predict(bestrep, indata, type='response')
  rep_vals[rep_vals<0] <- 0
  growth_vals <- predict(bestgrowth,indata, type='response')
  var_vals <- predict(bestvar,indata, type='response')
  var_vals[var_vals<0]=0.0001
  gmx <- matrix(data=NA, nrow=n.bin, ncol=n.bin)
  
  sdlgszcdf=pnorm(vec.bin,mnsldgsz,sdsldgsz)
  sdlgszprobs=sdlgszcdf[2:length(vec.bin)]-sdlgszcdf[1:(length(vec.bin)-1)]
  sdlgszprobs=sdlgszprobs/sum(sdlgszprobs)
  
  # growth probs using cdf fn
  for (ss in 1:(n.bin)) {
    growcdf <- pnorm(vec.bin,growth_vals[ss],sqrt(var_vals[ss]))
    grows <- growcdf[2:length(vec.bin)]-growcdf[1:(length(vec.bin)-1)]
    if(sum(grows)>0){grows <- grows/sum(grows)
    gmx[,ss] <- grows} else {gmx[,ss] <- NA} 
  } # end ss loop
  
  # make the surv*growth mx
  survgmx <- gmx*t(matrix( rep(sur_vals,(n.bin)),(n.bin))) 
  reprow <- rep_vals*reproyesnovals
  
  mx1 <- survgmx # growth and survival, without repro
  mx <- matrix(0, (n.bin+1), (n.bin+1))
  mx[2:(n.bin+1), 2:(n.bin+1)] = mx1
  mx[,1] = c(seedrates[1], seedrates[2]* sdlgszprobs)
  mx[1,2:(n.bin+1)] = reprow*seedsperfrt 
  
  lambdas_ipm_median[i] <- Re(eigen(mx)$values[1])
  dampratio_median[i]= damping.ratio(mx)
  lifespan_median[i] =lifespan(mx)
  print(i)
}# end IPM with medians


borderea.output.even <- as.data.frame(cbind(bin.num,lambdas_matrix, lambdas_ipm_mean,lambdas_ipm_median, dampratio_matrix,dampratio_mean,dampratio_median,lifespan_matrix,lifespan_mean,lifespan_median))
names(borderea.output.even) <- c("bin.num", 'lam.mx','lam,mn','lam.med','damp.mx','damp.mn','damp.med','life.mx','life.mn','life.med')


#Plotting########################################
# plotting lambda
dev.new(width=4, height = 8,noRStudioGD = TRUE)
par(mfrow=c(3,1))#,pty= "s")
plot(truebinsizes,lambdas_matrix, type = 'b', xlab = "# of bins in model", main="Borderea with even class breaks", col='blue', ylim=c(0.9,1.1))
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
plot(truebinsizes,dampratio_matrix, type = 'b', xlab = "# of bins in model", ylab = "damping ratio",col='blue', ylim=c(1,2))
points(truebinsizes,dampratio_median, type = 'b', xlab = "# of bins in model", col ='red', lty=3, pch=8)
linenames = expression(matrix,ipm)
points(truebinsizes,dampratio_mean, type = 'b', xlab = "# of bins in model", col ='red', lty=1, pch=1)





