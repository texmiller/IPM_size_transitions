###### This script illustrates how the mesh point method actually works when you think
###### a bit about how many bins are actually needed, instead of making uninformed guesses.  

graphics.off()
rm(list=ls())

# what bin numbers to try? 
bins =seq(100,500,by=50)
minsize=0
maxsize=101

means=seq(49,51,0.1)
sd = 0.5 # SD of growth

dev.new(width=6, height = 4)
# run through different mean sizes and bin numbers to generate estimates of summed probabilities of growth to any bin
alldats=NULL



system.time({
for (mn in means) {

  dats=NULL
  for (binnum in bins){
  vec.bin = c(minsize, minsize+1:binnum*(maxsize-minsize)*(1/binnum)) 
  binmids <- 0.5*(vec.bin[2:length(vec.bin)] + vec.bin[1:(length(vec.bin)-1)])
  n.bin = length(binmids)
  binwidth=(maxsize-minsize)/n.bin
  grows2=binwidth*dnorm(binmids,mn,sd) #mesh pt approach
   
  growcdf <- pnorm(vec.bin,mn,sd)
  grows1 <- growcdf[2:length(vec.bin)]-growcdf[1:(length(vec.bin)-1)] #CDF difference approach
  
  totgrow=pnorm(maxsize,mn,sd)-pnorm(minsize,mn,sd)

  dats=rbind(dats,c(binnum,sum(grows2)))
  alldats=rbind(alldats,cbind(mn,binnum,sum(grows2)))
} # end bin num loop

  if (mn==min(means)){
    plot(dats[,1],dats[,2], type='b', xlab='Class number', ylab='Summed growth probabilities',ylim=c(0.98,1.02))
    abline(1,0, col='red')
  }else {points(dats[,1],dats[,2], type='b')}
 
} # end mn loop 

})
