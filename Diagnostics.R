require(mgcv); require(moments); require(zoo); 


add_panel_label <- function(ltype="a",cex=1) {
    text <- paste(LETTERS[letters==ltype], ")", sep="")
    mtext(text=text, side=3, adj=0,cex=cex)
}

###########################################################################
### Spline scatterplot smoothing function
### Adjustable dof cost gamma, for use with rollapply outputs 
### show.quadratic adds a quadratic fit to the plot, for comparison 
###########################################################################
spline.scatter.smooth=function(x,y,gamma=2,show.quadratic=FALSE,...) {
  fit=gam(y~s(x,k=7),gamma=gamma,method="REML")
  plot(x,y,type="p",pch=16,cex=1.3,...);
  out=predict(fit,type="response"); 
  points(x,out,type="l",lwd=1)
  if(show.quadratic){
    fit2 = lm(y~x+I(x^2));
    points(x,fit2$fitted,type="l",col="blue",lwd=2,lty=2);
  }
}

#####################################################################################
## Rollapply moment diagnostics on mean, SD, skew, excess kurtosis 
## Inputs: 
##	  px is the covariate (vector), py is the response (vector of equal length)
##    windows (number) determines window size, by stating how many 
##       non-overlapping windows would include the entire data set. 
##		 The default windows=10 means each window includes 10% of the data
##    smooth=TRUE results in calling spline.scatter.smooth to make plots
##    scaled=TRUE adds benchmarks mean=0,sd=1,skew=excess kurtosis=0 to the plots
#####################################################################################
rollMoments=function(px,py,windows=10,smooth=TRUE,scaled=TRUE) {

  e = order(px); px=px[e]; py=py[e];  

  width=round(length(px)/windows); by=round(width/2); 
  rollx=rollapply(px,width=width,mean,by=by);
  rollmean = rollapply(py,width=width,mean,by=by); 
  rollsd=rollapply(py,width=width,sd,by=by); 
  rollkurt=rollapply(py,width=width,kurtosis,by=by);
  rollskew=rollapply(py,width=width,skewness,by=by);

  if(smooth) {
  par(mfrow=c(2,2),mar=c(4,4,2,1),cex.axis=1.3,cex.lab=1.4,bty="l"); 
  spline.scatter.smooth(rollx,rollmean,gamma=2,xlab="Fitted values",ylab="Mean", ylim=c(-1,1));
  if(scaled) abline(h=0,col="red",lty=2,lwd=2) 

  spline.scatter.smooth(rollx,rollsd,gamma=2,xlab="Fitted values",ylab="Std Dev",ylim=c(0,2)); 
  if(scaled) abline(h=1,col="red",lty=2,lwd=2) 

  spline.scatter.smooth(rollx,rollskew,gamma=2,xlab="Fitted values",ylab="Skew"); 
  if(scaled) abline(h=0,col="red",lty=2,lwd=2) 

  spline.scatter.smooth(rollx,rollkurt,gamma=2,xlab="Fitted values",ylab="Kurtosis"); 
  if(scaled) abline(h=3,col="red",lty=2,lwd=2)
}
return(list(rollx=rollx,rollmean=rollmean,rollsd=rollsd,rollkurt=rollkurt,rollskew=rollskew))
}

###############################################################################
# Nonparametric measures of scale, skew and excess kurtosis
###############################################################################
NPsd = function(x) {
	u = diff(quantile(x,c(0.25,0.75))/(qnorm(0.75)-qnorm(0.25)))
	as.numeric(u)
}	
	
NPskewness.Pearson = function(x) 3*(mean(x)-median(x))/sd(x) 

NPskewness = function(x,p=0.1) {
	q = quantile(x,c(p,0.5,1-p))
	u = (q[3]+q[1]-2*q[2])/(q[3]-q[1]);
	return(as.numeric(u)); 
	
}	

NPkurtosis=function(x,p=0.05) {
	q = quantile(x,c(p,0.25,0.75,1-p))
	qN = qnorm(c(p,0.25,0.75,1-p))
	u = (q[4]-q[1])/(q[3]-q[2]);
	uN = (qN[4]-qN[1])/(qN[3]-qN[2]);
	return (as.numeric(u/uN-1)) 
}
	
#####################################################################################
## Rollapply moment diagnostics on mean, SD, Nonparametric skew, excess Nonparametric kurtosis 
## Inputs: 
##	  px is the covariate (vector), py is the response (vector of equal length)
##    windows (number) determines window size, by stating how many 
##       non-overlapping windows would include the entire data set. 
##		 The default windows=10 means each window includes 10% of the data
##    smooth=TRUE results in calling spline.scatter.smooth to make plots
##    scaled=TRUE adds benchmarks mean=0,sd=1,skew=excess kurtosis=3 to the plots
#####################################################################################
rollMomentsNP=function(px,py,windows=10,smooth=TRUE,scaled=TRUE,xlab=NULL) {

  if(is.null(xlab)) xlab="Fitted values"; 
  e = order(px); px=px[e]; py=py[e];  

  width=round(length(px)/windows); by=round(width/2); 
  rollx=rollapply(px,width=width,mean,by=by);
  rollmean = rollapply(py,width=width,mean,by=by); 
  rollsd=rollapply(py,width=width,sd,by=by); 
  rollkurt=rollapply(py,width=width,NPkurtosis,by=by);
  rollskew=rollapply(py,width=width,NPskewness,by=by);

  if(smooth) {
  par(mfrow=c(2,2),mar=c(4,5,2,1),cex.axis=1.3,cex.lab=1.3,mgp=c(2.2,1,0),bty="l"); 
  spline.scatter.smooth(rollx,rollmean,gamma=2.8,xlab=xlab,ylab="Mean", ylim=c(-1,1));
  if(scaled) abline(h=0,col="red",lty=2,lwd=2) 
  add_panel_label("a"); 

  spline.scatter.smooth(rollx,rollsd,gamma=2.8,xlab=xlab,ylab="Std Dev",ylim=c(0,2)); 
  if(scaled) abline(h=1,col="red",lty=2,lwd=2) 
    add_panel_label("b"); 

  spline.scatter.smooth(rollx,rollskew,gamma=2.8,xlab=xlab,ylab="NP Skew", ylim=c(-1,1)); 
  if(scaled) abline(h=0,col="red",lty=2,lwd=2) 
  add_panel_label("c"); 

  spline.scatter.smooth(rollx,rollkurt,gamma=2.8,xlab=xlab,ylab="NP Kurtosis", ylim=c(-1,1)); 
  if(scaled) abline(h=0,col="red",lty=2,lwd=2)
  add_panel_label("d"); 
}
return(list(rollx=rollx,rollmean=rollmean,rollsd=rollsd,rollkurt=rollkurt,rollskew=rollskew))
}

#########################################################################################
## Comparison of conditional quantiles between simulated at actual data 
##
## The real data are a set of paired scalar observations (sortVariable_i, trueData_i) 
## Typically trueData_i is the i^th observation of future size, and sortVariable_i is
## the i^th observation of current size or fitted value. 
##
## simData is a matrix with nrow(simData)=length(trueData). Each column of simData is
## a simulated "data set" of the trueData variable (e.g., a simulated set of future sizes).
##
## nBins is the number of equal-size bins to use for the plotting  
## alpha_scale affects the colors used in plotting
##
## The tidyverse is required so that pipes and various other bits work. 
#########################################################################################
quantileComparePlot = function(sortVariable,trueData,simData,nBins,alpha_scale = 0.7) {
  xTrue = data.frame(x=sortVariable,y=trueData)
  qTrue <- xTrue %>% arrange(x) %>% 
    mutate(size_bin = cut_number(x,n=nBins)) %>% 
    group_by(size_bin) %>% 
    summarise(q1 = quantile(y,0.05),
              q2 = quantile(y,0.25),
              q3 = quantile(y,0.5),
              q4 = quantile(y,0.75),
              q5 = quantile(y,0.95),
              bin_mean = mean(x),
              bin_n = n()) 
  
  qSim = array(NA,dim=c(nBins,ncol(simData),5))
  qnames=c("q1","q2","q3","q4","q5"); 
  for(i in 1:ncol(simData)){
    xSim_i = data.frame(x=sortVariable,y=simData[,i])
    qSim_i <- xSim_i %>% arrange(x) %>% 
      mutate(size_bin = cut_number(x,n=nBins)) %>% 
      group_by(size_bin) %>% 
      summarise(q1 = quantile(y,0.05),
                q2 = quantile(y,0.25),
                q3 = quantile(y,0.5),
                q4 = quantile(y,0.75),
                q5 = quantile(y,0.95),
                bin_mean = mean(x),
                bin_n = n()) 
    for(j in 1:5) qSim[1:nBins,i,j]=unlist(qSim_i[1:nBins,qnames[j]])	
  }
  
  ylabs=c("5th","25th","50th","75th","95th");
  trueBinQ = simBinQ = matrix(NA,nBins,5); 
  par(mfrow=c(3,2),mar=c(4,4,2,1),cex.axis=1.3,cex.lab=1.3,mgp=c(2,1,0),bty="l");
  for(j in 1:5) {
    qSim_j = qSim[,,j]
    matplot(qTrue$bin_mean, qSim_j,col=alpha("gray",0.5),pch=16,xlab="Bin mean",
            ylab=paste(ylabs[j],"Percentile"),cex=1.4); 
    points(qTrue$bin_mean, unlist(qTrue[,qnames[j]]),pch=5,lwd=2,col=alpha("red",alpha_scale),cex=1.4)
    points(qTrue$bin_mean, apply(qSim_j,1,median),pch=1,lwd=2,col=alpha("black",alpha_scale),cex=1.4)
    trueBinQ[,j] = unlist(qTrue[,qnames[j]]);
    simBinQ[,j] = apply(qSim_j,1,median);
    if(j==1) legend("topleft",legend=c("Model simulations","Median of simulations","Data"),
                    col=c(alpha("gray",0.5),alpha("black",alpha_scale), alpha("red",alpha_scale)),
                    pch=c(1,1,5),lwd=2,cex=1.1,bty="n"); 
    add_panel_label(letters[j])
  }	
  
  Qcolors=c(rep(alpha("red",alpha_scale),5),rep(alpha("black",alpha_scale),5)); 
  matplot(qTrue$bin_mean,cbind(trueBinQ,simBinQ), col=Qcolors, pch=c(rep(5,5),rep(1,5)),lty=2, lwd=2, cex=1.4,type="o",
          xlab="Bin mean", ylab="All percentiles"); 
  add_panel_label("f"); 
}	

#########################################################################################
## Comparison of conditional moments between simulated at actual data 
##
## The real data are a set of paired scalar observations (sortVariable_i, trueData_i) 
## Typically trueData_i is the i^th observation of future size, and sortVariable_i is
## the i^th observation of current size or fitted value. 
##
## simData is a matrix with nrow(simData)=length(trueData). Each column of simData is
## a simulated "data set" of the trueData variable (e.g., a simulated set of future sizes).
##
## nBins is the number of equal-size bins to use for the plotting  
## alpha_scale affects the colors used in plotting
##
## The tidyverse is required so that pipes and various other bits work. 
#########################################################################################
momentsComparePlot = function(sortVariable,trueData,simData,normData=NULL,fittedMean=NULL,nBins,alpha_scale = 0.7) {

  # optional: skewness and kurtosis done on deviations from a fitted mean 
  if(is.null(fittedMean)) fittedMean=rep(0,length(sortVariable)); 
  
  xTrue = data.frame(x=sortVariable,y=trueData, z= fittedMean)
  qTrue <- xTrue %>% arrange(x) %>% 
    mutate(size_bin = cut_number(x,n=nBins)) %>% 
    group_by(size_bin) %>% 
    summarise(q1 = mean(y),
              q2 = sd(y),
              q3 = NPskewness(y-z),
              q4 = NPkurtosis(y-z),
              bin_mean = mean(x),
              bin_n = n()) 
              
       
  qSim = array(NA,dim=c(nBins,ncol(simData),4))
  qnames=c("q1","q2","q3","q4"); 
  for(i in 1:ncol(simData)){
    xSim_i = data.frame(x=sortVariable,y=simData[,i],z=fittedMean)
    qSim_i <- xSim_i %>% arrange(x) %>% 
      mutate(size_bin = cut_number(x,n=nBins)) %>% 
      group_by(size_bin) %>% 
      summarise(q1 = mean(y),
                q2 = sd(y),
                q3 = NPskewness(y-z),
                q4 = NPkurtosis(y-z),
                bin_mean = mean(x),
                bin_n = n()) 
    for(j in 1:4) qSim[1:nBins,i,j]=unlist(qSim_i[1:nBins,qnames[j]])	
  }
  
  # optional: compare with Gaussian model 
  if(!is.null(normData)) {  
  
  qNorm = array(NA,dim=c(nBins,ncol(normData),4))
  qnames=c("q1","q2","q3","q4"); 
  for(i in 1:ncol(simData)){
    xNorm_i = data.frame(x=sortVariable,y=normData[,i],z=fittedMean)
    qNorm_i <- xNorm_i %>% arrange(x) %>% 
      mutate(size_bin = cut_number(x,n=nBins)) %>% 
      group_by(size_bin) %>% 
      summarise(q1 = mean(y),
                q2 = sd(y),
                q3 = NPskewness(y-z),
                q4 = NPkurtosis(y-z),
                bin_mean = mean(x),
                bin_n = n()) 
    for(j in 1:4) qNorm[1:nBins,i,j]=unlist(qNorm_i[1:nBins,qnames[j]])	
  }
  
  } ## end of optional Gaussian 
  
  ylabs=c("Mean","Std. Deviation","NP skew","NP excess kurtosis");
  trueBinQ = simBinQ = matrix(NA,nBins,5); 
  par(mfrow=c(2,2),mar=c(4,4,2,1),cex.axis=1.3,cex.lab=1.3,mgp=c(2,1,0),bty="l");
  
  h = min(diff(qTrue$bin_mean)); 
  xlim = c(min(qTrue$bin_mean)-0.25*h, max(qTrue$bin_mean)+0.25*h); 
  
  for(j in 1:4) {
    qSim_j = qSim[,,j]; qNorm_j = qNorm[,,j]; 
    matplot(qTrue$bin_mean-0.2*h, qSim_j,col=alpha("gray",0.5),pch=16,xlab="Bin mean",
            ylab=ylabs[j],cex=1.4,xlim=xlim); 
    points(qTrue$bin_mean-0.2*h, apply(qSim_j,1,median),pch=1,lwd=2,col=alpha("black",alpha_scale),cex=1.4)        

    if(!is.null(normData)) {
    
    matpoints(qTrue$bin_mean+0.2*h, qNorm_j,col=alpha("cornflowerblue",0.5),pch=16,xlab="Bin mean",
            ylab=ylabs[j],cex=1.4,xlim=xlim); 
    points(qTrue$bin_mean+0.2*h, apply(qNorm_j,1,median),pch=1,lwd=2,col=alpha("blue",alpha_scale),cex=1.4)  
    
    }
    
    points(qTrue$bin_mean, unlist(qTrue[,qnames[j]]),pch=5,lwd=2,col=alpha("red",alpha_scale),cex=1.4)

    trueBinQ[,j] = unlist(qTrue[,qnames[j]]);
    simBinQ[,j] = apply(qSim_j,1,median);
    if(j==1 & is.null(normData)) {legend("topleft",legend=c("Model simulations","Median of simulations","Data"),
                    col=c(alpha("gray",0.5),alpha("black",alpha_scale), alpha("red",alpha_scale)),
                    pch=c(1,1,5),lwd=2,cex=1.1,bty="n"); }
     if(j==1 & !is.null(normData)) {legend("topleft",legend=c("Gaussian model","Non-Gaussian model","Data"),
                    col=c(alpha("blue",0.5),alpha("grey",alpha_scale), alpha("red",alpha_scale)),
                    pch=c(1,1,5),lwd=2,cex=1.1,bty="n"); }           
                    
    add_panel_label(letters[j])
  }	
  
}	

#########################################################################################
## Fitting parameters of a chosen distribution for a set of data bins  
##     y is the set of values to be fitted (unbinned)
##     sortVar is the variable used for sorting and binning 
##     DIST is the name of the distribution family (e.g., DIST = "JSU")
##
##     if rolling==TRUE, fitting is done on sliding pairs of bins (not yet implemented)  
##
##     For the moment this assumes a 4-parameter family. Will fix? 
##  
#########################################################################################
binnedPars <- function(y,sortVar,nBins,DIST,rolling=FALSE) {
    X = data.frame(sortVar=sortVar,response=y)
    X <- X %>% mutate(size_bin = cut_number(sortVar,n=nBins))
    bins = levels(X$size_bin); 
    mus = sigmas = nus = taus = bin_means = numeric(length(bins)); 
    for(j in 1:length(bins)){
        Xj=subset(X,size_bin==bins[j])
        fitj = gamlssMaxlik(Xj$response,DIST)
        pars=fitj$out[[1]]$estimate; 
        mus[j]=pars[1]; sigmas[j]=pars[2]; nus[j]=pars[3]; taus[j]=pars[4]; 
        bin_means[j]=mean(Xj$sortVar); 
    }	

    return(list(mus=mus,sigmas=sigmas,nus=nus,taus=taus,bin_means=bin_means)); 
}



