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
  fit=gam(y~s(x),gamma=gamma,method="REML")
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
  spline.scatter.smooth(rollx,rollmean,gamma=2,xlab="Fitted values",ylab="Mean");
  if(scaled) abline(h=0,col="red",lty=2,lwd=2) 

  spline.scatter.smooth(rollx,rollsd,gamma=2,xlab="Fitted values",ylab="Std Dev"); 
  if(scaled) abline(h=1,col="red",lty=2,lwd=2) 

  spline.scatter.smooth(rollx,rollskew,gamma=2,xlab="Fitted values",ylab="Skew"); 
  if(scaled) abline(h=0,col="red",lty=2,lwd=2) 

  spline.scatter.smooth(rollx,rollkurt,gamma=2,xlab="Fitted values",ylab="Kurtosis"); 
  if(scaled) abline(h=3,col="red",lty=2,lwd=2)
}
return(list(rollx=rollx,rollmean=rollmean,rollsd=rollsd,rollkurt=rollkurt,rollskew=rollskew))
}