rm(list=ls(all=TRUE)); 

setwd("c:/repos/IPM_size_transitions/gam practice"); 
require(minqa); require(fda); require(moments); require(zoo); 
source("JPfuns.R"); source("../Diagnostics.R"); source("../matrixImage.R"); 
graphics.off(); 

#######################################################
#  Display how skew and kurtosis depend on parameters
#######################################################

SkewMat = KurtMat = matrix(NA,150,151);
nu = seq(-3,3,length=150)
tau = seq(-1,1,length=151);  
for(i in 1:150) {
for(j in 1:151){
    epsilon=nu[i]; delta=exp(-tau[j]); 
    SkewMat[i,j]=JP_NPskewness(epsilon,delta);
    KurtMat[i,j]=JP_NPkurtosis(epsilon,delta);
}}

graphics.off(); 
dev.new(); 
image.plot(nu,tau,SkewMat,col=plasma(64));  title(main="NP Skewness"); 

dev.new(); 
image.plot(nu, tau,KurtMat,col=plasma(64)); title(main = "NP Kurtosis"); 
    

############################################################
#  Explore what diagnostic plots can reveal about residuals
#  250 data points is enough to find signs of trouble 
############################################################    
    
######### Create covariate for residuals 

z = rt(250,df=10); z=sort(z); hist(z); 

########### Create artificial "residuals" with known sgt parameters 
nu=exp(0.1*z); tau = -2+z;  
resids = rSJP(length(z), epsilon=nu, delta = exp(-tau)); 

jarque.test(resids)$p.value # normality test
agostino.test(resids)$p.value # skew 
anscombe.test(resids)$p.value # kurtosis 


graphics.off(); plot(z,resids); 
dev.new(); qqPlot(resids); 

dev.new(); 
out=rollMomentsNP(z,resids,windows=5);

fit = SJPMaxlik(resids,nstart=10)


###########################################################################
# Fit suitable distributions to binned data (FOR CORALS - must adapt!) 
###########################################################################
logResids <- data.frame(init=XH$logarea.t0,resids=scaledResids); 
logResids <- logResids %>% mutate(size_bin = cut_number(init,n=8))

source("../fitChosenDists.R"); 

tryDists=c("PE","GT","JSU", "SHASHo","SEP1","SEP3","SEP4"); 

bins = levels(logResids$size_bin); 
maxVals = matrix(NA,length(bins),length(tryDists)); 

colnames(maxVals) = tryDists; 
for(j in 1:length(bins)){
for(k in 1:length(tryDists)) {
	Xj=subset(logResids,size_bin==bins[j])
	fitj = gamlssMaxlik(y=Xj$resids,DIST=tryDists[k]); 
	maxVals[j,k] = fitj$aic;
	cat("Finished ", tryDists[k]," ",j,k, fitj$maximum,"\n"); 
}
}

 


if(FALSE) {

########## Create B-spline basis, and plot it 
B = create.bspline.basis(rangeval=range(z), norder=4, nbasis=6)
x = seq(min(z),max(z),length=200); 
out = eval.basis(x,B);   
matplot(x,out,type="l"); iprint=1000,maxfun=250000))

}








