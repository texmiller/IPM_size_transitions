########################################################
# Try out multiple_levene_test and multiple_bs_test 
########################################################

### move to the right local directory 
tom = "C:/Users/tm9/Dropbox/github/IPM_size_transitions"
steve = "c:/repos/IPM_size_transitions" 
home = ifelse(Sys.info()["user"] == "Ellner", steve, tom)
setwd(home); setwd("code"); 

add_panel_label <- function(ltype="a"){
    text <- paste(LETTERS[letters==ltype], ")", sep="")
    mtext(text=text, side=3, adj=0)
}

source("variance_diagnostics.R"); 

stopCluster(c1); 
c1<- makeCluster(8); 
registerDoParallel(c1);

N = 500; ## number of fitted values 
nreps = 250; ## number of replicate simulations 
R = 1000; ## number of randomizations for randomization test  

pbin = pspline = pbin2 = pspline2 = pbin3 = pspline3 = numeric(nreps); 
for(jrep in 1:nreps){
	cat("Rep ", jrep, "------------------------------", "\n")
	fitted_vals = sort(2*rbeta(N,3,3)); 
	
	sd_vals = rep(1,N); 
	scaled_resids = rJSU(N, mu=rep(0,N), sigma = sd_vals, nu = -2 + 2*fitted_vals, tau = 2); 
	scaled_resids = (scaled_resids-mean(scaled_resids))/sd(scaled_resids); 
	out = multiple_bartlett_test(fitted_vals, scaled_resids, 3, 10, R) 
	pbin[jrep]=out$p_value
	out = multiple_bs_test(fitted_vals, scaled_resids, 4, 10, R) 
	pspline[jrep]=out$p_value
	
	if(FALSE){
	sd_vals2 =  1  + exp(-2*fitted_vals) 
	scaled_resids2 = rJSU(N, mu=rep(0,N), sigma = sd_vals2, nu = -2 + 2*fitted_vals, tau = 2); 
	scaled_resids2 =(scaled_resids2-mean(scaled_resids2))/sd(scaled_resids2); 
	out = multiple_levene_test(fitted_vals, scaled_resids2, 3, 10, R) 
	pbin2[jrep]=out$p_value
	out = multiple_bs_test(fitted_vals, scaled_resids2, 4, 10, R) 
	pspline2[jrep]=out$p_value
	
	sd_vals3 = 1 + 0.4*dnorm(fitted_vals,mean=1,sd=0.3); 
	scaled_resids3 = rJSU(N, mu=rep(0,N), sigma = sd_vals3, nu = -2 + 2*fitted_vals, tau = 2); 
	scaled_resids3 = (scaled_resids3-mean(scaled_resids3))/sd(scaled_resids3); 
	out = multiple_levene_test(fitted_vals, scaled_resids3, 3, 10, R) 
	pbin3[jrep]=out$p_value
	out = multiple_bs_test(fitted_vals, scaled_resids3, 4, 10, R) 
	pspline3[jrep]=out$p_value
	}
	
	}
	
stopCluster(c1); 

graphics.off(); dev.new(width=14,height=11); 
par(mfcol=c(4,3),mar=c(4,4,2,1),cex.axis=1.3,cex.lab=1.5,mgp=c(2.1,1,0), bty="l");

hist(fitted_vals, 11, xlab = "Scaled residuals", main=""); add_panel_label("a"); 
plot(fitted_vals,scaled_resids,xlab="Fitted values", ylab="Scaled residuals"); add_panel_label("d"); 
matpoints(fitted_vals,cbind(-2*sd_vals,2*sd_vals),type="l",lty=2,col="black"); 
hist(pbin,xlab = "p-values", main=paste0("Multiple Levene test: type-I error rate ", mean(pbin<0.05)));  add_panel_label("g"); 
hist(pspline,xlab = "p-values", main=paste0("Multiple B-spline test: type-I error rate ", mean(pspline<0.05)));  add_panel_label("j"); 

hist(scaled_resids2, 11,xlab = "Scaled residuals", main=""); add_panel_label("b"); 
matplot(fitted_vals,cbind(scaled_resids2, -2*sd_vals2, 2*sd_vals2),type=c("p","l","l"), col="black", lty=2, pch=1,
		xlab="Fitted values", ylab="Scaled residuals"); add_panel_label("e"); 
hist(pbin2,20,xlab = "p-values", main=paste0("Multiple Levene test: power=", mean(pbin2<0.05)));  add_panel_label("h"); 
hist(pspline2,20,xlab = "p-values", main=paste0("Multiple B-spline test: power=", mean(pspline2<0.05)));  add_panel_label("k"); 

hist(scaled_resids3, 11, xlab = "Scaled residuals", main=""); add_panel_label("c"); 
matplot(fitted_vals,cbind(scaled_resids3, -2*sd_vals3, 2*sd_vals3),type=c("p","l","l"), col="black", lty=2, pch=1,
		xlab="Fitted values", ylab="Scaled residuals"); add_panel_label("f"); 
hist(pbin3,20,xlab = "p-values", main=paste0("Multiple Levene test: power=", mean(pbin3<0.05)));  add_panel_label("i"); 
hist(pspline3,20,xlab = "p-values", main=paste0("Multiple B-spline test: power=", mean(pspline3<0.05)));  add_panel_label("l"); 


## dev.copy2pdf(file="../manuscript/figures/test_variance_diagnostics_1000.pdf"); 

