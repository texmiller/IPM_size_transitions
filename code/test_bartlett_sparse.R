########################################################
# Try out multiple_levene_test and multiple_bs_test 
########################################################

### move to the right local directory 
tom = "C:/Users/tm9/Dropbox/github/IPM_size_transitions"
steve = "c:/repos/IPM_size_transitions" 
home = ifelse(Sys.info()["user"] == "Ellner", steve, tom)
setwd(home); setwd("code"); 

require(gamlss.dist); 

add_panel_label <- function(ltype="a"){
    text <- paste(LETTERS[letters==ltype], ")", sep="")
    mtext(text=text, side=3, adj=0)
}

source("variance_diagnostics.R"); 


stopCluster(c1); 
c1<- makeCluster(16); 
registerDoParallel(c1);

N = 500; ## number of fitted values 
nreps = 200; ## number of replicate simulations 
R = 1000; ## number of randomizations for randomization test  

pbin = pspline = pbin2 = pspline2 = pbin3 = pspline3 = numeric(nreps); 
for(jrep in 1:nreps){
	cat("Rep ", jrep, "------------------------------", "\n")
	fitted_vals = sort(exp(rnorm(N,0,0.25))); fitted_vals = 2*fitted_vals/max(fitted_vals);
	sd_vals = rep(1,N); 
	scaled_resids = rJSU(N, mu=rep(0,N), sigma = sd_vals, nu = -3 + fitted_vals, tau = 1.5); 
	scaled_resids = (scaled_resids-mean(scaled_resids))/sd(scaled_resids); 
	out = multiple_bartlett_test(fitted_vals, scaled_resids, 3, 10, R) 
	pbin[jrep]=out$p_value

	sd_vals2 =  1  + exp(-2*fitted_vals) 
	scaled_resids2 = rJSU(N, mu=rep(0,N), sigma = sd_vals2, nu = -3 + 1*fitted_vals, tau = 1.5); 
	scaled_resids2 =(scaled_resids2-mean(scaled_resids2))/sd(scaled_resids2); 
	out = multiple_bartlett_test(fitted_vals, scaled_resids2, 3, 10, R) 
	pbin2[jrep]=out$p_value
	
	}

	
stopCluster(c1); 

graphics.off(); dev.new(width=14,height=11); 
par(mfcol=c(3,2),mar=c(4,4,2,1),cex.axis=1.3,cex.lab=1.5,mgp=c(2.1,1,0), bty="l");

hist(scaled_resids, 11, xlab = "Scaled residuals", main=""); add_panel_label("a"); 
plot(fitted_vals,scaled_resids,xlab="Fitted values", ylab="Scaled residuals"); add_panel_label("d"); 
matpoints(fitted_vals,cbind(-2*sd_vals,2*sd_vals),type="l",lty=2,col="black"); 
hist(pbin,xlab = "p-values", main=paste0("Multiple Bartlett test: type-I error rate ", mean(pbin<0.05)));  add_panel_label("g");  

hist(scaled_resids2, 11,xlab = "Scaled residuals", main=""); add_panel_label("b"); 
matplot(fitted_vals,cbind(scaled_resids2, -2*sd_vals2, 2*sd_vals2),type=c("p","l","l"), col="black", lty=2, pch=1,
		xlab="Fitted values", ylab="Scaled residuals"); add_panel_label("e"); 
hist(pbin2,20,xlab = "p-values", main=paste0("Multiple Bartlett test: power=", mean(pbin2<0.05)));  add_panel_label("h"); 
 
