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
c1<- makeCluster(12); 
registerDoParallel(c1);

N = 500; ## number of fitted values 
nreps = 250; ## number of replicate simulations 
R = 1000; ## number of randomizations for randomization test  

pbin = pspline = pbin2 = pspline2 = numeric(nreps); 
for(jrep in 1:nreps){
	cat("Rep ", jrep, "------------------------------", "\n")
	fitted_vals = sort(2*rbeta(N,3,3)); 
	sd_vals = rep(1,N); 
	scaled_resids = rSST(N,mu=rep(0,N), sigma = sd_vals, nu = exp(-2 + 2*fitted_vals), tau = 5); 
	out = multiple_levene_test(fitted_vals, scaled_resids, 3, 8, R) 
	pbin[jrep]=out$p_value
	out = multiple_bs_test(fitted_vals, scaled_resids, 4, 8, R) 
	pspline[jrep]=out$p_value
	
	#sd_vals2 = 1 + 0.25*sin(3*pi*fitted_vals);  
	sd_vals2 =  1  + exp(-2*fitted_vals) 
	scaled_resids = rSST(N,mu=rep(0,N), sigma = sd_vals2, nu = exp(-2 + 2*fitted_vals), tau = 5); 
	out = multiple_levene_test(fitted_vals, scaled_resids, 3, 8, R) 
	pbin2[jrep]=out$p_value
	out = multiple_bs_test(fitted_vals, scaled_resids, 4, 8, R) 
	pspline2[jrep]=out$p_value
	
	
}

stopCluster(c1); 

graphics.off(); dev.new(width=8.5,height=8); 
par(mfcol=c(4,2),mar=c(4,4,2,1),cex.axis=1.3,cex.lab=1.4,mgp=c(2,1,0), bty="l");

hist(fitted_vals, xlab = "Fitted values", main=""); add_panel_label("a"); 
plot(fitted_vals,sd_vals,xlab="Initial size", ylab = "Standard deviation");   add_panel_label("c"); 
hist(pbin,xlab = "p-values", main=paste0("Multiple Levene test: type-I error rate ", mean(pbin<0.05)));  add_panel_label("d"); 
hist(pspline,xlab = "p-values", main=paste0("Multiple B-spline test: type-I error rate ", mean(pspline<0.05)));  add_panel_label("e"); 

scaled_resids = rSST(N,mu=rep(0,N), sigma = sd_vals, nu = exp(-2 + 2*fitted_vals), tau = 5); 
hist(scaled_resids, xlab = "Fitted values", main=""); add_panel_label("b"); 
plot(fitted_vals,sd_vals2,xlab="Initial size", ylab = "Standard deviation");   add_panel_label("f"); 
hist(pbin,xlab = "p-values", main=paste0("Multiple Levene test: power=", mean(pbin2<0.05)));  add_panel_label("g"); 
hist(pspline,xlab = "p-values", main=paste0("Multiple B-spline test: power=", mean(pspline2<0.05)));  add_panel_label("h"); 
