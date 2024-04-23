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

N = 1000; ## number of fitted values 
nreps = 500; ## number of replicate simulations 
R = 100; ## number of randomizations for randomization test  

pbin = pspline = pbin2 = pspline2 = pbin3 = pspline3 = numeric(nreps); 
for(jrep in 1:nreps){
	graphics.off()
	dev.new(width=9,height=4); par(mfrow=c(1,3))
	cat("Rep ", jrep, "------------------------------", "\n")
	fitted_vals = sort(2*rbeta(N,2,4)); 
	
	sd_vals = rep(1,N); 
	scaled_resids = rJSU(N, mu=rep(0,N), sigma = sd_vals, nu = -2 + 2*fitted_vals, tau = 2); 
	scaled_resids = (scaled_resids-mean(scaled_resids))/sd(scaled_resids); 
	out = multiple_levene_test(fitted_vals, scaled_resids, 3, 10, R) 
	pbin[jrep]=out$p_value
	cat(jrep,"case 1, best #bins=",out$bins_min_true,"\n")
	# out = multiple_bs_test(fitted_vals, scaled_resids, 4, 10, R) 
	# pspline[jrep]=out$p_value
	# plot(out$trend_x,out$trend_y,ylim=c(0,2*max(out$trend_y)))

	sd_vals2 =  1  + exp(-2*fitted_vals) 
	scaled_resids2 = rJSU(N, mu=rep(0,N), sigma = sd_vals2, nu = -2 + 2*fitted_vals, tau = 2); 
	scaled_resids2 =(scaled_resids2-mean(scaled_resids2))/sd(scaled_resids2); 
	out = multiple_levene_test(fitted_vals, scaled_resids2, 3, 10, R) 
	pbin2[jrep]=out$p_value
	cat(jrep,"case 2, best #bins=",out$bins_min_true,"\n")
	#out = multiple_bs_test(fitted_vals, scaled_resids2, 4, 10, R) 
	#pspline2[jrep]=out$p_value
	#cat(jrep,"case 2, best df=",out$df_min_true,"\n")
	#plot(out$trend_x,out$trend_y,ylim=c(0,2*max(out$trend_y)))
	
	sd_vals3 = 1 + 0.4*dnorm(fitted_vals,mean=1,sd=0.3); 
	sd_vals3 = 1 + 0.5*exp(-5*abs(fitted_vals-1)); 
	scaled_resids3 = rJSU(N, mu=rep(0,N), sigma = sd_vals3, nu = -2 + 2*fitted_vals, tau = 2); 
	scaled_resids3 = (scaled_resids3-mean(scaled_resids3))/sd(scaled_resids3); 
	out = multiple_levene_test(fitted_vals, scaled_resids3, 3, 10, R) 
	pbin3[jrep]=out$p_value
	cat(jrep,"case 3, best #bins=",out$bins_min_true,"\n")
	#out = multiple_bs_test(fitted_vals, scaled_resids3, 4, 10, R) 
	#pspline3[jrep]=out$p_value
	#cat(jrep," case 3, best df=",out$df_min_true,"\n")
	#plot(out$trend_x,out$trend_y,ylim=c(0,2*max(out$trend_y)))
	Sys.sleep(1); 
	}


