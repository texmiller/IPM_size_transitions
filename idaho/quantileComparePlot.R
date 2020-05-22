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
		points(qTrue$bin_mean, unlist(qTrue[,qnames[j]]),pch=1,lwd=2,col=alpha("red",alpha_scale),cex=1.4)
		points(qTrue$bin_mean, apply(qSim_j,1,median),pch=1,lwd=2,col=alpha("black",alpha_scale),cex=1.4)
		trueBinQ[,j] = unlist(qTrue[,qnames[j]]);
		simBinQ[,j] = apply(qSim_j,1,median);
		if(j==1) legend("topleft",legend=c("Model simulations","Median of simulations","Data"),
			col=c(alpha("gray",0.5),alpha("black",alpha_scale), alpha("red",alpha_scale)),
			pch=1,lwd=2,cex=1.1,bty="n"); 
		add_panel_label(letters[j])
	}	
	
	Qcolors=c(rep(alpha("red",alpha_scale),5),rep(alpha("black",alpha_scale),5)); 
	matplot(qTrue$bin_mean,cbind(trueBinQ,simBinQ), col=Qcolors, pch=1, lty=2, lwd=2, cex=1.4,type="o",
		xlab="Bin mean", ylab="All percentiles"); 
	add_panel_label("f"); 
}	
		


		
		