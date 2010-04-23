# Tämä tiedosto on peräisin Ilkalta 22.1.09.

## mu_intervals2
# This function has been edited from Ilkka's code. As an added functionality, empirical effects in the data are generated.
# 4.8.09 - for 3-way analysis

mu_intervals = function(mu,k,Ns,effN,stds=FALSE,mu.gen=NULL,plotQuantiles=FALSE) {

	N_iterations = dim(mu)[1]
	K = dim(mu)[2]
	N_simulations = dim(mu)[3]

	mrg = 2 # margin width of a figure
	#par(mar=c(2,1,1,1)*mrg)
	par(mar=c(2,3.5,0.3,0.3)*mrg)
	#par(mar=c(1,1,1,1)*mrg)
	par(lwd=8)
	
	#par(cex=2,cex.lab=1.5,cex.axis=1.25,cex.main=1.75)
	par(cex=3)
	par(cex.axis=2,cex.main=1.75)
	
	mu_q = array(dim=c(2,length(Ns)))
	if (dim(mu)[3]>1) {
		mu_means = apply(mu,c(2,3),mean,na.rm=T)
		mu_stds = apply(mu,c(2,3),sd,na.rm=T)
		mu_min = apply(mu,c(2,3),min,na.rm=T)
		mu_max = apply(mu,c(2,3),max,na.rm=T)
		if (K>1) {
			for (n in 1:length(Ns)) { # Unlike for other statistics, calculate quantile only for the cluster 'k' that is to be plotted.
				mu_q[,n] = quantile(mu[,k,n],c(0.025,0.975),na.rm=T)
			}
			ylim = c(-1,1)*1.1*max(abs(mu_q))
		} else {
			for (n in 1:length(Ns)) { # Unlike for other statistics, calculate quantile only for the cluster 'k' that is to be plotted.
				mu_q[,n] = quantile(mu[,k,n],c(0.025,0.975),na.rm=T) # CHANGE THIS - 15.5.09
			}
			ylim = c(-1,1)*1.1*max(abs(mu_q))
		}
	} else { # only one simulation
		mu_means = colMeans(mu,na.rm=T)
		mu_stds = apply(mu,2,sd,na.rm=T)
		mu_min = apply(mu,2,min,na.rm=T)
		mu_max = apply(mu,2,max,na.rm=T)
		if (K>1) {
# 			ylim = c(-1,1)*max(c(abs(mu_means[k,]-3*mu_stds[k,]),abs(mu_means[k,]+3*mu_stds[k,])))
			ylim = c(-1,1)*max(c(abs(mu_means[k]-3*mu_stds[k]),abs(mu_means[k]+3*mu_stds[k])))
		} else {
			ylim = c(-1,1)*max(c(abs(mu_means-3*mu_stds),abs(mu_means+3*mu_stds)))
		}
	}

	
	#plot(Ns,mu_means[k,],ylim=ylim,lty=1,main="",log='x',xlab="",ylab="")
	plot(Ns,mu_means[k,],ylim=ylim,lty=1,main="",log='x',xlab="",ylab="",cex=0.75,las=3,pch=19)
	if (plotQuantiles) {
		lines(Ns,mu_q[1,],lty=1)
		lines(Ns,mu_q[2,],lty=1)
	} else {
		lines(Ns,(mu_means[k,]-2*mu_stds[k,]),lty=3)
		lines(Ns,(mu_means[k,]+2*mu_stds[k,]),lty=3)
		lines(Ns,(mu_min[k,]),lty=5)
		lines(Ns,(mu_max[k,]),lty=5)
	}
	abline(h=0)
	axis(1,labels=F,lwd=8)
	axis(2,labels=F,lwd=8)
	
# 	if (k==1) {
# 		if (effN==1) {
# 			title(main=substitute(~alpha))
# 		} else if (effN==2) {
# 			title(main=substitute(~beta))
# 		} else if (effN==12) {
# 			title(main=substitute(~alpha~beta))
# 		}
# 	}
# 	if (effN==1) {
# 		if (k==1)
# 			title(ylab="Shared     ",line=2.75)
# 		else if (k==2)
# 			title(ylab="X-specific     ",line=2.75)
# 		else if (k==3)
# 			title(ylab="Y-specific     ",line=2.75)
# 	}
	par(cex.lab=3.25)
	if (k==1) {
		if (effN==1)
			title(ylab=substitute(~alpha),line=2.75)
		else if (effN==2)
			title(ylab=substitute(~beta),line=2.75)
		else if (effN==12)
			title(ylab=substitute(~"("~alpha~beta~")    "),line=2.75)
	} else if (k==2) {
		if (effN==1)
			title(ylab=substitute(~alpha^"x"),line=2.75)
		else if (effN==2)
			title(ylab=substitute(~beta^"x"),line=2.75)
		else if (effN==12)
			title(ylab=substitute(~"("~alpha~beta~")"^"x    "),line=2.75)
	} else if (k==3) {
		if (effN==1)
			title(ylab=substitute(~alpha^"y"),line=2.75)
		else if (effN==2)
			title(ylab=substitute(~beta^"y"),line=2.75)
		else if (effN==12)
			title(ylab=substitute(~"("~alpha~beta~")"^"y    "),line=2.75)
	}
	# Following part is added for 3-way analysis
	if (k==1) {
		if (effN==3)
			title(ylab=substitute(~gamma),line=2.75)
		else if (effN==13)
			title(ylab=substitute(~alpha~gamma),line=2.75)
		else if (effN==23)
			title(ylab=substitute(~beta~gamma),line=2.75)
		else if (effN==123)
			title(ylab=substitute(~alpha~beta~gamma),line=2.75)
	}
	# end of 3-way analysis
	
# 	par(cex.lab=2.5)
# 	if (k==3) {
# 		if (stds) {
# 			title(xlab=expression(~sigma))
# 		} else {
# 			title(xlab="n samples",line=4.5)
# 		}
# 	}
	#title(ylab="effect")
	
	# empirical mu
# 	if (!is.null(mu.gen)) {
# 		#points(Ns,mu.gen[,k],pch="x") # empirical 'mu_'
# 		points(Ns,mu.gen,pch="x") # empirical 'mu_'
# 	}

}
