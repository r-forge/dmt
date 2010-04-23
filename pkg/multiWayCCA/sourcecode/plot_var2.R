#####
# plot_var.R
# A function for visual inspection of posterior scale parameter ('lambda') and mu values. Lambda and mu are metabolite-specific parameters.
#
# Parameters:
# data - Rows correspond to samples (N), columns correspond to variables (metabolites, M). The size is NxM.
# case - Binary vector telling the case of each sample. Length is N.
# gender - Binary vector telling the gender of each sample. Length is N.
# state - Binary vector telling the temporal point of each sample. Length is N.
# vars - Matrix of posterior scale parameter (lambda) values. The size is MxT, where T is the number of posterior samples.
# sigmas - Matrix of posterior residual variance (sigma) values. The size is MxT.
# mus - Matrix of posterior metabolite-specific means (from function sample_mus). The size is MxT.
# path - A valid string of characters pointing the directory path where the figure will be saved.
# sorted - TRUE/FALSE: do/do not sort the values (in ascending order)
# varsGen - 
#
# Required libraries: 'gplots' for confidence interval plots.
#
# Comments:
# Some warnings occur because very small confidence intervals cannot be plotted into the figure.
#
# 26.1.09 - Tommi Suvitaival, tsuvitai@cc.hut.fi
#####

plot_var = function(data,case,gender,state,vars,sigmas,mus,path,sorted=TRUE,varsGen=NULL) {

	library('gplots') # library for confidence interval plots

	## Plot: metabolite-specific standard deviations

	# Empirical prior for the metabolite-specific mean and scale parameter is calculated using healthy gender 0 samples from the first week.
	mu0 = rowMeans(data[,(case==0 & gender==0 & state<2)])
	#print(mu0)
	var0 <- apply(data[,(case==0 & gender==0 & state<2)],1,sd)
	
	if (sorted) { # Sort the scale parameter vector into increasing order.
		var0_srt = sort(var0,index.return=T)
	} else { # Do not sort the scale parameter vector but create a similar list.
		var0_srt = list(x=var0,ix=1:length(var0))
	}
	
	sigmas = sqrt(sigmas) # variance to standard deviation

	# mean, max and mean values of the scale parameter over posterior samples, one for each dimension m
	# Length of the mean vector is M.
	vars_mean = rowMeans(vars)[var0_srt$ix]
	# Initialize a matrix for minimum and maximum values. The size is 2xM.
	vars_lim = matrix(nrow=2,ncol=length(var0_srt$x))
	# Find minimum and maximum values for each variable m.
	vars_lim[1,] = apply(vars,1,min)[var0_srt$ix]
	vars_lim[2,] = apply(vars,1,max)[var0_srt$ix]

	# mean, max and mean values of the residual variance over posterior samples, one for each dimension m
	# Length of the mean vector is M.
	sigmas_mean = rowMeans(sigmas)[var0_srt$ix]
	# Initialize a matrix for minimum and maximum values. The size is 2xM.
	sigmas_lim = matrix(nrow=2,ncol=length(var0_srt$x))
	# Find minimum and maximum values for each variable m.
	sigmas_lim[1,] = apply(sigmas,1,min)[var0_srt$ix]
	sigmas_lim[2,] = apply(sigmas,1,max)[var0_srt$ix]

	# vertical axis limits for the plot
	plot_ylim = c(min(c(vars_lim[1,],sigmas_lim[1,],min(var0),min(varsGen),0)),max(c(vars_lim[2,],sigmas_lim[2,],max(var0),max(varsGen))))

	# Open a png file.
	# Adjust the width and height parameters to be able to see details in the figure.
	if (sorted) {
		png(paste(path,"var-sorted.png",sep=""),width=840,height=840,type="cairo1")
	} else {
		png(paste(path,"var.png",sep=""),width=840,height=840,type="cairo1")
	}
	# confidence interval plot of the scale parameter
	plotCI(vars_mean,uiw=1,ui=vars_lim[2,],li=vars_lim[1,],col="blue",ylim=plot_ylim,xlab="",ylab="")
	# Following points are drawn into the same figure.
	par(new=T)
	# confidence interval plot of the scale parameter
	plotCI(sigmas_mean,uiw=1,ui=sigmas_lim[2,],li=sigmas_lim[1,],col="red",ylim=plot_ylim,xlab="",ylab="")
	# values of the empirical prior of the scale parameter
	points(1:length(var0_srt$x),var0_srt$x,col="black",pch="x",xlab="",ylab="")
	if (!is.null(varsGen)) {
		points(1:length(var0_srt$ix),varsGen[var0_srt$ix],col="green",pch="x",xlab="",ylab="")
		legend('topleft',legend=c("Posterior scale parameter (lambda)","Posterior residual standard deviation (sigma)","Empirical standard deviation","Actual scale parameter (generative lambda)"),col=c("blue","red","black","green"),lty=c(1,1,0,0),pch=c("o","o","x","x"))
	} else {
		legend('topleft',legend=c("Posterior scale parameter (lambda)","Posterior residual standard deviation (sigma)","Empirical standard deviation"),col=c("blue","red","black"),lty=c(1,1,0),pch=c("o","o","x"))
	}
	title(main="Confidence intervals of posterior standard deviations",xlab="Variables (metabolites) sorted by the ascending empirical scale parameter value.")
	dev.off() # Close the file.

	## Plot: metabolite-specific means

	# Initialize a matrix for minimum, mean and maximum values. The size is Mx3.
	mus_m = matrix(nrow=length(var0_srt$x),ncol=3)
	mus_m[,2] = rowMeans(mus)[var0_srt$ix] # metabolite-specific means
	mus_m[,1] = apply(mus,1,min)[var0_srt$ix]
	mus_m[,3] = apply(mus,1,max)[var0_srt$ix]

	# vertical axis limits for the plot
	plot_ylim = c(min(c(mus_m[,1],min(mu0))),max(c(mus_m[,3],max(mu0))))

	if (sorted) {
		png(paste(path,"mus_m-sorted.png",sep=""),width=840,height=840,type="cairo1")
	} else {
		png(paste(path,"mus_m.png",sep=""),width=840,height=840,type="cairo1")
	}
	# confidence interval of posterior parameters
	plotCI(mus_m[,2],uiw=1,ui=mus_m[,3],li=mus_m[,1],col="blue",barcol="red",ylim=plot_ylim)
	# empirical values
	points(1:length(var0_srt$ix),mu0[var0_srt$ix],col="black",pch="x",xlab="",ylab="")
	#matplot(mus_m,type="l")
	legend('topleft',legend=c("Posterior mean parameter (mu)","Empirical mean"),col=c("blue","black"),lty=c(1,0),pch=c("o","x"))
	title(main="Confidence intervals of metabolite-specific means",xlab="Variables (metabolites) sorted by the ascending empirical scale parameter value.")
	dev.off()

}