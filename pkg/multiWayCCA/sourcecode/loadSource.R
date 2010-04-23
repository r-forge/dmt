loadSource = function(path) {

	## Load required source code files


	source(paste(path,"sourcecode/cv04.R",sep=""))
	source(paste(path,"sourcecode/boxCI.R",sep=""))
	source(paste(path,"sourcecode/distributions.R",sep=""))
	source(paste(path,"sourcecode/functionsDR.R",sep=""))
	source(paste(path,"sourcecode/functionsManova-from1toS.R",sep="")) # Old 3-way
	source(paste(path,"sourcecode/GDRCCA3-manova-WwithZeros.R",sep=""))
	source(paste(path,"sourcecode/gener2wayCCAdata.R",sep=""))
	source(paste(path,"sourcecode/multiWayCCA.R",sep=""))
	source(paste(path,"sourcecode/mu_intervals8-3way.R",sep="")) # 4.8.09
	source(paste(path,"sourcecode/itsim.R",sep=""))
	source(paste(path,"sourcecode/plotCancors.R",sep=""))
	source(paste(path,"sourcecode/plotMatrixDistribution.R",sep=""))
	source(paste(path,"sourcecode/plotMatrixDistributionWithLabels.R",sep=""))
	source(paste(path,"sourcecode/plotMuIntervals-3way.R",sep="")) # 4.8.09
	source(paste(path,"sourcecode/plotSeries.R",sep=""))
	source(paste(path,"sourcecode/plotSurface.R",sep=""))
	source(paste(path,"sourcecode/plot_var2.R",sep=""))
	source(paste(path,"sourcecode/plotVectorDistribution.R",sep=""))
	source(paste(path,"sourcecode/plotVectorDistributionWithLabels.R",sep=""))
	source(paste(path,"sourcecode/plot_cor.R",sep=""))
	source(paste(path,"sourcecode/Rhat.R",sep=""))
	source(paste(path,"sourcecode/show_clusters.R",sep=""))
	source(paste(path,"sourcecode/write_clust.R",sep=""))

	require("mvtnorm") # Also other packages may be necessary.

}