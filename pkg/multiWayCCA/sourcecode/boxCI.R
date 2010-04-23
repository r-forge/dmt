## Tommi Suvitaival 4.11.09
# Plot a boxplot with normal quartiles but instead of whiskers of length 1.5x'interquantile range', plot whiskers according to the confidence intervals
# Arguments:
# x - position of the box on the x-axis
# y - data
# confidence - the confidence level; default is 0.95
# lwd - line width of the plot
# add - TRUE/FALSE: do/do not add the box to an existing plot

boxCI = function(x, confidence=0.95, add=F, at=1, lwd=1, width=NULL, xlim=NULL, ylim=NULL, xaxt="n", cex=0.5, ylab="") {

	library(gplots)
	
	uncertainty = 1-confidence
	if (add==F) { # Start a new figure.
		boxplot(x=x, at=at, range=1e-18, col="white", outline=F, lwd=lwd, add=F, width=width, xlim=xlim, ylim=ylim, xaxt=xaxt, cex=cex, ylab=ylab)
	}
	
	plotCI(x=at, y=median(x,na.rm=T), uiw=lwd, li=quantile(x,probs=uncertainty/2), ui=quantile(x,probs=(1-uncertainty/2)), add=T, lwd=lwd,lty=1,barcol="black",type="l", gap=0)
	boxplot(x=x, at=at, range=1e-18, col="white", outline=F, lwd=lwd, add=T, cex=cex, width=width, xaxt="n", ylab="", staplewex=0)

	return(TRUE)

}