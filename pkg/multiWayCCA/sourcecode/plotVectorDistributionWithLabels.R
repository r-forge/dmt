# plot a dim(X)[1] x dim(X)[2] matrix with a distribution of dim(X)[3] samples

plotVectorDistributionWithLabels = function(X,file,Xtrue=NULL,box=TRUE,effN=NULL,printYlab=TRUE) {

	nr = nrow(X)
	
	png(paste(file,".png",sep=""),width=200,height=200*nr,type="cairo1")
	layout(matrix(1:(nr),ncol=1,byrow=T))
	par(cex=1)
# 	par(cex.axis=1.5,cex.lab=1.75,cex.main=3)
	par(cex.axis=3,cex.lab=3.5,cex.main=3)
# 	if (printYlab)
	par(mar=c(0.25,3.4,0.25,0.1)*2) # margin width of a figure
# 	else
# 		par(mar=c(0.25,1,0.25,0.3)*2) # margin width of a figure
	par(lwd=3)
	for (i in 1:nr) {
		screen(i)
		if (box) {
# 			boxplot(X[i,],ylim=c(-1,1)*1.1*max(abs(X[i,])),pch=4,width=2)
			boxCI(X[i,],ylim=c(-1,1)*1.1*max(abs(X[i,])),width=2)
			abline(h=0)
			axis(2,labels=F,lwd=3)
			if (!is.null(Xtrue)) {
				abline(h=Xtrue[i],col="red")
			}
		} else {
			hist(X[i,],main="",xlab="")
		}
# 		if (i==1) {
# 			if (effN==1)
# 				title(main=substitute(~alpha))
# 			else if (effN==2)
# 				title(main=substitute(~beta))
# 			else if (effN==12)
# 				title(main=substitute(~alpha~beta))
# 		}
		if (printYlab) {
			if (i==1)
				title(ylab="Shared",line=4)
			if (i==2)
				title(ylab="X-specific",line=4)
			if (i==3)
				title(ylab="Y-specific",line=4)
		}
	}
	dev.off()

}