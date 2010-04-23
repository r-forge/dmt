# plot a dim(X)[1] x dim(X)[2] matrix with a distribution of dim(X)[3] samples

plotMatrixDistributionWithLabels = function(X,file,Xtrue=NULL) {

	nr = nrow(X)
	nc = ncol(X)
	
	if (nc>50) { # Only 50 samples fit into figure.
		nc = 50
		X = X[,1:nc,,drop=F]
	}
	png(paste(file,".png",sep=""),width=200*nc,height=200*nr,type="cairo1")
	layout(matrix(1:(nr*nc),ncol=nc,byrow=T))
	par(cex=1)
	par(cex.lab=2.25,cex.main=1.75,cex.axis=1.75)
	par(mar=c(0.1,2,1,0.2)*2) # margin width of a figure
	cnt = 1
	for (i in 1:nr) {
		for (j in 1:nc) {
			screen(cnt)
# 			boxplot(X[i,j,],ylim=c(-1,1)*1.1*max(abs(X[i,j,])))
			boxCI(X[i,j,],ylim=c(-1,1)*1.1*max(abs(X[i,j,])))
			abline(h=0)
			if (!is.null(Xtrue)) {
				abline(h=Xtrue[i,j],col="red")
			}
			# Labels
			if (j==1)
				title(ylab=paste("Cluster",i),line=2.5)
			if (i==1) {
				if (j==1)
					title(main="Shared loading    ",line=0.9)
				if (j==2)
					title(main="X-specific loading      ",line=0.9)
				if (j==3)
					title(main="Y-specific loading      ",line=0.9)
			}
			cnt = cnt+1
		}
	}
	dev.off()

}