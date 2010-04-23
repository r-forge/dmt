# plot a dim(X)[1] x dim(X)[2] matrix with a distribution of dim(X)[3] samples

plotMatrixDistribution = function(X,file,Xtrue=NULL) {

	nr = nrow(X)
	nc = ncol(X)
	
	if (nc>50) { # Only 50 samples fit into figure.
		nc = 50
		X = X[,1:nc,,drop=F]
	}
	png(paste(file,".png",sep=""),width=200*nc,height=200*nr,type="cairo1")
	layout(matrix(1:(nr*nc),ncol=nc,byrow=T))
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
			cnt = cnt+1
		}
	}
	dev.off()

}