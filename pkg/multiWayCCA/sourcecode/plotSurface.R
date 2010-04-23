plotSurface = function(X,file,Xtrue=NULL) {

	ymar = 1:nrow(X)
	if (!is.null(Xtrue)) {
		xmar = 1:(2*ncol(X))
		png(paste(file,".png",sep=""),width=50*ncol(X),height=20*nrow(X),type="cairo1")
		image(xmar,ymar,t(cbind(X,Xtrue)))
	} else {
		xmar = 1:ncol(X)
		png(paste(file,".png",sep=""),width=50*ncol(X),height=20*nrow(X),type="cairo1")
		image(xmar,ymar,t(X))
	}
	dev.off()

}