# Input:
# cors - nZ x Nps matrix, where nZ is the number of correlations and Nps the number of posterior samples

plotCancors = function(cors,file,trueCors=NULL) {

	Ncors = nrow(cors)
	png(paste(file,".png",sep=""),width=640,height=640,type="cairo1")
	layout(matrix(1:Ncors,nrow=1,byrow=T))
	for (i in 1:Ncors) {
		screen(i)
		boxplot(cors[i,],ylim=c(0,1))
		if (!is.null(trueCors)) {
			abline(h=trueCors[i],col="red")
		}
	}
	dev.off()

}