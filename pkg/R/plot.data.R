plot.data <- function(Xm, Ym) {

   data(chr17) # X, Y, locs
   
   X <- Xm; Y <- Ym # override X, Y, keep locs

   par(mfrow=c(3,1),mar=c(5.2,5.5,3,3))

        plot(locs, Y[,1], ylim = range(Y), col="gray", main=paste("Gene mutations (copy number) "),type='l',xlab="Chromosomal location (Mb)",ylab="Signal",yaxt='n',  ,cex.lab = 1.8,cex.main=2,las=1, cex.axis = 1.5)
        for (i in 2:ncol(Y)) {
          lines(locs, Y[,i],col="gray")
        }

        plot(locs, X[,1],ylim=range(X),col="gray",main=paste("Gene activity (expression)", sep=""), ylab="Signal", type='l',xlab="Chromosomal location (Mb)", yaxt = 'n', cex.lab=1.8,cex.main=2,las=1,mgp=c(3.5,1,0), cex.axis = 1.5)
        for (i in 2:ncol(X)) {
                lines(locs, X[,i],col="gray")
        }
}