plotSeries = function(x, fname, ylab=NULL, byrow=F) {

	nSamples = dim(x)[1]
	sp = 1:nSamples # Select which samples are to be plotted.
	if (nSamples>100) { # If the series is long, plot only every 5th sample.
		sp = sp[sp%%5==0]
	}
	
	if (length(dim(x))==2) {
		J = ncol(x)
		png(fname,height=200*J,width=200,type="cairo1")
		layout(matrix(1:J,ncol=1,byrow=byrow))
		if (is.null(ylab)) {
			ylab = 1:J
		}
		for (j in 1:J) {
			screen(j)
			plot(sp,x[sp,j],type="l",ylab=ylab[j])
		}
	} else if (length(dim(x))==3) {
		I = ncol(x)
		J = dim(x)[3]
		png(fname,height=200*J,width=200*I,type="cairo1")
		layout(matrix(1:(I*J),ncol=I,byrow=byrow))
		for (i in 1:I) {
			for (j in 1:J) {
				screen(i*j)
				plot(sp,x[sp,i,j],type="l",ylab=paste(i,",",j,sep=""))
			}
		}
	} else if (length(dim(x))==4) {
		I = ncol(x)
		J = dim(x)[3]
		png(fname,height=300*J,width=300*I,type="cairo1")
		layout(matrix(1:(I*J),ncol=I,byrow=byrow))
		for (i in 1:I) {
			for (j in 1:J) {
				screen(i*j)
				matplot(sp,x[sp,i,j,],type="l",ylab=paste(i,",",j,sep=""))
			}
		}
		
	}
	
	dev.off()
}