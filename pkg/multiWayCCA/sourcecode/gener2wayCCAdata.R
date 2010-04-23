# This function is from Ilkka, originally named "gener".
# Versions:
# 30.4.09 - treatments is a list containing MANOVA treatments. effects is a matrix containing component-specific effects of z.
# 22.4.09 - genParams is a list containing items 'Xnoise' and 'Ynoise'.

gener2wayCCAdata <- function(N_samples, nZ, xdim, ydim, xxdim, yydim, W=NULL, treatments=NULL, effects=NULL, params=NULL) {
	#xdim and ydim are number of latent variables, xx and yy the actual dimension of data

	require("mvtnorm") # Also other packages may be necessary.

	library(MCMCpack)
	z<-array(dim=c(N_samples,nZ))
	z2<-array(dim=c(N_samples,nZ))
	wx<-array(dim=c(xdim,nZ))
	wy<-array(dim=c(ydim,nZ))
	wx2<-array(dim=c(xdim,nZ))
	wy2<-array(dim=c(ydim,nZ))
	x<-array(dim=c(N_samples,xdim))
	y<-array(dim=c(N_samples,ydim))

	a0<-1
	b0<-0.1

	sigx<-diag(x=1,nrow=xdim,ncol=xdim)
	sigy<-diag(x=1,nrow=ydim,ncol=ydim)
	sigmin<-diag(x=1,nrow=nZ,ncol=nZ)
	mupri_x<-array(0,dim=c(xdim,1))
	mupri_y<-array(0,dim=c(ydim,1))
	mupri_min<-array(0,dim=c(nZ,1))


	###############
	#mupri_x[]<-0
	#mupri_y[]<-0
	#mupri_min[]<-0

	z<-mvrnorm(n = N_samples, mupri_min, sigmin, tol = 1e-6, empirical = FALSE)
	z2<-mvrnorm(n = N_samples, mupri_min, sigmin, tol = 1e-6, empirical = FALSE)
	
	# Treatment-specific effects
	if (!is.null(treatments)) {
		for (k in 1:nZ) {
			z[(treatments$a==1),k] = z[(treatments$a==1),k]+effects[1,k]
			z[(treatments$b==1),k] = z[(treatments$b==1),k]+effects[2,k]
			z[(treatments$a==1&treatments$b==1),k] = z[(treatments$a==1&treatments$b==1),k]+effects[3,k]
		}
	}

	mu_x<-mvrnorm(n = 1, mupri_x, sigx, tol = 1e-6, empirical = FALSE)
	mu_y<-mvrnorm(n = 1, mupri_y, sigy, tol = 1e-6, empirical = FALSE)
	beta_is<-rinvgamma(xdim, shape=a0, scale=b0)
	beta_is2<-rinvgamma(ydim, shape=a0, scale=b0)
	beta_is[]<-1
	beta_is2[]<-1

	psix<-diag(1,xdim)+4
	psiy<-diag(1,ydim)+4
	if (!is.null(W)) {
		if (nrow(W$X)!=xdim | ncol(W$X)!=nZ | nrow(W$Y)!=ydim | ncol(W$Y)!=nZ)
			print("ERROR: Wrong dimensions of W")
		wx = W$X
		wy = W$Y
	} else {
		for (i in 1:nZ){
		bi<-mvrnorm(n = 1,mupri_x , sigx*beta_is[i], tol = 1e-6, empirical = FALSE)
		wx[,i]<-bi
		}

		for (i in 1:nZ){
		bi<-mvrnorm(n = 1,mupri_y , sigy*beta_is2[i], tol = 1e-6, empirical = FALSE)
		wy[,i]<-bi
		}
	}


	for (i in 1:nZ) {
		bi<-mvrnorm(n = 1,mupri_x , sigx*beta_is[i], tol = 1e-6, empirical = FALSE)
		wx2[,i]<-bi
	}

	for (i in 1:nZ) {
		bi<-mvrnorm(n = 1,mupri_y , sigy*beta_is2[i], tol = 1e-6, empirical = FALSE)
		wy2[,i]<-bi
	}


	for(i in 1:N_samples){
		x[i,]<-mvrnorm(n = 1,mu_x+wx %*% z[i,],psix , tol = 1e-6, empirical = FALSE)
		y[i,]<-mvrnorm(n = 1,mu_y+wy %*% z[i,],psiy , tol = 1e-6, empirical = FALSE)
	}
	print(psix)
	print(psiy)
	init_clus_x<-floor(runif(xxdim, min=1, max=xdim+0.99))
	init_clus_y<-floor(runif(yydim, min=1, max=ydim+0.99))
	V_xx<-array(0,dim=c(xxdim,xdim))
	V_yy<-array(0,dim=c(yydim,ydim))
	for(i in 1:xxdim){
		if(init_clus_x[i]>0) {
			V_xx[i,init_clus_x[i]]=1;
		}
	}
	for(i in 1:yydim){
		if(init_clus_y[i]>0){
			V_yy[i,init_clus_y[i]]=1;
		}
	}
	xx<-V_xx%*%t(x)+params$Xnoise*array(rnorm(xxdim*N_samples),dim=c(xxdim,N_samples))
	yy<-V_yy%*%t(y)+params$Ynoise*array(rnorm(yydim*N_samples),dim=c(yydim,N_samples))
	print(params$Xnoise)

	gener <- list(xlat=x, ylat=y, WX=wx, WY=wy, X=xx, Y=yy, covariates=list(), Vx=V_xx, Vy=V_yy)
	gener$covariates$a = treatments$a
	gener$covariates$b = treatments$b
	gener$covariates$c = rep(0,length(treatments$a))
	
	return(gener)
	
}
 
