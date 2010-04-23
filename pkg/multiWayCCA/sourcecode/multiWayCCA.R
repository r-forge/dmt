## Multi-way, multi-view analysis
# This script is for the analysis of 3-way 2-view data.
# Part of the "multi-way, multi-view analysis" package.
# You will have to write your own procedure for loading your data.
# Tommi Suvitaival and Ilkka Huopaniemi 2009

# ## Set path for the root folder of the package
# 
# path = "/share/mi/workspace/tsuvitai/mi/CCA/package/package-100317/" # This is the root folder of the package. The code should be located in folder "sourcecode" and the results will be saved into sub-folders of the folder "results".
# 
# 
# ## PARAMETERS
# 
# ## Set unique name for the experiment
# 
# runId = "GDRCCA-100325-Kiina-plasmaAndSpleen-outliersRemoved-short" # Keep this ID unique for each simulation. The program will save the results into a folder named after this.
# 
# ## Simulation parameters
# 
# maxBurnIn = 8 # Maximum length of the burn-in if the chains have not converged.
# Niterfinal <- 10 # Number of iterations
# 
# dataset = 1 # 0/1: use generated data (0) or real data (1)
# doPlotting = TRUE # TRUE/FALSE: do/do not perform plotting
# realDataFiles = c(5,4) # which two Kiina data files to select: 1~cerebra, 2~heart, 3~liver, 4~lung, 5~plasma, 6~spleen, 7~tumor. This parameter selects which views from the complete data set to use in the analysis.
# takeLog = TRUE # TRUE/FALSE: do/do not perform a log-transform on real data
# 
# nXlat <- 6 # dimensionality of the latent variable of the first data set in learning
# nYlat <- 5 # dimensionality of the latent variable of the second data set in learning
# 
# ## DATA
# 
# ## Generate or load data
# if (dataset==1) { ## Load and process real data. Write your own data loading procedure here!
# 	data = processKiinaData(realDataFiles) # Imports the two views and covariates 'case', 'gender' and 'state'.
# 	N = ncol(data$X) # Number of samples (columns) in the data set. Must be equal in both data sets 'X' and 'Y'. Also covariate vectors are required to be of this length.
# 	
# 	# If performing less-than-3-way analysis, set some covariates (e.g. 'state') to zero. The vectors are anyhow required by the sampler.
# # 	data$case = rep(0,N)
# # 	data$gender = rep(0,N)
# 	data$state = rep(0,N)
# 	
# 	# Number of features in each view: they do not have to be equal.
# 	nX = nrow(data$X) # the number of features in view X
# 	nY = nrow(data$Y) # the number of features in view Y
# 	
# 	if (takeLog) { # Optionally, take logarithm of the data to make it more normally distributed. Watch out for zero values when taking the logarithm!
# 		data$X = log(data$X+0.01)
# 		data$Y = log(data$Y+0.01)
# 	}
# 	if (sum(realDataFiles==7)>0) { # If tumor used, only one-way analysis possible.
# 		data$case = data$case*0
# 	}
# 	
# } else if (dataset==0) { ## Generate synthetic data from the model
# 	## Generated data parameters
# 	N <- 100 # number of samples in the generated data, if model behavior is not analysed as a function of the sample-size
# 	nX = 25 # number of features in view X
# 	nY = 25 # number of features in view Y
# 	nXlatGen = nXlat # the number of view-specific latent components from which the view X is generated from
# 	nYlatGen = nYlat # the number of view-specific latent components from which the view Y is generated from
# 	nZgen = 3 # dimensionality of the "higher level" latent variable in learning. Should be 3: one shared latent variable plus one view-specific variable for each view.
# 	fixedW = TRUE # TRUE/FALSE: do/do not use fixed projection from z to xlat and ylat on generated data, defined by user
# 	genParams = list()
# 	genParams$Xnoise = 1 # noise of the data set X
# 	genParams$Ynoise = genParams$Xnoise # noise of the data set Y
# 	## Projection matrix 'W'
# 	# Proper CCA data generation: data sets X and Y have a shared latent variable 'z'.
# 	if (fixedW) { # Use fixed projection 'W' from 'z' to 'xlat' and 'ylat'.
# 		Wgen = list()
# 		Wgen$X = matrix(c(1,0,-1, 1,0,0, 0,0,0),ncol=nZgen,byrow=F) # projection of 'z' to 'xlat'
# 		Wgen$Y = matrix(c(1,0,-1, 0,0,0, 0,1,0),ncol=nZgen,byrow=F) # projection of 'z' to 'ylat'
# 	} else { # Use random projection.
# 		Wgen = NULL
# 	}
# 	## Covariates
# 	treatments = list()
# 	#data$case = rep(0,N)
# 	treatments$case = rep(c(0,1),length.out=N) # the first covariate
# 	treatments$gender = rep(c(0,0,1,1),length.out=N) # the second covariate
# 	#treatments$gender = rep(0,length.out=N)
# 	treatments$state = rep(c(0,0,0,0,1,1,1,1),length.out=N) # the third covariate
# 	#treatments$state = rep(0,N)
# 	## Covariate effects
# 	# Three-way effects of latent variable 'z'. Number of columns is 'nZ'. In each component: 'mu_c', 'mu_g', 'mu_s', 'mu_cg', 'mu_sc', 'mu_sg', 'mu_scg'. Warning: effects here are in wrong order compared to 'sampleEffects'!
# 	effects = array(c(0,0,0,0,0,0,0, 0,0,0,2,0,0,0, 0,0,0,0,0,0,0),dim=c(7,nZgen))
# 	## Generate data
# 	data = generCCAdata3way(N,nZgen,nXlatGen,nYlatGen,nX,nY,W=Wgen,treatments=treatments,effects=effects,params=genParams)
# } ## End of data loading/generation

multiWayCCA = function(data=data, nXlat=nXlat, nYlat=nYlat, takeLog=takeLog, maxBurnIn=NburnIn, Niterfinal=Niterfinal, path=path, runId=runId) {

	time.start = Sys.time()

	## Convert covariates to be of the same form as in the implementation.

	data$case = data$covariates$a
	data$gender = data$covariates$b
	data$state = data$covariates$c
	data$covariates = NULL

	## Create directories

	pathRun = paste(path,"results/",runId,"/",sep="")
	dir.create(pathRun) # Create a directory for plottings.
	file.copy(paste(path,"scripts/",runId,".R",sep=""),paste(pathRun,runId,".R",sep="")) # Save this script.
	
	nZ = 3 # dimensionality of the common latent variable in learning. Should be 3: one shared latent variable plus one data set-specific variable for each data set.

	## More parameters. Do not change these unless you know what you are doing.

	Niter <- 5 # number of Gibbs samples to draw between printing
	Nseq <- 1 # How many parallel sampling sequences. Currently, only samples from the first sequence are saved and analyzed.
	unitScale = TRUE # TRUE/FALSE: do/do not scale variables to be of unit scale

	# Number of features in each view: they do not have to be equal.
	nX = nrow(data$X) # the number of features in view X
	nY = nrow(data$Y) # the number of features in view Y
	
	if (takeLog) { # Optionally, take logarithm of the data to make it more normally distributed. Watch out for zero values when taking the logarithm!
		data$X = log(data$X+0.01)
		data$Y = log(data$Y+0.01)
	}

	set.seed(as.numeric(Sys.time()))
	seed = .Random.seed # Save the random generator seed!
	# set.seed(seed) # Recover previously used seed.

	## Define prior parameters
	setParameters <- function() {
		betaFirst <<- 1e-1
		betaPrior <<- 1e-1
		sigFirst <<- 1e-1
		sigPrior <<- 1e-1
		psi_prior_type <<- 1
		psi_prior <<- 1
		mu_vec <<- NULL # 'mu' no longer is a parameter of data but of latent variables.
		mu_prior <<- 1
		
		arvonta <<- TRUE # FALSE/TRUE: do/do not allow empty clusters in 'xLat' and 'yLat'
		jama <<- FALSE # TRUE/FALSE: do/do not use left-over cluster in 'xLat' and 'yLat'
		sampling <<- list()
		sampling$MuX <<- FALSE # TRUE/FALSE: do/do not estimate the feature-specific mean parameter
		sampling$Mu <<- FALSE
		if (unitScale) { # If the data is scaled, do not estimate scale parameter.
			sampling$scale <<- FALSE
		} else {
			sampling$scale <<- TRUE
		}
		## Covariate effects to be estimated
		# Decide which covariate effects are to be sampled. The order is (1) 'mu_c', (2) 'mu_g', (3) 'mu_cg', (4) 'mu_s', (5) 'mu_sc', (6) 'mu_sg', (7) 'mu_scg'. Some of these can be set to false to allow non-standard analyses.
		sampleEffects <<- rep(FALSE,7)
		sampleEffects[1] <<- (sum(data$case)>0) # mu_c
		sampleEffects[2] <<- (sum(data$gender)>0) # mu_g
		sampleEffects[3] <<- (sum(data$case&data$gender)>0) # mu_cg
		sampleEffects[4] <<- (sum(data$state>1)>0) # mu_s
		sampleEffects[5] <<- (sum(data$state>1 & data$case==1)>0) # mu_sc
		sampleEffects[6] <<- (sum(data$state>1 & data$gender==1)>0) # mu_sg
		sampleEffects[7] <<- (sum(data$state>1 & data$case==1 & data$gender==1)>0) # mu_scg

		## Zero values of projection matrix 'W'
		# Limit some values of projection matrix 'W' to zero, leading to data set-specific components of latent variable 'z'.
		setWtoZero <<- array(0,dim=c(2,nZ))
		setWtoZero[1,3] <<- 1 # Third z does not interact with xlat.
		setWtoZero[2,2] <<- 1 # Second z does not interact with ylat.
	}
	setParameters()

	## END OF PARAMETERS


	S = max(data$state) # last time-point or state

	## The data pre-processing
	# Transform the features to be zero-meaned.
	data$X = data$X-rowMeans(data$X[,(data$case==0 & data$gender==0 & data$state==0)]) # 2-level 3-way 6.8.09 
	data$Y = data$Y-rowMeans(data$Y[,(data$case==0 & data$gender==0 & data$state==0)]) # 2-level 3-way 6.8.09 
	# Transform the features to be unit-scale.
	if (unitScale) {
		data$X = data$X/apply(data$X[,(data$case==0 & data$gender==0 & data$state==0)],1,sd) # 2-level 3-way 6.8.09 
		data$Y = data$Y/apply(data$Y[,(data$case==0 & data$gender==0 & data$state==0)],1,sd) # 2-level 3-way 6.8.09 
	}

	## Compute empirical priors
	priors = list()
	priors$muX0 = rowMeans(data$X[,(data$case==0 & data$gender==0 & data$state==0)]) # 2-level 3-way 6.8.09 
	priors$muY0 = rowMeans(data$Y[,(data$case==0 & data$gender==0 & data$state==0)]) # 2-level 3-way 6.8.09 
	priors$varsX0 = apply(data$X[,(data$case==0 & data$gender==0 & data$state==0)],1,sd) # 2-level 3-way 6.8.09 
	priors$varsY0 = apply(data$Y[,(data$case==0 & data$gender==0 & data$state==0)],1,sd) # 2-level 3-way 6.8.09 
	priors$N0 = N #sum(data$case==0) # weight of the prior
	#priors$N0 = sum(data$case==0 & data$gender==0 & data$state<2) # weight of the prior


	## SAMPLING

	## Initialize sampling
	first <- TRUE
	lists <- array(list(),Nseq)
	models <- array(list(),Nseq)
	rnd <- 0

	print(paste("Sampling with ",N," samples",sep=""))

	temp <- sampleMPCCA(data,priors,Niter,nXlat,nYlat,nZ)

	while (rnd*Niter<maxBurnIn) { ## Burn-in
		rnd <- rnd + 1
		for (i in 1:Nseq) {
			print(runId)
			print(Sys.time())
			print(paste("Sampling sequence",i,"/ samples from",(rnd-1)*Niter+1,"to",rnd*Niter))
			print("====================")
			if(first) { # Create a new model and start sampling
				temp <- sampleMPCCA(data,priors,Niter,nXlat,nYlat,nZ)
	# 			lists[[i]] <- stripList(temp$storage)
				first <- FALSE
			} else { # Continue sampling an existing model
				temp <- sampleMPCCAcont(data,temp$model,priors,Niter)
	# 			lists[[i]] <- concatenateLists(lists[[i]],temp$storage,dropN=Niter/2)
			}
	# 		model <- temp$model # Save the Gibbs sample.
		}
	} ## End of burn-in

	modelBi = temp$model

	for (i in 1:Nseq) { ## Sampling after burn-in
		print(runId)
		print(Sys.time())
		print(paste("Sampling sequence ",i,"/ samples from",rnd*Niter+1,"to",rnd*Niter+Niterfinal))
		print("====================")
		samples <- sampleMPCCAcont(data,modelBi,priors,Niterfinal)
	# 	lists[[i]] <- concatenateLists(lists[[i]],temp$storage)
	# 	models[[i]] <- temp$model
	} ## End of sampling after burn-in

	## END OF SAMPLING

	## Save the Gibbs samples into file
	save(list=ls(),file=paste(pathRun,"simulation.RData",sep="")) # Save the results.

	## Collect posterior samples

	# The lists are decomposed into matrices corresponding to each variable.
	Nps = length(samples$storage) # number of posterior samples
	posterior = list() # a list for the posterior samples

	# The dimensionality reduction part
	posterior$Vx = array(dim=c(nX,nXlat,Nps)) # clusterings of the data X
	posterior$Vy = array(dim=c(nY,nYlat,Nps)) # clusterings of the data Y
	posterior$xlat = array(dim=c(nXlat,N,Nps)) # 'xlats'
	posterior$ylat = array(dim=c(nYlat,N,Nps)) # 'ylats'
	posterior$SigmaX = array(dim=c(nX,Nps))
	posterior$SigmaY = array(dim=c(nY,Nps))
	posterior$varsX = array(dim=c(nX,Nps)) # variable-specific variance for data X
	posterior$varsY = array(dim=c(nY,Nps)) # variable-specific variance for data Y
	posterior$muX = array(dim=c(nX,Nps)) # variable-specific mean of data X
	posterior$muY = array(dim=c(nY,Nps)) # variable-specific mean of data Y
	# The bCCA part
	posterior$z = array(dim=c(nZ,N,Nps)) # latent variable of Xlat and Ylat
	posterior$mu = array(dim=c(nXlat+nYlat,Nps)) # latent variable-specific mean
	posterior$PsiX = array(dim=c(nXlat,nXlat,Nps))
	posterior$PsiY = array(dim=c(nYlat,nYlat,Nps))
	posterior$Ux = array(dim=c(nXlat,nZ,Nps))
	posterior$Uy = array(dim=c(nYlat,nZ,Nps))
	posterior$Wx = array(dim=c(nXlat,nZ,Nps))
	posterior$Wy = array(dim=c(nYlat,nZ,Nps))
	posterior$WxRSS = array(dim=c(nXlat,Nps))
	posterior$WyRSS = array(dim=c(nYlat,Nps))
	posterior$beta = array(dim=c(nZ,Nps))
	posterior$like = list()
	posterior$like$z = rep(NA,Nps)
	posterior$like$xLat = rep(NA,Nps)
	posterior$like$yLat = rep(NA,Nps)
	posterior$like$X = rep(NA,Nps)
	posterior$like$Y = rep(NA,Nps)
	# The MANOVA part
	if (sampleEffects[1])
		posterior$mu_c = array(dim=c(nZ,Nps))
	if (sampleEffects[2])
		posterior$mu_g = array(dim=c(nZ,Nps))
	if (sampleEffects[3])
		posterior$mu_cg = array(dim=c(nZ,Nps))
	if (sampleEffects[4]) {
		if (S==1)
			posterior$mu_s = array(dim=c(nZ,Nps))
		else
			posterior$mu_s = array(dim=c(nZ,S,Nps))
	}
	if (sampleEffects[5]) {
		if (S==1)
			posterior$mu_sc = array(dim=c(nZ,Nps))
		else
			posterior$mu_sc = array(dim=c(nZ,S,Nps))
	}
	if (sampleEffects[6]) {
		if (S==1)
			posterior$mu_sg = array(dim=c(nZ,Nps))
		else
			posterior$mu_sg = array(dim=c(nZ,S,Nps))
	}
	if (sampleEffects[7]) {
		if (S==1)
			posterior$mu_scg = array(dim=c(nZ,Nps))
		else
			posterior$mu_scg = array(dim=c(nZ,S,Nps))
	}

	# Canonical correlations
	corBCCA = array(dim=c(nZ,Nps))
	corCCA = array(dim=c(min(nXlat,nYlat),Nps))
	tmp = list()
	# Also effects need to be saved, later.

	# Take posterior samples from the first chain only.
	for (j in 1:Nps) {
		# The dimensionality reduction part
		posterior$Vx[,,j] = samples$storage[[j]]$Vx
		posterior$Vy[,,j] = samples$storage[[j]]$Vy
		posterior$xlat[,,j] = samples$storage[[j]]$xlat
		posterior$ylat[,,j] = samples$storage[[j]]$ylat
		posterior$SigmaX[,j] = samples$storage[[j]]$SigmaX
		posterior$SigmaY[,j] = samples$storage[[j]]$SigmaY
		posterior$varsX[,j] = samples$storage[[j]]$varsX
		posterior$varsY[,j] = samples$storage[[j]]$varsY
		posterior$muX[,j] = samples$storage[[j]]$muX
		posterior$muY[,j] = samples$storage[[j]]$muY
		# The bCCA part
		posterior$z[,,j] = samples$storage[[j]]$z
		posterior$mu[,j] = samples$storage[[j]]$mu
		posterior$PsiX[,,j] = samples$storage[[j]]$PsiX
		posterior$PsiY[,,j] = samples$storage[[j]]$PsiY
		posterior$Wx[,,j] = samples$storage[[j]]$W[1:nXlat,]
		posterior$Wy[,,j] = samples$storage[[j]]$W[(nXlat+(1:nYlat)),]
		posterior$WxRSS[,j] = rowSums((posterior$Wx[,,j,drop=F])^2)
		posterior$WyRSS[,j] = rowSums((posterior$Wy[,,j,drop=F])^2)
		posterior$beta[,j] = samples$storage[[j]]$beta
		# Feed function 'findRotation' with posterior samples 'j'.
		tmp = findRotation(posterior$xlat[,,j,drop=F],posterior$ylat[,,j,drop=F],samples$storage[[j]])
		corBCCA[,j] = tmp$cor
		posterior$Ux[,,j] = tmp$Ux
		posterior$Uy[,,j] = tmp$Uy
		if (sampleEffects[1])
			posterior$mu_c[,j] = samples$storage[[j]]$mu_c
		if (sampleEffects[2])
			posterior$mu_g[,j] = samples$storage[[j]]$mu_g
		if (sampleEffects[3])
			posterior$mu_cg[,j] = samples$storage[[j]]$mu_cg
		if (sampleEffects[4])
			posterior$mu_s[,j] = samples$storage[[j]]$mu_s
		if (sampleEffects[5])
			posterior$mu_sc[,j] = samples$storage[[j]]$mu_sc
		if (sampleEffects[6])
			posterior$mu_sg[,j] = samples$storage[[j]]$mu_sg
		if (sampleEffects[7])
			posterior$mu_scg[,j] = samples$storage[[j]]$mu_scg
		if (nZ>1) {
			corCCA[,j] = cancor(t(posterior$xlat[,,j]),t(posterior$ylat[,,j]))$cor
		} else if (nXlat==1 & nYlat==1) {
			corCCA[,j] = cancor(as.matrix(posterior$xlat[,,j]),as.matrix(posterior$ylat[,,j]))$cor
		}
		# likelihoods
		posterior$like$z[j] = samples$storage[[j]]$like$z
		posterior$like$xLat[j] = samples$storage[[j]]$like$xlat
		posterior$like$yLat[j] = samples$storage[[j]]$like$ylat
		posterior$like$X[j] = samples$storage[[j]]$like$X
		posterior$like$Y[j] = samples$storage[[j]]$like$Y
	} ## End of collect posterior samples

		## For generated data, arrange the components of 'ylat' and 'xlat' into order of generation.
		Vxsum = apply(posterior$Vx,c(1,2),sum)
		Vysum = apply(posterior$Vy,c(1,2),sum)
		Vx.avg = Vxsum/max(Vxsum) # set values to be from zero to one
		Vy.avg = Vysum/max(Vysum) # set values to be from zero to one

		Vx.mode = (Vx.avg==apply(Vx.avg,1,max))*1
		Vy.mode = (Vy.avg==apply(Vy.avg,1,max))*1

		if (dataset==0 & nXlat==nYlat) { # For generated data the correct ordering of clusters is known.
			Vx.ind = compare_clusters(data$Vx,Vx.avg)
			Vy.ind = compare_clusters(data$Vy,Vy.avg)
		} else { # real data - no correct ordering
			Vx.ind = 1:nXlat
			Vy.ind = 1:nYlat
		}

		## Save the Gibbs samples into file
		save(list=ls(),file=paste(pathRun,"simulation.RData",sep="")) # Save the results.

	## Printing and plotting
	if (doPlotting) {
		print("Plotting")
		plotSeries(cbind(posterior$like$z,posterior$like$xLat,posterior$like$yLat,posterior$like$X,posterior$like$Y), fname=paste(pathRun,"like.png",sep=""), ylab=c("z", "xLat", "yLat", "X", "Y")) # likelihoods
		if (dataset!=0) { # If real data, write metabolite names under each cluster.
			write_clust(NULL,data$Xnames,Vx.mode,Vx.avg,paste(pathRun,"Vx.txt",sep=""))
			write_clust(NULL,data$Ynames,Vy.mode,Vy.avg,paste(pathRun,"Vy.txt",sep=""))
		}
		textFile = paste(pathRun,"cors.txt",sep="")
		write(paste("\n",N," common samples.\nDataset X: ",nX," variables\nDataset Y: ",nY," variables",sep=""),file=textFile,append=T)
	# 		write(paste("\nSampler converged in ",rnd*Niter," iterations.",sep=""),file=textFile,append=T)
		write(paste("\nNumber of posterior samples: ",Nps,sep=""),file=textFile,append=T)
		# Plot the posterior distributions of variables
		plotVectorDistribution(posterior$mu,file=paste(pathRun,"mu",sep="")) # Wx
		plotMatrixDistribution(posterior$PsiX[Vx.ind,Vx.ind,],file=paste(pathRun,"PsiX",sep="")) # PsiX
		plotMatrixDistribution(posterior$PsiY[Vy.ind,Vy.ind,],file=paste(pathRun,"PsiY",sep="")) # PsiY
	# 		plotMatrixDistribution(posterior$Wx[Vx.ind,,,drop=F],file=paste(pathRun,"Wx",sep="")) # Wx
	# 		plotMatrixDistribution(posterior$Wy[Vy.ind,,,drop=F],file=paste(pathRun,"Wy",sep="")) # Wy
		plotMatrixDistributionWithLabels(posterior$Wx[Vx.ind,,,drop=F],file=paste(pathRun,"Wx",sep="")) # Wx
		plotMatrixDistributionWithLabels(posterior$Wy[Vy.ind,,,drop=F],file=paste(pathRun,"Wy",sep="")) # Wy
		plotVectorDistribution(posterior$WxRSS[Vx.ind,],file=paste(pathRun,"WxRSS",sep=""),box=FALSE) # Wx row square sums
		plotVectorDistribution(posterior$WyRSS[Vy.ind,],file=paste(pathRun,"WyRSS",sep=""),box=FALSE) # Wy row square sums
		if (sum(is.na(posterior$Ux))==0 & sum(is.nan(posterior$Ux))==0)
			plotMatrixDistribution(posterior$Ux[Vx.ind,,,drop=F],file=paste(pathRun,"Ux",sep="")) # Ux
		if (sum(is.na(posterior$Uy))==0 & sum(is.nan(posterior$Uy))==0)
			plotMatrixDistribution(posterior$Uy[Vy.ind,,,drop=F],file=paste(pathRun,"Uy",sep="")) # Uy
		plotMatrixDistribution(posterior$z,file=paste(pathRun,"z",sep="")) # z
		plotMatrixDistribution(posterior$xlat,file=paste(pathRun,"xlat",sep="")) # xlat
		plotMatrixDistribution(posterior$ylat,file=paste(pathRun,"ylat",sep="")) # ylat
		plotVectorDistribution(posterior$beta,file=paste(pathRun,"beta",sep=""),box=FALSE) # beta
		if (dataset==0) {
			plotSurface(Vx.avg[,Vx.ind],file=paste(pathRun,"Vx",sep=""),Xtrue=data$Vx) # Vx
			plotSurface(Vy.avg[,Vy.ind],file=paste(pathRun,"Vy",sep=""),Xtrue=data$Vy) # Vy
		} else {
			plotSurface(Vx.avg[,Vx.ind],file=paste(pathRun,"Vx",sep="")) # Vx
			plotSurface(Vy.avg[,Vy.ind],file=paste(pathRun,"Vy",sep="")) # Vy
		}
		# MANOVA part
		if (sampleEffects[1])
	# 			plotVectorDistribution(posterior$mu_c,file=paste(pathRun,"mu_c",sep="")) # mu_c
			plotVectorDistributionWithLabels(posterior$mu_c,file=paste(pathRun,"effect_a",sep=""),effN=1) # mu_c
		if (sampleEffects[2])
	# 			plotVectorDistribution(posterior$mu_g,file=paste(pathRun,"mu_g",sep="")) # mu_g
			plotVectorDistributionWithLabels(posterior$mu_g,file=paste(pathRun,"effect_b",sep=""),effN=2) # mu_g
		if (sampleEffects[3])
	# 			plotVectorDistribution(posterior$mu_cg,file=paste(pathRun,"mu_cg",sep="")) # mu_cg
			plotVectorDistributionWithLabels(posterior$mu_cg,file=paste(pathRun,"effect_ab",sep=""),effN=12) # mu_cg
		if (sampleEffects[4])
			plotVectorDistribution(posterior$mu_s,file=paste(pathRun,"effect_c",sep="")) # mu_s
		if (sampleEffects[5])
			plotVectorDistribution(posterior$mu_sc,file=paste(pathRun,"effect_ac",sep="")) # mu_sc
		if (sampleEffects[6])
			plotVectorDistribution(posterior$mu_sg,file=paste(pathRun,"effect_bc",sep="")) # mu_sg
		if (sampleEffects[7])
			plotVectorDistribution(posterior$mu_scg,file=paste(pathRun,"effect_abc",sep="")) # mu_scg

		plot_var(data$X,data$case,data$gender,data$state,posterior$varsX,posterior$SigmaX,posterior$muX,paste(pathRun,"varsX",sep=""))
		plot_var(data$Y,data$case,data$gender,data$state,posterior$varsY,posterior$SigmaY,posterior$muY,paste(pathRun,"varsY",sep=""))

		# Correlation series
		plot_cor(data$X,data$case,data$gender,posterior$Vx,Vxsum,NULL,paste(pathRun,"X-",sep=""))
		plot_cor(data$Y,data$case,data$gender,posterior$Vy,Vysum,NULL,paste(pathRun,"Y-",sep=""))
	} else {
		print("No plotting performed!")
	}

	## End the script
	print(runId)
	print(paste("Started",time.start))
	print(paste("Complete at",Sys.time()))

	return(posterior)

}