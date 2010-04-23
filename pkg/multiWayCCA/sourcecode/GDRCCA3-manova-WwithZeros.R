
#
# Bayesian CCA, (C) Arto Klami, arto.klami@tkk.fi
# Implements the Gibbs sampler for CCA as described in
# Klami, Kaski: Local dependent components, ICML 2007.
# http://www.machinelearning.org/proceedings/icml2007/papers/278.pdf
#

# =================================================================
#
# Some general functions needed for the sampler
#
# =================================================================


eye <- function(N) {
  diag(array(1,N))
}

invSqrt <- function(mat) {
  x <- as.matrix(mat)

  fac <- svd(x)
  mat_res <- (fac$v %*% sqrt(diag(1/fac$d,nrow=length(fac$d))) %*% t(fac$u))
  mat_res
}

sqrtMat <- function(mat) {
  x <- as.matrix(mat)

  fac <- svd(x)
  # 28.4.09 - Parameter "nrow" passed to function 'diag' to enable latent variable dimensionality of one.
  mat_res <- (fac$v %*% sqrt(diag(fac$d,nrow=length(fac$d))) %*% t(fac$u))
  mat_res
}

debug <- FALSE

# =================================================================
#
# Functions for sampling all the necessary parameters with Gibbs.
#
# =================================================================

#
# SampleZ
# - if noSample is true, then z is not sampled but instead the
#   function returns the parameters of the conditional posterior
#   distribution
#

sampleZ <- function(x,W,z,Psi,beta,noSample=FALSE) {
  if(debug) print("Entering sampleZ")

  b <- t(W) %*% solve(crossprod(t(W)) + Psi)
  dz <- ncol(W)
  mu <- b %*% x
  cv <- eye(dz) - b %*% W

  randsample <- mvrnorm(dim(x)[2],array(0,c(1,dz)),cv)

  # Transpose the result for one-dimensional latent space
  if(ncol(mu)==1)
    randsample <- t(randsample)

  z <- mu + t(randsample)

  like <- sum(dmvnorm(randsample,array(0,c(1,dz)),cv,log=TRUE))

  if(!noSample)
    list(z=z, like=like)
  else
    list(mu=mu, cov=cv)
}

#
# SamplePsi
#

samplePsi <- function(x,W,z) {
  dfPrior <- nrow(W)+1

  if(psi_prior_type==1) {
    Sprior <- psi_prior*eye(nrow(W))
  } else {
    Sprior <- array(psi_prior,c(nrow(W),nrow(W))) + psi_prior/2*eye(nrow(W))
  }
  df <- dfPrior + psi_prior + dim(x)[2]

  xtemp <- x - W %*% z
  S <- Sprior + crossprod(t(xtemp))

  Psi <- rinvwishart(df,S)
  like <- dinvwishart(Psi,df,S)

  list(Psi=Psi, like=like)
}

#
# SampleW
# - samples each column separately; the new ones are immediately
#   used when sampling the following ones
#
# Argument 'nXlat' added. -28.5.09

sampleW <- function(x,W,z,Psi,beta,nXlat=NULL,newW=NULL) {
  if(debug) print("Entering sampleW")

  iPsi <- solve(Psi)

  xtemp <- x

  like <- 0

  for (i in 1:ncol(W)) {
    if(beta[i]>0) {
      Wrest <- as.matrix(W[,-i,drop=FALSE])

      zrest <- as.matrix(z[-i,,drop=FALSE])
      temp1 <- z[i,,drop=TRUE]
      temp2 <- xtemp - Wrest %*% zrest
 
      a <- apply(temp2 * (rep(1,nrow(W)) %o% temp1),1,sum)
      b <- sum(temp1^2) * iPsi + 1/beta[i] * eye(nrow(W))
      cv <- solve(b)

      mu <- iPsi %*% cv %*% a

      W[,i] <- mvrnorm(1,mu,cv)

		# Modify W for dataset-specific latent components. -28.5.09
		if (setWtoZero[1,i] & !is.null(nXlat))
			W[1:nXlat,i] = 0
		if (setWtoZero[2,i] & !is.null(nXlat))
			W[(nXlat+1):nrow(W),i] = 0

      like <- like + dmvnorm(W[,i],mu,cv,log=TRUE)
    } else { # beta[i] was zero
      W[,i] <- array(0,c(1,nrow(W)))
    }
  }

  list(W=W, like=like)
}

#
# SampleBeta
# - ARD prior parameter for W
#

sampleBeta <- function(x,W,beta) {
  if(debug) print("Entering sampleBeta")

  a_prior <- betaFirst
  b_prior <- betaPrior

  a_list <- array(0,ncol(W))
  b_list <- array(0,ncol(W))

  like <- 0
  for (i in 1:ncol(W)) {
    a <- nrow(W)/2 + a_prior
    b <- 0.5*sum(W[,i]^2) + b_prior

    beta[i] <- 1 / rgamma(1,a,b)

    like <- like + dinvgamma(beta[i],a,b,log=TRUE)
  }

  list(beta=beta, like=like)
}

#
# sampleMu
#

sampleMu <- function(x,model) {
  prior <- mu_prior

  if(debug) print("SampleMu")

  PsiX <- model$PsiX
  PsiY <- model$PsiY
  Psi <- rbind(cbind(PsiX,array(0,c(nrow(PsiX),ncol(PsiY)))),
         cbind(array(0,c(nrow(PsiY),ncol(PsiX))),PsiY))
  Psi <- Psi + crossprod(t(model$W))
  iPsi <- solve(Psi)

  if(length(mu_vec)==0) {
    me <- array(0,c(nrow(x),1))
  } else {
    me <- mu_vec
  }
  for (j in 1:ncol(x)) {
    me <- me + x[,j]
  }

  cv <- solve(ncol(x) * iPsi + prior*eye(ncol(Psi)))

  me <- cv %*% (iPsi %*% me)

  cen <- mvrnorm(1,me,cv)

  like <- dmvnorm(cen,me,cv,log=TRUE)

  list(mu=cen, like=like)
}

#
# updateModel
# - runs one iteration of Gibbs for the model
# 'data' consists of matrices 'X' and 'Y', and of vectors 'case', 'gender' and 'state'.
#

updateModel <- function(data,model,priors) {
	
	## Update covariate effects.
	
	if (sampleEffects[1])
		model$mu_c = sampleMu_c(model$z,data$case,data$gender,data$state,model$mu_g,model$mu_s,model$mu_cg,model$mu_sc,model$mu_sg,model$mu_scg)
	if (sampleEffects[2])
		model$mu_g = sampleMu_g(model$z,data$case,data$gender,data$state,model$mu_c,model$mu_s,model$mu_cg,model$mu_sc,model$mu_sg,model$mu_scg)
	if (sampleEffects[3])
		model$mu_cg = sampleMu_cg(model$z,data$case,data$gender,data$state,model$mu_c,model$mu_g,model$mu_s,model$mu_sc,model$mu_sg,model$mu_scg)
	if (sampleEffects[4])
		model$mu_s = sampleMu_s(model$z,data$case,data$gender,data$state,model$mu_c,model$mu_g,model$mu_cg,model$mu_sc,model$mu_sg,model$mu_scg)
	if (sampleEffects[5])
		model$mu_sc = sampleMu_sc(model$z,data$case,data$gender,data$state,model$mu_c,model$mu_g,model$mu_s,model$mu_cg,model$mu_sg,model$mu_scg)
	if (sampleEffects[6])
		model$mu_sg = sampleMu_sg(model$z,data$case,data$gender,data$state,model$mu_c,model$mu_g,model$mu_s,model$mu_cg,model$mu_sc,model$mu_scg)
	if (sampleEffects[7])
		model$mu_scg = sampleMu_scg(model$z,data$case,data$gender,data$state,model$mu_c,model$mu_g,model$mu_s,model$mu_cg,model$mu_sc,model$mu_sg)
	
	## Update CCA-related variables.
	
	vlat <- rbind(model$xlat,model$ylat)

	#
	# Update the mean of the cluster
	#
	if (sampling$Mu) {
		temp <- sampleMu(vlat,model)
		model$mu <- temp$mu
		model$like$mu <- temp$like
	}

	#
	# Create temporary data sets with mean removed; everything else
	# is hence conditioned on current values of mu
	#
	xlattemp <- model$xlat - rep(model$mu[1:nrow(model$xlat)],ncol(model$xlat))
	ylattemp <- model$ylat - rep(model$mu[(nrow(model$xlat)+1):nrow(vlat)],ncol(model$ylat))
	vlattemp <- vlat - rep(model$mu,ncol(vlat))

	#
	# Update the marginal model covariances (Psi-parameters)
	#
	Wx <- model$W[1:nrow(model$xlat),,drop=FALSE]
	temp <- samplePsi(xlattemp,Wx,model$z)
	model$PsiX <- temp$Psi
	model$like$PsiX <- temp$like

	Wy <- model$W[(nrow(model$xlat)+1):nrow(vlat),,drop=FALSE]
	temp <- samplePsi(ylattemp,Wy,model$z)
	model$PsiY <- temp$Psi
	model$like$PsiY <- temp$like

	#
	# Update the priors for W (beta).
	#
	temp <- sampleBeta(vlattemp,model$W,model$beta)
	model$beta <- temp$beta
	model$like$Beta <- temp$like

	#
	# Update the projection matrix W
	#
	PsiX <- model$PsiX
	PsiY <- model$PsiY
	Psi <- rbind(cbind(PsiX,array(0,c(nrow(PsiX),ncol(PsiY)))),
				cbind(array(0,c(nrow(PsiY),ncol(PsiX))),PsiY))
	# Some of cells of W may be set to zero in the function. -28.5.09
	temp <- sampleW(vlattemp,model$W,model$z,Psi,model$beta,nrow(model$xlat))
	model$W <- temp$W
	model$like$W <- temp$like

	#
	# Sample z
	#
# 	temp <- sampleZ(vlattemp,model$W,model$z,Psi,model$beta)
# 	model$z <- temp$z
	# 30.4.09 - Use new sampler on 'z'. This function does not have likelihood calculation yet.
	tmp = sampleXlats2(vlat,data$case,data$gender,data$state,model$mu,model$mu_c,model$mu_g,model$mu_s,model$mu_cg,model$mu_sc,model$mu_sg,model$mu_scg,Psi,model$W,NULL)
	model$z = tmp$z
	model$like$z <- tmp$like # likelihood computation -24.3.10

	# Arrange the columns of W, Bx and By in a nice order; does not
	# affect the model in any way.
	# 27.4.09 - ordering commented out
# 	ord <- order(-apply(model$W^2,2,sum))
# 	model$beta <- model$beta[ord]
# 	model$W <- model$W[,ord,drop=FALSE]
# 	model$z <- model$z[ord,,drop=FALSE]

	## Update FA-related variables.
	
	# Update the dataset-specific clusterings 'Vx' and 'Vy'.
	model$Vx = sampleV(data$X,model$xlat,model$muX,model$SigmaX,model$thetaX,ccasta=FALSE,model$varsX)
	model$Vy = sampleV(data$Y,model$ylat,model$muY,model$SigmaY,model$thetaY,ccasta=FALSE,model$varsY)
	
	# Update the dataset-specific residual covariances 'SigmaX' and 'SigmaY'.
	model$SigmaX = sampleSigma(data$X,model$muX,model$Vx,model$xlat,model$varsX)
	model$SigmaY = sampleSigma(data$Y,model$muY,model$Vy,model$ylat,model$varsY)
	
	# Update the dataset-specific latent variables 'xlat' and 'ylat'.
	# Effects still need to be added here.
	Wx <- model$W[1:nrow(model$xlat),,drop=FALSE]
	Elat = Wx%*%model$z # expectation value vector of the latent variables
	tmp <- sampleXlatsWithPriorDR(data$X,model$muX,model$SigmaX,model$Vx,model$xlat,Elat,model$PsiX,NULL,model$varsX)
	model$xlat = tmp$xLat
	model$like$xlat = tmp$like
	Wy <- model$W[(nrow(model$xlat)+1):nrow(vlat),,drop=FALSE]
	Elat = Wy%*%model$z # expectation value vector of the latent variables
	tmp <- sampleXlatsWithPriorDR(data$Y,model$muY,model$SigmaY,model$Vy,model$ylat,Elat,model$PsiY,NULL,model$varsY)
	model$ylat = tmp$xLat
	model$like$ylat = tmp$like

	# Update the dataset-specific scale parameters 'varsX' and 'varsY'.
	if (sampling$scale) {
		model$varsX = sampleVars4(data$X,model$muX,model$Vx,model$xlat,model$SigmaX,priors$varsX0,priors$N0)
		model$varsY = sampleVars4(data$Y,model$muY,model$Vy,model$ylat,model$SigmaY,priors$varsY0,priors$N0)
	}
	
	# Update the dataset-specific mean variables.
	if (sampling$MuX) { # Parameter decides whether metabolite-specific mu is sampled.
		model$muX = sampleMuDR(data$X,model$SigmaX,model$Vx,model$varsX,model$xlat,priors$muX0,priors$varsX0)
		model$muY = sampleMuDR(data$Y,model$SigmaY,model$Vy,model$varsY,model$ylat,priors$muY0,priors$varsY0)
	}

	## Compute the data likelihood -24.3.10
	
	model$like$X = 0
	Vx.tmp = matrix((model$varsX),nrow(data$X),nrow(model$xlat))*model$Vx
	for (n in 1:ncol(data$X)) {
		model$like$X = model$like$X + sum(dnorm(x=data$X[,n], mean=model$muX+Vx.tmp%*%model$xlat[,n], sd=sqrt(model$SigmaX), log=T))
	}
	model$like$Y = 0
	Vy.tmp = matrix((model$varsY),nrow(data$Y),nrow(model$ylat))*model$Vy
	for (n in 1:ncol(data$Y)) {
		model$like$Y = model$like$Y + sum(dnorm(x=data$Y[,n], mean=model$muY+Vy.tmp%*%model$ylat[,n], sd=sqrt(model$SigmaY), log=T))
	}
	
	return(model)
	
}

# =================================================================
# 
# Functions to be called for sampling a model
# - first initialize a few global parameters with setParameters()
# - then run the sampler with:
#   result <- sampleGibbs(x,y,nIter)
#   model <- result$model
#
# =================================================================

#
# A prototype of a function that sets all the global variables.
# Replace this with a version that sets the desired values and
# call once before creating a model.
#

setParameters <- function() {
  betaFirst <<- 1e-1
  betaPrior <<- 1e-1

  sigFirst <<- 1e-1
  sigPrior <<- 1e-1

  # Prior type 1 favors diagonal covariance. Prior type 2 favors
  # full covariance.
  psi_prior_type <<- 1
  psi_prior <<- 1

  # Mean vector, NULL means that it is a vector of zeroes
  mu_vec <<- NULL
  mu_prior <<- 1
}

#
# sampleMPCCA
# - The main function for sampling from the model
# - Input:
#   - x,y: matrices where each column is a sample
#   - nIter: number of posterior samples (optional)
#   - zDim: dimensionality of the latent space (optional)
# - Output:
#   - model: list describing the model, contains the following elements:
#     - xDim, yDim: data dimensionalities
#     - zDim: dimensionality of the latent space
#     - mu: mean of data
#     - W: projection matrix; first xDim rows are for x-space, the rest
#       for y-space
#     - beta: ARD-variances for columns of W
#     - z: latent variables
#     - PsiX, psiY: covariance matrices
#     - like: current likelihood
#   - storage: list of earlier models, one for each posterior
#     sample. Otherwise identical with the final model but do not
#     contain the latent variables
#

#sampleMPCCA <- function(X,Y,nIter=100,zDim=0) {
# Also give dimensions of latent variables.
sampleMPCCA <- function(data,priors,nIter=100,xlatDim=4,ylatDim=5,zDim=0) {
	#PsiXDim <- nrow(xlat)
	PsiXDim = xlatDim
	#PsiYDim <- nrow(ylat)
	PsiYDim <- ylatDim
	if(zDim==0) {
#     zDim <- min(nrow(xlat),nrow(ylat))
		zDim = min(xlatDim,ylatDim)
	}

	# Store the dimensionalities of the data
# 	xlatDim <- nrow(xlat)
# 	ylatDim <- nrow(ylat)

	# Initialize the variables of the dimension reduction part.
	modelDR = sampleFromPriorDR(data,priors,xlatDim,ylatDim,N)
	# Combine the dataset-specific latent variables into one array.
	vlat <- rbind(modelDR$xlat,modelDR$ylat)

	# Temporary model for prior sampling
	model <- list(xDim=xlatDim, yDim=ylatDim, zDim=zDim, like=NULL)
	model <- sampleFromPrior(case=data$case, gender=data$gender, state=data$state, x=vlat, N=ncol(vlat), model=model)

	# Include the variables of the dimension reduction part.
	model = c(model,modelDR)
	rm(modelDR)

	like <- NULL
	storage <- array(list(),nIter)

	if(nIter>0) {
		for (i in 1:nIter) {
			model <- updateModel(data,model,priors)
			
			storage[[i]] <- collectData(storage,model)
		}
	}
	
	list(model=model, storage=storage)
}

#
# sampleMPCCAcont:
# - continue sampling from a given model
# - enables checking convergence between batches of iterations
#

sampleMPCCAcont <- function(data,model,priors,nIter) {
  storage <- array(list(),nIter)

  if(nIter>0) {
    for (i in 1:nIter) {
      model <- updateModel(data,model,priors)

      storage[[i]] <- collectData(storage,model)
    }
  }

  list(model=model, storage=storage)
}

#
# sampleFromPrior
# - creates a new model by sampling from the prior
# - not supposed to be called directly
#

sampleFromPrior <- function(case, gender, state, x, N, model) {
	# Do not allow very large beta-values that would result in
	# near zero variance; we don't want to have components turned
	# of in the initialization
	beta <- 1 / rgamma(model$zDim,betaFirst,betaPrior)
	beta[which(beta>10)] <- 10
	beta <- beta[order(-beta)]

	z <- array(rnorm(model$zDim*N),c(model$zDim,N))
	W <- mvrnorm(nrow(x),array(0,c(1,model$zDim)),diag(beta,nrow=model$zDim))

	dfPrior <- model$xDim+1
	if(psi_prior_type==1) {
		Sprior <- psi_prior*eye(model$xDim)
	} else {
		Sprior <- array(psi_prior,c(model$xDim,model$xDim)) + psi_prior/2*eye(model$xDim)
	}
	df <- dfPrior + psi_prior

	PsiX <- rinvwishart(df,Sprior)

	dfPrior <- model$yDim+1
	if(psi_prior_type==1) {
		Sprior <- psi_prior*eye(model$yDim)
	} else {
		Sprior <- array(psi_prior,c(model$yDim,model$yDim)) + psi_prior/2*eye(model$yDim)
	}
	df <- dfPrior + psi_prior
	PsiY <- rinvwishart(dfPrior,Sprior)

	if(length(mu_vec)==0) {
		me <- array(0,c(nrow(x),1))
	} else {
		me <- mu_vec
	}
	mu <- mvrnorm(1,me,mu_prior*eye(nrow(x)))
	if (!sampling$Mu) { # If no sampling of mu; set it to zero. -090504
		mu = mu*0
	}

	## Treatment effects

	if (sum(case>0)!=0)
		mu_c = rnorm(n=model$zDim)
	else
		mu_c = rep(0,model$zDim)
	if (sum(gender>0)!=0)
		mu_g = rnorm(n=model$zDim)
	else
		mu_g = rep(0,model$zDim)
	if (sum(case>0)!=0 & sum(gender>0)!=0)
		mu_cg = rep(0,model$zDim)
	else
		mu_cg = rep(0,model$zDim)
	
	if (sum(state==0)==0) { # state information included
		S = max(state)
		mu_s = matrix(rnorm(n=(model$zDim*S)),nrow=model$zDim,ncol=S)
		if (sum(case>0)!=0)
			mu_sc = matrix(rnorm(n=(model$zDim*S)),nrow=model$zDim,ncol=S)
		else
			mu_sc = matrix(0,nrow=model$zDim,ncol=S)
		if (sum(gender>0)!=0)
			mu_sg = matrix(0,nrow=model$zDim,ncol=S)
		else
			mu_sg = matrix(0,nrow=model$zDim,ncol=S)
		if (sum(case>0)!=0 & sum(gender>0)!=0)
			mu_scg = matrix(0,nrow=model$zDim,ncol=S)
		else
			mu_scg = matrix(0,nrow=model$zDim,ncol=S)
		rm(S)
	} else {
		mu_s = matrix(0,nrow=model$zDim,ncol=1)
		mu_sc = matrix(0,nrow=model$zDim,ncol=1)
		mu_sg = matrix(0,nrow=model$zDim,ncol=1)
		mu_scg = matrix(0,nrow=model$zDim,ncol=1)
	}

	##

	model$W <- W
	model$z <- z
	model$beta <- beta
	model$PsiX <- PsiX
	model$PsiY <- PsiY
	model$mu <- mu

	model$mu_c = mu_c
	model$mu_g = mu_g
	model$mu_cg = mu_cg
	model$mu_s = mu_s
	model$mu_sc = mu_sc
	model$mu_sg = mu_sg
	model$mu_scg = mu_scg

	model
}

# =================================================================
#
# Functions for extracting information from the posterior samples
#
# =================================================================

#
# collectData
# - Store the current state of parameters, leaving the latent
#   variables out as they take quite a bit of space
#

collectData <- function(list, model) {
  temp <- model

  # Do not store the latent variables that take a lot of space
  #temp$z <- NULL # Do store them. -20.4.09

  temp
}

#
# Load/save a set of posterior samples
#

saveSamples <- function(filename, samples) {
  sampleList <- samples

  save(sampleList, file=filename, compress=TRUE)
}

loadSamples <- function(filename) {
  load(filename)

  samples <- sampleList

  samples
}

#
# Concatenates two lists so that part of the first is dropped.
# This creates a process that always retains a given proportion
# of samples in the result list.
# Used just to save space; we are anyway going to ignore some
# samples as burn-in period, so do not store them either
#

concatenateLists <- function(list1, list2, dropN=0) {
  l1 <- length(list1)
  l2 <- length(list2)

  dropN <- round(dropN)

  newlist <- array(list(),l1-dropN+l2)

  if(l1>dropN)
    for(i in 1:(l1-dropN)) {
      newlist[[i]] <- list1[[dropN + i]]
    }
  for(i in 1:l2) {
    newlist[[l1-dropN+i]] <- list2[[i]]
  }

  newlist
}

stripList <- function(list) {
  l <- length(list)
  newList <- array(list(),round(l/2))

  j <- 1
  for(i in (round(l/2)+1):l) {
    newList[[j]] <- list[[i]]
    j <- j+1
  }

  newList
}

#
# Go though the posterior samples and collect some parameters
# as a vector. Used for fetching data for convergence check.
#

concludeData <- function(list, mode, skip=1) {
  l <- length(list)

  vec <- array(0,l/skip)

  j <- 1
  for (i in seq(1,l,skip)) {
    if(mode==1) {           # likelihood
      vec[j] <- list[[i]]$like
    } else if (mode==2) {   # largest beta
      vec[j] <- list[[i]]$beta[1]
    } else if (mode==3) {   # magnitude of the first column of W
      vec[j] <- sum(list[[i]]$W[,1]^2)
    } else if (mode==4) {   # determinant of PsiX
      vec[j] <- determinant(list[[i]]$PsiX,log=TRUE)$modulus[1]
    } else if (mode==5) {   # determinant of PsiY
      vec[j] <- determinant(list[[i]]$PsiY,log=TRUE)$modulus[1]
    } else if (mode==6) {   # mu
      vec[j] <- sum(list[[i]]$mu^2)
    }
    j <- j + 1
  }

  vec
}

#
# Check whether the sampler has converged, using the
# potential scale reduction factor by Gelman et al.
# - see Rhat.R and itsim.R for more details
#

testConvergence <- function(lists, modes=1:6) {
  mat <- array(0,c(length(lists[[1]]),length(lists),length(modes)))

  for(m in 1:length(modes)) {
    for(i in 1:length(lists)) {
      vec <- concludeData(lists[[i]],modes[m])
      mat[,i,m] <- vec
    }

    val <- Rhat(mat,keep.all=TRUE)
  }

  list(val=val, mat=mat)
}

#
# projectData:
# Find the "CCA-projection" of the data, meaning the expected
# values of the latent variables z given the data samples.
# Returns E[z|x,y], E[z|x] and E[z|y]
#

projectData <- function(x,y,model) {
  v <- rbind(x,y)
  vtemp <- v - rep(model$mu,ncol(v))
 
  PsiX <- model$PsiX
  PsiY <- model$PsiY
  Psi <- rbind(cbind(PsiX,array(0,c(nrow(PsiX),ncol(PsiY)))),
               cbind(array(0,c(nrow(PsiY),ncol(PsiX))),PsiY))
  temp <- sampleZ(vtemp,model$W,model$z,Psi,model$beta,noSample=TRUE)

  tempx <- sampleZ(vtemp[1:model$xDim,],model$W[1:model$xDim,],model$z,PsiX,model$beta,noSample=TRUE)
  tempy <- sampleZ(vtemp[(model$xDim+1):nrow(vtemp),],model$W[(model$xDim+1):nrow(model$W),],model$z,PsiY,model$beta,noSample=TRUE)

  list(xy=temp$mu,x=tempx$mu,y=tempy$mu)
}

#
# findRotation:
# - find the real CCA components by solving the rotational ambiquity
# - uses the method given by Archambeau et al. "Robust probabilistic
#    projections. The correct version is in the errata:
#    http://www.cs.ucl.ac.uk/staff/c.archambeau/publ/icml_ca06_errata.pdf
# - notation also follows that paper
#

findRotation <- function(x,y,model) {
  dx <- nrow(x)
  dy <- nrow(y)

  Wx <- model$W[1:dx,,drop=FALSE]
  Wy <- model$W[(dx+1):(dx+dy),,drop=FALSE]
  Bx <- t(Wx) %*% solve(model$PsiX) %*% Wx + eye(ncol(Wx))
  By <- t(Wy) %*% solve(model$PsiY) %*% Wy + eye(ncol(Wy))

  Rt1 <- sqrtMat(eye(ncol(Wx))-solve(Bx)) %*% (eye(ncol(Wy))-solve(By)) %*% sqrtMat(eye(ncol(Wx))-solve(Bx))
  Rt2 <- sqrtMat(eye(ncol(Wy))-solve(By)) %*% (eye(ncol(Wx))-solve(Bx)) %*% sqrtMat(eye(ncol(Wy))-solve(By))

  Rte1 <- eigen(Rt1)
  Rte2 <- eigen(Rt2)
  R1 <- Rte1$vectors
  R2 <- Rte2$vectors

  covX <- crossprod(t(Wx)) + model$PsiX
  covY <- crossprod(t(Wy)) + model$PsiY

  isqX <- invSqrt(eye(ncol(Wx))-solve(Bx))
  isqY <- invSqrt(eye(ncol(Wy))-solve(By))

  Ux <- solve(covX) %*% Wx %*% isqX %*% R1
  Uy <- solve(covY) %*% Wy %*% isqY %*% R2

  Vx <- invSqrt(covX) %*% Wx %*% isqX %*% R1
  Vy <- invSqrt(covY) %*% Wy %*% isqY %*% R2

  list(Ux = Ux, Uy = Uy, Vx = Vx, Vy = Vy, cor = sqrt(Rte1$values))
}


