
# W update
#   - free
#   - positive
#   - Wx ~ Wy
#   - Wx = Wy


# phi update
#     - full
#     - diagonal
#     - isotropic

########################################

library(dmt)

#source("../R/pcca.R")
#source("../R/internals.R")
#source("../R/simCCA.optimize3.R")
#source("../R/optimize.parameters.R")
#source("../R/costs.R")
#source("../R/get.W.R")
#source("../R/M.set.R")
#source("../R/phi.updates.R")
#source("../R/dependency.score.R")

#####################################

# Two data sets, different dimensions
priors <- NULL
xdim <- 8
ydim <- 5
cors <- c()
#for (zdim in seq(1, min(xdim, ydim), 1)) {
zDimension <- zdim <- 2

  toy <- generate.toydata(N = 100, zDim = zdim, xDim = xdim, yDim = ydim, marginal.covariances = "full")

  covX.true <- toy$Wx%*%t(toy$Wx)
  covY.true <- toy$Wy%*%t(toy$Wy)  
  phiX.true <- toy$Bx%*%t(toy$Bx)
  phiY.true <- toy$By%*%t(toy$By)

  ##############################
  
  # Initialize  
  # samples are always matched i.e. ncol(X) = ncol(Y)
  Nsamples <- ncol(X)	
  if ( length(priors) == 0 ) { priors <- list() }
  if ( is.null(priors$Nm.wxwy.sigma) ) { priors$Nm.wxwy.sigma <- Inf } # tune similarity constraint Wx ~ Wy

  Dim <- list()
  Dim$X <- nrow(X)
  Dim$Y <- nrow(Y)
  Dim$Z <- zDimension
	
  Dcov <- list()
  Dcov$X <- cov(t(X))
  Dcov$Y <- cov(t(Y))
  Dcov$total <- cov(t(rbind(X, Y)))
		
  # initialize with scalar diagonal noise on the marginals (shared by all features)
  phi.init <- list(X = diag(var(as.vector(X)), Dim$X), Y = diag(var(as.vector(Y)), Dim$Y)) 

  ############################################  

  nullmat  <- matrix(0, nrow = Dim$X, ncol = Dim$Y)
 

  # Update phi (full covariances)		        			
  phi <- phi.init
  phi.inv <- list()
  phi.inv$X <- solve(phi$X)
  phi.inv$Y <- solve(phi$Y)    
  phi.inv$total <- rbind(cbind(phi.inv$X, nullmat), cbind(t(nullmat), phi.inv$Y))    
  W.old <- W.true <- list(X = toy$Wx, Y = toy$Wy, total = rbind(toy$Wx, toy$Wy))
  M <- set.M.full2(W.old, phi.inv, dz = Dim$Z) 
  phi <- update.phi.EM.fullcov(Dcov, W.true, phi.inv, W.old, M, nullmat)

  ##############################


cor(as.vector(phi$X), as.vector(phiX.true))
cor(as.vector(phi$Y), as.vector(phiY.true))


#############  
#  # Initialize W's
#  W.init   <- list()
#  W.init$X <- as.matrix(eigen(Dcov$X)$vectors[, 1:Dim$Z])
#  W.init$Y <- as.matrix(eigen(Dcov$Y)$vectors[, 1:Dim$Z])
#  W.init$total <- rbind(W.init$X, W.init$Y) # can be possibly removedin some special cases#
#
#  # This optimizes Wx and Wy assuming they are independent
#  opt <- optim(c(as.vector(W$X), as.vector(W$Y)), cost.W.exponential, method = "L-BFGS-B", phi = phi, priors = priors, Dim = Dim, Dcov = Dcov, control = list(maxit = 1e6), lower = -10*max(Dcov$total), upper = 10*max(Dcov$total))
#  W <- get.W2(opt$par, Dim)		#

  ##################################################

