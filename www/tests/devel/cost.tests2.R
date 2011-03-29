
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
xdim <- 5
ydim <- 5
cors <- c()
#for (zdim in seq(1, min(xdim, ydim), 1)) {
zDimension <- zdim <- 2

  toy <- generate.toydata(N = 10000, zDim = zdim, xDim = xdim, yDim = ydim, marginal.covariances = "full")

  covX.true <- toy$Wx%*%t(toy$Wx)
  covY.true <- toy$Wy%*%t(toy$Wy)  
  phiX.true <- toy$Bx%*%t(toy$Bx)
  phiY.true <- toy$By%*%t(toy$By)

  ##############################

model <- fit.dependency.model(toy$X, toy$Y, zDimension,
	             marginalCovariances = "full",
		     covLimit = 1e-3,
		     priors = list(), 
		     matched = TRUE)

covX.estimated <- model@W$X%*%t(model@W$X)
covX.true <- toy$Wx%*%t(toy$Wx)   
print(cor(as.vector(covX.estimated), as.vector(covX.true)))

covY.estimated <- model@W$Y%*%t(model@W$Y)
covY.true <- toy$Wx%*%t(toy$Wx)   
print(cor(as.vector(covY.estimated), as.vector(covY.true)))

phiX.estimated <- model@phi$X%*%t(model@phi$X)
phiX.true <- toy$Bx%*%t(toy$Bx)   
print(cor(as.vector(phiX.estimated), as.vector(phiX.true)))

phiY.estimated <- model@phi$Y%*%t(model@phi$Y)
phiY.true <- toy$Bx%*%t(toy$Bx)   
print(cor(as.vector(phiY.estimated), as.vector(phiY.true)))

