
# Wx = Wy
# W >= 0
# NOTE: only implemented with marginalCovariances = "full"
# TODO: isotropic cases need considerable speedup


library(dmt)

#source("~/local/Rpackages/dmt/SVN/dmt/pkg/R/fit.dependency.model.R")
#source("~/local/Rpackages/dmt/SVN/dmt/pkg/R/internals.R")
#source("~/local/Rpackages/dmt/SVN/dmt/pkg/R/initialize2.R")
#source("~/local/Rpackages/dmt/SVN/dmt/pkg/R/costs.R")
#source("~/local/Rpackages/dmt/SVN/dmt/pkg/R/get.W.R")
#source("~/local/Rpackages/dmt/SVN/dmt/pkg/R/phi.updates.R")
#source("~/local/Rpackages/dmt/SVN/dmt/pkg/R/dependency.score.R")
#source("~/local/Rpackages/dmt/SVN/dmt/pkg/R/optimize.parameters.R")


priors <- list(W = 1e-3, Nm.wx.wy.sigma = 0) # W>=0; Wx = Wy

N <- 50           
zDim <- 1
xDim <- 8
yDim <- 8                                                        

###############################################################

# For each marginal covariance structure, 
# test model performance

cors.list <- list() 
zdims <- seq(1, min(xDim, yDim), 2) # test with various latent variable dimensionalities

for (marginalCovariances in c("isotropic", "identical isotropic", "diagonal", "full")) {
  # TODO: other marginal covariance structures for Wx = Wy; W>=0
  #for (marginalCovariances in c("full")) {

  cors <- c()
  for (zDim in zdims) {
 
    print(paste(marginalCovariances, " / zDim: ", zDim))

    toy <- generate.toydata(N = N, zDim = zDim, xDim = xDim, yDim = yDim, 
      	 		    marginal.covariances = marginalCovariances, 
      	 		    priors = priors)

    res <- fit.dependency.model(toy$X, toy$Y, zDimension = zDim,
	            marginalCovariances = marginalCovariances,
		    priors = priors, matched = TRUE, verbose = TRUE)
  
    vec <- compare.estimate.and.truth(res, toy) 

    cors <- rbind(cors, vec)

  }

  colnames(cors) <- c("wx", "wy", "phix", "phiy")
  rownames(cors) <- as.character(zdims)
  cors.list[[marginalCovariances]] <- cors

}


print(cors.list)













