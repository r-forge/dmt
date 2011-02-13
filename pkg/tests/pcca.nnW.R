# Wx != Wy
# W >= 0
# OK, sufficiently good

library(dmt)

priors <- list(W = 1e-3, Nm.wx.wy.sigma = Inf)

N <- 50           
zDim <- 1
xDim <- 8
yDim <- 6      

###############################################################

# For each marginal covariance structure, 
# test model performance

cors.list <- list() 
zdims <- seq(1, min(xDim, yDim), 2) # test with various latent variable dimensionalities

# TODO: other marginal covariance structures for Wx = Wy; W>=0
#for (marginalCovariances in c("full")) {
for (marginalCovariances in c("diagonal", "identical isotropic", "isotropic", "full")) {

  cors <- c()
  for (zDim in zdims) {
 
    print(paste(marginalCovariances, " / zDim: ", zDim))

    toy <- generate.toydata(N = N, zDim = zDim, xDim = xDim, yDim = yDim, 
      	 		    marginal.covariances = marginalCovariances, 
      	 		    priors = priors)

    res <- fit.dependency.model(toy$X, toy$Y, zDimension = zDim,
	            marginalCovariances = marginalCovariances,
		    priors = priors, matched = FALSE, verbose = TRUE)
  
    vec <- compare.estimate.and.truth(res, toy) 

    cors <- rbind(cors, vec)

  }

  colnames(cors) <- c("wx", "wy", "phix", "phiy")
  rownames(cors) <- as.character(zdims)
  cors.list[[marginalCovariances]] <- cors

}


print(cors.list)

