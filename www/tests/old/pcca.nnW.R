# Wx != Wy
# W >= 0

library(dmt)

priors <- list(W = 1e-3)
#priors <- list(W = 0.1, Nm.wx.wy.sigma = NULL)
marginalCovariances = "full"

####################################################

xdim <- 6
ydim <- 4
N <- 1e3
cors <- c()
#for (zdim in seq(1, min(xdim, ydim), 1)) {
#for (zdim in 1:2) {
 
 zdim <- 1
 
  print(zdim)

  print("Generating toy data")
  toy <- generate.toydata(N = N, zDim = zdim, xDim = xdim, yDim = ydim, marginal.covariances = marginalCovariances, 
  priors = priors)

  print("Fitting the model")
  res <- fit.dependency.model(toy$X, toy$Y, zDimension = zdim,
	            marginalCovariances = marginalCovariances,
		    priors = priors, matched = FALSE, verbose = TRUE)
  
  vec <- compare.estimate.and.truth(res, toy) 
  cors <- rbind(cors, vec)

#}

#colnames(cors) <- c("wx", "wy", "phix", "phiy")
print(cors)

