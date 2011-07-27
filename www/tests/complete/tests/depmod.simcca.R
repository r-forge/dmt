# Wx = Wy
# W free (no regularization)
# mostly OK, could be improved with full covariances

library(dmt)

#fs <- list.files("~/local/Rpackages/dmt/SVN/dmt/pkg/R", full.names = TRUE)
#lapply(fs, source)

priors <- list( Nm.wx.wy.sigma = 0 )

N <- 100
yDim <- xDim <- 8 # simcca

###############################################################

# For each marginal covariance structure, 
# test model performance

cors.list <- list() 
zdims <- seq(1, min(xDim, yDim), 2) # test with various latent variable dimensionalities

# TODO: other marginal covariance structures for Wx = Wy; W>=0
#for (marginalCovariances in c("full")) {
#for (marginalCovariances in c("isotropic")) {
#for (marginalCovariances in c("diagonal")) {
#for (marginalCovariances in c("identical isotropic")) {
for (marginalCovariances in c("isotropic", "diagonal", "identical isotropic", "full")) {
#for (marginalCovariances in c("isotropic", "full")) {

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

#pdf("pic/depmod.simcca.pdf")
#plot(as.vector(covX.estimated), as.vector(covX.true)); abline(0,1)
#dev.off()
