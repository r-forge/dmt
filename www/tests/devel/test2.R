# Wx = Wy
# W free (no regularization)
# mostly OK, could be improved with full covariances

library(dmt)

#source("~/local/Rpackages/dmt/SVN/dmt/pkg/R/generate.toydata.R")
#source("~/local/Rpackages/dmt/SVN/dmt/pkg/R/fit.dependency.model.R")  
#source("~/local/Rpackages/dmt/SVN/dmt/pkg/R/internals.R")            
#source("~/local/Rpackages/dmt/SVN/dmt/pkg/R/initialize2.R")   
#source("~/local/Rpackages/dmt/SVN/dmt/pkg/R/costs.R")    
#source("~/local/Rpackages/dmt/SVN/dmt/pkg/R/get.W.R")    
#source("~/local/Rpackages/dmt/SVN/dmt/pkg/R/phi.updates.R")   
#source("~/local/Rpackages/dmt/SVN/dmt/pkg/R/dependency.score.R")  
#source("~/local/Rpackages/dmt/SVN/dmt/pkg/R/optimize.parameters.R")  
#source("~/local/Rpackages/dmt/SVN/dmt/pkg/R/pfa.R")  
#source("~/local/Rpackages/dmt/SVN/dmt/pkg/R/ppca.R")  
#source("~/local/Rpackages/dmt/SVN/dmt/pkg/R/pcca.R")  
#source("~/local/Rpackages/dmt/SVN/dmt/pkg/R/pcca.with.isotropic.margins.R")  
#source("~/local/Rpackages/dmt/SVN/dmt/pkg/R/M.set.R")  
#source("~/local/Rpackages/dmt/SVN/dmt/pkg/R/set.beta.R")  
#source("~/local/Rpackages/dmt/SVN/dmt/pkg/R/W.cca.EM.R")  

priors <- list(Nm.wxwy.mean = NA,  Nm.wx.wy.sigma = 0 )
N <- 50
yDim <- xDim <- 8 # simcca

###############################################################

# For each marginal covariance structure, 
# test model performance

cors.list <- list() 
zdims <- seq(1, min(xDim, yDim), 2) # test with various latent variable dimensionalities

# TODO: other marginal covariance structures for Wx = Wy; W>=0
#for (marginalCovariances in c("full")) {
for (marginalCovariances in c("isotropic", "diagonal", "identical isotropic", "full")) {

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
