# Wx = Wy
# W >= 0
# NOTE: only implemented with marginalCovariances = "full"
# TODO: isotropic cases need considerable speedup

set.seed(21)

library(dmt)

source("~/local/Rpackages/dmt/SVN/dmt/pkg/R/optimize.parameters.R")
source("~/local/Rpackages/dmt/SVN/dmt/pkg/R/utilities.R")
#source("~/local/Rpackages/dmt/SVN/dmt/pkg/R/fit.dependency.model.R")
source("~/local/Rpackages/dmt/SVN/dmt/pkg/R/internals.R")
#source("~/local/Rpackages/dmt/SVN/dmt/pkg/R/initialize2.R")
source("~/local/Rpackages/dmt/SVN/dmt/pkg/R/update.W.R")
source("~/local/Rpackages/dmt/SVN/dmt/pkg/R/get.W.R")
source("~/local/Rpackages/dmt/SVN/dmt/pkg/R/update.phi.R")
source("~/local/Rpackages/dmt/SVN/dmt/pkg/R/dependency.score.R")
source("~/local/Rpackages/dmt/SVN/dmt/pkg/R/set.M.R")


priors <- list(W = 1e-3, Nm.wx.wy.sigma = 0) # W>=0; Wx = Wy

N <- 50           
zDim <- 1
xDim <- 8
yDim <- 8                                                        

cors.list <- list() 
zdims <- seq(1, min(xDim, yDim), 2) # test with various latent variable dimensionalities

marginalCovariances = "isotropic"
zDim <- 1

 toy <- generate.toydata(N = N, zDim = zDim, xDim = xDim, yDim = yDim, 
      	 		    marginal.covariances = marginalCovariances, 
      	 		    priors = priors)



# res <- fit.dependency.model(toy$X, toy$Y, zDimension = zDim,#
#	            marginalCovariances = marginalCovariances,#
#		    priors = priors, matched = TRUE, verbose = TRUE)
  

#fit.dependency.model <- function (
X <- toy$X
Y <- toy$Y
matched = TRUE
zDimension = zDim
          marginalCovariances = "isotropic"
          covLimit = 1e-3
          includeData = TRUE 
	  calculateZ = TRUE
	  verbose = FALSE

  if ( verbose ) { cat("Checking data\n") }
  dat <- check.data(X, Y, zDimension)
  X <- dat$X
  Y <- dat$Y
  zDimension <- dat$zDimension
	
  # FIXME store/return intercepts as well; further dependency models including intercepts
  if ( nrow(X) < nrow(Y) ) {stop("If the two data matrices do not have equal dimensionality, place the smaller one in Y.")} # FIXME automate
  if ( !nrow(X) == nrow(Y) ) {
    #message("The data sets have unequal dimensionality.")
    if ( matched ) { stop("Cannot use matched methods for nonmatched data.") }
  }

  if ( verbose ) { cat("Checking inputs\n") }
  if ( !is.null(priors$Nm.wxwy.sigma) && priors$Nm.wxwy.sigma == Inf ) { matched <- FALSE; message("priors$Nm.wxwy.sigma == Inf; Wx ~ Wy independendent i.e. matched = FALSE") }  
  if ( covLimit == 0 )  { covLimit <- 1e-3 } # avoid numerical overflows
  res <- NA; method <- ""
  
  if (matched) {

    if ( verbose ) { cat("Matched features case\n") }
      
    # Matrix normal distribution variance not specified
    if (is.null(priors$Nm.wxwy.sigma)) {
      message("Matched variables but priors$Nm.wxwy.sigma not given, using strong matching with Wx = Wy.")
      priors$Nm.wxwy.sigma <- 0
    }
      
    # Matrix normal distribution mean matrix not specified
    if (is.null(priors$Nm.wxwy.mean)) {
      message("The matrix Nm.wxwy.mean is not specified. Using identify matrix.")
      priors$Nm.wxwy.mean <- 1
    }    

    method <- "pSimCCA"
        
    # Case IIa: fully constrained case Wx = Wy
    if (priors$Nm.wxwy.sigma == 0) { #Wx = Wy        
        
      if ( verbose ) { cat("Assuming Wx = Wy\n") }
	
      #  SimCCA with full covariances with constraint Wx = Wy
      #  "probsimcca.full" = aucs.simcca.full      
      #  Denoting Wy = T*Wx = TW; W = Wx this fits the case T = I with
      #  full-rank Wx, Wy, Sigmax, Sigmay: (dxd-matrices where d equals to
      #  number of features in X and Y)

      # If prior for W is given, we must optimize W (no analytical
      # solution to EM)
          
      # Regularization for W (W > 0 implemented)
      if (!is.null(priors$W)) {
            
	if ( verbose ) { cat("Wx = Wy with regularized W (W>=0)\n") }
	if ( verbose ) { cat(marginalCovariances); cat("\n") }	

        res <- optimize.parameters(X, Y, zDim = zDimension, priors = priors, marginalCovariances = marginalCovariances, epsilon = covLimit, par.change = 1e6, verbose = verbose)

	method <- "pCCA with W prior"
       	
      } else if (is.null(priors$W)) {
        
	if ( verbose ) {cat("Wx = Wy; free W.\n")}

          # mlsp'09 simcca
          # message("Case Wx = Wy. No regularization for W.")
	  
	  # use this for full W (EM algorithm, unstable for n ~ p)
	  res <- optimize.parameters(X, Y, zDimension, priors, 
	      	 		     marginalCovariances, epsilon = covLimit,
				     par.change = 1e6)
				     
          method <- "matched case Wx = Wy with unconstrained W. Check covariances from parameters."
          # FIXME: speeups possible here when Wx = Wy but not yet implemented with other than full covs
      }
      
    } else if ( priors$Nm.wxwy.sigma > 0 ) {
      # Case IIb: partially constrained Wx ~ Wy
                
      if ( verbose ) { cat("partially constrained Wx ~ Wy.\n") }
		
      if ( !is.null(priors$W) ) {
        if ( verbose ) {cat("regularized W.\n")}
        # FIXME: consider adding fast option with simply nonnegative W but no distributional priors
        stop("Not implemented regularization for W with partially constrained Wx ~ Wy.")
      } else if (is.null(priors$W)) {
        if ( verbose ) { cat("Partially constrained Wx ~ Wy. No regularization for W.\n") }
        if ( verbose ) { cat(marginalCovariances); cat("\n") }            		  
			  
        if ( marginalCovariances == 'isotropic' ) {
	  # message("SimCCA with isotropic covariances and regularized H (through sigmas).")
	
          # FIXME: consider later adding other covariance structures if needed?
	  # note that the computation is slow then          		

          res <- optimize.parameters(X, Y, zDimension, priors, marginalCovariances, epsilon = covLimit)
          method <- "constrained Wx~Wy with matrix normal distribution prior"

        } else if ( !marginalCovariances == 'isotropic' ) {
          stop("Only isotropic marginal covariances implemented with constrained Wx ~ Wy in the general case.")
        }
      } 
    }
  }

  ##################################################################

  if ( verbose ) { cat("Checking the model..\n") }

  # Test whether model exists for given arguments
  if (any(is.na(unlist(res)))) {
    stop("Error with model parameters.")
  } else {
    params <- list(marginalCovariances = marginalCovariances, Nm.wxwy.mean = priors$Nm.wxwy.mean, Nm.wxwy.sigma = priors$Nm.wxwy.sigma, zDimension = zDimension, covLimit = covLimit)
    score <- dependency.score(res)
  }
  
  if ( verbose ) {cat("Creating DependencyModel object..\n")}
  
  model <- new("DependencyModel", W = res$W, phi = res$phi, score = score, method = method, params = params)	
  
