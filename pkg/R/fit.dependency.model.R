fit.dependency.model <-
function (X, Y,
          zDimension = 1,
          marginalCovariances = "full",
          covLimit = 1e-3,
          priors = list(), matched = TRUE,
          includeData = TRUE, calculateZ = TRUE)
{

  # (C) 2008-2011 Olli-Pekka Huovilainen and Leo Lahti 
  # FreeBSD license (keep this notice).

  # zDimension = 1; marginalCovariances = "full"; H = 1; sigmas = 0; covLimit = 1e-3; mySeed = 123; priors = NULL
  
  # Fits the generative model
  # X = Wx * z + epsx
  # Y = Wy * z + epsy
  # with various modeling assumptions

  dat <- check.data(X, Y, zDimension)
  X <- dat$X
  Y <- dat$Y
  zDimension <- dat$zDimension
	
  # FIXME store/return intercepts as well; further dependency models including intercepts
  if (nrow(X) < nrow(Y)) {stop("If the two data matrices do not have equal dimensionality, place the smaller one in Y.")} # FIXME automate
  if (!nrow(X) == nrow(Y)) {
    #message("The data sets have unequal dimensionality.")
    if (matched) {stop("Cannot use matched methods for nonmatched data.")}
  }

  if (!is.null(priors$Nm.wxwy.sigma) && priors$Nm.wxwy.sigma == Inf) { matched <- FALSE; message("priors$Nm.wxwy.sigma == Inf; Wx ~ Wy independendent i.e. matched = FALSE") }
  
  if ( covLimit == 0 )  { covLimit <- 1e-3 } # avoid numerical overflows
  res <- NA; method <- ""
  
  if (!matched) {

    if ( !is.null(priors$Nm.wxwy.sigma )) {
      warning("priors$Nm.wxwy.sigma not implemented for non-matched variables. Setting priors$Nm.wxwy.sigma = NULL.")
      priors$Nm.wxwy.sigma <- NULL
    }

    if ( is.null(priors$W) ) { 

      # message("Wx ~ Wy free. No regularization for W.")

      if (marginalCovariances == "full") { # standard pCCA

        method <- "pCCA"
        res <- calc.pcca(X, Y, zDimension)

      } else if (marginalCovariances == "diagonal") { 

        # pFA for two data sets corresponds to 
        # standard pFA for concatenated data

        res <- calc.pfa(X, Y, zDimension)    
        method <- "pFA"

      } else if (marginalCovariances == "isotropic") {
      
        # phiX != phiY in general
        # FIXME: add tests
        res <- calc.pcca.with.isotropic.margins(X, Y, zDimension)
        method <- "pCCA with isotropic margins"
	
      } else if (marginalCovariances == "identical isotropic") {
      
        # FIXME: add tests for this
        # pPCA for two data sets corresponds to 
        # standard pPCA for concatenated data

        res <- calc.ppca(X, Y, zDimension)
        method <- "pPCA"
	
      } else { stop("Erroneous marginalCovariances parameter provided!") }

     } else if ( !is.null(priors$W) ) { # for some reason stating priors$W caused crash before   
      # Wx ~ Wy free; priors for W given    

      # Prior for W is given -> no analytical solution to EM
      # Exponential prior for W,
      # Used to force positive W with exponential distribution.
      # priors$W is the rate parameter of the exponential. 
      # The smaller, the flatter tail.

      # FIXME: in addition to exponential prior W ~ exp(alpha)
      # implement sparsity prior W ~ N(0, sd*I)

      res <- optimize.parameters(X, Y, zDimension, priors = priors, marginalCovariances, epsilon = covLimit)
      method <- "Unconstrained Wx ~ Wy with exponential priors for Wx and Wy. Check marginal covariances from parameters."
	
    }
    
  } else if (matched) {

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

    # FIXME siirra sisemmas
    if (priors$Nm.wxwy.sigma == 0 && !marginalCovariances == "full") {
      stop("With priors$sigma.w = 0 the matched simcca model is implemented only with full marginal covariances.")
    }    

    method <- "pSimCCA"
        
    # Case IIa: fully constrained case Wx = Wy
    if (priors$Nm.wxwy.sigma == 0) { #Wx = Wy        
        
      #  SimCCA with full covariances with constraint Wx = Wy
      #  "probsimcca.full" = aucs.simcca.full      
      #  Denoting Wy = T*Wx = TW; W = Wx this fits the case T = I with
      #  full-rank Wx, Wy, Sigmax, Sigmay: (dxd-matrices where d equals to
      #  number of features in X and Y)

      # If prior for W is given, we must optimize W (no analytical
      # solution to EM)
          
      # Regularization for W (W > 0 implemented)
      if (!is.null(priors$W)) {
            
        # SimCCA Wx = Wy with regularized W (W>=0)
        #message("Case Wx = Wy and regularized W.")

        if (marginalCovariances == "full") {
          # use this for full W (EM algorithm, unstable for n ~ p) (?)
          res <- optimize.parameters(X, Y, zDimension, epsilon = covLimit,                                  
                                  par.change = 1e6, priors = priors)                                 
          method <- "pCCA with W prior"
        } else if (marginalCovariances == "diagonal") {	
          stop("Matched case with regularized W implemented only with full marginalCovariances")
	  # FIXME: add diagonal covs
        } else if (marginalCovariances == "isotropic") {
          stop("Matched case with regularized W implemented only with full marginalCovariances")
	  # FIXME: add this and also identical isotropic
        } else {
          stop("Matched case with regularized W implemented only with full marginalCovariances")
	}      
       	
      } else if (is.null(priors$W)) {
        
          # mlsp'09 simcca
          # message("Case Wx = Wy. No regularization for W.")
	  
	  # use this for full W (EM algorithm, unstable for n ~ p)
	  res <- optimize.parameters(X, Y, zDimension, priors, 
	      	 		     marginalCovariances, epsilon = covLimit,
				     par.change = 1e6)
				     
          method <- "matched case Wx = Wy with unconstrained W. Check covariances from parameters."
          # FIXME: speeups possible here when Wx = Wy but not yet implemented with other than full covs
      }
      
    } else if (priors$Nm.wxwy.sigma > 0) {
      # Case IIb: partially constrained Wx ~ Wy
                
      if (!is.null(priors$W)) {
        # FIXME: consider adding fast option with simply nonnegative W but no distributional priors
        stop("Not implemented regularization for W with partially constrained Wx ~ Wy.")
      } else if (is.null(priors$W)) {
        # message("Partially constrained Wx ~ Wy. No regularization for W.")
            		    
        if (marginalCovariances == 'isotropic') {
	  # message("SimCCA with isotropic covariances and regularized H (through sigmas).")
	
          # FIXME: consider later adding other covariance structures if needed?
	  # note that the computation is slow then          		

          res <- optimize.parameters(X, Y, zDimension, priors, marginalCovariances, epsilon = covLimit)
          method <- "constrained Wx~Wy with matrix normal distribution prior"

        } else if (!marginalCovariances == 'isotropic') {
          stop("Only isotropic marginal covariances implemented with constrained Wx ~ Wy in the general case.")
        }
      } 
    }
  }

  ##################################################################

  # Test whether model exists for given arguments
  if (any(is.na(unlist(res)))) {
    stop("Error with model parameters.")
  } else {
    params <- list(marginalCovariances = marginalCovariances, Nm.wxwy.mean = priors$Nm.wxwy.mean, Nm.wxwy.sigma = priors$Nm.wxwy.sigma, zDimension = zDimension, covLimit = covLimit)
    score <- dependency.score(res)
  }
  
  model <- new("DependencyModel", W = res$W, phi = res$phi, score = score, method = method, params = params)	
  if (includeData) model@data <- list(X = X, Y = Y)
  if (calculateZ) model@z <- z.expectation(model, X, Y)
  model
}



# FIME: consider whether we should keep this or remove			    
#pcca.isotropic <- function(X, Y, zDimension = NULL, matched = FALSE, covLimit = 1e-6, includeData = TRUE, calculateZ = TRUE){
#  if (is.null(zDimension)) { zDimension = min(nrow(X), nrow(Y)) }#
#
#  fit.dependency.model(X,Y,zDimension,marginalCovariances = "isotropic", covLimit = 1e-6,
#                       includeData = includeData, calculateZ = calculateZ)          
#}                                                                      
