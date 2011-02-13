fit.dependency.model <- function (X, Y,
          zDimension = 1,
          marginalCovariances = "full",
          covLimit = 1e-3,
          priors = list(), matched = TRUE,
          includeData = TRUE, calculateZ = TRUE, verbose = FALSE)
{

  # (C) 2008-2011 Olli-Pekka Huovilainen and Leo Lahti 
  # FreeBSD license (keep this notice).

  # zDimension = 1; marginalCovariances = "full"; H = 1; sigmas = 0; covLimit = 1e-3; mySeed = 123; priors = NULL
  
  # Fits the generative model
  # X = Wx * z + epsx
  # Y = Wy * z + epsy
  # with various modeling assumptions

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

  if (verbose) { cat("Checking inputs\n") }
  if (!is.null(priors$Nm.wxwy.sigma) && priors$Nm.wxwy.sigma == Inf) { matched <- FALSE; message("priors$Nm.wxwy.sigma == Inf; Wx ~ Wy independendent i.e. matched = FALSE") }  
  if ( covLimit == 0 )  { covLimit <- 1e-3 } # avoid numerical overflows
  res <- NA; method <- ""
  
  if (!matched) {

    if (verbose) {cat("Model for non-matched case\n")}

    if ( is.null(priors$Nm.wxwy.sigma )) {
      #warning("priors$Nm.wxwy.sigma not implemented for non-matched variables. Setting priors$Nm.wxwy.sigma = Inf.")
      priors$Nm.wxwy.sigma <- Inf
    }

    if ( is.null(priors$W) ) { 

      if ( verbose ) {cat("Wx ~ Wy free. No regularization for W.\n")}
      if ( verbose ) {cat(marginalCovariances); cat("\n")}
      

      if ( marginalCovariances == "full" ) { # standard pCCA

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
      
      if ( verbose ) { cat("Wx ~ Wy free; exponential (nonnegative) prior for W.\n") }

      # Prior for W is given -> no analytical solution to EM
      # Exponential prior for W,
      # Used to force positive W with exponential distribution.
      # priors$W is the rate parameter of the exponential. 
      # The smaller, the flatter tail.

      # TODO: implement also sparsity prior W ~ N(0, sd*I)

      res <- optimize.parameters(X, Y, zDimension, priors = priors, marginalCovariances, epsilon = covLimit, verbose = verbose)

      method <- "Free Wx ~ Wy with exponential priors for Wx and Wy. Check marginal covariances from parameters."
	
    }
    
  } else if (matched) {

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
	
        res <- optimize.parameters(X, Y, zDimension, epsilon = covLimit, par.change = 1e6, priors = priors, verbose = verbose)
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
  if (includeData) model@data <- list(X = X, Y = Y)
  if (calculateZ) model@z <- z.expectation(model, X, Y)
  
  if ( verbose ) {cat("fit.dependency.model OK.\n")}
  
  model
}



# FIME: consider whether we should keep this or remove			    
#pcca.isotropic <- function(X, Y, zDimension = NULL, matched = FALSE, covLimit = 1e-6, includeData = TRUE, calculateZ = TRUE){
#  if (is.null(zDimension)) { zDimension = min(nrow(X), nrow(Y)) }#
#
#  fit.dependency.model(X,Y,zDimension,marginalCovariances = "isotropic", covLimit = 1e-6,
#                       includeData = includeData, calculateZ = calculateZ)          
#}                                                                      
