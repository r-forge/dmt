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

  # Check that data is given as a matrix
  if (!is.matrix(X) || !is.matrix(Y)){
    stop("Data needs to be given as a matrix")
  }
  
  if (covLimit == 0)  { covLimit <- 1e-3 } # avoid numerical overflows

  # Center data
  X <- t(centerData(t(X), rm.na = TRUE))
  Y <- t(centerData(t(Y), rm.na = TRUE))

  # Check dimensionality
  if(zDimension > nrow(X) || zDimension > nrow(Y)) {
    warning("zDimension cannot exceed data dimensionality; using full dimensionality min(ncol(X), ncol(Y)).")
    zDimension <- min(nrow(X), nrow(Y))
  }       

  if (nrow(X) < nrow(Y)) {stop("If the two data matrices do not have equal dimensionality, then place the smaller one in Y.")} # FIXME automate
  
  # Storage for results
  res <- NA; method <- ""

  if (!nrow(X) == nrow(Y)) {
    #message("The data sets have unequal dimensionality.")
    if (matched) {stop("Cannot use matched methods for nonmatched data.")}
  }

  
  # Exponential prior for elements in W. Forces nonnegative solutions.
  if (!is.null(priors$W)) {
    nonnegative.w <- TRUE
    #if (!priors$W > 0) { stop("Nonnegative prior for priors$W required!") }
  } else {
    nonnegative.w <- FALSE
    # warning("No priors for W (priors$W) provided with nonmatched variables; no constraints applied.")
    # FIXME: implement/apply here the fast non-prior case a.k.a pCCA when no
    # priors are provided (for nonmatched case at least)?
  }
  
  #####################################

  if (!matched) {

    #message("Assuming non-matched variables.")
    if (!is.null(priors$Nm.wxwy.sigma)) { 
      warning("priors$Nm.wxwy.sigma not implemented for non-matched variables. Setting priors$Nm.wxwy.sigma = NULL.")
      priors$Nm.wxwy.sigma <- NULL
    }
  
    if (marginalCovariances == "full") {

      if (!is.null(priors$W)){
        # Prior for W is given -> need to optimize W (no analytical solution to EM)
        # Wx ~ Wy free; priors for W given    
	# Used to force positive W with exponential distribution.
        # priors$W is the rate parameter of the exponential. 
        # The smaller, the flatter tail.
        res <- simCCA.optimize3(X, Y, zDimension, epsilon = covLimit, priors = priors)
        method <- "pCCA with W priori"
      } else {
        # Using normal pcca when no W prior is given
	# FIXME: do we really need two separate functions: pcca and calc.pcca; and likewise with pca and pfa?
        method <- "pCCA"
        res <- calc.pcca(X, Y, zDimension)
      }
    } else if (marginalCovariances == "diagonal"){          
      # Probabilistic factor analysis model as proposed in          
      # EM Algorithms for ML Factoral Analysis, Rubin D. and      
      # Thayer D. 1982    
      res <- calc.pfa(X, Y, zDimension)    
      method <- "pFA"

    } else if(marginalCovariances == "identical isotropic"){        
      res <- calc.ppca(X, Y, zDimension)
      method <- "pPCA"
    } 
    
  } else if (matched) {

    # Matrix normal distribution variance not specified
    if (is.null(priors$Nm.wxwy.sigma)) {
      #message("priors$Nm.wxwy.sigma not given, using noninformative priors$Nm.wxwy.sigma = Inf, which corresponds to uncoupled Wx~Wy (pCCA).")
      priors$Nm.wxwy.sigma <- Inf
    }
      
    # Matrix normal distribution mean matrix not specified
    if (is.null(priors$Nm.wxwy.mean) && priors$Nm.wxwy.sigma != Inf) {
      warning("The matrix Nm.wxwy.mean is not specified. Setting H = 1.")
      priors$Nm.wxwy.mean <- 1
    }    
      
    if (priors$Nm.wxwy.sigma == 0 && !marginalCovariances == "full") {
      stop("With priors$sigma.w = 0 the matched simcca model is implemented only with full marginal covariances.")
    }    

    # Case I : Relation Wx, Wy not constrained.

    # Wx ~ Wy not constrained
    if ( priors$Nm.wxwy.sigma == Inf ) {

      #message("Relation between Wx, Wy is unconstrained.")

      if (!is.null(priors$Nm.wxwy.mean)){
        warning("Given prior Nm.wxwy.mean is ignored, beacuse Nm.wxwy.sigma == Inf")
        priors$Nm.wxwy.mean <- NULL
      }
      
      if (!nonnegative.w) {
       # message("Wx ~ Wy free. No regularization for W.")

        if (marginalCovariances == "full") {           
          # Analytical ML solution for CCA 
          res <- calc.pcca(X, Y, zDimension)            
          method <- "pCCA"
          #message("Full marginal covariances.")
        } else if (marginalCovariances == "diagonal"){          
          # Probabilistic factor analysis model as proposed in          
          # EM Algorithms for ML Factoral Analysis, Rubin D. and      
          # Thayer D. 1982    
          res <- calc.pfa(X, Y, zDimension)    
          method <- "pFA"
          #message("Diagonal marginal covariances.")
        } else if (marginalCovariances == "isotropic") {
          # pCCA assuming isotropic margins with phiX != phiY 
          res <- calc.pcca.with.isotropic.margins(X, Y, zDimension, epsilon = covLimit)  
          method <- "pCCA with isotropic margins"
          #message("Isotropic marginal covariances.")          
        } else if(marginalCovariances == "identical isotropic"){        
          res <- calc.ppca(X, Y, zDimension)
          #message("Identical isotropic marginal covariances.")
          method <- "pPCA"
        }                         
      } else if (nonnegative.w) {

        #message("Wx ~ Wy free. Regularized W (W>=0).")        

        if (!marginalCovariances == "full") {           
          warning("Only full marginal covariances implemented with uncontrained Wx/Wy and regularized W. Using full marginal covariances.")
          marginalCovariances <- "full"           
        }
        
        # Currently implemented exponential prior for W,
        # priors$W is the rate parameter of the exponential.
        # The smaller, the flatter. 
        # By default, this function does not constrain Wx~Wy
        # (priors$sigma.w = Inf) but imposes prior on W.  FIXME:
        # this should work also with constrained Wx ~ Wy, test and
        # compare FIXME: check and message how full/diag/isotropic
        # cov go here

        # any zDimension should work here
        res <- simCCA.optimize2(X, Y, zDimension, 
                                 epsilon = covLimit, priors = priors)
        method <- "pCCA with W priori"
      }

    } else if ( !priors$Nm.wxwy.sigma == Inf ){

      method <- "pSimCCA"
      #message("Regulating the relationship Wx ~ Wy.")
        
      # Case IIa: fully constrained case Wx = Wy
      if (priors$Nm.wxwy.sigma == 0) { #Wx = Wy        
        
        #  SimCCA with full covariances with constraint Wx = Wy
        #  "probsimcca.full" = aucs.simcca.full      
        #  Denoting Wy = T*Wx = TW; W = Wx this fits the case T = I with
        #  full-rank Wx, Wy, Sigmax, Sigmay: (dxd-matrices where d equals to
        #  number of features in X and Y)

        # If prior for W is given, we must optimize W (no analytical
        # solution to EM)
          
        # Regularization (W > 0)
        if (nonnegative.w) {
            
          # SimCCA Wx = Wy with regularized W (W>=0)
          #message("Case Wx = Wy and regularized W.")
          # Initialize (FIXME: make initialization as in the other options)
          inits <- initialize2(X, Y)
          phi.init <- inits$phi
          W.init <- inits$W
          Dcov <- inits$Dcov
          Dim <- inits$Dim
          Dim$Z <- zDimension
          nullmat <- inits$nullmat
          Nsamples <- inits$Nsamples
 
          # use this for full W (EM algorithm, unstable for n ~ p)
          res <- optimize.simCCA.W2(W.init$X,
                                      phi.init,
                                      Dim = Dim,
                                      Dcov = Dcov,
                                      nullmat = nullmat,
                                      epsilon = covLimit,
                                      par.change = 1e6,
                                      dz = zDimension,
                                      priors = priors)
          method <- "pCCA with W priori"

        } else if (!nonnegative.w) {
          # mlsp'09 simcca
          #message("Case Wx = Wy. No regularization for W.")
          res <- simCCA.optimize.fullcov.EM(X, Y, zDimension, epsilon = covLimit)
          method <- "pSimCCA"
          # FIXME: priors for phi todo?            
            
        }
      } else if (priors$Nm.wxwy.sigma != 0) {
        # Case IIb: partially constrained Wx ~ Wy
                
        if (nonnegative.w) {
          stop("Not implemented regularized (nonnegative) W with partially constrained Wx ~ Wy. Only special cases priors$sigma.w = 0 and priors$sigma.w = Inf are available.")            
        } else if (!nonnegative.w) {
          #message("Partially constrained Wx ~ Wy. No regularization for W.")
            
          if (marginalCovariances == 'isotropic') {

            # Make H identity matrix if scalar is given
            if(length(priors$Nm.wxwy.mean) == 1){ priors$Nm.wxwy.mean <- diag(1, nrow(X), nrow(Y)) }
            if(ncol(priors$Nm.wxwy.mean) != nrow(X)){ stop("columns of H must match rows of X") }
            if(nrow(priors$Nm.wxwy.mean) != nrow(Y)){ stop("rows of H must match rows of Y") }
                 
            #message("SimCCA with isotropic covariances and regularized H (through sigmas).")
            res <- simCCA.optimize(X, Y, zDimension, priors$Nm.wxwy.mean, sigma2.T = priors$Nm.wxwy.sigma, sigma2.W = 1e12, epsilon = covLimit)
            method <- "pSimCCA with T priori"
          } else if (!marginalCovariances == 'isotropic') {
            stop("With priors$sigma.w !=0, only isotropic marginal covariances implemented.")
          }
        }
      } 
    }
  }

  ##################################################################

  # Test whether model exists for given arguments
  if (any(is.na(res))) {
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


ppca <- function(X, Y = NULL, zDimension = NULL, matched = FALSE, includeData = TRUE, calculateZ = TRUE){                

  if (is.null(zDimension)) { zDimension = min(nrow(X), nrow(Y)) }

  if (!is.null(Y)) {                 
    model <- fit.dependency.model(X,Y,zDimension,marginalCovariances = "identical isotropic", matched = matched)        
  } else {         

    # Check that data is given as a matrix
    if ( !is.matrix(X) ){ stop("Data needs to be given as a matrix") }

    #if (ncol(X) > 1)                         
    X <- t(centerData(t(X), rm.na = TRUE))
           
    # Check if dimensionality is too big               
    if(zDimension > nrow(X)) {warning("Latent variable dimensionality cannot exceed data dimensionality. Setting zDimension = nrow(X)")}
                      
    res <- calc.ppca(X, Y, zDimension) 
    method <- "pPCA"    
                               
    params <- list(marginalCovariances = "isotropic", zDimension = zDimension, covLimit = 0)      
    score <- dependency.score( res ) 
    model <- new("DependencyModel", W = res$W, phi = res$phi, score = score, method = method, params = params)   
    if (includeData) model@data <- list(X = X)
    if (calculateZ) model@z <- z.expectation(model, X)

  }
  
  return(model)        
}    

               	      
pfa <- function(X, Y = NULL, zDimension = NULL, matched = FALSE, includeData = TRUE, calculateZ = TRUE){       

  if (is.null(zDimension)) { zDimension = min(nrow(X), nrow(Y)) }

  if (!is.null(Y)) { 
    model <- fit.dependency.model(X, Y, zDimension, marginalCovariances = "diagonal", matched = matched,
                         includeData = includeData, calculateZ = calculateZ)       
  } else { 
    # Check that data is given as a matrix
    if (!is.matrix(X)){
      stop("Data needs to be given as a matrix")
    }
    #if (ncol(X) > 1)                         
    X <- t(centerData(t(X), rm.na = TRUE))               
           
    # Check if dimensionality is too big               
    if(zDimension > nrow(X)) {warning("Latent variable dimensionality cannot exceed data dimensionality. Setting zDimension = nrow(X)")}
                      
    res <- calc.pfa(X, Y, zDimension)        
    method <- "pFA"
                               
    params <- list(marginalCovariances = "diagonal", zDimension = zDimension, covLimit = 0)      
    score <- dependency.score( res ) 

    model <- new("DependencyModel", W = res$W, phi = res$phi, score = score, method = method, params = params)                                                                    
    if (includeData) model@data <- list(X = X, Y = Y)
    if (calculateZ) model@z <- z.expectation(model, X, Y)

  }                     

  return(model)
}                                                                     
                            
# FIME: consider whether we should keep this or remove			    
pcca.isotropic <- function(X, Y, zDimension = NULL, matched = FALSE, covLimit = 1e-6, includeData = TRUE, calculateZ = TRUE){
  if (is.null(zDimension)) { zDimension = min(nrow(X), nrow(Y)) }

  fit.dependency.model(X,Y,zDimension,marginalCovariances = "isotropic", covLimit = 1e-6,
                       includeData = includeData, calculateZ = calculateZ)          
}                                                                      
                                                                   
pcca <- function(X, Y, zDimension = NULL, matched = FALSE, includeData = TRUE, calculateZ = TRUE){  

  if (is.null(zDimension)) { zDimension = min(nrow(X), nrow(Y)) }
  
  fit.dependency.model(X, Y, zDimension, marginalCovariances = "full", matched = matched,
                       includeData = includeData, calculateZ = calculateZ)                                                     
}                                                                      
          
