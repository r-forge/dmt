optimize.parameters <- function (X, Y, zDimension = 1, priors = NULL, 
                                 marginalCovariances = "full", 
				 epsilon = 1e-6, par.change = 1e6, verbose = FALSE) {

  # Suitable for at least:
  # nonmatched, prior$W, full marginals

  # Different from simCCA.optimize.R in that T is not optimized here
  # (not included in the model) but there is option to set prior on W
  # (W.prior)

  if ( verbose ) {cat("Initialize\n")}
  inits <- initialize2(X, Y, zDimension, marginalCovariances)
  phi <- inits$phi
  phi.inv <- inits$phi.inv
  W <- inits$W
  Dcov <- inits$Dcov
  Dim <- inits$Dim
  nullmat <- inits$nullmat
  Nsamples <- inits$Nsamples
	      
  # FIXME: handle priors completely outside this function later!
  
  if ( length(priors) == 0 ) { priors <- list() }

  ###  Wx ~ Wy prior inits  ###

  if ( verbose ) {cat("Checking the priors\n")}
  
  if ( !is.null(priors$Nm.wxwy.mean) ) {
    if ( length(priors$Nm.wxwy.mean) == 1 ){ 
      priors$Nm.wxwy.mean <- priors$Nm.wxwy.mean*diag(1, nrow(X), nrow(Y)) 
    }
    if ( ncol(priors$Nm.wxwy.mean) != nrow(X)){ stop("columns of priors$Nm.wxwy.mean must match rows of X") }
    if ( nrow(priors$Nm.wxwy.mean) != nrow(Y)){ stop("rows of priors$Nm.wxwy.mean must match rows of Y") }  
  }

  if ( is.null(priors$Nm.wxwy.sigma) || priors$Nm.wxwy.sigma == Inf ) { 
    # Wx, Wy relation not constrained
        if ( verbose ) {cat("Wx ~ Wy free\n")}
    
    priors$Nm.wxwy.sigma <- Inf 
    # cost.W.exponential accepts also priors$W = NULL i.e. no W prior
    cost.new <- cost.W.exponential(c(as.vector(W$X), as.vector(W$Y)), phi, priors, Dim, Dcov)

  } else if (priors$Nm.wxwy.sigma > 0) { # Wx ~ Wy constrained

    if ( verbose ) {cat("Wx ~ Wy constrained\n")}

    priors$T.tmp <- 1/(2 * Nsamples * priors$Nm.wxwy.sigma)
    # We assume here that Wy = T%*%Wx. Optimizing also T.
    T <- array(rnorm(Dim$Y*Dim$X,0,sd(W$X)), dim = c(Dim$Y, Dim$X))
    # Ensure that Wy = T * Wx:
    W$Y <- T%*%W$X
    # W not necessarily positive here

    cost.new <- cost.W(c(as.vector(W$X),as.vector(T)), phi, priors, Dim, Dcov)

  } else if (priors$Nm.wxwy.sigma == 0) { # Wx = Wy

        if ( verbose ) {cat("Wx = Wy \n")}

    # Ensure that the dimensionality of given w matches with given zDimension
    w <- as.matrix(inits$W$X[, 1:zDimension], ncol = zDimension)
    W <- list(X = w, Y = w, total = rbind(w, w))
    if ( !is.null(priors$W) ) {
      if ( verbose ) {cat(paste("prior for W: ", priors$W, "\n"))}
      cost.new <- cost7(abs(as.vector(W$X)), phi, Dcov, Dim, priors)    
    } else {
       if ( verbose ) {cat(paste("no prior for W. \n"))}
      cost.new <- cost7(as.vector(W$X), phi, Dcov, Dim, priors)        
    }  
  }
  
  
  if ( verbose ) {cat(paste("Starting iterations \n"))}
  while (par.change > epsilon || par.change < 0) {

    cost.old <- cost.new

    ###################################################

    # Update W: initialize with previous W	

    W.old <- W

    if (!is.null(priors$W)) {
      # Regularized W

      if (is.null(priors$Nm.wxwy.sigma) || priors$Nm.wxwy.sigma == Inf) {
        # optimizes Wx and Wy assuming they are independent
        opt <- optim(c(as.vector(W$X), as.vector(W$Y)), cost.W.exponential, 
	             method = "L-BFGS-B", phi = phi, priors = priors, 
		     Dim = Dim, Dcov = Dcov, control = list(maxit = 1e6), 
		     lower = -10*max(Dcov$total), upper = 10*max(Dcov$total))

        # Convert optimized W parameter vector to actual matrices
        # Note that here we always assume that W is positive
        W <- get.W2(opt$par, Dim)
      } else if (priors$Nm.wxwy.sigma == 0) {
      
      	# assuming Wx = Wy, we can speed up (FIXME; analytical alternatives?)

        # SimCCA Wx = Wy with regularized W (W>=0)
        # message("Case Wx = Wy and regularized W.")

	opt <- optim(as.vector(W$X), cost7, method = "L-BFGS-B", 
	             phi = phi, priors = priors, Dim = Dim, Dcov = Dcov, 
		     control = list(maxit = 1e6), 
                     lower = -10*max(Dcov$total), upper = 10*max(Dcov$total))

        w <- W$X <- W$Y <- get.W4(abs(opt$par), Dim)$X
	W$total <- rbind(w, w)
		
      } else {
        stop("W regularization implemented only for identical or independent Wx, Wy ie. priors$Nm.wxwy.sigma = 0 and priors$Nm.wxwy.sigma = Inf")
        # FIXME add the intermediates, should be straightforward by combining penalized optimizations
      }

    } else { # Unconstrained W
        
      if (priors$Nm.wxwy.sigma == 0) { # Wx = Wy

        # assuming Wx = Wy
        # see Bach-Jordan 2005, sec. 4.1 for details
	# equations modified from there to match Wx = Wy case

        W <- W.simcca.EM(W, phi, Dim, Dcov)
		    		    
      } else if (priors$Nm.wxwy.sigma > 0 && priors$Nm.wxwy.sigma < Inf) { # Wx ~ Wy constrained

        # Update W: initialize with previous W			       
        opt <- optim(c(as.vector(W$X), as.vector(T)), cost.W, method = "L-BFGS-B", phi = phi, priors = priors, Dim = Dim, Dcov = Dcov,
               control = list(maxit = 1e6),lower = -10*max(Dcov$total), upper = 10*max(Dcov$total))
      
        # Convert optimized W parameter vector to actual matrices
        wt <- get.W(opt$par, Dim)
        W <- wt$W
        T <- wt$T		
      } else { # Wx, Wy, independent a priori -> pCCA
        stop("Special case, corresponding to pCCA. No need to loop over W, phi in optimization. Use pCCA function for direct solution.")
      }
    }
    
    W.new <- W # redundant?
    
    ##################################################

    # Update phi

    phi.inv$X <- solve(phi$X)
    phi.inv$Y <- solve(phi$Y)    
    phi.inv$total <- rbind(cbind(phi.inv$X, nullmat),
                           cbind(t(nullmat), phi.inv$Y))    

    if (marginalCovariances == "full") {

      if (priors$Nm.wxwy.sigma == 0) { # Wx = Wy
 
        # FIXME: implement this also for other covariance structures

        # see Bach-Jordan 2005, sec. 4.1 for details
        M <- solve(t(W.old$X)%*%(phi.inv$X + phi.inv$Y)%*%W.old$X + diag(zDimension))
        # beta <- M%*%t(W$X)%*%phi.inv.sum # ct. set.beta.fullcov(M, W$total, phi.inv$total)
        # FIXME: replace this with the general M.set functions
        # FIXME: speedup by sharing M/beta with W.simcca.EM?

        # Update phi
        phi <- phi.EM.simcca(Dcov, W.new, phi.inv, W.old, M)
	     
      } else {  # assuming in general Wx != Wy
      
        # also check from optimize.fullcov.R
        M <- set.M.full2(W.old, phi.inv, dz = Dim$Z) 
        phi <- phi.EM.cca(Dcov, W.new, phi.inv, W.old, M, nullmat)
      }

     } else if (marginalCovariances == "isotropic") {

        # FIXME: speedups possible when Wx = Wy. Implement.

        # set.M is for isotropic X or Y
        M <- list()
        M$X <- set.M(W$X, phi$X)
        M$Y <- set.M(W$Y, phi$Y)
	       
        beta <- list()
        beta$X <- set.beta(M$X, W$X, phi$X)
        beta$Y <- set.beta(M$Y, W$Y, phi$Y)

        phi <- update.phi(Dcov, M, beta, W, phi)

     } else if ( marginalCovariances == "diagonal" ) {
     
       # FIXME: speedups possible when Wx = Wy. Implement.     
       # FIXME: compare M, beta to isotropic/full cases and join common parts
       # FIXME needs to be checked!
       phi <- phi.diagonal.double(W$total, phi.inv$total, Dcov$total, Dim)
       #phi$X <- phi.diagonal.single(W$total, phi.inv$total, Dcov$X, Dim)       
       #phi$Y <- phi.diagonal.single(W$total, phi.inv$total, Dcov$Y, Dim)            
     
     } else {
       stop("Unknown marginalCovariances parameter!")
     }

    ##########################################################################

    # MONITORING CONVERGENCE

    if (!is.null(priors$W)) { # W regularized

      # FIXME: cost7 and cost.W.exponential should both optimize nonneg W; combine?
      if (priors$Nm.wxwy.sigma == 0) {
        cost.new <- cost7(abs(as.vector(W$X)), phi, Dcov, Dim, priors)
      } else if (priors$Nm.wxwy.sigma == Inf) {
        cost.new <- cost.W.exponential(c(as.vector(W$X), as.vector(W$Y)), phi, priors, Dim, Dcov)      
      } else {
        stop("W regularization implemented only for independent and identical Wx, Wy cases.")
        # FIXME add the intermediates also
      }
    
    } else { # W not regularized; Wx ~ Wy is regularized

      if (priors$Nm.wxwy.sigma == 0) { # Extreme case Wx = Wy
        #message("Corresponds to simcca: W free; Wx = Wy.")
        cost.new <- cost7(as.vector(W$X), phi, Dcov, Dim, priors)    	
      } else if (priors$Nm.wxwy.sigma > 0 && priors$Nm.wxwy.sigma < Inf) {
        cost.new <- cost.W(c(as.vector(W$X), as.vector(T)), phi, priors, Dim, Dcov)
      } else if  (priors$Nm.wxwy.sigma == Inf)  {
      	stop("Corresponds to pCCA. No need for (slow) iterative optimization. Use pcca function instead.")
      } else {
        stop("Provide proper (nonneg. real) value for priors$Nm.wxwy.sigma")
      }
    }
    
    par.change <- (cost.old - cost.new)

  }

    if ( verbose ) {cat(paste(" Iterations OK. \n"))}

  if (marginalCovariances == "isotropic") {
    # force these scalars into diagonal matrices
    phi$X <- diag(phi$X, nrow(X))
    phi$Y <- diag(phi$Y, nrow(Y))
    phi$total <- diag(c(diag(phi$X),diag(phi$Y)),(nrow(X)+nrow(Y)))
  }

  W$total <- rbind(W$X, W$Y)  
  rownames(W$X) <- rownames(X)
  rownames(W$Y) <- rownames(Y)
  rownames(W$total) <- c(rownames(X), rownames(Y))
  
  rownames(phi$X) <- colnames(phi$X) <- rownames(X)
  rownames(phi$Y) <- colnames(phi$Y) <- rownames(Y)
  rownames(phi$total) <- colnames(phi$total) <- c(rownames(X), rownames(Y))

  return( list(W = W, phi = phi) )

}

