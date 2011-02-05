optimize.parameters <- function (X, Y, zDimension = 1, priors = NULL, 
                                 marginalCovariances = "full", 
				 epsilon = 1e-6, par.change = 1e6, cost.old = 1e6) {

  # Suitable for at least:
  # nonmatched, prior$W, full marginals

  # Different from simCCA.optimize.R in that T is not optimized here
  # (not included in the model) but there is option to set prior on W
  # (W.prior)
  
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
  # tune similarity constraint Wx ~ Wy
  if ( is.null(priors$Nm.wxwy.sigma) ) { priors$Nm.wxwy.sigma <- Inf } 
  if ( length(priors$Nm.wxwy.mean) == 1 ){ priors$Nm.wxwy.mean <- diag(1, nrow(X), nrow(Y)) }
  if ( ncol(priors$Nm.wxwy.mean) != nrow(X)){ stop("columns of H must match rows of X") }
  if ( nrow(priors$Nm.wxwy.mean) != nrow(Y)){ stop("rows of H must match rows of Y") }  
  if ( is.null(priors$W) ) { sigma2.W <- Inf }
  priors$W.tmp <- 1/(2 * sigma2.W * Nsamples)
  priors$T.tmp <- 1/(2 * Nsamples * priors$Nm.wxwy.sigma)
  if (!is.null(priors$W)) {
    # cost.W.exponential accepts also priors$W = NULL i.e. no W prior
    cost.new <- cost.W.exponential(c(as.vector(W$X), as.vector(W$Y)), phi, priors, Dim, Dcov)
  } else {
     # We assume here that Wy = T%*%Wx. Optimizing also T.
     T <- T.init <- array(rnorm(Dim$Y*Dim$X,0,sd(W$X)), dim = c(Dim$Y, Dim$X))
     # Ensure that Wy = T * Wx:
     W$Y <- T.init%*%W$X
     # here W not necessarily positive
     cost.new <- cost.W(c(as.vector(W$X),as.vector(T)), phi, priors, Dim, Dcov)
  }
  
  if (priors$Nm.wxwy.sigma == 0) { # Wx = Wy
    w <- inits$W$X  
    # Ensure that the dimensionality of given w matches with given zDimension
    w <- w[, 1:zDimension]
    W <- list()
    W$X <- W$Y <- w
    W$total <- rbind(w, w)
 
    cost.new <- initcost <- cost7(abs(as.vector(W$X)), phi, Dcov, Dim, priors)    
    
  }

  while (par.change > epsilon || par.change < 0) {

    cost.old <- cost.new

    ###################################################

    # Update W: initialize with previous W	

    W.old <- W

    if (!is.null(priors$W)) {

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

    } else {
        			   
       # Update W: initialize with previous W			       
       opt <- optim(c(as.vector(W$X), as.vector(T)), cost.W, method = "L-BFGS-B", phi = phi, priors = priors, Dim = Dim, Dcov = Dcov,
                control = list(maxit = 1e6),lower = -10*max(Dcov$total), upper = 10*max(Dcov$total))
       
      # Convert optimized W parameter vector to actual matrices
      wt <- get.W(opt$par, Dim)
      W <- wt$W
      T <- wt$T		

    }
    
    W.new <- W # redundant?
    
    ##################################################

    # Update phi

    if (marginalCovariances == "full") {

      phi.inv$X <- solve(phi$X)
      phi.inv$Y <- solve(phi$Y)    
      phi.inv$total <- rbind(cbind(phi.inv$X, nullmat),
                           cbind(t(nullmat), phi.inv$Y))    


      if (priors$Nm.wxwy.sigma == 0) {
     
        # Wx = Wy
        # FIXME: implement this also for other covariance structures

        # see Bach-Jordan 2005, sec. 4.1 for details
        M <- solve(t(W.old$X)%*%(phi.inv$X + phi.inv$Y)%*%W.old$X + diag(zDimension))
        # FIXME: replace this with the general M.set functions

        # Update phi
        phi <- phi.EM.simcca(Dcov, W.new, phi.inv, W.old, M)
	     
      } else {
      
        # also check from optimize.fullcov.R
        M <- set.M.full2(W.old, phi.inv, dz = Dim$Z) 
        #M <- set.M.full2(W.new, phi.inv, dz = Dim$Z) # new or old here? check

        # assuming in general Wx != Wy
        phi <- phi.EM.cca(Dcov, W.new, phi.inv, W.old, M, nullmat)

     }
     
    } else if (marginalCovariances == "isotropic") {

       # set.M is for isotropic X or Y
       M <- list()
       M$X <- set.M(W$X, phi$X)
       M$Y <- set.M(W$Y, phi$Y)
	       
       beta <- list()
       beta$X <- set.beta(M$X, W$X, phi$X)
       beta$Y <- set.beta(M$Y, W$Y, phi$Y)

       if (is.null(priors$W)) {
         phi <- update.phi(Dcov, M, beta, W, phi)
       }     

    }

    # Check and print marginal likelihood (-logP) for the data
    # the smaller, the better are the parameters
    if (!is.null(priors$W)) {

      # FIXME: cost7 and cost.W.exponential should both optimize nonneg W; combine
      if (priors$Nm.wxwy.sigma == 0) {
        cost.new <- cost7(abs(as.vector(W$X)), phi, Dcov, Dim, priors)
      } else if (priors$Nm.wxwy.sigma == Inf) {
        cost.new <- cost.W.exponential(c(as.vector(W$X), as.vector(W$Y)), phi, priors, Dim, Dcov)      
      } else {
        stop("W regularization implemented only for independent and identical Wx, Wy cases.")
        # FIXME add the intermediates also
      }
    
    } else {
      cost.new <- cost.W(c(as.vector(W$X), as.vector(T)), phi, priors, Dim, Dcov)
    }
    
    par.change <- (cost.old - cost.new)

  }

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

