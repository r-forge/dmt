optimize.parameters <- function (X, Y, zDimension = 1, priors = NULL, marginalCovariances = "full", epsilon = 1e-6, par.change = 1e6, cost.old = 1e6) {

  # Suitable for at least:
  # nonmatched, prior$W, full marginals

  # Different from simCCA.optimize.R in that T is not optimized here
  # (not included in the model) but there is option to set prior on W
  # (W.prior)
  
  # FIXME: make this universal optimization function which combines also T as optional thing

  ###############################################

  # samples are always matched i.e. ncol(X) = ncol(Y)
  Nsamples <- ncol(X)
  
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
  
  ##############################################

  Dim <- list()
  Dim$X <- nrow(X)
  Dim$Y <- nrow(Y)
  Dim$Z <- zDimension
  nullmat  <- matrix(0, nrow = Dim$X, ncol = Dim$Y)
  
  Dcov <- list()
  Dcov$X <- cov(t(X))
  Dcov$Y <- cov(t(Y))
  Dcov$total <- cov(t(rbind(X, Y)))

  # initialize with scalar diagonal noise on the marginals (shared by all features)
  phi <- phi.init <- list(X = diag(var(as.vector(X)), Dim$X), Y = diag(var(as.vector(Y)), Dim$Y)) 

  if (marginalCovariances == "isotropic") {
    phi$X <- var(as.vector(X))
    phi$Y <- var(as.vector(Y))   
    phi.init <- phi # FIXME can phi.init be replaced solely by phi everywhere??
  }

  # FIXME: if phi$Y is scalar (as in segmented/mir case) we can speed up here. Do later.
  phi.inv  <- list()
  phi.inv$X <- solve(phi$X)
  phi.inv$Y <- solve(phi$Y)  
  phi.inv$total <- rbind(cbind(phi.inv$X, nullmat), cbind(t(nullmat), phi.inv$Y))

  # Initialize W's
  W.init   <- list()
  W.init$X <- as.matrix(eigen(Dcov$X)$vectors[, 1:Dim$Z])
  W.init$Y <- as.matrix(eigen(Dcov$Y)$vectors[, 1:Dim$Z])
  W.init$total <- rbind(W.init$X, W.init$Y) # can be possibly removed in some special cases
  W <- W.init

  if (!is.null(priors$W)) {
    # cost.W.exponential accepts also priors$W = NULL i.e. no W prior
    cost.new <- cost.W.exponential(c(as.vector(W$X), as.vector(W$Y)), phi, priors, Dim, Dcov)
  } else {
     # We assume here that Wy = T%*%Wx. Optimizing also T.
     T <- T.init <- array(rnorm(Dim$Y*Dim$X,0,sd(W.init$X)), dim = c(Dim$Y, Dim$X))
     # Ensure that Wy = T * Wx:
     W.init$Y <- T.init%*%W.init$X
     W <- W.init
      
     # here W not necessarily positive
     cost.new <- cost.W(c(as.vector(W$X),as.vector(T)), phi, priors, Dim, Dcov)
  }
    

  while (par.change > epsilon) {

    cost.old <- cost.new

    ###################################################

    # Update W: initialize with previous W	

    if (!is.null(priors$W)) {

      # optimizes Wx and Wy assuming they are independent
      opt <- optim(c(as.vector(W$X), as.vector(W$Y)), cost.W.exponential, method = "L-BFGS-B", phi = phi, priors = priors, Dim = Dim, Dcov = Dcov, control = list(maxit = 1e6), lower = -10*max(Dcov$total), upper = 10*max(Dcov$total))

      # Convert optimized W parameter vector to actual matrices
      # Note that here we always assume that W is positive
      W.old <- W
      W.new <- W <- get.W2(opt$par, Dim)

    } else {
        			   
       # Update W: initialize with previous W			       
       opt <- optim(c(as.vector(W$X), as.vector(T)), cost.W, method = "L-BFGS-B", phi = phi, priors = priors, Dim = Dim, Dcov = Dcov,
                control = list(maxit = 1e6),lower = -10*max(Dcov$total), upper = 10*max(Dcov$total))
       
      # Convert optimized W parameter vector to actual matrices
      wt <- get.W(opt$par, Dim)
      W <- wt$W
      T <- wt$T		

    }
    
    ##################################################

    # Update phi

    if (marginalCovariances == "full") {

      phi.inv$X <- solve(phi$X)
      phi.inv$Y <- solve(phi$Y)    
      phi.inv$total <- rbind(cbind(phi.inv$X, nullmat),
                           cbind(t(nullmat), phi.inv$Y))    

      # also check from optimize.fullcov.R
      M <- set.M.full2(W.old, phi.inv, dz = Dim$Z) 
      #M <- set.M.full2(W.new, phi.inv, dz = Dim$Z) # new or old here? check

      # assuming in general Wx != Wy
      phi <- phi.EM.cca(Dcov, W.new, phi.inv, W.old, M, nullmat)

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
      cost.new <- cost.W.exponential(c(as.vector(W$X), as.vector(W$Y)), phi, priors, Dim, Dcov)
    } else{
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

