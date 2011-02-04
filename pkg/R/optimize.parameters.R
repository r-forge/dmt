optimize.parameters <-
function (W, phi, Dim, Dcov, priors, marginalCovariances = "full", epsilon = 1e-6, par.change = 1e6, cost.old = 1e6) {

  # NOTE:
  # 1) cost.W.exponential accepts also priors$W = NULL i.e. no W prior
  # 2) FIXME: implement other marginalCov structures and sparse W prior

  # initializing
  cost.new <- cost.W.exponential(c(as.vector(W$X), as.vector(W$Y)), phi, priors, Dim, Dcov)
  costs    <- c(cost.new)
  nullmat  <- matrix(0, nrow = Dim$X, ncol = Dim$Y)

  # FIXME: if phi$Y is scalar (as in segmented/mir case) we can speed up here. Do later.
  phi.inv  <- list()
  phi.inv$X <- solve(phi$X)
  phi.inv$Y <- solve(phi$Y)  
  phi.inv$total <- rbind(cbind(phi.inv$X, nullmat), cbind(t(nullmat), phi.inv$Y))

  cnt <- 1
  while (par.change > epsilon) {

    cost.old <- cost.new

    ##############################################

    # Update W: initialize with previous W	
    
    # This optimizes Wx and Wy assuming they are independent
    opt <- optim(c(as.vector(W$X), as.vector(W$Y)), cost.W.exponential, method = "L-BFGS-B", phi = phi, priors = priors, Dim = Dim, Dcov = Dcov, control = list(maxit = 1e6), lower = -10*max(Dcov$total), upper = 10*max(Dcov$total))

    # Convert optimized W parameter vector to actual matrices
    # Note that here we always assume that W is positive
    W.old <- W
    W.new <- W <- get.W2(opt$par, Dim)

    ##################################################

    # Update phi

    if (marginalCovariances = "full") {

      phi.inv$X <- solve(phi$X)
      phi.inv$Y <- solve(phi$Y)    
      phi.inv$total <- rbind(cbind(phi.inv$X, nullmat),
                           cbind(t(nullmat), phi.inv$Y))    

      # also check from optimize.fullcov.R
      M <- set.M.full2(W.old, phi.inv, dz = Dim$Z) 
      #M <- set.M.full2(W.new, phi.inv, dz = Dim$Z) # new or old here? check

      # assuming in general Wx != Wy
      phi <- phi.EM.cca(Dcov, W.new, phi.inv, W.old, M, nullmat)

    } #else if (marginalCovariances = "isotropic") {
      #   
      #phi <- update.phi(Dcov, M, beta, W, phi)
   # 
   #   stop("With priors for W, only marginalCovariances = full has been implemented.")
   #   # FIXME: add other cov structures
   # }

    ###################################################

    # Check and print marginal likelihood (-logP) for the data
    # the smaller, the better are the parameters
    cost.new <- cost.W.exponential(c(as.vector(W$X), as.vector(W$Y)), phi, priors, Dim, Dcov)
    
    par.change <- (cost.old - cost.new)
    cnt <- cnt + 1
    costs[[cnt]] <- cost.new

  }

  list(W = W, phi = phi, costs = costs, score = costs[[length(costs)]])

}
