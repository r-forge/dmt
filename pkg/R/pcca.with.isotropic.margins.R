pcca.with.isotropic.margins <- function (X, Y, zDimension = 1, epsilon = 1e-6, delta = 1e6) {

  # epsilon and delta are convergence parameters
  # zDimension determines the dimensionality of the shared latent variable Z

  #  Dependency model
  #  X ~ N(Wx*z, sigmax*I)
  #  y ~ N(Wy*z, sigmay*I)
  #  i.e. isotropic marginals but in general  sigmax != sigmay
  # This is a special case of probabilistic factor analysis model

  # FIXME: ensure that X, Y have zero-mean (shift if necessary);
  # alternatively add mean parameter in the model

  res <- calc.pcca.with.isotropic.margins(X, Y, zDimension, epsilon = 1e-6, delta = 1e6)

  phi <- res$phi  
    W <- res$W   


  colnames(phi$X) <- rownames(phi$X) <- rownames(X)
  colnames(phi$Y) <- rownames(phi$Y) <- rownames(Y)
  colnames(phi$total) <- rownames(phi$total) <- c(rownames(X), rownames(Y))
  
  list(phi = phi, W = W)

  # FIXME provide here proper DependencyModel object as in pcca, pfa and ppca
}


calc.pcca.with.isotropic.margins <- function (X, Y, zDimension, epsilon = 1e-6, delta = 1e6) {

  dat <- check.data(X, Y, zDimension)
  X <- dat$X
  Y <- dat$Y
  zDimension <- dat$zDimension
  
  # initialize
     inits <- initialize2(X, Y, zDimension, marginalCovariances = "isotropic")
      Dcov <- inits$Dcov
       Dim <- inits$Dim
         W <- inits$W  

  # FIXME: ensure that X, Y have zero-mean (shift if necessary);
  # alternatively add mean parameter in the model
  phi <- list(X = 1, Y = 1)
  
  # iterate until convergence:
  while (delta > epsilon) {

    W.old <- W
          
    ##########################################

    # Update Phi
    phi$X <- update.phi.isotropic(Dcov$X, W$X, phi$X, Dim$X) 
    phi$Y <- update.phi.isotropic(Dcov$Y, W$Y, phi$Y, Dim$Y)
          
    #######################################

    # Full CCA update for W

    phi.inv.full <- diag(c(rep(1/phi$X, Dim$X), rep(1/phi$Y, Dim$Y)))
          M <- set.M.full(W$total, phi.inv.full, Dim$Z) # corresponds to G in Bishop's book
       beta <- set.beta.fullcov(M, W$total, phi.inv.full)
    W$total <- W.cca.EM(Dcov, M, beta)
        W$X <- matrix(W$total[1:Dim$X,], nrow = Dim$X)
        W$Y <- matrix(W$total[-(1:Dim$X),], nrow = Dim$Y)

    ########################################
          
    # check convergence (enough to check W)
    delta <- max(abs(as.vector(W$total - W.old$total)))
          
  }

  # Format scalars into matrices
  #phi$X <- diag(phi$X, Dim$X)
  #phi$Y <- diag(phi$Y, Dim$Y)
  #phi$total <- diag(c(diag(phi$X), diag(phi$Y)))

  list(W = W, phi = phi)

}

