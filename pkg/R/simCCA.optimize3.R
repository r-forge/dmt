simCCA.optimize3 <- function (X, Y, zDimension = 1, epsilon = 1e-6, priors = NULL, marginalCovariances = "full") {

  # Suitable for at least:
  # nonmatched, prior$W, full marginals

  #################################################

  # Different from simCCA.optimize.R in that T is not optimized here
  # (not included in the model) but there is option to set prior on W
  # (W.prior)
  # FIXME: make this universal optimization function which combines also T as optional thing

  #################################################

  # Initialize

  # samples are always matched i.e. ncol(X) = ncol(Y)
  Nsamples <- ncol(X)

  if ( length(priors) == 0 ) { priors <- list() }
  if ( is.null(priors$Nm.wxwy.sigma) ) { priors$Nm.wxwy.sigma <- Inf } # tune similarity constraint Wx ~ Wy

  Dim <- list()
  Dim$X <- nrow(X)
  Dim$Y <- nrow(Y)
  Dim$Z <- zDimension

  Dcov <- list()
  Dcov$X <- cov(t(X))
  Dcov$Y <- cov(t(Y))
  Dcov$total <- cov(t(rbind(X, Y)))

  # initialize with scalar diagonal noise on the marginals (shared by all features)
  phi.init <- list(X = diag(var(as.vector(X)), Dim$X), Y = diag(var(as.vector(Y)), Dim$Y)) 
  
  # Initialize W's
  W.init   <- list()
  W.init$X <- as.matrix(eigen(Dcov$X)$vectors[, 1:Dim$Z])
  W.init$Y <- as.matrix(eigen(Dcov$Y)$vectors[, 1:Dim$Z])
  W.init$total <- rbind(W.init$X, W.init$Y) # can be possibly removed in some special cases
  
  ##################################################

  # optimize until convergence
  res <- optimize.parameters(W.init, phi.init, Dim, Dcov, priors, marginalCovariances = "full", epsilon)

  W <- list(X = res$W$X, Y = res$W$Y, total = rbind(res$W$X, res$W$Y))
  phi <- res$phi
  
  rownames(W$X) <- rownames(X)
  rownames(W$Y) <- rownames(Y)
  rownames(W$total) <- c(rownames(X), rownames(Y))
      
  #phiX <- diag(res$phi$X, nrow(X))
  rownames(phi$X) <- colnames(phi$X) <- rownames(X)

  #phiY <- diag(res$phi$Y, nrow(Y))
  rownames(phi$Y) <- colnames(phi$Y) <- rownames(Y)

  #phitotal <- diag(c(diag(phiX), diag(phiY)), nrow(X) + nrow(Y))
  rownames(phi$total) <- colnames(phi$total) <- c(rownames(X), rownames(Y))

  return( list(W = W, phi = phi) )

}

