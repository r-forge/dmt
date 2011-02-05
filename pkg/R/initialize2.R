initialize2 <- function (X, Y, zDimension = NULL, marginalCovariances) {

  zDimension <- ifelse(is.null(zDimension), min(nrow(X), nrow(Y)), zDimension)



  Nsamples <- ncol(X)
  Dim <- list(X = nrow(X), Y = nrow(Y), Z = zDimension)
  nullmat  <- matrix(0, nrow = Dim$X, ncol = Dim$Y)

  if (marginalCovariances == "isotropic") {
    # Scalar values
    phi$X <- var(as.vector(X))
    phi$Y <- var(as.vector(Y))
    phi$total <- c(phi$X, phi$Y)
  } else if (marginalCovariances == "identical isotropic") {
    # Scalar values phix = phiy
    phi$X <- phi$Y <- var(c(as.vector(X), as.vector(Y)))
    phi$total <- c(phi$X, phi$Y)
  } else {
    # diagonal matrices
    # initialize with scalar diagonal noise on the marginals (shared by all features)
    phi <- list(X = diag(var(as.vector(X)), Dim$X), 
                Y = diag(var(as.vector(Y)), Dim$Y))  
    phi$total <- rbind(cbind(phi$X,nullmat), cbind(nullmat, phi$Y))

  }

  # FIXME: if phi$Y is scalar (as in segmented/mir case) we can speed up here. Do later.
  phi.inv  <- list()
  phi.inv$X <- solve(phi$X)
  phi.inv$Y <- solve(phi$Y)
  phi.inv$total <- rbind(cbind(phi.inv$X, nullmat), cbind(t(nullmat), phi.inv$Y))

  Dcov <- list()
  Dcov$X <- cov(t(X), use = "pairwise.complete.obs")
  Dcov$Y <- cov(t(Y), use = "pairwise.complete.obs")
  Dcov$xy <- cov(t(X), t(Y), use = "pairwise.complete.obs")
  Dcov$yx <- t(Dcov$xy)
  Dcov$total <- rbind(cbind(Dcov$X, Dcov$xy), cbind(Dcov$yx, Dcov$Y))
  Dcov$sum   <- Dcov$X + Dcov$Y + Dcov$xy + Dcov$yx
  Dcov$sum <- cov(t(X + Y), use = "pairwise.complete.obs")
  
  # It is possible that covariances calculated with pairwise complete
  # observations are not positive semi-definite.
  # Check this. If not pos.sem.def, then replace with the closest
  # pos. semidefinite matrix
  if (any(eigen(Dcov$total)$values < 0)) {
    #message("Covariance approximation used to avoid numerical instability.")
    Dcov$X  <- as.matrix(nearPD(Dcov$X)$mat)
    Dcov$Y  <- as.matrix(nearPD(Dcov$Y)$mat)
    Dcov$total <- rbind(cbind(Dcov$X, Dcov$xy), cbind(Dcov$yx, Dcov$Y))
    Dcov$sum   <- Dcov$X + Dcov$Y + Dcov$xy + Dcov$yx
  }

  # Initialize W's
  W <- list()
  W$X <- as.matrix(eigen(Dcov$X)$vectors[, 1:Dim$Z])
  W$Y <- as.matrix(eigen(Dcov$Y)$vectors[, 1:Dim$Z])	
  W$total <- rbind(W$X, W$Y) 
  
  list(phi = phi, W = W, Dcov = Dcov, Dim = Dim, nullmat = nullmat, Nsamples = Nsamples)

}

