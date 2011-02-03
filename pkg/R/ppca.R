ppca <- function (X, Y = NULL, zDimension = NULL, includeData = TRUE, calculateZ = TRUE) {

  dat <- check.data(X, Y, zDimension)
  X <- dat$X  
  Y <- dat$Y
  zDimension <- dat$zDimension

  res <- calc.ppca(X, Y, zDimension)
  
  method <- "pPCA"	
  params <- list(marginalCovariances = "isotropic", zDimension = zDimension)
  score <- dependency.score( res )
  model <- new("DependencyModel", W = res$W, phi = res$phi, score = score, method = method, params = params)
  if (includeData) { model@data <- list(X = X, Y = Y) }
  if (calculateZ) { model@z <- z.expectation(model, X, Y) }
  model
 	    
}


# FIXME: now calc.ppca used only in one-data case by function ppca
# either remove two-data case, or test and compare with fit.dependency.model and take into use
calc.ppca <- function (X, Y, zDimension) {

  # Replaces function solve.CCA

  # if zDimension = NULL then full-rank solution is computed
  
  # Probabilistic PCA
  # (See Tipping and Bishop 1999)

  # ML estimates W, sigma for probability model
  # X ~ N(Wz, sigma*I)
  # i.e. latent variable model with isotropic noise

  # If only X is given in the argument, compute
  # pPCA for X

  # If both X and Y are given in the argument, compute
  # pPCA for concatenated [X; Y]
  # Assuming isotropic and identical marginal noise, the
  # principal subspace will capture the dependencies between X and Y.
  # This corresponds to the model (sigmax = sigmay = sigma)
  # X ~ N(Wx*z, sigma*I); Y ~ N(Wy*z, sigma*I)
  # This provides a simple comparison method for more
  # detailed dependency models.

  if (is.null(Y)) {
    res <- ppca.calculate(X, zDimension)
    phi <- list(total = diag(res$phi, nrow(X)))
    rownames(res$W) <- rownames(X)
    colnames(phi$total) <- rownames(phi$total) <- rownames(X)
    #W <- list(X = res$W, total = res$W)
    W <- list(total = res$W)    
  } else {
    # If second argument (Y) given, compute
    # pPCA with two concatenated data sets
    res <- ppca.calculate(rbind(X,Y), zDimension)

    # Make phi diagonal matrix
    phitotal <- diag(res$phi,(nrow(X) + nrow(Y)))

    # Variable names to W and phi
    rownames(res$W) <- c(rownames(X),rownames(Y))
    rownames(phitotal) <- c(rownames(X),rownames(Y))
    colnames(phitotal) <- rownames(phitotal)

    # Divide W and phi to X and Y parts
    phi <- list(X = phitotal[(1:nrow(X)),(1:nrow(X))], 
                Y = phitotal[-(1:nrow(X)),-(1:nrow(X))],
    	        total = phitotal)
    W <- list(X = as.matrix(res$W[(1:nrow(X)),]), 
              Y = as.matrix(res$W[-(1:nrow(X)),]), 
	      total = res$W)

  }
  # Note that if X, Y given then phi$X = phi$Y in the pCCA model
  # Here W corresponds to W$total of other functions when X, Y both given
  list(W=W, phi=phi)
}


ppca.calculate <- function (X, zDimension) {

 # FIXME: ensure that X is zero-mean

  # Probabilistic PCA
  # (See Tipping and Bishop 1999 / 3.2)

  # ML estimates W, sigma for probability model
  # X ~ N(Wz, sigma*I)
  # i.e. latent variable model with isotropic noise

  # Use full-rank if dimensionality is not specified
  zDimension <- ifelse(is.null(zDimension), nrow(X), zDimension)

  # eigenvalues D and eigenvectors U
    duv <- svd(X)
    U <- duv$u
    D <- sort(duv$d, decreasing=TRUE)

    # ML estimate of the variance given dimensionality zDimension for the latent
    # variable z in model X ~ Wz + N(0,sigma) where sigma is isotropic
    d <- nrow(U)
    phi <- sum(D[-seq(zDimension)])/(d-zDimension) 

    # ML estimate for W given variance
    # Here set R <- I (R is an arbitrary orthogonal rotation matrix)
    
    W <- as.matrix(U[,1:zDimension])%*%sqrt(diag(D)[seq(zDimension),seq(zDimension)]-phi*diag(1,zDimension,zDimension))


  if (zDimension == nrow(X)) {

      # If W is full-rank then the isotropic error term disappears assuming data X is gaussian
      # then X%*%t(X) i.e. cov(t(X)) approximates W%*%t(W) (since X ~ Wz and z ~ N(0,I))

    # Note rotational ambiguity for W, Z 
    cat("Full-rank PCA calculated, isotropic error term is zero. Consider checking the principal components.\n")
    W <- matrix.sqrt(cov(t(X)))
    phi <- 0
  }
  
    list(W = W, phi = phi)
}
