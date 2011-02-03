pcca <- function (X, Y, zDimension = NULL, includeData = TRUE, calculateZ = TRUE) {

  # (C) 2008-2011 Olli-Pekka Huovilainen and Leo Lahti
  # License: FreeBSD (keep this notice)

  # replaces: solve.CCA.full
  
  # If zDimension given, then
  # only pick zDimension first principal components
  # and estimate marginals accordingly
  # relies on the fact that the principal components
  # can be estimated consecutively in pCCA
  
  # Add here centering of the data matrices X, Y
  # (center the dimensions to 0)

  dat <- check.data(X, Y, zDimension)
  X <- dat$X
  Y <- dat$Y
  zDimension <- dat$zDimension

  res <- calc.pcca(X, Y, zDimension)

  method <- "pCCA"
  params <- list(marginalCovariances = "full", zDimension = zDimension)
  score <- dependency.score( res )
  model <- new("DependencyModel", W = res$W, phi = res$phi, score = score, method = method, params = params)
  if ( includeData ) model@data <- list(X = X, Y = Y)
  if ( calculateZ )  model@z <- z.expectation(model, X, Y) 
  model
  
}


calc.pcca <- function (X, Y, zDimension) {

  Dcov <- list()
  Dcov$X <- cov(t(X))
  Dcov$Y <- cov(t(Y))

  # Solve W (solve.w utilizes Archambeau06 equations from PCA for shortcut)
  # FIXME: compare accuracy and speed to direct EM update scheme?
  W <- solve.w(t(X), t(Y), Dcov$X, Dcov$Y, zDimension)

  # Then take only the zDimension first components if defined
  W$X <- as.matrix(W$X)
  W$Y <- as.matrix(W$Y)
  W$total <- rbind(W$X, W$Y)
        
  # estimate
  phi <- list()
  phi$X <- Dcov$X - W$X%*%t(W$X)
  phi$Y <- Dcov$Y - W$Y%*%t(W$Y)

  # Retrieve principal canonical components from the prob.CCA model
  # assuming that W and phi are known (at least in the full-rank case)
  # U <- solve.archambeau(X, Y, W$X, W$Y, phi$X, phi$Y)

  list(W = W, phi = phi)
  
}









