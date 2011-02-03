pfa <- function (X, Y = NULL, zDimension = NULL, includeData = TRUE, calculateZ = TRUE) {

  # Probabilistic factorial analysis model as proposed in
  # EM Algorithms for ML Factoral Analysis, Rubin D. and 
  # Thayer D. 1982

  # Assumption: R = I

  dat <- check.data(X, Y, zDimension)
  X <- dat$X
  Y <- dat$Y
  zDimension <- dat$zDimension

  res <- calc.pfa(X, Y, zDimension)

  method <- "pFA"
  params <- list(marginalCovariances = "diagonal", zDimension = zDimension)
  score <- dependency.score( res )  
  model <- new("DependencyModel", W = res$W, phi = res$phi, score = score, method = method, params = params)
  if ( includeData ) model@data <- list(X = X, Y = Y)
  if ( calculateZ )  model@z <- z.expectation(model, X, Y)
  model

}

calc.pfa <- function (X, Y, zDimension) {

  ##################################################
  
  # Y.rubin is Y in (Rubin & Thayer, 1982)
  # Variables on columns and samples on rows
  if (is.null(Y)){
    Y.rubin <- t(X)
    # Factor loading matrix
    beta <- t(eigen(cov(t(X)))$vectors[, 1:zDimension])
  }
  else {
    Y.rubin <- cbind(t(X),t(Y))
    # Use different initialization for beta when data has inequal dimensionalities
    if (nrow(X) != nrow(Y)) {
      beta <- t(eigen(cov(Y.rubin))$vectors[,1:zDimension])
    }
    else {
      init <- initialize2(X,Y)
      # Factor loading matrix
      beta <- t(init$W$total[,1:zDimension])
    }
  }
  epsilon <- 1e-3
  colnames(beta) <- colnames(Y.rubin)
  tau2 <- diag(ncol(Y.rubin))

  Cyy <- cov(Y.rubin)
  delta <- 1e12
  # EM
  while(delta > epsilon){

    beta.old <- beta
    tau2.old <- tau2

    # E step
    invtau2 <- solve(tau2)
    tbb <- invtau2 - (invtau2%*%t(beta))%*%solve(diag(zDimension) + beta%*%invtau2%*%t(beta))%*%(beta%*%invtau2)

    d <- tbb%*%t(beta)
    D <- diag(zDimension) - beta%*%d
	
    # M step
    beta <- solve(t(d)%*%Cyy%*%d+D)%*%t(Cyy%*%d)
    tau2 <- diag(diag(Cyy - Cyy%*%d%*%solve(t(d)%*%Cyy%*%d + D)%*%t(Cyy%*%d)))
    delta <- max(sum(abs(tau2-tau2.old)),sum(abs(beta-beta.old)))
  }

  # Convert names as same in other methods
  if (is.null(Y)){
    W <- list(total = t(beta))
    phi <- list(total = tau2)
  }
  else {
    W <- list(X = as.matrix(t(beta)[(1:nrow(X)),]), Y = as.matrix(t(beta)[-(1:nrow(X)),]), total = t(beta))
    phi <- list(X = tau2[1:nrow(X),1:nrow(X)], Y = tau2[-(1:nrow(X)),-(1:nrow(X))], total = tau2)                
  }
  
  res <- list(W = W, phi = phi)

  res
}

