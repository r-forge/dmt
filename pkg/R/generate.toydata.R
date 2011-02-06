generate.toydata <- function (N = 100, zDim = 2, xDim = 3, yDim = 3, marginal.covariances = "full", priors = NULL) {

  # FIXME: add Wx = Wy prior option and Wx ~ Wy or Wy = TWx option

  Z <- matrix(rnorm(N*zDim), nrow = zDim)

  # So far, have marginal noise dimensionality equal 
  # the whole span of X/Y data
  zxDim <- xDim
  Zx <- matrix(rnorm(N*zxDim), nrow = zxDim)

  zyDim <- yDim
  Zy <- matrix(rnorm(N*zyDim), nrow = zyDim)  

  if (is.null(priors$W)) {
    Wx <- matrix(rnorm(zDim*xDim), nrow = xDim)
    Wy <- matrix(rnorm(zDim*yDim), nrow = yDim)
  } else if (priors$W > 0) {
    Wx <- matrix(rexp(zDim*xDim, rate = priors$W), nrow = xDim)
    Wy <- matrix(rexp(zDim*yDim, rate = priors$W), nrow = yDim)
  }
  
  if (marginal.covariances == "full") {

    Bx <- matrix(rnorm(zxDim*xDim), nrow = xDim)
    By <- matrix(rnorm(zyDim*yDim), nrow = yDim)

  } else if (marginal.covariances == "diagonal") {

    Bx <- diag(rnorm(xDim))
    By <- diag(rnorm(yDim))

  } else if (marginal.covariances == "isotropic") {

    Bx <- rnorm(1)*diag(xDim)
    By <- rnorm(1)*diag(yDim)

  } else if (marginal.covariances == "identical isotropic") {

    const <- rnorm(1)
    Bx <- const*diag(xDim)
    By <- const*diag(yDim)

  }  

  # Marginal noise. Note: full marginal noise assumed for the time being
  nx <- Bx%*%Zx
  ny <- By%*%Zy    

  X <- as.matrix(Wx%*%Z + nx, nrow = xDim)
  Y <- as.matrix(Wy%*%Z + ny, nrow = yDim)
		      
  list(Z = Z, X = X, Y = Y, Wx = Wx, Wy = Wy, Bx = Bx, By = By, Zx = Zx, Zy = Zy)

}

  
compare.estimate.and.truth <- function (res, toy) {
  
  # res: output from fit.dependency.model
  # toy: toydata object containing the toydata used to train the model 
  #      (including true parameters used to generate the toydata)
  
  wtw.x.estimated <- res@W$X%*%t(res@W$X)
  wtw.x.true <- toy$Wx%*%t(toy$Wx)
  
  wtw.y.estimated <- res@W$Y%*%t(res@W$Y)
  wtw.y.true <- toy$Wy%*%t(toy$Wy)  
  
  phiX.estimated <- res@phi$X
  phiX.true <- toy$Bx%*%t(toy$Bx)

  phiY.estimated <- res@phi$Y
  phiY.true <- toy$By%*%t(toy$By)

  corsx <- cor(as.vector(wtw.x.estimated), as.vector(wtw.x.true))
  corsy <- cor(as.vector(wtw.y.estimated), as.vector(wtw.y.true))
  
  cormx <- cor(as.vector(phiX.estimated), as.vector(phiX.true))
  cormy <- cor(as.vector(phiY.estimated), as.vector(phiY.true))  

  unlist(list(wtw.x = corsx, wtw.y = corsy, phi.x = cormx, phi.y = cormy))

}

