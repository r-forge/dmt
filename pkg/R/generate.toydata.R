generate.toydata <- function (N = 100, zDim = 2, xDim = 3, yDim = 3,
marginal.covariances = "full") {

  Z <- matrix(rnorm(N*zDim), nrow = zDim)

  # So far, have marginal noise dimensionality equal 
  # the whole span of X/Y data
  zxDim <- xDim
  Zx <- matrix(rnorm(N*zxDim), nrow = zxDim)

  zyDim <- yDim
  Zy <- matrix(rnorm(N*zyDim), nrow = zyDim)  

  Wx <- matrix(rnorm(zDim*xDim), nrow = xDim)
  Wy <- matrix(rnorm(zDim*yDim), nrow = yDim)

  if (marginal.covariances == "full") {
    Bx <- matrix(rnorm(zxDim*xDim), nrow = xDim)
    By <- matrix(rnorm(zyDim*yDim), nrow = yDim)
  }
  
  if (marginal.covariances == "diagonal") {
    Bx <- diag(rnorm(xDim))
    By <- diag(rnorm(yDim))
  }  
  
  if (marginal.covariances == "isotropic") {
    Bx <- rnorm(1)*diag(xDim)
    By <- rnorm(1)*diag(yDim)
  }
  
  X <- as.matrix(Wx%*%Z + Bx%*%Zx, nrow = xDim)
  Y <- as.matrix(Wy%*%Z + By%*%Zy, nrow = xDim)
		      
  list(Z = Z, X = X, Y = Y, Wx = Wx, Wy = Wy, Bx = Bx, By = By, Zx = Zx, Zy = Zy)

}

  # Center data? (should be built into the functions)
  #X <- t(apply(X,1,function(x){(x-mean(x))}))
  #Y<-t(apply(Y,1,function(x){(x-mean(x))}))