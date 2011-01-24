z.expectation <- function (model, X, Y = NULL) {

  W <- getW(model)
  phi <- getPhi(model)

  # Center data
  X <- t(centerData(t(X), rm.na = TRUE))
  
  if (!is.null(Y)){

    Y <- t(centerData(t(Y), rm.na = TRUE))
    
    S <- solve(t(W$X)%*%solve(phi$X)%*%W$X + t(W$Y)%*%solve(phi$Y)%*%W$Y + diag(ncol(W$X)))

    return(S%*%(t(W$X)%*%solve(phi$X)%*%X + t(W$Y)%*%solve(phi$Y)%*%Y))
  }
  else {
    W <- W$total
    phi <- phi$total
    S <- solve(t(W) %*% solve(phi) %*% W + diag(ncol(W)))
    return(S %*% (t(W) %*% solve(phi) %*% X))
  }
}