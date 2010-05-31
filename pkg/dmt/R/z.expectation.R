z.expectation <- function (model, X, Y = NULL) {

  W <- getW(model)
  phi <- getPhi(model)

    if (!is.null(Y)) {
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
