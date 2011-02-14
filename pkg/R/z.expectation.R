
# (C) 2008-2011 Leo Lahti and Olli-Pekka Huovilainen          
# All rights reserved. 
# FreeBSD License (keep this notice)     

# "I am among those who think that science has great beauty. A scientist
# in his laboratory is not only a technician: he is also a child placed
# before natural phenomena which impress him like a fairy tale."
# - Marie Curie


z.expectation <- function (model, X, Y = NULL) {

  W <- getW(model)
  phi <- getPhi(model)

  # Center data
  X <- t(centerData(t(X), rm.na = TRUE))
  
  if (!is.null(Y)){

    Y <- t(centerData(t(Y), rm.na = TRUE))

    # FIXME: use the set.M functions here?
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