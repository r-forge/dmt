# (C) 2008-2011 Leo Lahti and Olli-Pekka Huovilainen          
# All rights reserved. 
# FreeBSD License (keep this notice)     


# "Reality is merely an illusion, albeit a very persistent one."
# - Albert Einstein
        

get.W <- function (vec, Dim) {

  # Convert parameter vector into matrices Wx and Wy
  W <- list()
  W$X <- array(vec[1:(Dim$X*Dim$Z)], dim = c(Dim$X,Dim$Z))
  T <- array(vec[-seq(Dim$X*Dim$Z)], dim = c(Dim$Y,Dim$X))
  W$Y <- T%*%W$X

  list(W = W, T = T)
}

get.W2 <- function (vec, Dim) {

  # Convert parameter vector into matrices Wx and Wy

  # Note that we always assume that W is positive Therefore remove the
  #sign here to speed up optimization
  vec <- abs(vec)
  
  # NOTE: assumes dx = dy  
  W <- list()
  nx <- (Dim$X*Dim$Z)
  W$X <- array(vec[1:nx], dim = c(Dim$X, Dim$Z))
  W$Y <- array(vec[(nx + 1):length(vec)], dim = c(Dim$Y, Dim$Z))
  W$total <- rbind(W$X, W$Y)
  
  W

}

get.W4 <- function (vec, Dim) {

  # Convert parameter vector into matrices Wx and Wy
  # assumes Wx = Wy

  # NOTE: assumes dx = dy  
  W <- list()
  W$X <- W$Y <- array(abs(vec), dim = c(Dim$X, Dim$Z))
  # Note that we always assume that W is positive
  # Therefore remove the sign here to speed up optimization
  # FIXME: later, allow individual scales for Wx, Wy?
  
  W

}

