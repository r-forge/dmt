library('mvtnorm')
library(MASS)

rdirichlet <- function(alpha) {
  dim <- length(alpha)
  y <- rep(0,dim)
  for(i in 1:dim) y[i] = rgamma(1,alpha[i])
  y <- y/sum(y)
}

dinvgamma <- function(x, shape, rate=1, log=TRUE) {
  val <- shape*log(rate) - lgamma(shape) - (shape+1)*log(x) - rate / x

  if(!log)
    val <- exp(val)

  val
}

pinvgamma <- function(x, shape, rate) {
  val <- pgamma(1/x,shape,rate,lower.tail=FALSE,log.p=TRUE)

  val
}

qinvgamma <- function(p, shape, rate, maxIter=20) {
  logp <- log(p)

  low <- 0
  high <- 1

  for(i in 1:maxIter) {
    mid <- (low + high)/2
    if(pinvgamma(mid,shape,rate)<logp) {
      low <- mid
    } else {
      high <- mid
    }
  }
  mid <- (low + high) / 2
}

#
# Changed to mvrnorm because rmvnorm seemed to freeze in some
# cases with certain S matrices
#

rwishart <- function(df, S) {
  #sam <- rmvnorm(df, array(0,nrow(S)), S)
  sam <- mvrnorm(df, array(0,c(1,nrow(S))), S)

  theta <- crossprod(sam)#/df)

  #solve(theta)
}

rinvwishart <- function(df, S) {
  solve(rwishart(df,solve(S)))
}

dinvwishart <- function(W, v, S) {
  k <- nrow(W)

  i <- 1:k

  val1 <- - (v*k/2*log(2) + k*(k-1)/4*log(pi) + sum(lgamma((v+1-i)/2)))
  val2 <- v/2*determinant(S,logarithm=TRUE)$modulus[1] - (v+k+1)/2*determinant(W,logarithm=TRUE)$modulus[1]
  val3 <- -1/2*sum(diag(S %*% solve(W)))

  #print(v)
  #print(c(val1,val2,val3))

  val <- val1+val2+val3
}