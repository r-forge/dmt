# (C) 2008-2011 Leo Lahti and Olli-Pekka Huovilainen          
# All rights reserved. 
# FreeBSD License (keep this notice)     


# "The important thing in science is not so much to obtain new facts as
#  to discover new ways of thinking about them."
#  - Sir William Bragg 


dependency.score <- function (model) {

  W <- model$W
  phi <- model$phi

  if (!is.null(W$X)){
    # this equals to the trace of the full Phi
    noise <- sum(diag(as.matrix(phi$X))) + sum(diag(as.matrix(phi$Y))) 

    wtw <- rbind(W$X,W$Y)%*%t(rbind(W$X,W$Y))
  }
  else {
    # For single data case, check proportion between 
    # latent covariance and expected diagonal noise
    noise <- sum(diag(as.matrix(phi$total)))
    wtw <- W$total %*% t(W$total)
  }

  signal <- sum(diag(wtw)) # trace of full WWt covariance 

  cost <- signal/noise
  cost
}

