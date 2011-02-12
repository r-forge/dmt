cost.W <- function (vec, phi, priors, Dim, Dcov) {

  # Wx ~ Wy constrained
  # no W prior

  H <- priors$Nm.wxwy.mean

  # Retrieve the actual W and T from the parameter vector
  wt <- get.W(vec, Dim)
  W <- wt$W
  T <- wt$T

  # Marginal cost for the whole data set
  # integrated over z
  # given parameters W, phi
  # P(X,Y | W, phi) = integral N(X|Wx*z,phix)*N(Y|Wy*z,phiy)*N(z|0,I)
  # We report -logP here  
  # Data prob. Taken from probCCA paper, section 4, l1

  wtw.xy <- W$X%*%t(W$Y)
  Sigma <- rbind(cbind(W$X%*%t(W$X) + phi$X*diag(Dim$X), wtw.xy),
          cbind(t(wtw.xy),W$Y%*%t(W$Y) + phi$Y*diag(Dim$Y)))

  # -logP for the data
  cost.data <- log(det(Sigma)) + sum(diag(solve(Sigma)%*%Dcov$total))

  # -logP for T prior
  tcost <- sum((T - H)^2) * priors$T.tmp

  # -logP for W prior - skip since not used now
  # priors$W.tmp <- 1/(2 * Nsamples * priors$W)
  # NOTE considerable speed increase in optimize iteration if 
  # this is calculated outside this function!
  #wcost <- sum((W$X)^2) * priors$W.tmp

  cost.data + tcost #+ wcost

}


cost.W.exponential <- function (vec, phi, priors = NULL, Dim, Dcov) {

  # allows exponential prior for W
  # in general, Wx != Wy

  # remove sign as we assume W always positive here
  vec <- abs(vec)
  
  # Retrieve W from the parameter vector
  W <- get.W2(vec, Dim)

  # Marginal cost for the whole data set
  # integrated over z
  # given parameters W, phi
  # P(X,Y | W, phi) = integral N(X|Wx*z,phix)*N(Y|Wy*z,phiy)*N(z|0,I)
  # report -logP here
  
  # Data prob. Taken from probCCA paper, section 4, l1
  wtw.xy <- W$X%*%t(W$Y)

  Sigma <- rbind(cbind(W$X%*%t(W$X) + phi$X, wtw.xy),
                 cbind(t(wtw.xy), W$Y%*%t(W$Y) + phi$Y))
  
  # -logP for the data
print("HERHE")
print("PHI")
print(phi)
print("WYW")
print(W$Y%*%t(W$Y))
print("WXW")
print(W$X%*%t(W$X))
print("WXYW")
print(wtw.xy)
print("Sigma")
print(Sigma)
print("detSigma")
print(det(Sigma))
  cost.data <- log(det(Sigma)) + sum(diag(solve(Sigma)%*%Dcov$total))
print("HERHE2")  
  
  # -logP for W prior
  # wcost <- sum((W$X)^2) * priors$W
  # Assuming exponential prior distribution with rate parameter priors$W
  wcost <- 0 # no effect
  if (!is.null(priors$W)){
    if (priors$W > 0) {
      wcost <- -sum(dexp(vec, rate = priors$W, log = TRUE))
    }
  }
  
  cost.data + wcost

}

cost7 <- function (Wvec, phi, Dcov, Dim, priors) {

  # SimCCA: identical Wx = Wy
  # allows W prior

  # NOTE: possible to optimize quite much by removing W matrix conversions?

  if ( !is.null(priors$W) ) { Wvec <- abs(Wvec) }    

  W <- get.W4(Wvec, Dim)$X
  wtw <- W%*%t(W)

  Sigma <- rbind(cbind(wtw + phi$X, wtw),
                 cbind(wtw, wtw + phi$Y))

  # Marginal cost for the whole data set
  # integrated over z
  # given parameters W, phi
  # P(X,Y | W, phi) = integral N(X|Wx*z,phix)*N(Y|Wy*z,phiy)*N(z|0,I)
  # We report -logP here

  # restrict solutions to cases where det(Sigma)>=0

  # -logP for the data
  detsigma <- det( Sigma )

  if (detsigma > 0) {
    cost.data <- log(detsigma) + sum(diag(solve(Sigma)%*%Dcov$total))
  } else { cost.data <- Inf }
    
  # -logP for W prior
  # wcost <- sum((W$X)^2) * priors$W
  # Assuming exponential prior distribution with rate parameter priors$W
  wcost <- 0 # no effect by efault
  if ( !is.null(priors$W) ) {
    #multiply by 2 to count for both wx and wy
    wcost <- -2*sum(dexp(Wvec, rate = priors$W, log = TRUE))
  } 

  #print(paste("wcost", wcost))
  #print(paste("cost.data", cost.data))
  
  cost.data + wcost

}
