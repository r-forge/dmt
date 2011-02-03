cost.W <-
function (vec, phi, priors, Dim, Dcov, H) {

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
        tcost <- sum((T - H)^2) * priors$T
    
  # -logP for W prior - skip since not used now
  #wcost <- sum((W$X)^2) * priors$W

  cost.data + tcost #+ wcost

}



cost.W.exponential <-
function (vec, phi, priors = NULL, Dim, Dcov) {

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
  cost.data <- log(det(Sigma)) + sum(diag(solve(Sigma)%*%Dcov$total))
  
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



cost5 <-
function (W, phi, Dcov) {

  # SimCCA - just use identical Wx = Wy
  wtw <- W%*%t(W)

  Sigma <- rbind(cbind(wtw + phi$X, wtw),
          cbind(wtw,wtw + phi$Y))

        # Marginal cost for the whole data set
        # integrated over z
        # given parameters W, phi
        # P(X,Y | W, phi) = integral N(X|Wx*z,phix)*N(Y|Wy*z,phiy)*N(z|0,I)
        # We report -logP here

        # restrict solutions to cases where det(Sigma)>=0
  detsigma <- det( Sigma )
  if (detsigma > 0) {
                #print(detsigma)
    log(detsigma) + sum(diag(solve(Sigma)%*%Dcov$total))
  } else { Inf }

}



cost7 <-
function (Wvec, phi, Dcov, Dim, priors) {

  if (!is.null(priors$W)) {
    Wvec <- abs(Wvec)
  }
    
  W <- get.W4(Wvec, Dim)$X
  
  # SimCCA - just use identical Wx = Wy
  wtw <- W%*%t(W)

  Sigma <- rbind(cbind(wtw + phi$X, wtw),
          cbind(wtw,wtw + phi$Y))

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
  vec <- abs(as.vector(W))
  if (!is.null(priors)) {
    #multiply by 2 to count for both wx and wy
    wcost <- -2*sum(dexp(vec, rate = priors$W, log = TRUE))

  } else {
    wcost <- 0 # no effect
  }

  #print(paste("wcost", wcost))
  #print(paste("cost.data", cost.data))
  
  cost.data + wcost

}

