# (C) 2008-2011 Leo Lahti and Olli-Pekka Huovilainen          
# All rights reserved. 
# FreeBSD License (keep this notice)     


# "In science the credit goes to the man who convinces the world, not the
#  man to whom the idea first occurs."
# - Sir Francis Darwin



update.phi <- function (Dcov, M, beta, W, phi) {

  # Empirical estimate. Should give similar results than
  # update.phi.isotropic. Compare and combine.
  # FIXME: remove M, not neeed here

  # simple marginal covariance: sigma2 * I
  # i.e. just one parameter; in general phix != phiy

  # Mean equals to dividing with dimension
  # add small constant to avoid singularities
  phix <- max(mean(diag(Dcov$X - Dcov$X%*%t(beta$X)%*%t(W$X))), 0) + 1e-3
  phiy <- max(mean(diag(Dcov$Y - Dcov$Y%*%t(beta$Y)%*%t(W$Y))), 0) + 1e-3

  # Return phi
  list(X = phix, Y = phiy)

}

update.phi.isotropic <- function (Xcov, W, epsilon, dx) {

  # used to update phix and phiy, one at a time

  # ML estimate of the variance given W
  # See Roweis, 'EM Algorithms for PCA and SPCA'
  # section 4 (SPCA; M-step for epsilon)

  # epsilon is the variance from previous
  # iteration step
  
  # auxiliary variables
  wtw <- W%*%t(W)
  M <- set.M.isotropic(wtw, epsilon, dx)

  # Calculate updated phi (= epsilon) and return
  # add small constant to avoid singularity
  max(sum(diag(Xcov - wtw%*%M%*%Xcov))/dx, 0) + 1e-3

}


phi.EM.simcca <- function (Dcov, W.new, phi.inv, W.old, M) {

  # assuming Wx = W

  # modifiying BachJordan sec. 4.1
  dx <- ncol(Dcov$X)
 
  # Reduces to this when Wx = Wy:
  mat <- W.old$X%*%M%*%t(W.new$X)
  nullmat <- matrix(0, nrow = dx, ncol = dx)
  mat2 <- Dcov$total - Dcov$total%*%rbind(cbind(phi.inv$X%*%mat,nullmat), cbind(nullmat,phi.inv$Y%*%mat))

  # Diagonal is regularized to avoid singluar matrix
  phi <- list()
  phi$total <- mat2 + (1e-2)*diag(nrow(mat2))  
  phi$X <- matrix(phi$total[1:dx,1:dx], dx)
  phi$Y <- matrix(phi$total[(dx+1):(2*dx), (dx+1):(2*dx)], dx)
  #phi$Y <- matrix(phi$total[-(1:dx), -(1:dx)], dy)

  phi

}

phi.EM.cca <- function (Dcov, W.new, phi.inv, W.old, M, nullmat) {

  # for general Wx != Wy
  # with full covariances

  # From BachJordan sec. 4.1
  dx <- ncol(Dcov$X)
  dy <- ncol(Dcov$Y)  

  # From BachJordan sec. 4.1	
  mat <- W.old$total%*%M%*%t(W.new$total)  

  ## FIXME: check below for speedup; confirm that gives similar results and add
  ## (should work if M = [M0; M0]?)
  ## From BachJordan sec. 4.1	
  ## vrt. beta.fullcov: M%*%t(W)%*%phi.inv
  ## phi.inv$X%*%W.old$X%*%M%*%t(W.new$X) = t(beta.fullcov)%*%t(W.new$X) tms.
  #mat.x <- phi.inv$X%*%mat
  #mat.y <- phi.inv$Y%*%mat
  #mat2  <- Dcov$total - Dcov$total%*%rbind(cbind(mat.x, mat.x), cbind(mat.y, mat.y))
   mat2  <- Dcov$total - Dcov$total%*%phi.inv$total%*%mat
  
  # Diagonal is regularized to avoid singluar matrix	  
  phi <- list()
  phi$total <- mat2 + (1e-2)*diag(nrow(mat2))
  phi$X <- matrix(phi$total[1:dx,1:dx], dx)
  phi$Y <- matrix(phi$total[-(1:dx), -(1:dx)], dy)
  
  phi
   
}
    


