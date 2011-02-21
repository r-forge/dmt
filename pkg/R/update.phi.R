# (C) 2008-2011 Leo Lahti and Olli-Pekka Huovilainen          
# All rights reserved. 
# FreeBSD License (keep this notice)     


# "In science the credit goes to the man who convinces the world, not the
#  man to whom the idea first occurs."
# - Sir Francis Darwin

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

  # Calculate updated phi (= epsilon);
  # add small constant to avoid singularity
  max(sum(diag(Xcov - wtw%*%M%*%Xcov))/dx, 0) + 1e-3

}

phi.pca <- function (dat, zDim) {

  # isotropic phi for probabilistic PCA
  # with given latent variable dimensionality;
  # note that the isotropic phi explains everything that
  # cannot be explained by the latent variable:
  # XXt = WWt + phi

  # With concatenated data:
  #dat <- rbind(X, Y)

  # eigenvalues D and eigenvectors U             
  duv <- svd(dat)       
  U <- duv$u
  D <- sort(duv$d, decreasing = TRUE)
 
  sum(D[-seq(zDim)])/(nrow(duv$u) - zDim)
	             
}

phi.EM.simcca <- function (Dcov, W.new, phi.inv, W.old, M) {

  # assuming Wx = W

  # modifiying BachJordan sec. 4.1
  dx <- ncol(Dcov$X)
 
  # Reduces to this when Wx = Wy:
  mat <- W.old$X%*%M%*%t(W.new$X)
  nullmat <- matrix(0, nrow = dx, ncol = dx)
  mat2 <- Dcov$total - Dcov$total%*%rbind(cbind(phi.inv$X%*%mat, nullmat), cbind(nullmat, phi.inv$Y%*%mat))

  # Regularize diagonal to avoid singularity
  phi <- list()
  phi$total <- mat2 + (1e-2)*diag(nrow(mat2))  
  phi$X <- matrix(phi$total[1:dx,1:dx], dx)
  phi$Y <- matrix(phi$total[(dx+1):(2*dx), (dx+1):(2*dx)], dx)

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
    




