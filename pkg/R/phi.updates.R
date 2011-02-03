update.phi <-
function (Dcov, M, beta, W, phi) {

  # This assumes simple marginal covariance: sigma2 * I
  # i.e. just one parameter; in general phix != phiy

  # Mean equals to dividing with dimension
  phix <- mean(diag(Dcov$X - Dcov$X%*%t(beta$X)%*%t(W$X)))
  phiy <- mean(diag(Dcov$Y - Dcov$Y%*%t(beta$Y)%*%t(W$Y)))

  # Return phi
  list(X = phix, Y = phiy)

}


update.phi.isotropic <- function (Xcov, W, epsilon, dx) {

  # apparently assuming phix = phiy

  # ML estimate of the variance given W
  # See Roweis, 'EM Algorithms for PCA and SPCA'
  # section 4 (SPCA; M-step for epsilon)

  # epsilon is the variance from previous
  # iteration step
  
  # auxiliary variables
  wtw <- W%*%t(W)
  M <- set.M.isotropic(wtw, epsilon, dx)

  # Calculate updated phi (= epsilon) and return
  sum(diag(Xcov - wtw%*%M%*%Xcov))/dx

}


phi.EM.simcca <- function (Dcov, W.new, phi.inv, W.old, M) {

  # assuming Wx = Wy
  
  # From BachJordan sec. 4.1
  dx <- ncol(Dcov$X)

  # Reduces to this when Wx = Wy:
  mat <- W.old$X%*%M%*%t(W.new$X)
  nullmat <- matrix(0, nrow = dx, ncol = dx)
  mat2 <- Dcov$total - Dcov$total%*%rbind(cbind(phi.inv$X%*%mat,nullmat), cbind(nullmat,phi.inv$Y%*%mat))

  # ensure that diagonals are nonnegative by adding a small positive constant
  dd <- diag(mat2)
  dd[dd < 1e-6] <- 1e-6
  diag(mat2) <- dd

  # Diagonal is regularized to avoid singluar matrix
  mat2 <- mat2 + (1e-2)*diag(ncol(Dcov$total))

  phi <- list()
  phi$total <- mat2
  phi$X <- matrix(phi$total[1:dx,1:dx], dx)
  phi$Y <- matrix(phi$total[(dx+1):(2*dx), (dx+1):(2*dx)], dx)

  phi

}

update.phi.EM.fullcov <- function (Dcov, W.new, phi.inv, W.old, M, nullmat) {

  # From BachJordan sec. 4.1	
  mat2  <- Dcov$total - Dcov$total%*%phi.inv$total%*%W.old$total%*%M%*%t(W.new$total)

  ## FIXME: check below for speedup; confirm that gives similar results and add
  ## From BachJordan sec. 4.1	
  ## vrt. beta.fullcov: M%*%t(W)%*%phi.inv
  ## phi.inv$X%*%W.old$X%*%M%*%t(W.new$X) = t(beta.fullcov)%*%t(W.new$X) tms.
  #mat   <- W.old$total%*%M%*%t(W.new$total)
  #mat.x <- phi.inv$X%*%mat
  #mat.y <- phi.inv$Y%*%mat
  #foo   <- rbind(cbind(mat.x, mat.x), cbind(mat.y, mat.y))
  #mat2  <- Dcov$total - Dcov$total%*%foo

  # ensure that diagonals are nonnegative by adding a small positive constant
  dd <- diag(mat2)
  dd[dd < 1e-6] <- 1e-6
  diag(mat2) <- dd

  # for large dimensionality, regularize more the diagonal
  # the more dimensions, more regularization
  if (ncol(Dcov$X) > 5) {mat2 <- mat2 + (1e-2)*diag(ncol(Dcov$total))}

  phi <- list()	
  phi$total <- mat2
  phi$X <- as.matrix(phi$total[1:nrow(Dcov$X),1:nrow(Dcov$X)], nrow = nrow(Dcov$X))
  phi$Y <- as.matrix(phi$total[-seq(nrow(Dcov$X)), -seq(nrow(Dcov$X))], nrow = nrow(Dcov$Y))

  phi

}


