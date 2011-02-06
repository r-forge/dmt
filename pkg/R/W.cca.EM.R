################################

# Functions for solving W, given phi and data
# FIXME: Compare performance, merge?

################################

W.cca.EM <- function (Dcov, M, beta) {
  
  # EM update for W when phi is assumed known
  # Return total W, i.e. [Wx; Wy]

  #beta <-M%*%t(W$total)%*%phi.inv$total
  #M <- solve(t(W)%*%phi.inv%*%W + I)
  Dcov$total%*%t(beta)%*%solve(M + beta%*%Dcov$total%*%t(beta))

}


#######################################

W.simcca.EM <- function (W, phi, Dim, Dcov) {

  # CCA update for W, assuming Wx = Wy
  # see Bach-Jordan 2005, sec. 4.1 for details
  # equations modified from there to match Wx = Wy case
  # FIXME: speedup by sharnig M/beta with phi updates as in in W.cca.EM?
  # (phi.EM.simcca)
  
  what <- 2*W$X
  phihat.inv <- solve(phi$X + phi$Y)
  M.w <- solve(t(what)%*%phihat.inv%*%what + diag(Dim$Z))
  beta.w <- M.w%*%t(what)%*%phihat.inv
  what <- Dcov$sum%*%t(beta.w)%*%solve(M.w + beta.w%*%Dcov$sum%*%t(beta.w))
  w <- what/2
  W$X <- W$Y <- w 
  W$total <- rbind(w, w)		        
			
  ## alternative update; test, compare
  #beta <-M%*%t(W$total)%*%phi.inv$total
  #M <- solve(t(W)%*%phi.inv%*%W + I)  
  #W$total <- W.cca.EM(Dcov, M, beta)
  #W$X = W$total[1:Dim$X,]
  #W$Y = W$total[-c(1:Dim$X),]
  ## robust solution: take average of the two estimates
  #W$X = W$Y = (W$X+W$Y)/2
  #W$total = rbind(W$X,W$Y)

  W

}

					
#######################################
	
solve.w <- function (Xc, Yc, Cxx, Cyy, dz = NULL) {

  # FIXME: compare with the other W updates, e.g. W.cca.EM

  # assumes Xc, Yc : samples x features, zero-mean features
  # Cxx and Cyy are covariances from cov(Xc) and cov(Yc)
  # dz shows the desired rank of latent Z

  # NOTE: here the dimensions of Xc and Yc do not need to match
  # Note: in previous solve.w the input data was features x samples

  # Traditional CCA solution (modified from cancor function):
  nr <- nrow(Xc)

  qx <- qr(Xc)
  qy <- qr(Yc)
  dx <- qx$rank
  dy <- qy$rank

  if (dx < ncol(Xc) || dy < ncol(Yc)) {
    stop("Unable to calculate the pCCA model; the sample covariance matrix is not invertible.")
  }
  
  z <- svd(qr.qty(qx, qr.qy(qy, diag(1, nr, dy)))[1L:dx, , drop = FALSE], dx, dy)

  xcoef <- backsolve((qx$qr)[1L:dx, 1L:dx, drop = FALSE], z$u)
  ycoef <- backsolve((qy$qr)[1L:dy, 1L:dy, drop = FALSE], z$v)

  # Solve W using Archambeau06 equations
  
  # Note: only requirement for Q is that Qx%*%t(Qy) = canonical correlations
  # Q corresponds to M (dz x dz) in Bach-Jordan 2005, p.8 (before sec 4.1)
  Qx <- diag(z$d[1:dz],dz,dz) # dz x dz matrix
  #Qy <- diag(1, nrow(Qx)) # also a dz x dz matrix: identity matrix -> omit

  # ML estimates for the prob. model W:
  dz <- ifelse(is.null(dz), length(z$d), dz)
  Wx <- Cxx%*%xcoef[,1:dz]%*%Qx
  Wy <- Cyy%*%ycoef[,1:dz]#%*%Qy # Qy is identity matrix -> omit

  list(X = Wx, Y = Wy)
}

solve.archambeau <- function (X, Y, Wx, Wy, btb.x, btb.y) {

  # Use the trick introduced in Archambeau et al., ICML 2006. Robust
  # probabilistic projections, and the later correction appendix.

  # This is used to retrieve the original principal components from
  # the probabilistic CCA model solution i.e. get round the rotational
  # invariance problem

  # NOTE: Archambeau uses B in a different meaning than Bach-Jordan.                                  
  # Here denote his B with B.arch  
	
  # btb = B%*%t(B)   

  Nd <- nrow(Wx)

  Bx.arch <- Wx%*%solve(btb.x)%*%t(Wx) + diag(Nd)
  By.arch <- Wy%*%solve(btb.y)%*%t(Wy) + diag(Nd)
				
  # R is a rotation matrix given by                                                                     
  R <- eigen((diag(Nd) - solve(Bx.arch))%*%(diag(Nd) - solve(By.arch)))$vector

  Ux <- solve(cov(t(X)))%*%Wx%*%solve(matrix.sqrt((diag(Nd)-solve(Bx.arch))))%*%R
  Uy <- solve(cov(t(Y)))%*%Wy%*%solve(matrix.sqrt((diag(Nd)-solve(By.arch))))%*%R

  list(X = Ux, Y = Uy)

}