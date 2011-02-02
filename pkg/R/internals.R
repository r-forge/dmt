



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
	
solve.w <- function (Xc, Yc, Cxx, Cyy, dz = NULL) {

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

  if (dx < ncol(Xc) || dy < ncol(Yc))
    stop("Unable to calculate the pCCA model; the sample covariance matrix is not invertible.\nMake sure that you are not using segmented data")

  z <- svd(qr.qty(qx, qr.qy(qy, diag(1, nr, dy)))[1L:dx, , drop = FALSE], dx, dy)

  xcoef <- backsolve((qx$qr)[1L:dx, 1L:dx, drop = FALSE], z$u)
  ycoef <- backsolve((qy$qr)[1L:dy, 1L:dy, drop = FALSE], z$v)

  #rownames(xcoef) <- colnames(Xc)[qx$pivot][1L:dx]
  #rownames(ycoef) <- colnames(Yc)[qy$pivot][1L:dy]
  #cca <- list(cor = z$d, xcoef = xcoef, ycoef = ycoef)

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


matrix.sqrt <-
function (A) {
	#Solve square root (A^(1/2)) of a diagonalizable nxn-matrix A
	#First diagonalize A, then compute sqrt for the diagonal elements
	#of the diagonalized matrix Adot. Then compute matrix sqrt from these
	#Eigenvectors of A

	#Author: Leo Lahti
	#Date: 13.2.2007
	#Comment: I did not easily find suitable function in R although there might well be a better optimized one. 

	e<-eigen(A) 
	V<-e$vector
	#try inverse the eigvec matrix. 
	#if not diagonalizable then some error will occur from try I think
	Vinv<-try(solve(V)) 
	D <- diag(e$values)

	V%*%sqrt(D)%*%Vinv
}



 # the inverse of square root of a square matrix

"SqrtInvMat" <- function(matrix)
{

 x <- as.matrix(matrix)

 fac <- svd(x)

 res <-(fac$v %*% (diag(1/sqrt(fac$d),nrow = length(fac$d))) %*% t(fac$u))

 return(res)

}

