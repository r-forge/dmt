
# (C) 2008-2011 Leo Lahti and Olli-Pekka Huovilainen          
# All rights reserved. 
# FreeBSD License (keep this notice)     






check.data <- function (X, Y, zDimension = NULL) {

  if (is.null(zDimension)) { 
    # No latent dimensionality specified: using full latent dim
    zDimension <- min(nrow(X), nrow(Y))   
  } else if (zDimension > min(nrow(X), nrow(Y))) {
    warning("Latent variable dimensionality cannot exceed data dimensionality. Setting zDimension = min(nrow(X), nrow(Y))")
    zDimension <- min(nrow(X), nrow(Y))
  }
        
  # Check that data is given as a matrix
  if ( !is.matrix(X) ){
    stop("X data needs to be a matrix")
  }
	    
  # Check that data is given as a matrix
  if ( !is.null(Y) && !is.matrix(Y) ){
    stop("Y data needs to be a matrix")
  }

  #cat("Centering the data..\n")
  X <- t(centerData(t(X), rm.na = TRUE))		        
  if (!is.null(Y)) {			  
    Y <- t(centerData(t(Y), rm.na = TRUE))			      
  }			        				
      
  list(X = X, Y = Y, zDimension = zDimension)

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



SqrtInvMat <- function( matrix ) {

  # Calculates square root of inverse of an (invertible) square matrix.

  x <- as.matrix(matrix)

  fac <- svd(x)

  fac$v %*% (diag(1/sqrt(fac$d),nrow = length(fac$d))) %*% t(fac$u)

}




#internal subroutine extra

concatenate <- function(datasets)
{

   mat <- datasets #list of data sets

  m <- length(mat)

  com <- mat[[1]]

  if(m > 1)
  {
       for(i in 2:m)
       {
         com <- cbind(com,mat[[i]])
       }
  }


  return(com)
}
