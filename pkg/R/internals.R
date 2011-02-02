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
	



