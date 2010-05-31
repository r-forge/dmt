z.effects <- function(model,X,Y = NULL){

  W <- getW(model)
  W <- W$total
  
  if (!is.null(Y)) {
    z <- z.expectation(model,X,Y)

    # Calculate first component of PCA for W*z
    pca <- princomp(t(W%*%z))
    projvec <- pca$loadings[,1]

    # Project data to this component
    data <- rbind(X,Y)
    proj <- t(data)%*%projvec
   
    return(proj)
  }
  
  # for models with one data set
  else {
    W <- W$total
    z <- z.expectation(model,X)
    
    # Calculate first component of PCA for W*z
    pca <- princomp(t(W%*%z))
    projvec <- pca$loadings[,1]

    # Project data to this component
    data <- X
    proj <- t(data)%*%projvec
   
    return(proj)
  
  }
}

