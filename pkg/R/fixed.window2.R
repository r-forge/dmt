fixed.window2 <-
function (X, Y, middleIndex, windowSize){
		
  # Indices for the window
  mini <- (middleIndex - (trunc((windowSize - 1)/2)))
  maxi <- (middleIndex + (trunc(windowSize/2)))	
  inds <- mini:maxi	

  if (length(inds) >= windowSize) {
  
    if (min(inds) < 1) {mini <- 1; maxi <- mini + windowSize - 1}
    if (max(inds) > nrow(X)) {maxi <- nrow(X); mini <- maxi - windowSize + 1}

  } else {stop(paste("windowSize should not exceed the number of features (", nrow(X), ")" ))}

  inds <- mini:maxi

  # Check that indices don't get out of bounds
  indsOutBounds <- (min(inds) < 1 || max(inds) > nrow(X))

  if(!indsOutBounds){
    Xm <- t(centerData(t(X[inds,]), rm.na = TRUE))
    Ym <- t(centerData(t(Y[inds,]), rm.na = TRUE))
    res <- list(X = Xm, Y = Ym, fail = FALSE)
  }
  else {
    res <- list(fail = TRUE)
  }
  res
}

