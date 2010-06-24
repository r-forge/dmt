

 # the inverse of square root of a square matrix

"SqrtInvMat" <-
 function(matrix)
{
 x <- as.matrix(matrix)

 fac <- svd(x)


 res <-(fac$v %*% (diag(1/sqrt(fac$d),nrow = length(fac$d))) %*% t(fac$u))

 return(res)

}

