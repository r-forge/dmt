set.M.full2 <- function (W, phi.inv, dz) {

  # This corresponds to 
  # M <- set.M.full(W$total, phi.inv$total, dz = Dim$Z) # for non-matched case
  # but should be faster

  # for full marginal covariance
  solve(t(W$X)%*%phi.inv$X%*%W$X + t(W$Y)%*%phi.inv$Y%*%W$Y + diag(dz))

}


 