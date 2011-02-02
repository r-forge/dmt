
set.M.full <- function (W, phi.inv, dz) {
  # for full marginal covariance
  solve(t(W)%*%phi.inv%*%W + diag(dz))
}

#  solve(t(W$X)%*%phi.inv$X%*%W$X + t(W$Y)%*%phi.inv$Y%*%W$Y + diag(dz))



