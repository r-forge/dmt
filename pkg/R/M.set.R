

set.M.full <- function (W, phi.inv, dz) {
  # for full marginal covariance
  solve(t(W)%*%phi.inv%*%W + diag(dz))
}


set.M.full2 <- function (W, phi.inv, dz) {

  # This corresponds to 
  # M <- set.M.full(W$total, phi.inv$total, dz = Dim$Z) # for non-matched case
  # but should be faster

  # for full marginal covariance
  solve(t(W$X)%*%phi.inv$X%*%W$X + t(W$Y)%*%phi.inv$Y%*%W$Y + diag(dz))

}


 
set.M.isotropic <- function (wtw, sigma.sq, dx) {
  solve(wtw + diag(sigma.sq, dx, dx))
}


set.M <-
function (W, phi) {solve(t(W)%*%W/phi + diag(ncol(W)))}

