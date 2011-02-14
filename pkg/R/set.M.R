# (C) 2008-2011 Leo Lahti and Olli-Pekka Huovilainen          
# All rights reserved. 
# FreeBSD License (keep this notice)     



# "Do not worry about your difficulties in Mathematics. I can assure you
# mine are still greater."
# - Albert Einstein


set.M.isotropic <- function (wtw, sigma.sq, dx) {solve(wtw + diag(sigma.sq, dx, dx)) }


set.M.full <- function (W, phi.inv, dz) {
  # for full marginal covariance
  solve(t(W)%*%phi.inv%*%W + diag(dz))
}


set.M <- function (W, phi) {solve(t(W)%*%W/phi + diag(ncol(W)))}

set.M.full2 <- function (W, phi.inv, dz) {

  # This corresponds to 
  # M <- set.M.full(W$total, phi.inv$total, dz = Dim$Z) # for non-matched case
  # but should be faster

  # for full marginal covariance
  solve(t(W$X)%*%phi.inv$X%*%W$X + t(W$Y)%*%phi.inv$Y%*%W$Y + diag(dz))

}
 
