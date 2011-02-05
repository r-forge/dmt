simCCA.optimize.fullcov.EM <-
function (X, Y, zDimension = NULL, epsilon = 1e-6) {

  # use this for full W (EM algorithm, unstable for n ~ p)
  res <- optimize.simCCA.W(W.init$X, phi.init, Dim, Dcov,
                           nullmat, epsilon, par.change = 1e6,
                           mySeed = mySeed + 1, dz = zDimension)

  res

}

