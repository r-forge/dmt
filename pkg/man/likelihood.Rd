%\Rdversion{1.1}
\name{pfa.log.likelihood}
\alias{pfa.log.likelihood}
\title{Likelihood for the models.}
\description{Likelihood functions for the implemented models.}
\usage{pfa.log.likelihood(X, wtw, phi)}
\arguments{
  \item{X}{Data: features x samples matrix.}
  \item{wtw}{Latent covariance structure used in PFA, correspond to W*t(W) in the model X = Wz + epsilon.}
  \item{phi}{Marginal covariance in the model X = Wz + epilon with epsilon ~ N(0, phi).}
}
\details{Other likelihood functions will be added later.}
\value{Log-likelihood of the data, given the model parameters.}
\references{See citation("dmt").}
\author{Leo Lahti \email{leo.lahti@iki.fi}}
\seealso{fit.dependency.model, pfa}
\examples{
library(dmt)
# Generate toydata
N <- 100
xdim <- 10
zdim <- 3
toy <- generate.toydata(N = N, zDim = zdim, xDim = xdim, yDim = xdim, marginal.covariances = "diagonal")
# Estimate model parameters
res <- pfa(toy$X, zDimension = zdim)
W <- res@W$total
phi <- res@phi$total
wtw <- crossprod(t(W)) # is the same as W * t(W)
# Calculate log-likelihood for the model
L <- pfa.log.likelihood(toy$X, wtw, phi)
}
\keyword{ utilities }

