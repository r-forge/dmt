\name{z.expectation}
%\Rdversion{1.1}
\alias{z.expectation}
\title{The model parameter z}
\description{
  Expectation of the latent variable z.
}
\usage{
z.expectation(model, X, Y = NULL)
}

\arguments{
    \item{model}{
		The fitted dependency model.
	}
	\item{X, Y}{

	  Data sets used in fitting the dependency modeling functions
	  %(\code{\link{screen.cgh.mrna}} or
	  %\code{link{fit.dependency.model}}). 
	  Note: Arguments must be given
	  in the same order as in \code{\link{fit.dependency.model}} %or \code{\link{screen.cgh.mrna}}.
          Only X is needed for dependency model for one data set.
	}
	
}
\details{

  \code{z.expectation} gives ML estimate of the shared latent variable
  Z, given data X, Y and the model parameters in \code{model}.


}

\value{
\code{z.expectation} gives the matrix z.
}
\references{
Dependency Detection with Similarity Constraints,
Lahti et al., 2009 Proc. MLSP'09 IEEE International Workshop on Machine Learning for Signal Processing,
See \url{http://www.cis.hut.fi/lmlahti/publications/mlsp09_preprint.pdf}

A Probabilistic Interpretation of Canonical Correlation Analysis,
Bach Francis R. and Jordan Michael I. 2005 Technical Report 688. Department of
Statistics, University of California, Berkley.
\url{http://www.di.ens.fr/~fbach/probacca.pdf}

Probabilistic Principal Component Analysis,
Tipping Michael E. and Bishop Christopher M. 1999. \emph{Journal of the Royal Statistical Society},
Series B, \bold{61}, Part 3, pp. 611--622.
\url{http://research.microsoft.com/en-us/um/people/cmbishop/downloads/Bishop-PPCA-JRSS.pdf}
}
\author{
Olli-Pekka Huovilainen \email{ohuovila@gmail.com} and Leo Lahti \email{leo.lahti@iki.fi}
}


\seealso{
\code{\link{DependencyModel-class}}, 
%\code{\link{screen.cgh.mrna}}
}
\examples{
# Load example data X, Y
data(modelData) 

## Fit default model for the two data sets X, Y
model <- fit.dependency.model(X, Y)

## Retrieve the shared latent variable Z between data sets
z.expectation(model, X, Y)


}
\keyword{math}
