\name{calculate}
%\Rdversion{1.1}
\alias{dependency.screen}
\title{Screening dependency in matched matrices.}
\description{The function assumes two matrices with matched and
ordered features. Dependency model is fitted within sliding windows of
fixed size. The dependency scores between the windows can be compared
to identify windows with the highest dependency. This function can be
used for instance to detect local chromosomal regions where
neighboring genes show high dependency between two measurement
sources, such as gene copy number and expression data.}
\usage{dependency.screen(X, Y, windowSize, method = "", params = list(), max.dist = 1e7)}
\arguments{
  \item{X,Y}{Data sets in a matrix form (samples x features)}
  \item{windowSize}{
    Determine the window size. This specifies the number of nearest
    features to be included in the sliding window, and
    therefore the scale of the investigated region.
  }
  \item{method}{
    Specify the dependency model: 
    \describe{
      \item{'pCCA'}{probabilistic canonical correlation analysis \cite{Bach & Jordan 2005}} 
      \item{'pPCA'}{probabilistic principal component analysis \cite{Tipping & Bishop 1999}}
      \item{'pFA'}{probabilistic factor analysis \cite{Rubin & Thayer 1982}}
      \item{'pSimCCA'}{probabilistic similarity constrained canonical correlation analysis \cite{Lahti et al. 2009}}
      \item{'TPriorpSimCCA'}{probabilistic similarity constarined canonical correlation analysis
	with possibility to tune T prior (Lahti et al. 2009)}
    }
    If anything else, the model is specified by the given parameters.
  }
  \item{params}{
    List of parameters for the dependency model.
    \describe{
      \item{sigmas}{Variance parameter for the matrix normal prior
	distribution of the transformation matrix T. This describes the
	deviation of T from H}
      \item{H}{Mean parameter for the matrix normal prior distribution
	prior of transformation matrix T}
      \item{zDimension}{Dimensionality of the latent variable}
      \item{mySeed}{Random seed.}
      \item{covLimit}{Convergence limit. Default depends on the selected
	method: 1e-3 for pSimCCA with full marginal covariances and 1e-6
	for pSimCCA in other cases.}
    }
  }

  \item{max.dist}{Maximum allowed distance between probes. Used in
    automated matching of the probes between the two data sets based on
    chromosomal location information.}
  
} 

\details{Function \code{dependency.screen} assumes that data is
already matched. This can be done with \code{\link{pint.match}}. It
takes sliding gene windows with \code{\link{fixed.window}} and fits
dependency models to each window with
\code{\link{fit.dependency.model}} function.  If method name is
specified, this overrides the corresponding model parameters,
corresponding to the modeling assumptions of the specified
model. Otherwise method for dependency models is determined by
parameters. Dependency scores are plotted with \link{dependency score
plotting}.}

\value{ Depending on the arguments, returns a
\linkS4class{DependencyScreenModels}. This contains the models for
each sliding window.}

\references{

Dependency Detection with Similarity Constraints, Lahti et al., 2009
Proc. MLSP'09 IEEE International Workshop on Machine Learning for Signal
Processing, See
\url{http://www.cis.hut.fi/lmlahti/publications/mlsp09_preprint.pdf}

A Probabilistic Interpretation of Canonical Correlation Analysis, Bach
Francis R. and Jordan Michael I. 2005 Technical Report 688. Department
of Statistics, University of California, Berkley.
\url{http://www.di.ens.fr/~fbach/probacca.pdf}

Probabilistic Principal Component Analysis, Tipping Michael E. and
Bishop Christopher M. 1999. \emph{Journal of the Royal Statistical
Society}, Series B, \bold{61}, Part 3, pp. 611--622.
\url{http://research.microsoft.com/en-us/um/people/cmbishop/downloads/Bishop-PPCA-JRSS.pdf}

EM Algorithms for ML Factoral Analysis, Rubin D. and Thayer
D. 1982. \emph{Psychometrika}, \bold{vol. 47}, no. 1.

}

\author{
Leo Lahti \email{leo.lahti@iki.fi} and Olli-Pekka Huovilainen \email{ohuovila@gmail.com}
}

\seealso{
To fit a dependency model: \code{\link{fit.dependency.model}}.
\linkS4class{ScreenDependencyModels} holds dependency models for chromosomal arm,
for plotting see \link{dependency score plotting}
}
\examples{
data(chromosome17)
X <- geneExp
Y <- geneCopyNum

## pSimCCA model
models <- dependency.screen(X, Y)
#plot(models)

## pCCA model with 3-dimensional latent variable z, and other parameters explicitly specified
models <- dependency.screen(X, Y, windowSize = 10, method="pCCA", params = list(zDimension = 3))
#plot(models) FIXME

}
\keyword{math}
\keyword{iteration}
