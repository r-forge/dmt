\name{fit.dependency.model}
%\Rdversion{1.1}
\alias{fit.dependency.model}
\alias{ppca}
\alias{pcca}
\alias{pcca.isotropic}
\alias{pfa}
\title{Fit dependency model between two data sets.}

\description{Fits the generative latent variable model (see vignette
for model specification) between two data sets with selected prior
constraints on marginal covariance structures, the structure of W,
latent dimensionality etc. Probabilistic versions of PCA, factor
analysis and CCA are available as special cases.}

\usage{
fit.dependency.model(X, Y, zDimension = 1, marginalCovariances = "full",
                     covLimit = 1e-3,
                     priors = list(), matched = TRUE,
                     includeData = TRUE, calculateZ = TRUE)

ppca(X, Y, zDimension = NULL, matched = FALSE, includeData = TRUE, calculateZ = TRUE)
pfa(X, Y = NULL, zDimension = NULL, matched = FALSE, includeData = TRUE, calculateZ = TRUE)
pcca.isotropic(X, Y, zDimension = NULL, matched = FALSE, covLimit = 1e-6, includeData = TRUE, calculateZ = TRUE)
pcca(X, Y, zDimension = NULL, matched = FALSE, includeData = TRUE, calculateZ = TRUE)

}
\arguments{
  \item{X, Y}{
The two data sets. 'Variables x samples'. The second data set
(\code{Y}) is optional.}
  \item{zDimension}{Dimensionality of the shared latent variable.}
  \item{marginalCovariances}{
    Structure of marginal covariances. Options: \code{"identical isotropic"},
    \code{"isotropic"}, \code{"diagonal"} and \code{"full"}
  }
  
%  \item{H}{ Mean of the matrix normal prior distribution for the
%    transformation matrix T. Must be a matrix of size (variables in
%    first data set) x (variables in second data set). If value is
%    \code{1}, H will be made identity matrix of appropriate size.}

%  \item{sigmas}{ Variance parameter for the matrix normal prior
%    distribution of the transformation matrix T. Described the allowed
%    deviation scale of the transformation matrix T from the mean matrix
%    H.}

  \item{covLimit}{ Convergence limit. The default value depends on the
selected model. }

  \item{priors}{Prior parameters for the model.}
    \describe{
  \item{Nm.wxwy.mean}{ Mean of the matrix normal prior distribution for the
    transformation matrix T. Must be a matrix of size (variables in
    first data set) x (variables in second data set). If value is
    \code{1}, \code{Nm.wxwy.mean} will be made identity matrix of appropriate size.}

  \item{Nm.wxwy.sigma}{ Variance parameter for the matrix normal prior
    distribution of the transformation matrix \code{T}. Described the allowed
    deviation scale of the transformation matrix \code{T} from the mean matrix
    \code{Nm.wxwy.mean}.}
    
    }

  \item{matched}{Logical indicating if the variables (dimensions) are
matched between X and Y. Applicable only when dimX = dimY. Affects the
results only when prior on the relationship Wx ~ Wy is set.}
  
  \item{includeData}{Logical indicating whether the original data is
included to the model output. Using \code{FALSE} can be used to
    save memory.}
  
  \item{calculateZ}{Logical indicating whether an expectation of the latent variable Z
is included in the model output. Otherwise the expectation can be
calculated with \code{getZ} or \code{z.expectation}. Using \code{FALSE} speeds up the calculation of the dependency model.}
}

\details{ Various choices for the model structure and parameters can
  be set through the general function \code{fit.dependency.model}.
  Special cases of the model, obtained with particular prior
  assumptions, include probabilistic canonical correlation analysis
  (\code{pcca}; \cite{Bach & Jordan 2005}), probabilistic principal
  component analysis (\code{ppca}; \cite{Tipping & Bishop 1999}),
  probabilistic factor analysis (\code{pfa}; \cite{Rubin & Thayer
  1982}), and regularized canonical correlation analysis (pSimCCA;
  \cite{Lahti et al. 2009}).
 
  The standard probabilistic PCA and factor analysis are methods for
  modeling a single data set as X ~ N(WZ, Bx*t(Bx)), where the
  covariance Bx*t(Bx) is isotropic for pPCA and diagonal for pFA.
  However, these methods can here be additionally applied to model
  dependencies between two data sets (X, Y), assuming the model X ~
  N(Wx*Z, Bx*t(Bx)); Y ~ N(Wy*Z, By*t(By)) i.e. shared latent variable
  Z plus dataset specific variation.
  
 The special cases are obtained with the following choices in the
 \code{fit.dependency.model} function:

\describe{

  \item{pPCA}{ \code{marginalCovariances = "identical isotropic"}
(\cite{Tipping & Bishop 1999}) }

\item{pFA}{ \code{marginalCovariances = "diagonal"} (\cite{Rubin &
Thayer 1982}) }

\item{pCCA}{ \code{marginalCovariances = "full"} or
\code{"isotropic"} (\cite{Bach & Jordan 2005}) }

\item{pSimCCA}{ \code{marginaCovariances = "full", priors = list(Nm.wxwy.mean = I, Nm.wxwy.sigma = 0)}. This
is the default method.  (\cite{Lahti et al. 2009}) }

\item{pSimCCA with T prior}{ \code{H = I, marginalCovariances =
"isotropic", priors = list(Nm.wxwy.mean=1} (\cite{Lahti et al. 2009})
}

}

%Resulting \linkS4class{DependencyModel} object does not have location or z variable. 
%Location can be set with \code{setLoc} method (see examples)
%Expectation of the latent variable z can be calculated with
%\code{link{z.expectation}}.

To avoid computational singularities, the covariance matrix phi is
regularised by adding a small constant to the diagonal.
}

\value{
\linkS4class{DependencyModel}
}

\references{
Dependency Detection with Similarity Constraints, Lahti et al., 2009
Proc. MLSP'09 IEEE International Workshop on Machine Learning for Signal
Processing, \url{http://arxiv.org/abs/1101.5919}

A Probabilistic Interpretation of Canonical Correlation Analysis, Bach
Francis R. and Jordan Michael I. 2005 Technical Report 688. Department
of Statistics, University of California, Berkley.
\url{http://www.di.ens.fr/~fbach/probacca.pdf}

Probabilistic Principal Component Analysis, Tipping Michael E. and
Bishop Christopher M. 1999. \emph{Journal of the Royal Statistical
Society}, Series B, \bold{61}, Part 3, pp. 611--622.
\url{http://research.microsoft.com/en-us/um/people/cmbishop/downloads/Bishop-PPCA-JRSS.pdf}

EM Algorithms for ML Factorial Analysis, Rubin D. and Thayer
D. 1982. \emph{Psychometrika}, \bold{vol. 47}, no. 1.

}
\author{ Olli-Pekka Huovilainen \email{ohuovila@gmail.com} and Leo Lahti
\email{leo.lahti@iki.fi} }

\seealso{ 
%For windowing data: \code{\link{fixed.window}}. 
Reults from
this function: \linkS4class{DependencyModel}. 
%Calculating dependency
%models to chromosomal arm, chromosome or genome
%\code{\link{screen.cgh.mrna}}. For calculation of latent variable z:
%\code{link{z.expectation}}.  
}


\examples{
data(modelData) # Load example data X, Y

# pSimCCA model
model <- fit.dependency.model(X, Y, zDimension = 1)

# pPCA model for two data sets
model <- ppca(X, Y)

# Getting the latent variable Z when it has been calculated with the model
#getZ(model)

#model2 <- fit.dependency.model(X, Y, includeData = FALSE, calculateZ = FALSE)
# Calculating the latent variable Z when its not calculated and the original data is not included
#getZ(model, X, Y)

}
\keyword{math}
\keyword{iteration}
