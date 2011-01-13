\name{plot}
%\Rdversion{1.1}
\alias{plot.DependencyModel}
\alias{dependency score plotting}
\title{
Dependency score plotting
}
\description{
Plot the contribution of the samples and variables to the dependency
model 
%or dependency model fitting scores of chromosomal arm, chromosome
%or genome.
}
\usage{
plot.DependencyModel(x, X, Y, ann.types = NULL, ann.cols = NULL, legend.x = 0, 
        legend.y = 1, legend.xjust = 0, legend.yjust = 1, order = FALSE, 
        cex.z = 0.6, cex.WX = 0.6, cex.WY = 0.6, ...)

}

\arguments{
	\item{x}{
		\code{\link{DependencyModel-class}} models to be plotted. 
	}
	\item{X, Y}{
		data sets used in dependency modeling.
	}
	\item{ann.types}{
		a factor for annotation types for samples. Each value corresponds one sample in 
		datasets. Colors are used to indicate different types.
	}
	\item{ann.cols}{
		colors used to indicate different annotation types. Gray scale is used if 'NULL' given.
	}
	\item{legend.x, legend.y}{
		the x and y co-ordinates to be used to position the legend for annotation types.
	}
	\item{legend.xjust, legend.yjust}{
		how the legend is to be justified relative to the legend x and y location.
		A value of 0 means left or top justified, 0.5 means centered and 1 means 
  		right or bottom justified.
	}
	\item{order}{
		logical; if 'TRUE', values for sample contributions are ordered according
 		to their values.
	}
	\item{cex.z, cex.WX, cex.WY}{
		Text size for variable names.
	}

    \item{...}{
	optional plotting parameters
      }
}
%\details{
%Function plots scores of each dependency model of a gene for the whole
%chromosomal arm, chromosome or genome according to used
%method. \code{plot(x, cancerGenes = NULL, showDensity = FALSE, ...)} is
%also usable and chosen according to class of \code{models}.
%}
\references{
Dependency Detection with Similarity Constraints Lahti et al.,
MLSP'09. See
\url{http://www.cis.hut.fi/lmlahti/publications/mlsp09_preprint.pdf}
}
\author{
Olli-Pekka Huovilainen \email{ohuovila@gmail.com}
}


\seealso{
\code{\link{DependencyModel-class}}, 
}
\examples{
data(modelData) # load example data X, Y
model <- fit.dependency.model(X, Y) # fit default model
plot(model, X, Y)
}
\keyword{hplot}
