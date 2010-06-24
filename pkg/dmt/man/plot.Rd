\name{plot}
%\Rdversion{1.1}
\alias{plot.DependencyModel}
\alias{plot.ChromosomeArmModels}
\alias{plot.ChromosomeModels}
\alias{plot.GenomeModels}
\alias{dependency score plotting}
\title{
Dependency score plotting
}
\description{
Plot the contribution of the samples and variables to the dependency
model or dependency model fitting scores of chromosomal arm, chromosome
or genome.
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
	\item{hilightGenes}{
		vector of strings; Name of genes to be hilighted with dots.
	}
	\item{showDensity}{
		logical; if 'TRUE' small vertical lines are drwan in the bottom of 
		the plot under each gene.
	}
	\item{showTop}{
	  numeric; Number of models with highest dependencies to be
	hilighted. A horizontal dashed line is drawn to show threshold
	value. With \code{0} no line is drawn.
      }
      \item{topName}{
	logical; If \code{TRUE}, gene names are printed to hilighted
	models with highest dependecies. Otherwise hilighted models are
	numbered according to their rank in dependency score.
      }

      \item{type, xlab, ylab, main}{
	plot type and labels. See \code{\link{plot}} for details. In
	\code{plot.GenomeModels} these affet only if \code{onePlot} is
	\code{TRUE}.
      }
      
      \item{onePlot}{
	If \code{TRUE}, all dependency scores are plotted in one plot
	window. Otherwise one plot window is used for each chromosome.
      }

      \item{pch, cex}{
	symbol type and size for hilightGenes. See \code{\link{points}}
	for details.
      }
      
      \item{tpch, tcex}{
	symbol type and size for genes with highest scores. See
	\code{\link{points}} for details.
      }
      
      \item{ylim, xlim}{
	axis limits. Default values are calculated from
	data. Lower limit for y is 0 and upper limit is either 1
	or maximum score value. X limits are gene location
	range. See \code{\link{plot}} for details.
      }
      \item{mfrow, mar, ps, mgp}{
	chromosome plots' layout, marginals, text size and margin
	line. See \code{\link{par}} for details.
      }
      \item{...}{
	optional plotting parameters
      }
}
\details{
Function plots scores of each dependency model. 
}
\references{
Dependency Detection with Similarity Constraints Lahti et al.,
MLSP'09. See
\url{http://www.cis.hut.fi/lmlahti/publications/mlsp09_preprint.pdf}
}
\author{
Leo Lahti \email{leo.lahti@iki.fi} and Olli-Pekka Huovilainen \email{ohuovila@gmail.com}
}


\seealso{
\code{\link{DependencyModel-class}}, 
\code{\link{DependencyScreenModels-class}}, 
\code{\link{dependency.screen}}
}
\examples{
## NOT RUN
#data(chromosome17)
#window <- fixed.window(geneExp, geneCopyNum, 150, 10)
#models <- dependency.screen(geneExp$data, geneCopyNum$data)
#plot(models) FIXME
}
\keyword{hplot}
