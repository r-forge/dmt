\name{window}
%\Rdversion{1.1}
\alias{fixed.window}
\alias{sparse.window}
\title{
Form data with a selected window size for the model fitting
}
\description{

Forms a chosen window of two data matrices to use for
\code{fit.dependency.model} picking same number of features from both directions. \code{sparse.window}
forms a window around one sample in the first data set with a number of
samples from the second data set.

}
\usage{
fixed.window(X, Y, middleIndex, windowSize)
sparse.window(X, Y, xIndex, windowSize)
}
\arguments{
  \item{X}{
First data set. In \code{sparse.window} windows will be formed around each sample in this data set.
}
  \item{Y}{
Second data set.
}
  \item{middleIndex}{
Index of middle position for window.
}
  \item{xIndex}{
Index of middle position in \code{X} for window.
}
  \item{windowSize}{
Number of genes in window. In \code{sparse.window} \code{X} has always one sample in window.
}
}
\details{
Window contains windowSize nearest genes. Warning is given if windowSize genes is not found in the same chromosomal arm.
Data of both data sets is normalised so that each genes data has zero mean.
}

\value{
List of window data:
\item{X}{window of the first data set}
\item{Y}{window of the second data set}
\item{loc}{location of gene}
\item{geneName}{name of the gene}
\item{edge}{logical; TRUE if iteration to one direction has stopped because edge of data in chromosomal arm has been found.}
\item{fail}{logical; TRUE if chromosomal arm contains less than windowSize genes.}
}
\author{
Leo Lahti (leo.lahti@iki.fi and Olli-Pekka Huovilainen \email{ohuovila@gmail.com}
}

\seealso{
Dependency model fitting: \code{\link{fit.dependency.model}}
}
\examples{
#data(chr7)
#window <- fixed.window2(X, Y, 10, 10)
#model <- fit.dependency.model(window$X, window$Y, H = NA)
#setFeatureName(model) <- window$geneName # FIXME: add this functionality
}

