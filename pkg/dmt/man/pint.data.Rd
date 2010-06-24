\name{pint.data}
%\Rdversion{1.1}
\alias{pint.data}
\alias{pint.match}
\title{
Forms a data set and pairs samples in two data sets.
}
\description{
Forms a data set for use in functions in 'pint' package (e.g. \code{\link{screen.cgh.mrna}}). 
Pairs samples in two data sets.
}
\usage{
pint.data(data, info)
pint.match(X, Y, max.dist)
}
\arguments{
  \item{data}{
    Probe-level data in a matrix or data frame. 
  }
  \item{info}{
    Additional information of the features in a data frame.
  }
  \item{X, Y}{
    Data sets to be paired.
  }
  \item{max.dist}{maximum distance between paired genes in base pairs.
  }
}
\details{
Function \code{pint.match} goes through every sample in \code{X} and
finds the nearest sample in \code{Y} which is in the same chromosome
arm. If more than one sample in \code{X} has same nearest sample in
\code{Y}, all but one is discarded. Samples with longer distance
than \code{max.dist} are discarded.
}
\value{
 \code{pint.data} returns a list with a matrix with sample data and a
data frame with \code{chr} (chromosome), \code{arm} (chromosome arm) and
\code{loc} (location).

\code{pint.match} return a list with two data sets. These can be used in
\code{\link{dependency.screen}} function. Experimental.

}
\author{
Leo Lahti \email{leo.lahti@iki.fi} and Olli-Pekka Huovilainen \email{ohuovila@gmail.com}
}

\seealso{
\code{\link{dependency.screen}}
}

\examples{
#data(chromosome17)
#newData <- pint.match(geneExp, geneCopyNum)
}

