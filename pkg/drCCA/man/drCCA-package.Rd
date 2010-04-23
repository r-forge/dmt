\name{drCCA}
\alias{drCCA}
\docType{package}
\title{
    drCCA - dimensionality reduction with CCA
}
\description{
    
   A package for data fusion task in bioinformatics applications.
   The package can be used to combine two or more data sets with
   co-occurring samples (such as genes) into a common representation
   of low dimensionality. The package implements the drCCA method
   that uses linear transformation to create a representation
   focusing on effects shared by several data sources, while
   reducing the dimensionality by ignoring variation that is
   specific to any of the data sources alone. The underlying assumption
   is that such data-specific variation is often noise or less
   interesting due to being related to a phenomenom not visible in
   the other sources. Check the software's webpage for the latest version.


}
\details{

\tabular{ll}{
Package: \tab drCCA\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2007-12-20\cr
License: \tab LGPL\cr
URL: \tab http://www.cis.hut.fi/projects/mi/software/drCCA \cr
}

 This package contains functions implementing the drCCA method
 for data fusion.

 For a complete list of functions, use 'library(help="drCCA")'

}
\author{
 Abhishek Tripathi, Arto Klami

Maintainer:Abhishek Tripathi<abhishek.tripathi@cs.helsinki.fi>

}

\keyword{ package }

\examples{
     data(expdata1)
     data(expdata2)
     drCCAcombine(list(expdata1,expdata2)) #use default values

}
