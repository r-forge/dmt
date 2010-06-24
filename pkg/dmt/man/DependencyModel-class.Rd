\name{DependencyModel-class}
%\Rdversion{1.1}
\docType{class}
\alias{DependencyModel-class}
\alias{getW}
\alias{getPhi}
\alias{getScore}
\alias{getW,DependencyModel-method}
\alias{getPhi,DependencyModel-method}
\alias{getScore,DependencyModel-method}
\alias{getParams,DependencyModel-method}
\alias{getWindowSize,DependencyModel-method}
\alias{getModelMethod,DependencyModel-method}
\alias{getFeatureName,DependencyModel-method}
\title{Class "DependencyModel"}
\description{A Dependency model for one or two data sets}
\section{Objects from the Class}{
  Returned by \code{\link{fit.dependency.model}}, \code{\link{ppca}}, \code{\link{pfa}}, \code{\link{pcca}} and \code{\link{pcca.isotropic}} functions.
}
\section{Slots}{
  \describe{
    \item{W}{a list of X, Y and total components containing the relationship between two data sets; for dependency model
             for one dataset, only total is given}	    
    \item{phi}{a list of X, Y and total components containing the data set specific covariances; for dependency model
               for one dataset, only total is given}	 
    \item{score}{score for fitness of model}
    \item{loc}{middle location of the window in base pairs}   
    \item{windowSize}{size of the window}
    \item{method}{name of the used method}
    \item{params}{list of parameters used in dependency model}
  }
}

\section{Methods}{
  \describe{
    \item{getW}{\code{signature(model = "DependencyModel")}: Returns a list of model variable \code{W}s \code{X} 
      , \code{Y} and \code{total} component}
    \item{getPhi}{\code{signature(model = "DependencyModel")}:  Returns a list of model variable \code{phi}s \code{X}
      and \code{Y} and \code{total} component}
    \item{getScore}{\code{signature(model = "DependencyModel")}: Returns the dependency score of model}
    \item{getParams}{\code{signature(model = "DependencyModel")}: Returns a list of used parameters for the method}
    \item{getModelMethod}{\code{signature(model = "DependencyModel")}: Returns the name of the used method}
    \item{getWindowSize}{\code{signature(model = "DependencyModel")}: Returns the size of window}
  }
}

\author{Leo Lahti \email{leo.lahti@iki.fi} and Olli-Pekka Huovilainen \email{ohuovila@gmail.com}}

\seealso{Calculation of latent variable z: \code{\link{z.expectation}}.
For calculation of dependency models for chromosomal arm, chromosome or
genome: \code{\link{dependency.screen}}. For plotting dependency
scores see \link{dependency score plotting}.}
\examples{
## NOT RUN
#data(chromosome17)
#window <- fixed.window(geneExp, geneCopyNum, 10, 10)
#model <- fit.dependency.model(window$X, window$Y)
# Contributions of samples and variables to model
#plot(model,geneExp,geneCopyNum) FIXME
}
\keyword{classes}
