
\name{SqrtInvMat}
\alias{SqrtInvMat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Squareroot of inverse of a square matrix}
\description{
    The function calculates the  squareroot of inverse of a square matrix, 
    provided the given matrix is invertible. The function is based on
    Singular Value Decomposition.
}
\usage{
SqrtInvMat(matrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{matrix}{A square and invertible matrix}
}
\details{
         The function uses \code{\link{svd}} function. If the matrix is not 
  invertible, it gives an error.

}
\value{

  The function returns the squareroot of inverse of given matrix.
}

\references{Gene H. Golub, Charles F. Van Loan,1989 \emph{The john Hopkins University Press}, Matrix computations
}
\author{Abhishek Tripathi, Arto Klami }
\seealso{\code{\link{svd}},\code{\link{eigen}} }

\examples{
 
   # create a randome matrix
#     x <- matrix(runif(50),10)
#     cov_x <- cov(x) # to assure invertibility
#     SqrtInvMat(cov_x)

}
\keyword{internal}

