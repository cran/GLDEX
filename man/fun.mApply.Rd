\name{fun.mApply}
\alias{fun.mApply}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Applying functions based on an index for a matrix. }
\description{
  This is a generic function that can be used to find mean, variance, sum or 
  other operations according to some
index imposed on the matrix or vector.
}
\usage{
fun.mApply(X, INDEX, FUN = NULL, ..., simplify = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{ Matrix with n rows. }
  \item{INDEX}{ Vector or list of vectors of length n. }
  \item{FUN}{ Function to operate on submatrices of \code{X} by \code{INDEX}}
  \item{\dots}{ Arguments to function. }
  \item{simplify}{ Set as \code{TRUE} by default, see \code{\link{sapply}} fo 
  details. }
}
\value{
 If \code{FUN} returns more than one number, \code{fun.mApply} returns a matrix 
 with rows corresponding to unique values of \code{INDEX}.
}

\author{ Tony Plate }
\examples{ 
# Finding the row medians of a matrix (matrix(1:20,nrow=5))
fun.mApply(matrix(1:20,nrow=5),list(1:5),median)
}
\keyword{manip}