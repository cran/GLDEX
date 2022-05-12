\name{fun.gen.qrn}
\alias{fun.gen.qrn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Finds the low discrepancy quasi random numbers }
\description{
  This function calls the \code{runif.sobol}, \code{runif.sobol.owen} and 
  \code{runif.halton} essentially from the \pkg{spacefillr} package.
}
\usage{
fun.gen.qrn(n, dimension, scrambling, FUN = "runif.sobol")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{ Number to generate. }
  \item{dimension}{ Number of dimensions. }
  \item{scrambling}{ seed used, or leap as in the case of 
  \code{\link{QUnif}}. }
  \item{FUN}{ This can be \code{"runif.sobol"} (default), \code{"runif.sobol.owen"}, \code{"runif.halton"} 
  or \code{"QUnif"}. }
}

\value{
 A vector of values if dimension=1, otherwise a matrix of values between 0 
 and 1.
}

\author{ Steve Su }

\examples{
fun.gen.qrn(1000,5,3,"runif.sobol")

fun.gen.qrn(1000,5,409,"QUnif")
}
\keyword{datagen}
