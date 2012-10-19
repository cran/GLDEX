\name{fun.gen.qrn}
\alias{fun.gen.qrn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Finds the low discrepancy quasi random numbers }
\description{
  This function calls the \code{link{runif.sobol}} and 
  \code{\link{runif.halton}} from the \pkg{fOptions} package and 
  \code{\link{QUnif}} from \pkg{sfsmisc} package.
}
\usage{
fun.gen.qrn(n, dimension, scrambling, FUN = "runif.sobol")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{ Number to generate. }
  \item{dimension}{ Number of dimensions. }
  \item{scrambling}{ Scrambling method used, or leap as in the case of 
  \code{\link{QUnif}}. }
  \item{FUN}{ This can be \code{"runif.sobol"} (default), \code{"runif.halton"} 
  or \code{"QUnif"}. }
}

\value{
 A vector of values if dimension=1, otherwise a matrix of values between 0 
 and 1.
}
\references{ Bratley P., Fox B.L. (1988); Algorithm 659: Implementing Sobol's
     Quasirandom Sequence Generator, ACM Transactions on Mathematical
     Software 14, 88-100.

     Joe S., Kuo F.Y. (1998); Remark on Algorithm 659: Implementing
     Sobol's Quaisrandom Seqence Generator.
 }
\author{ Steve Su }
\seealso{ \code{link{runif.sobol}}, \code{\link{runif.halton}}, 
\code{\link{QUnif}} }
\examples{
fun.gen.qrn(1000,5,3,"runif.sobol")

fun.gen.qrn(1000,5,409,"QUnif")
}
\keyword{datagen}
