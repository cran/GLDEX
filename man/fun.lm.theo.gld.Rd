\name{fun.lm.theo.gld}
\alias{fun.lm.theo.gld}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Find the theoretical first four L moments of the generalised lambda 
distribution. }
\description{
  This function computes the first four L moments for both RS and FMKL  
  generalised lambda distributions.
}
\usage{
fun.lm.theo.gld(L1, L2, L3, L4, param)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{L1}{ Lambda 1. Or c(Lambda 1,Lambda 2,Lambda 3,Lambda 4). }
  \item{L2}{ Lambda 2. }
  \item{L3}{ Lambda 3. }
  \item{L4}{ Lambda 4. }
  \item{param}{ \code{rs} or \code{fmkl} distribution. }
}
\value{
A vector listing the first four L moments
}
\references{ Asquith, W. (2007). 
"L-moments and TL-moments of the generalized lambda distribution." 
Computational Statistics and Data Analysis 51(9): 4484-4496.

Karvanen, J. and A. Nuutinen (2008). 
"Characterizing the generalized lambda distribution by L-moments." 
Computational Statistics and Data Analysis 52(4): 1971-1983. }
\author{ Steve Su }
\note{ Sometimes the theoretical moments may not exist, in those cases, 
\code{NA} is returned.
}
\seealso{ \code{\link{fun.theo.mv.gld}}  }
\examples{
fun.lm.theo.gld(1, 2, 3, 4, "rs")
fun.lm.theo.gld(1, 2, 3, 4, "fmkl")
}
\keyword{univar} 