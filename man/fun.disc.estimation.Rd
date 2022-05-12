\name{fun.disc.estimation}
\alias{fun.disc.estimation}

\title{ Estimates the mean and variance after cutting up a vector of variable 
into evenly spaced categories. }
\description{
This function supplements \code{\link{fun.nclass.e}} and it is not intended to 
be used by the users directly.
}
\usage{
fun.disc.estimation(x, nint)
}
\arguments{
  \item{x}{ A vector of observations. }
  \item{nint}{ Number of intervals to cut the vectors into.}
}
\details{
The function cuts the vector into evenly spaced categories and estimate 
the mean and variance of
the actual data based on the categorisation.
}
\value{
Two numerical values, the first being the mean and the second being the 
variance.
}
\author{ Steve Su }
\seealso{ \code{\link{fun.nclass.e}} }
\examples{
\donttest{
## Cut up a randomly normally distributed observations into 5 evenly spaced 
## categories and estimate the mean and variance based on this cateogorisation.
junk<-rnorm(1000,3,2)
fun.disc.estimation(junk,5)
}}
\keyword{univar}