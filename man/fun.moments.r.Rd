\name{fun.moments.r}
\alias{fun.moments.r}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculate mean, variance, skewness and kurtosis of a numerical vector }
\description{
  This function evaluates the mean, variance, skewness and kurtosis of a 
  numerical vector. Missing values are automatically removed.
}
\usage{
fun.moments.r(x, normalise = "N")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A numeric vector }
  \item{normalise}{ "Y" if you want kurtosis to be calculated with reference
  to kurtosis = 0 under Normal distribution. Default is "N".
}
}

\value{
A vector of mean, variance, skewness and kurtosis.
}
\author{ Steve Su }
\seealso{ \code{\link{fun.theo.mv.gld}} }
\note{ Please contact the author directly if you find a bug! }
\examples{

fun.moments.r(rnorm(1000))
fun.moments.r(rnorm(1000),normalise="Y")

}
\keyword{ univar }

