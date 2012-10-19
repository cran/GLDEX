\name{fun.nclass.e}
\alias{fun.nclass.e}
\title{ Estimates the number of classes or bins to smooth over in the 
discretised method of fitting generalised lambda distribution to data. }
\description{
  Support function for discretised method of fitting distribution to data.
}
\usage{
fun.nclass.e(x)
}
\arguments{
  \item{x}{ Vector of data. }
}
\details{
  This function calculates the mean and variance of the discretised data from 1 
  to the very last observation and chooses the best number
of categories that represent the mean and variance of the actual data set 
through the criterion of squared deviations.
}
\value{
  A numerical value suggesting the best number of class that can be used to 
  represent the mean and variane of the original data set.
}
\author{ Steve Su }
\note{ This is not designed to be called directly by end user.
}
\seealso{ \code{\link{fun.disc.estimation}} }
\examples{
fun.nclass.e(rnorm(100,3,2))
}
\keyword{univar}