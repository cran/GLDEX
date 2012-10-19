\name{fun.which.zero}
\alias{fun.which.zero}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Determine which values are zero. }
\description{
  Returns an integer vector showing the position of zero values in the data.
}
\usage{
fun.which.zero(data)
}
\arguments{
  \item{data}{ A vector of data. }
}
\value{
An integer vector showing the position of zero values in the data.
}

\author{ Steve Su }
\note{ Any missing values will be returned as missing. }
\seealso{ \code{\link{fun.zero.omit}} }
\examples{
# Finding where the zeros are in this vector: c(0,1,2,3,4,0,2)
fun.which.zero(c(0,1,2,3,4,0,2))
# Finding where the zeros are in this vector: c(0,1,2,3,NA,0,2)
fun.which.zero(c(0,1,2,3,NA,0,2))
}
\keyword{manip} 