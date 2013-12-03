\name{fun.zero.omit}
\alias{fun.zero.omit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Returns a vector after removing all the zeros. }
\description{
  This function returns a vector after removing all the zeros.
}
\usage{
fun.zero.omit(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ A vector of data. }
}

\value{
  Returns a vector after removing zeros and also give information on the number 
  of zeros in the data removed.
}
\author{ Steve Su }
\note{ Missing value and Inf values are not removed in this zero removing 
process.
}
\seealso{ \code{\link{fun.which.zero}} }
\examples{
# Removing zero entries from the vector c(0,1,2,3,4,0,2)
fun.zero.omit(c(0,1,2,3,4,0,2))
}
\keyword{manip} 