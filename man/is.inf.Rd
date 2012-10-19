\name{is.inf}
\alias{is.inf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Returns a logical vecto, TRUE if the value is Inf or -Inf. }
\description{
  This function works in similar fashion as in \code{\link{is.na}} and 
  \code{\link{is.notinf}}.
}
\usage{
is.inf(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A numerical value or a vector of data. }
}

\value{
  A logical vector, \code{T} if the value is Inf or -Inf.
}

\author{ Steve Su }
\note{ In the presence of missing value, the function will return a missing 
value.
}
\seealso{ \code{\link{is.na}},\code{\link{is.notinf}} }
\examples{
is.inf(c(Inf,2,2,1,-Inf))
}
\keyword{manip}