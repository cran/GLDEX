\name{pretty.su}
\alias{pretty.su}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ An alternative to the normal pretty function in R. }
\description{
Divide a range of values into equally spaced divisions. End points are given as 
output.
}
\usage{
pretty.su(x, nint = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A vector of values. }
  \item{nint}{ Number of intervals required. }
}
\details{
This is also used for the plotting of histogram in the \code{\link{histsu}} 
function.
}
\value{
A vector of endpoints dividing the data into equally spaced regions.
}

\author{ Steve Su }

\seealso{ \code{\link{pretty}} }
\examples{
# Generate random numbers from normal distribution:
junk<-rnorm(1000,2,3)

# Cut them into 7 regions, 8 endpoints.
pretty.su(junk,7)

}
\keyword{manip} 