\name{which.na}
\alias{which.na}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Determine Missing Values }
\description{
Returns a vector showing the position of 
missing values in a vector.
}
\usage{
which.na(x)
}
\arguments{
  \item{x}{ An object which should be of \code{logical}, \code{numeric}, 
  or \code{complex} }
}
\value{
This returns the indices of values in \code{x} which are missing or 
"Not a Number". 
}
\examples{
# A non-zero number divided by zero creates infinity,
# zero over zero creates a NaN 
weird.values <- c(1/0, -2.9/0, 0/0, NA) 
which.na(weird.values) 
# Produces:
#  [1] 3 4 
}
\keyword{manip} 