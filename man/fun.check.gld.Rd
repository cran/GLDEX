\name{fun.check.gld}
\alias{fun.check.gld}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Check whether the RS or FMKL/FKML GLD is a valid GLD for single values of L1, 
L2, L3 and L4
}
\description{
This function will return a single logical value showing whether a combination 
of
L1, L2, L3 and L4 is a valid GLD.
}
\usage{
fun.check.gld(lambda1, lambda2, lambda3, lambda4, param)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lambda1}{
A numerical vector for L1 of GLD
}
  \item{lambda2}{
A numerical vector for L2 of GLD
}
  \item{lambda3}{
A numerical vector for L3 of GLD
}
  \item{lambda4}{
A numerical vector for L4 of GLD
}
  \item{param}{
Can be "rs", "fmkl", or "fkml"
}
}
\value{
A single logical value indicating whether the specified GLD is a valid 
probability density function
}
\author{
Steve Su
}
\seealso{
\code{\link{fun.check.gld.multi}}
}
\examples{

fun.check.gld(1,4,3,2,"rs")

fun.check.gld(1,4,3,2,"fkml")

fun.check.gld(1,4,3,-2,"rs")
}

\keyword{distribution}

