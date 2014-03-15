\name{fun.check.gld.multi}
\alias{fun.check.gld.multi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Check whether the RS or FMKL/FKML GLD is a valid GLD for vectors of L1, L2, L3 
and L4
}
\description{
This function will return a logical vector showing whether vector combinations 
of L1, L2, L3 and L4 are valid GLDs.
}
\usage{
fun.check.gld.multi(l1, l2, l3, l4, param)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{l1}{
A numerical vector for L1 of GLD
}
  \item{l2}{
A numerical vector for L2 of GLD
}
  \item{l3}{
A numerical vector for L3 of GLD
}
  \item{l4}{
A numerical vector for L4 of GLD
}
  \item{param}{
Can be "rs", "fmkl", or "fkml"
}
}
\details{
This is an extension to \code{\link{fun.check.gld}}
}
\value{
A logical vector indicating whether the specified parameters are valid GLDs
}
\author{
Steve Su
}
\seealso{
\code{\link{fun.check.gld}}}

\examples{

fun.check.gld.multi(c(1,2,3),c(4,5,6),c(7,8,9),c(10,11,12),param="rs")

fun.check.gld.multi(c(1,2,3),c(4,5,6),c(7,8,9),c(10,11,-12),param="rs")

}

\keyword{distribution}

