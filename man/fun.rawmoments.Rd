\name{fun.rawmoments}
\alias{fun.rawmoments}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Computes the raw moments of the generalised lambda distribution up to 
4th order. }
\description{ This function is of theoretical interest only.
}
\usage{
fun.rawmoments(L1, L2, L3, L4, param = "fmkl")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{L1}{ Location parameter of the generalised lambda distribution. }
  \item{L2}{ Scale parameter of the generalised lambda distribution. }
  \item{L3}{ First shape parameter of the generalised lambda distribution. }
  \item{L4}{ Second shape parameter of the generalised lambda distribution. }
  \item{param}{ \code{"rs"} or \code{"fmkl"} specifying the type of 
  the generalised lambda distribution. }
}
\details{
  This function is the building block for \code{\link{fun.theo.bi.mv.gld}}.
}
\value{
A vector showing the raw moments of the specified generalised lambda 
distribution up to the fourth order.
}
\references{  Freimer, M., Mudholkar, G. S., Kollia, G. & Lin, C. T. (1988), A
     study of the generalized tukey lambda family, Communications in
     Statistics - Theory and Methods *17*, 3547-3567.

     Karian, Zaven A. and Dudewicz, Edward J. (2000), Fitting
     statistical distributions: the Generalized Lambda Distribution and
     Generalized Bootstrap methods, Chapman & Hall

     Ramberg, J. S. & Schmeiser, B. W. (1974), An approximate method
     for generating asymmetric random variables, Communications of the
     ACM *17*, 78-82. }
\author{ Steve Su }
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{

## Generate some random numbers using FMKL and RS generalised lambda 
## distributions and then compute the empirical and theoretical 
## E(X), E(X^2), E(X^3), E(X^4) 

junk<-rgl(100000,1,2,3,4)
mean(junk)
mean(junk^2)
mean(junk^3)
mean(junk^4)

junk<-rgl(100000,1,2,3,4,"rs")
mean(junk)
mean(junk^2)
mean(junk^3)
mean(junk^4)

fun.rawmoments(1,2,3,4)
fun.rawmoments(1,2,3,4,"rs")
}
\keyword{univar}