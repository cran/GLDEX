\name{fun.RPRS.lm}
\alias{fun.RPRS.lm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fit RS generalised lambda distribution to data set using L moment 
matching }
\description{
  This function fits RS generalised lambda distribution to data set 
  using L moment matching
}
\usage{
fun.RPRS.lm(data, rs.init = c(-1.5, 1.5), leap = 3, FUN = "runif.sobol", 
no = 10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ Dataset to be fitted }
  \item{rs.init}{ Initial values for RS distribution optimization, 
  \code{c(-1.5,1.5)} tends to work well. }
  \item{leap}{ See scrambling argument in \code{\link{fun.gen.qrn}}. }
  \item{FUN}{ A character string of either  \code{"runif.sobol"} (default), \code{"runif.sobol.owen"}, 
  \code{"runif.halton"} or \code{"QUnif"}. }
  \item{no}{ Number of initial random values to find the best initial values for 
  optimisation. }
}

\details{
This function provides method of L moment fitting scheme for RS GLD. Note this 
function can fail if there are no defined percentiles from the data set or if 
the initial values do not lead to a valid RS generalised lambda distribution.

This function is based on scheme detailed in the literature below but it has
been modified by the author (Steve Su).
}
\value{
A vector representing four parameters of the RS generalised lambda distribution.
}
\references{
Asquith, W. (2007). 
"L-moments and TL-moments of the generalized lambda distribution." 
Computational Statistics and Data Analysis 51(9): 4484-4496.

Karvanen, J. and A. Nuutinen (2008). 
"Characterizing the generalized lambda distribution by L-moments." 
Computational Statistics and Data Analysis 52(4): 1971-1983. }
\author{ Steve Su }
\seealso{ \code{\link{fun.RPRS.ml}}, \code{\link{fun.RPRS.mm}}, 
\code{\link{fun.RPRS.qs}}, \code{\link{fun.data.fit.ml}} 
\code{\link{fun.data.fit.lm}}, \code{\link{fun.data.fit.qs}},  
\code{\link{fun.data.fit.mm}} }
\examples{
\donttest{
# Fitting the normal distribution
 fun.RPRS.lm(data=rnorm(1000,2,3),rs.init=c(-1.5,1.5),leap=3)
}}
\keyword{smooth} 