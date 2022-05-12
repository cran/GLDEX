\name{fun.RPRS.mm}
\alias{fun.RPRS.mm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fit RS generalised lambda distribution to data set using moment 
matching }
\description{
  This function fits RS generalised lambda distribution to data set 
  using moment matching
}
\usage{
fun.RPRS.mm(data, rs.init = c(-1.5, 1.5), leap = 3, FUN = "runif.sobol", 
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
This function provides method of moment fitting scheme for RS GLD. Note this 
function can fail if there are no defined percentiles from the data set or if 
the initial values do not lead to a valid RS generalised lambda distribution.

This function is based on scheme detailed in the literature below but it has
been modified by the author (Steve Su).
}
\value{
A vector representing four parameters of the RS generalised lambda distribution.
}
\references{ Karian, Z. and E. Dudewicz (2000). Fitting Statistical 
Distributions: The Generalized Lambda Distribution and Generalised Bootstrap 
Methods. New York, Chapman and Hall. }
\author{ Steve Su }
\seealso{ \code{\link{fun.RPRS.ml}}, \code{\link{fun.RPRS.lm}}, 
\code{\link{fun.RPRS.qs}}, \code{\link{fun.data.fit.ml}} 
\code{\link{fun.data.fit.lm}}, \code{\link{fun.data.fit.qs}},  
\code{\link{fun.data.fit.mm}} }
\examples{
\donttest{
# Fitting the normal distribution
 fun.RPRS.mm(data=rnorm(1000,2,3),rs.init=c(-1.5,1.5),leap=3)
}}
\keyword{smooth} 