\name{fun.RPRS.qs}
\alias{fun.RPRS.qs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fit RS generalised lambda distribution to data set using quantile 
matching }
\description{
  This function fits RS generalised lambda distribution to data set 
  using quantile matching
}
\usage{
fun.RPRS.qs(data, rs.init = c(-1.5, 1.5), leap = 3, FUN = "runif.sobol", 
trial.n = 100, len = 1000, type = 7, no = 10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ Dataset to be fitted }
  \item{rs.init}{ Initial values for RS distribution optimization, 
  \code{c(-1.5,1.5)} tends to work well. }
  \item{leap}{ See scrambling argument in \code{\link{fun.gen.qrn}}. }
  \item{FUN}{ A character string of either  \code{"runif.sobol"} (default), \code{"runif.sobol.owen"}, 
  \code{"runif.halton"} or \code{"QUnif"}. }
  \item{trial.n}{ Number of evenly spaced quantile ranging from 0 to 1 to be
  used in the checking phase, to find the best set of initial values for 
  optimisation, this is intended to be lower than \code{len} to speed up the 
  fitting algorithm. Default is 100. }
  \item{len}{ Number of evenly spaced quantile ranging from 0 to 1 to be used,
  default is 1000 }
  \item{type}{ Type of quantile to be used, default is 7, see \code{quantile} }
  \item{no}{ Number of initial random values to find the best initial values for 
  optimisation. }
}

\details{
This function provides quantile matching fitting scheme for RS GLD. Note this 
function can fail if there are no defined percentiles from the data set or if 
the initial values do not lead to a valid RS generalised lambda distribution.
}
\value{
A vector representing four parameters of the RS generalised lambda distribution.
}
\references{ 
Su (2008). Fitting GLD to data via quantile matching method. (Book chapter to
appear) }
\author{ Steve Su }
\seealso{ \code{\link{fun.RPRS.ml}}, \code{\link{fun.RPRS.lm}}, 
\code{\link{fun.RPRS.mm}}, \code{\link{fun.data.fit.ml}} 
\code{\link{fun.data.fit.lm}}, \code{\link{fun.data.fit.qs}},  
\code{\link{fun.data.fit.mm}} }
\examples{
\donttest{
# Fitting the normal distribution
 fun.RPRS.qs(data=rnorm(1000,2,3),rs.init=c(-1.5,1.5),leap=3)
}}
\keyword{smooth} 