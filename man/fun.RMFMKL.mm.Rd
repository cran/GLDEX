\name{fun.RMFMKL.mm}
\alias{fun.RMFMKL.mm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fit FMKL generalised lambda distribution to data set using moment 
matching }
\description{
  This function fits FMKL generalised lambda distribution to data set 
  using moment matching
}
\usage{
fun.RMFMKL.mm(data, fmkl.init = c(-0.25, 1.5), leap = 3, FUN = "runif.sobol", 
no = 10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ Dataset to be fitted }
  \item{fmkl.init}{ Initial values for FMKL distribution optimization, 
  \code{c(-0.25,1.5)} tends to work well. }
  \item{leap}{ Scrambling (0,1,2,3) for the Sobol sequence for the distribution 
  fit. See scrambling/leap argument for \code{\link{runif.sobol}}, 
  \code{\link{runif.halton}} or \code{\link{QUnif}}. }
  \item{FUN}{ A character string of either  \code{"runif.sobol"} (default), 
  \code{"runif.halton"} or \code{"QUnif"}. }
  \item{no}{ Number of initial random values to find the best initial values for 
  optimisation. }
}
\details{
This function provides method of moment fitting scheme for FMKL GLD. Note this 
function can fail if there are no defined percentiles from the data set or if 
the initial values do not lead to a valid FMKL generalised lambda distribution.

This function is based on scheme detailed in the literature below but it has
been modified by the author (Steve Su).
}
\value{
A vector representing four parametefmkl of the FMKL generalised lambda distribution.
}
\references{ Karian, Z. and E. Dudewicz (2000). Fitting Statistical 
Distributions: The Generalized Lambda Distribution and Generalised Bootstrap 
Methods. New York, Chapman and Hall. }
\author{ Steve Su }
\seealso{ \code{\link{fun.RMFMKL.ml}}, \code{\link{fun.RMFMKL.lm}}, 
\code{\link{fun.RMFMKL.qs}}, \code{\link{fun.data.fit.ml}} 
\code{\link{fun.data.fit.lm}}, \code{\link{fun.data.fit.qs}},  
\code{\link{fun.data.fit.mm}} }
\examples{
## Fitting the normal distribution
# fun.RMFMKL.mm(data=rnorm(1000,2,3),fmkl.init=c(-0.25,1.5),leap=3)
}
\keyword{smooth} 