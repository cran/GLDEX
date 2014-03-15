\name{fun.RMFMKL.qs}
\alias{fun.RMFMKL.qs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fit FMKL generalised lambda distribution to data set using quantile 
matching }
\description{
  This function fits FMKL generalised lambda distribution to data set 
  using quantile matching
}
\usage{
fun.RMFMKL.qs(data, fmkl.init = c(-0.25, 1.5), leap = 3, FUN = "runif.sobol", 
trial.n = 100, len = 1000, type = 7, no = 10000)
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
This function provides quantile matching fitting scheme for FMKL GLD. Note this 
function can fail if there are no defined percentiles from the data set or if 
the initial values do not lead to a valid FMKL generalised lambda distribution.
}
\value{
A vector representing four parametefmkl of the FMKL generalised lambda 
distribution.
}
\references{ 
Su (2008). Fitting GLD to data via quantile matching method. (Book chapter to
appear) }
\author{ Steve Su }
\seealso{ \code{\link{fun.RMFMKL.ml}}, \code{\link{fun.RMFMKL.lm}}, 
\code{\link{fun.RMFMKL.mm}}, \code{\link{fun.data.fit.ml}} 
\code{\link{fun.data.fit.lm}}, \code{\link{fun.data.fit.qs}},  
\code{\link{fun.data.fit.mm}} }
\examples{
## Fitting the normal distribution
# fun.RMFMKL.qs(data=rnorm(1000,2,3),fmkl.init=c(-0.25,1.5),leap=3)
}
\keyword{smooth} 