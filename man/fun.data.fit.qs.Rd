\name{fun.data.fit.qs}
\alias{fun.data.fit.qs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fit data using quantile matching estimation for RS and FMKL GLD }
\description{
This function fits generalised lambda distributions to data using quantile
matching method
}
\usage{
fun.data.fit.qs(data, rs.leap = 3, fmkl.leap = 3, rs.init = c(-1.5, 1.5), 
fmkl.init = c(-0.25, 1.5), FUN = "runif.sobol", trial.n = 100, len = 1000, 
type = 7, no = 10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ Dataset to be fitted. }
  \item{rs.leap}{ Scrambling (0,1,2,3) for the sobol sequence for the RPRS 
  distribution fit. See scrambling/leap argument for \code{\link{runif.sobol}}, 
  \code{\link{runif.halton}} or \code{\link{QUnif}}. }
  \item{fmkl.leap}{ Scrambling (0,1,2,3) for the sobol sequence for the RMFMKL 
  distribution fit. See scrambling/leap argument for \code{\link{runif.sobol}}, 
  \code{\link{runif.halton}} or \code{\link{QUnif}}. }
  \item{rs.init}{ Inititial values (lambda3 and lambda4) for the RS generalised 
  lambda distribution. }
  \item{fmkl.init}{ Inititial values (lambda3 and lambda4) for the FMKL 
  generalised lambda distribution. }
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
This function consolidates \code{\link{fun.RPRS.qs}} and  
\code{\link{fun.RMFMKL.qs}} and gives all the fits in 
one output.
}
\value{
A matrix showing the parameters of RS and FMKL generalised lambda distributions.
}
\references{ Su (2008). Fitting GLD to data via quantile matching method. 
(Book chapter to appear) }
\author{ Steve Su }
\note{ RPRS can sometimes fail if it is not possible to calculate the 
percentiles of the data set. This usually happens
when the number of data point is small. }
\seealso{ \code{\link{fun.RPRS.qs}}, \code{\link{fun.RMFMKL.qs}}, 
\code{\link{fun.data.fit.ml}}, \code{\link{fun.data.fit.hs}}, 
\code{\link{fun.data.fit.hs.nw}} , \code{\link{fun.data.fit.mm}}, 
\code{\link{fun.data.fit.lm}} }
\examples{
## Fitting normal(3,2) distriution using the default setting
# junk<-rnorm(50,3,2)
# fun.data.fit.qs(junk)
}
\keyword{smooth}
