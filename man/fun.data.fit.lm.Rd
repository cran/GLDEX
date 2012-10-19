\name{fun.data.fit.lm}
\alias{fun.data.fit.lm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fit data using L moment matching estimation for RS and FMKL GLD }
\description{
This function fits generalised lambda distributions to data using L moment
matching method
}
\usage{
fun.data.fit.lm(data, rs.leap = 3, fmkl.leap = 3, rs.init = c(-1.5, 1.5), 
fmkl.init = c(-0.25, 1.5), FUN = "runif.sobol", no = 10000)
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
  \item{no}{ Number of initial random values to find the best initial values for 
  optimisation. }
}
\details{
This function consolidates \code{\link{fun.RPRS.lm}} and  
\code{\link{fun.RMFMKL.lm}} and gives all the fits in 
one output.
}
\value{
A matrix showing the parameters of RS and FMKL generalised lambda distributions.
}
\references{
Asquith, W. (2007). 
"L-moments and TL-moments of the generalized lambda distribution." 
Computational Statistics and Data Analysis 51(9): 4484-4496.

Karvanen, J. and A. Nuutinen (2008). 
"Characterizing the generalized lambda distribution by L-moments." 
Computational Statistics and Data Analysis 52(4): 1971-1983. }
\author{ Steve Su }
\seealso{ \code{\link{fun.RPRS.qs}}, \code{\link{fun.RMFMKL.qs}}, 
\code{\link{fun.data.fit.ml}}, \code{\link{fun.data.fit.hs}}, 
\code{\link{fun.data.fit.hs.nw}} , \code{\link{fun.data.fit.qs}}, 
\code{\link{fun.data.fit.mm}} }
\examples{
## Fitting normal(3,2) distriution using the default setting
# junk<-rnorm(50,3,2)
# fun.data.fit.lm(junk)
}