\name{fun.data.fit.ml}
\alias{fun.data.fit.ml}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fit data using RS, FMKL maximum likelihood estimation and the FMKL 
starship method. }
\description{
This function fits generalised lambda distributions to data using RPRS, RMFMKL 
and starship methods.
}
\usage{
fun.data.fit.ml(data, rs.leap = 3, fmkl.leap = 3, rs.init = c(-1.5, 1.5), 
fmkl.init = c(-0.25, 1.5),FUN="runif.sobol",no=10000)
}
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
This function consolidates \code{\link{fun.RPRS.ml}}, 
\code{\link{fun.RMFMKL.ml}} and \code{\link{starship}} and gives all the fits in 
one output.
}
\value{
A matrix showing the parameters of generalised lambda distribution for RPRS, 
FMFKL and STAR methods.
}
\references{  
King, R.A.R. & MacGillivray, H. L. (1999), A starship method for fitting the 
generalised lambda distributions, Australian and New
Zealand Journal of Statistics, 41, 353-374

Su, S. (2007). Numerical Maximum Log Likelihood Estimation for Generalized 
Lambda Distributions. Computational statistics and data analysis 51(8) 
3983-3998. 

Su (2007). Fitting Single and Mixture of Generalized Lambda Distributions to 
Data via Discretized and Maximum Likelihood Methods: GLDEX in R. 
Journal of Statistical Software: *21* 9.
}
\author{ Steve Su }
\note{ RPRS can sometimes fail if it is not possible to calculate the 
percentiles of the data set. This usually happens
when the number of data point is small. }
\seealso{ \code{\link{fun.RPRS.ml}}, \code{\link{fun.RMFMKL.ml}}, 
\code{\link{starship}}, \code{\link{fun.data.fit.hs}}, 
\code{\link{fun.data.fit.hs.nw}} ,
\code{\link{fun.data.fit.qs}} ,
\code{\link{fun.data.fit.mm}} ,
\code{\link{fun.data.fit.lm}} }
\examples{
## Fitting normal(3,2) distriution using the default setting
# junk<-rnorm(50,3,2)
# fun.data.fit.ml(junk)
}
\keyword{smooth}
