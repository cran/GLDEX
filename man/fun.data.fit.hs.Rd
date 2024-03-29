\name{fun.data.fit.hs}
\alias{fun.data.fit.hs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fit RS and FMKL generalised distributions to data using discretised 
approach with weights. }
\description{
This function fits RS and FMKL generalised distribution to data using 
discretised approach with weights. It is designed to act as a smoother device 
rather than a definitive fit. 
}
\usage{
fun.data.fit.hs(data, rs.default = "Y", fmkl.default = "Y", rs.leap = 3, 
fmkl.leap = 3, rs.init = c(-1.5, 1.5), fmkl.init = c(-0.25, 1.5), no.c.rs = 50, 
no.c.fmkl = 50,FUN="runif.sobol", no=10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ Dataset to be fitted }
  \item{rs.default}{ If yes, this function uses the default method 
  \code{\link{fun.nclass.e}} to calculate number of classes required for the RS 
  distribution fits. }
  \item{fmkl.default}{ If yes, this function uses the default method 
  \code{\link{fun.nclass.e}} to calculate number of classes required for the 
  FMKL distribution fits. }
  \item{rs.leap}{ See scrambling argument in \code{\link{fun.gen.qrn}}. }
  \item{fmkl.leap}{ See scrambling argument in \code{\link{fun.gen.qrn}}. }
  \item{rs.init}{ Initial values for RS distribution optimization, 
  \code{c(-1.5,1.5)} tends to work well. }
  \item{fmkl.init}{ Initial values for FMKL distribution optimization,
   \code{c(-0.25,1.5)} tends to work well. }
  \item{no.c.rs}{ Number of classes or bins of histogram to be optimized over 
  for the RS GLD. This argument is ineffective if \code{default="Y"}. }
  \item{no.c.fmkl}{ Number of classes or bins of histogram to be optimized over 
  for the FMKL GLD. This argument is ineffective if \code{default="Y"}. }
  \item{FUN}{ A character string of either  \code{"runif.sobol"} (default), \code{"runif.sobol.owen"}, 
  \code{"runif.halton"} or \code{"QUnif"}. }
  \item{no}{ Number of initial random values to find the best initial values for 
  optimisation. }
}
\details{
  This function optimises the deviations of frequency of the bins to that of the 
  theoretical so it has the effect of "fitting clothes" onto the data set. The 
  user can decide the frequency of the bins they want the distribution to smooth 
  over. The resulting fit may or may not be an adequate fit from a formal 
  statistical point of view such as satisfying the goodness of fit for example, 
  but it can be useful to suggest the range of different distributions exhibited 
  by the data set. The default number of classes calculates the mean and 
  variance after categorising the data into different bins and uses the number 
  of classes that best matches the mean and variance of the original, ungrouped 
  data.The weighting is designed to accentuate the peak or the dense part of the 
  distribution and suppress the tails. 
}
\value{
 A matrix showing the four parameters of the RS and FMKL distribution fit.
}
\references{ Su, S. (2005). A Discretized Approach to Flexibly Fit Generalized 
Lambda Distributions to Data. Journal of Modern Applied Statistical 
Methods (November): 408-424. 

Su (2007). Fitting Single and Mixture of Generalized Lambda Distributions to 
Data via Discretized and Maximum Likelihood Methods: GLDEX in R. 
Journal of Statistical Software: *21* 9. }
\author{ Steve Su }
\note{ In some cases, the resulting fit may not converge, there are currently no 
checking mechanism in place to ensure global convergence.
The RPRS method can sometimes fail if there are no valid percentiles in the data 
set or if initial values do not give a valid distribution.
}
\seealso{ \code{\link{fun.RPRS.hs.nw}}, \code{\link{fun.RMFMKL.hs.nw}}, 
\code{\link{fun.RMFMKL.hs}}, \code{\link{fun.RPRS.hs}}, 
\code{\link{fun.data.fit.hs.nw}}, \code{\link{fun.data.fit.ml}} }
\examples{
\donttest{
# Fitting normal(3,2) distriution using the default setting
 junk<-rnorm(1000,3,2)
 fun.data.fit.hs(junk)
}}
\keyword{smooth}
