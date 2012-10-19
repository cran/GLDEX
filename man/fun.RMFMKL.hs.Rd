\name{fun.RMFMKL.hs}
\alias{fun.RMFMKL.hs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fit FMKL generalised distribution to data using discretised approach 
with weights. }
\description{
  This function fits FMKL generalised distribution to data using discretised 
  approach with weights. It is designed to act as a smoother device rather than 
  as a definitive fit. 
}
\usage{
fun.RMFMKL.hs(data, default = "Y", fmkl.init = c(-0.25, 1.5), no.c.fmkl = 50, 
leap = 3,FUN="runif.sobol",no=10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ Dataset to be fitted }
  \item{default}{ If yes, this function uses the default method 
  \code{\link{fun.nclass.e}} to calculate number of classes required. }
  \item{fmkl.init}{ Initial values for FMKL distribution optimization, 
  \code{c(-0.25,1.5)} tends to work well. }
  \item{no.c.fmkl}{ Number of classes or bins of histogram to be optimized over. 
  This argument is ineffective if \code{default="Y"}. }
  \item{leap}{ Scrambling (0,1,2,3) for the Sobol sequence for the distribution 
  fit. See scrambling/leap argument for \code{\link{runif.sobol}}, 
  \code{\link{runif.halton}} or \code{\link{QUnif}}. }
  \item{FUN}{ A character string of either  \code{"runif.sobol"} (default), 
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
by the data set. The default number of classes calculates the mean and variance 
after categorising the data into different bins and uses the number of classes 
that best matches the mean and variance of the original, ungrouped data. The 
weighting is designed to accentuate the peak or the dense part of the 
distribution and suppress the tails. 
}
\value{
A vector representing four parametefmkl of the FMKL generalised lambda 
distribution.
}
\references{ Su, S. (2005). A Discretized Approach to Flexibly Fit Generalized 
Lambda Distributions to Data. Journal of Modern Applied Statistical Methods 
(November): 408-424. 

Su (2007). Fitting Single and Mixture of Generalized Lambda Distributions to 
Data via Discretized and Maximum Likelihood Methods: GLDEX in R. 
Journal of Statistical Software: *21* 9.
}
\author{ Steve Su }
\note{ In some cases, the resulting fit may not converge, there are currently no 
checking mechanism in place to ensure global convergence.
}
\seealso{ \code{\link{fun.RMFMKL.hs.nw}}, \code{\link{fun.RPRS.hs.nw}}, 
\code{\link{fun.RPRS.hs}}, \code{\link{fun.data.fit.hs}}, 
\code{\link{fun.data.fit.hs.nw}} }
\examples{
## Using the default number of classes
# fun.RMFMKL.hs(data=rnorm(1000,3,2),default="Y",fmkl.init=c(-0.25,1.5),leap=3)
## Using 20 classes
# fun.RMFMKL.hs(data=rnorm(1000,3,2),default="N",fmkl.init=c(-0.25,1.5),
# no.c.fmkl=20,leap=3)
}
\keyword{smooth}
