\name{fun.RMFMKL.ml.m}
\alias{fun.RMFMKL.ml.m}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit RS generalised lambda distribution to data set using maximum 
likelihood estimation
}
\description{
This function fits FMKL generalised lambda distribution to data set using 
maximum likelihood estimation using faster implementation through C programming
}
\usage{
fun.RMFMKL.ml.m(data, fmkl.init = c(-0.25, 1.5), leap = 3, FUN = "runif.sobol", 
no = 10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ Dataset to be fitted }
  \item{fmkl.init}{ Initial values for FMKL distribution optimization, 
  \code{c(-0.25,1.5)} tends to work well. }
  \item{leap}{ See scrambling argument in \code{\link{fun.gen.qrn}}. }
  \item{FUN}{ A character string of either  \code{"runif.sobol"} (default), \code{"runif.sobol.owen"}, 
  \code{"runif.halton"} or \code{"QUnif"}. }
  \item{no}{ Number of initial random values to find the best initial values for 
  optimisation. }
}
\details{
This function provides one of the definitive fit to data set using generalised 
lambda distributions.
}
\value{
A vector representing four parameters of the FMKL generalised lambda 
distribution.
}
\references{ Su, S. (2007). Numerical Maximum Log Likelihood Estimation for 
Generalized Lambda Distributions. Journal of Computational statistics and 
data analysis 51(8) 3983-3998. 

Su (2007). Fitting Single and Mixture of Generalized Lambda Distributions to 
Data via Discretized and Maximum Likelihood Methods: GLDEX in R. 
Journal of Statistical Software: *21* 9.
}
\author{ Steve Su }
\seealso{ \code{\link{fun.RMFMKL.ml}} }
\examples{
\donttest{
# Fitting the normal distribution
 fun.RMFMKL.ml.m(data=rnorm(1000,2,3),fmkl.init=c(-0.25,1.5),leap=3)
}}
\keyword{smooth} 
