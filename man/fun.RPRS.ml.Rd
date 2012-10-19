\name{fun.RPRS.ml}
\alias{fun.RPRS.ml}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fit RS generalised lambda distribution to data set using maximum 
likelihood estimation }
\description{
  This function fits RS generalised lambda distribution to data set using 
  maximum likelihood estimation.
}
\usage{
fun.RPRS.ml(data, rs.init = c(-1.5, 1.5), leap = 3,FUN="runif.sobol",no=10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ Dataset to be fitted }
  \item{rs.init}{ Initial values for RS distribution optimization, 
  \code{c(-1.5,1.5)} tends to work well. }
  \item{leap}{ Scrambling (0,1,2,3) for the Sobol sequence for the distribution 
  fit. See scrambling/leap argument for \code{\link{runif.sobol}}, 
  \code{\link{runif.halton}} or \code{\link{QUnif}}. }
  \item{FUN}{ A character string of either  \code{"runif.sobol"} (default), 
  \code{"runif.halton"} or \code{"QUnif"}. }
  \item{no}{ Number of initial random values to find the best initial values for 
  optimisation. }
}
\details{
This function provides one of the definitive fit to data set using generalised 
lambda distributions. Note this function can fail if there are no defined 
percentiles from the data set or if the initial values do not lead to a valid RS 
generalised lambda distribution.
}
\value{
A vector representing four parameters of the RS generalised lambda distribution.
}
\references{ Su, S. (2007). Numerical Maximum Log Likelihood Estimation for 
Generalized Lambda Distributions. Computational statistics and data analysis 
51(8) 3983-3998. 

Su (2007). Fitting Single and Mixture of Generalized Lambda Distributions to 
Data via Discretized and Maximum Likelihood Methods: GLDEX in R. 
Journal of Statistical Software: *21* 9.
}
\author{ Steve Su }
\seealso{ \code{\link{fun.RMFMKL.ml}}, \code{\link{fun.data.fit.ml}},
\code{\link{fun.data.fit.qs}} }
\examples{
## Fitting the normal distribution
# fun.RPRS.ml(data=rnorm(1000,2,3),rs.init=c(-1.5,1.5),leap=3)
}
\keyword{smooth} 
