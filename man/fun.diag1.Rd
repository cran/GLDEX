\name{fun.diag1}
\alias{fun.diag1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Diagnostic function for theoretical distribution fits through the 
resample Kolmogorov-Smirnoff tests }
\description{
This function is primarily designed to be used for testing the fitted 
distribution with reference to a theoretical distribution. It is also
tailored for output obtained from the \code{\link{fun.data.fit.ml}} function.
}
\usage{
fun.diag1(result, test, no.test = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{result}{ Output from \code{\link{fun.data.fit.ml}} function. }
  \item{test}{ Simulated observations from theoretical distribution, the length 
  should be no.test\eqn{\mbox{\textasciitilde}}{^}2. }
  \item{no.test}{ Number of times to do the KS tests. }
}
\value{
A vector showing the number of times the KS p-value is greater than 0.05 for 
each of the distribution fit strategy.
}
\references{ Su, S. (2005). A Discretized Approach to Flexibly Fit Generalized 
Lambda Distributions to Data. Journal of Modern Applied Statistical 
Methods (November): 408-424. 

Su, S. (2007). Numerical Maximum Log Likelihood Estimation for Generalized 
Lambda Distributions. Journal of Computational statistics and data analysis 
51(8) 3983-3998. 

Su (2007). Fitting Single and Mixture of Generalized Lambda Distributions to 
Data via Discretized and Maximum Likelihood Methods: GLDEX in R. 
Journal of Statistical Software: *21* 9. }
\author{ Steve Su }
\note{ If there are ties, jittering is used in \code{\link{ks.gof}}.
}
\seealso{ \code{\link{fun.diag2}}, \code{\link{fun.diag.ks.g}}, 
\code{\link{fun.diag.ks.g.bimodal}} }
\examples{
## Fits a Weibull 5,2 distribution:
# weibull.approx.ml<-fun.data.fit.ml(rweibull(1000,5,2))

## Compute the resample K-S test results.
# fun.diag1(weibull.approx.ml, rweibull(100000, 5, 2))
}
\keyword{htest}