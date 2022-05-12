\name{fun.diag2}
\alias{fun.diag2}

\title{ Diagnostic function for empirical data distribution fits through the 
resample Kolmogorov-Smirnoff tests }
\description{
This function is primarily designed to be used for testing the fitted 
distribution with reference to an empirical data. It is also
tailored for output obtained from the \code{\link{fun.data.fit.ml}} function.
}
\usage{
fun.diag2(result, data, no.test = 1000, len=100, alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{result}{ Output from \code{\link{fun.data.fit.ml}} function. }
  \item{data}{ Observations in which the distribution was fitted upon. }
  \item{no.test}{ Number of times to do the KS tests. }
  \item{len}{ Number of observations to sample from the data. This is also the 
  number of observations sampled from the fitted distribution in each KS test. }
  \item{alpha}{ Significance level of KS test. }
}
\value{
A vector showing the number of times the KS p-value is greater than alpha for 
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
\seealso{ \code{\link{fun.diag1}}, \code{\link{fun.diag.ks.g}}, 
\code{\link{fun.diag.ks.g.bimodal}} }
\examples{
\donttest{
# Fits a Normal 3,2 distribution:
 junk<-rnorm(1000,3,2)
 fit<-fun.data.fit.ml(junk)

# Compute the resample K-S test results.
 fun.diag2(fit,junk)
}}
\keyword{htest}