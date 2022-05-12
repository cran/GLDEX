\name{fun.diag.ks.g}
\alias{fun.diag.ks.g}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Compute the simulated Kolmogorov-Smirnov tests for the unimodal 
dataset }
\description{
  This function counts the number of times the p-value exceed 0.05 for the null 
  hypothesis that the observations simulated from the fitted distribution is the same 
  as the observations simulated from the unimodal data set. 

}
\usage{
fun.diag.ks.g(result, data, no.test = 1000, len = floor(0.9 * length(data)), 
param, alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{result}{ A vector representing the four parameters of the generalised 
  lambda distribution.  }
  \item{data}{ The unimodal dataset. }
  \item{no.test}{ Total number of tests required. }
  \item{len}{ Number of data to sample. }
  \item{param}{ Type of the generalised lambda distribution, \code{"rs"} or 
  \code{"fmkl"}. }
  \item{alpha}{ Significance level of KS test. }
}

\value{
A numerical value representing number of times the p-value exceeds alpha.
}
\references{ Stephens, M. A. (1986). Tests based on EDF statistics. 
In Goodness-of-Fit Techniques. D'Agostino, R. B. and Stevens, M. A., eds. 
New York: Marcel Dekker.

Su, S. (2005). A Discretized Approach to Flexibly Fit Generalized Lambda 
Distributions to Data. Journal of Modern Applied Statistical Methods (November): 
408-424.

Su (2007). Nmerical Maximum Log Likelihood Estimation for Generalized Lambda 
Distributions.  Computational Statistics and Data Analysis: *51*, 8, 3983-3998.

Su (2007). Fitting Single and Mixture of Generalized Lambda Distributions to 
Data via Discretized and Maximum Likelihood Methods: GLDEX in R. 
Journal of Statistical Software: *21* 9.

}
\author{ Steve Su }
\note{
If there are ties, jittering is used in \code{\link{ks.gof}}.
}

\seealso{ \code{\link{fun.diag.ks.g.bimodal}} }

\examples{
\donttest{
# Generate 1000 random observations from Normal distribution with mean=100, 
# standard deviation=10. Save this as junk
 junk<-rnorm(1000,100,10)

# Fit junk using RPRS method via the maxmum likelihood.
 fit1<-fun.RPRS.ml(junk, c(-1.5, 1.5), leap = 3)

# Calculate the simulated KS test result:
 fun.diag.ks.g(fit1,junk,param="rs")
}}
\keyword{htest}