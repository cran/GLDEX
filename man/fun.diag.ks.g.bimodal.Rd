\name{fun.diag.ks.g.bimodal}
\alias{fun.diag.ks.g.bimodal}
\title{ Compute the simulated Kolmogorov-Smirnov tests for the bimodal dataset }
\description{
This function counts the number of times the p-value exceed 0.05 for the null 
hypothesis that the observations simulated tom fitted distribution is the same 
as the observations simulated from the bimodal data set. 
}
\usage{
fun.diag.ks.g.bimodal(result1, result2, prop1, prop2, data, no.test = 1000, 
len = floor(0.9 * length(data)), param1, param2)
}
\arguments{
  \item{result1}{ A vector representing the four parameters of the first 
  generalised lambda distribution. }
  \item{result2}{ A vector representing the four parameters of the second 
  generalised lambda distribution. }
  \item{prop1}{ Proportion of the first distribution fitted to the bimodal 
  dataset. }
  \item{prop2}{ Proportion of the second distribution fitted to the bimodal 
  dataset. }
  \item{data}{ The bimodal dataset. }
  \item{no.test}{ Total number of tests required. }
  \item{len}{ Number of data to sample. }
  \item{param1}{ Type of first generalised lambda distribution, can be 
  \code{"rs"} or \code{"fmkl"}. }
  \item{param2}{ Type of second generalised lambda distribution, can be 
  \code{"rs"} or \code{"fmkl"}. }
}
\value{
A numerical value representing number of times the p-value exceeds 0.05.
}
\references{ Stephens, M. A. (1986). Tests based on EDF statistics. In Goodness-
of-Fit Techniques. D'Agostino, R. B. and Stevens, M. A., eds. New York: Marcel 
Dekker.

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
\seealso{ \code{\link{fun.diag.ks.g}} }
\examples{
## Fit the faithful[,1] data from the MASS library
# fit1<-fun.auto.bimodal.ml(faithful[,1],init1.sel="rprs",init2.sel="rmfmkl",
# init1=c(-1.5,1,5),init2=c(-0.25,1.5),leap1=3,leap2=3)
## Run diagnostic KS tests
# fun.diag.ks.g.bimodal(fit1$par[1:4],fit1$par[5:8],prop1=fit1$par[9],
# data=faithful[,1],param1="rs",param2="fmkl")
}
\keyword{htest}