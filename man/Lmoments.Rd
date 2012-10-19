\name{Lmoments}
\alias{Lmoments}
\alias{Lcoefs}
\alias{Lmomcov}
\alias{Lmoments_calc}
\alias{Lmomcov_calc}
\concept{moments}
\concept{skewness}
\concept{kurtosis}
\concept{descriptive statistics}
\title{L-moments}
\description{
  Calculates sample L-moments, L-coefficients and covariance matrix of L-moments. 
}
\usage{
Lmoments(data,rmax=4,na.rm=FALSE,returnobject=FALSE,trim=c(0,0))
Lcoefs(data,rmax=4,na.rm=FALSE,trim=c(0,0))
Lmomcov(data,rmax=4,na.rm=FALSE)
Lmoments_calc(data,rmax=4)
Lmomcov_calc(data,rmax=4)
}
\arguments{
  \item{data}{matrix or data frame.}
  \item{rmax}{maximum order of L-moments.}
  \item{na.rm}{a logical value indicating whether 'NA' values should be
          removed before the computation proceeds.}
  \item{returnobject}{a logical value indicating whether a list object should be 
  returned instead of an array of L-moments.}
  \item{trim}{c(0,0) for ordinary L-moments and c(1,1) for trimmed (t=1) 
  L-moments}
}
\value{
\item{Lmoments}{ returns an array of L-moments containing a row for each 
variable in data, or if returnobject=TRUE, a list containing the following: }
\item{lambdas}{an array of L-moments} 
\item{ratios}{an array of mean, L-scale and L-moment ratios}
\item{trim}{the value of the parameter 'trim'}
\item{source}{a string with value "Lmoments" or "t1lmoments"}
\item{Lcoefs}{returns an array of L-coefficients (mean, L-scale, L-skewness, 
L-kurtosis, ...)} 
\item{Lmomcov}{returns the covariance matrix of L-moments or a list of 
covariance matrices if the input has multiple columns.}
\item{Lmoments_calc}{is internal function.}
\item{Lmomcov_calc}{is internal function.}}

\note{Functions \code{Lmoments} and \code{Lcoefs} calculate trimmed L-moments 
if you specify \code{trim=c(1,1)}.}

\references{
Karvanen, J. and A. Nuutinen (2008). 
"Characterizing the generalized lambda distribution by L-moments." 
Computational Statistics and Data Analysis 52(4): 1971-1983.

Asquith, W. (2007). 
"L-moments and TL-moments of the generalized lambda distribution." 
Computational Statistics and Data Analysis 51(9): 4484-4496.

Elamir, E. A., Seheult, A. H. 2004. 
"Exact variance structure of sample L-moments" 
Journal of Statistical Planning and Inference 124 (2) 337-359.

Hosking, J. 1990. "L-moments: Analysis and estimation distributions using linear 
combinations of order statistics", 
Journal of Royal Statistical Society B 52, 105-124.
}
\seealso{ \code{\link{t1lmoments}} for trimmed L-moments }
\examples{

x<-rnorm(500)
Lmoments(x)

}
\author{Juha Karvanen <\email{juha.karvanen@ktl.fi}>}
\keyword{robust}
\keyword{univar}

