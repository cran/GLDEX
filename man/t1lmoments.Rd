\name{t1lmoments}
\alias{t1lmoments}
\concept{moments}
\concept{skewness}
\concept{kurtosis}
\concept{descriptive statistics}
\concept{heavy tails}
\title{Trimmed L-moments}
\description{
Calculates sample trimmed L-moments with trimming parameter 1. }
\usage{
t1lmoments(data,rmax=4)
}
\arguments{
  \item{data}{matrix or data frame.}
    \item{rmax}{maximum order of trimmed L-moments.}
}
\value{array of trimmed L-moments (trimming parameter = 1) up to order 4 
containing a row for each variable in data.}

\note{Functions \code{link{Lmoments}} and \code{link{Lcoefs}} calculate 
trimmed L-moments if you specify \code{trim=c(1,1)}.}

\references{Karvanen, J. and A. Nuutinen (2008). 
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
\seealso{\code{\link{Lmoments}} for L-moments }
\examples{

x<-rnorm(500)
t1lmoments(x);

}
\author{Juha Karvanen \email{juha.karvanen@ktl.fi}}
\keyword{robust}
\keyword{univar}


