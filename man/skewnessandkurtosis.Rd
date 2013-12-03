\name{skewness and kurtosis}
\alias{skewness}
\alias{kurtosis}
\title{ Compute skewness and kurtosis statistics }
\description{ This uses the S+ version directly. }
\usage{
skewness(x, na.rm = FALSE, method = "fisher")
kurtosis(x, na.rm = FALSE, method = "fisher")
}
\arguments{
  \item{x}{ Any numerical object. Missing values \code{NA} are allowed. }
  \item{na.rm}{ Logical flag: if na.rm=TRUE, missing values are removed from x 
  before doing the computations. If \code{na.rm=FALSE} and \code{x} contains 
  missing values, then the return value is NA. }
  \item{method}{ Character string specifying the computation method. The two 
  possible values are \code{fisher} for Fisher's g1 (skewness) and g2 (kurtosis) 
  versions, and \code{moment} for the functional forms of the statistics. Only 
  the first character of the string needs to be supplied. 
}
}
\details{
The \code{moment} forms are based on the definitions of skewness and kurtosis 
for distributions; these forms should be used when resampling (bootstrap or 
jackknife). The "fisher" forms correspond to the usual "unbiased" definition of 
sample variance, though in the case of skewness and kurtosis exact unbiasedness 
is not possible. 
}
\value{
A single value of skewness or kurtotis. 

If y = x - mean(x), then the "moment" method computes the skewness value as 
mean(y\eqn{\mbox{\textasciitilde}}{^}3)/mean(y\eqn{\mbox{\textasciitilde}}{^}2)
\eqn{\mbox{\textasciitilde}}{^}1.5 
and the kurtosis value as mean(y\eqn{\mbox{\textasciitilde}}{^}4)/mean(y
\eqn{\mbox{\textasciitilde}}{^}2)\eqn{\mbox{\textasciitilde}}{^}2 - 3. 
To see the "fisher" calculations, print out the functions. 
}
\author{ Splus }
\seealso{ \code{var} }
\examples{
x <- runif(30) 
skewness(x) 
skewness(x, method="moment") 
kurtosis(x) 
kurtosis(x, method="moment") 
}
\keyword{univar} 