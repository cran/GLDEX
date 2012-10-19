\name{ks.gof}
\alias{ks.gof}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Kolmogorov-Smirnov test }
\description{
Performs one or two sample Kolmogorov-Smirnov tests.
}
\usage{
ks.gof(x, y, ..., alternative = c("two.sided", "less", "greater"), 
exact = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A numeric vector of data values. }
  \item{y}{ Either a numeric vector of data values, or a character string naming 
  a distribution function. }
  \item{\dots}{ Parameters of the distribution specified (as a character string) 
  by 'y'. }
  \item{alternative}{ Indicates the alternative hypothesis and must be one of 
  \code{"two.sided"} (default),\code{"less"}, or \code{"greater"}.  
  You can specify just the initial letter. }
  \item{exact}{ \code{NULL} or a logical indicating whether an exact p-value 
  should be computed.  See Details for the meaning of \code{NULL}. 
  Not used for the one-sided two-sample case. }
}
\details{

If \code{y} is numeric, a two-sample test of the null hypothesis that
\code{x} and \code{y} were drawn from the same continuous distribution is
performed.

Alternatively, \code{y} can be a character string naming a continuous
distribution function.  In this case, a one-sample test is carried
out of the null that the distribution function which generated \code{x}
is distribution \code{y} with parameters specified by "...".

The possible values \code{"two.sided"} (default),\code{"less"}
and \code{"greater"} of \code{alternative} specify the null hypothesis that the 
true distribution function of \code{x} is equal to, not less than or not
greater than the hypothesized distribution function (one-sample case) or the 
distribution function of \code{y} (two-sample case),
respectively.

Exact p-values are not available for the one-sided two-sample case, or in the 
case of ties.  \code{exact = NULL} (the default), an exact p-value is computed 
if the sample size if less than 100 in the one-sample case, and if the product 
of the sample sizes is less than 10000 in the two-sample case.  Otherwise, 
asymptotic distributions are used whose approximations may be inaccurate in
small samples.  In the one-sample two-sided case, exact p-values are obtained as 
described in Marsaglia, Tsang & Wang (2003).  The formula of Birnbaum & Tingey 
(1951) is used for the one-sample one-sided case.

If a single-sample test is used, the parameters specified in "..." must be 
pre-specified and not estimated from the data. There is some more refined 
distribution theory for the KS test with estimated parameters 
(see Durbin, 1973), but that is not implemented in \code{ks.gof}.

}
\value{
   A list with class \code{"htest"} containing the following components: 
  \item{statistic}{ Value of test statistics. }
  \item{p.value}{ P-value. }
  \item{alternative}{Character string describing the alternative hypothesis. }
  \item{method}{ Character string indicating what type of test was performed. }  
  \item{data.name}{ Character string giving the name(s) of the data. }
}
\references{ 
Z. W. Birnbaum & Fred H. Tingey (1951), One-sided confidence
contours for probability distribution functions. The Annals of
Mathematical Statistics, *22*/4, 592-596.

William J. Conover (1971), Practical nonparametric statistics.
New York: John Wiley & Sons. Pages 295-301 (one-sample
"Kolmogorov" test), 309-314 (two-sample "Smirnov" test).

Durbin, J. (1973) Distribution theory for tests based on the
sample distribution function.  SIAM.

George Marsaglia, Wai Wan Tsang & Jingbo Wang (2003), Evaluating
Kolmogorov's distribution. Journal of Statistical Software,
*8*/18. <URL: http://www.jstatsoft.org/v08/i18/>.
}
\author{ R }
\note{ This function handle ties by jittering, adding a very small uniform 
random number generated from the minimal value of the data set divided by 1e+08 
to minimal value divided by 1e+07.
}
\seealso{ \code{\link{ks.test}} }
\examples{
x <- rnorm(50)
y <- runif(30)
# Do x and y come from the same distribution?
ks.gof(x, y)
# Does x come from a shifted gamma distribution with shape 3 and rate 2?
ks.gof(x+2, "pgamma", 3, 2) # two-sided, exact
ks.gof(x+2, "pgamma", 3, 2, exact = FALSE)
ks.gof(x+2, "pgamma", 3, 2, alternative = "gr")
}
\keyword{htest} 
