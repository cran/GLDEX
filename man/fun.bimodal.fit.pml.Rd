\name{fun.bimodal.fit.pml}
\alias{fun.bimodal.fit.pml}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Finds the final fits using partition maximum likelihood estimation for 
the bimodal dataset. }
\description{
This is the secondary optimization procedure to evaluate the final bimodal 
distribution fits using the partition
maximum likelihood. It usually relies on initial values found by 
\code{fun.bimodal.init} function.
}
\usage{
fun.bimodal.fit.pml(data1, data2, first.fit, second.fit, prop, param1, param2, 
selc1, selc2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data1}{ First data set, usually obtained by 
  \code{\link{fun.class.regime.bi}}. }
  \item{data2}{ Second data set, usually obtained by 
  \code{\link{fun.class.regime.bi}}. }
  \item{first.fit}{ The distribution parameters or the initial values of the 
  first distribution fit. }
  \item{second.fit}{ The distribution parameters or the initial values of the 
  second distribution fit. }
  \item{prop}{ The proportion of the data set, usually obtained from 
  \code{fun.bimodal.init}. }
  \item{param1}{ Can be either \code{rs} or \code{fmkl}, depending on the type 
  of first distribution used. }
  \item{param2}{ Can be either \code{rs} or \code{fmkl}, depending on the type 
  of second distribution used. }
  \item{selc1}{ Selection of initial values for the first distribution, can be 
  either \code{"rs"}, \code{"fmkl"} or \code{"star"}. Choose initial values from 
  RPRS (ML), RMFMKL (ML) or STAR method. }
  \item{selc2}{ Selection of initial values for the second distribution, can be 
  either \code{"rs"}, \code{"fmkl"} or \code{"star"}. Choose initial values from 
  RPRS (ML), RMFMKL (ML) or STAR method. }
} 
\details{
This function should be used in tandem with \code{\link{fun.bimodal.init}} 
function.
}
\value{
  \item{par}{ The first four numbers are the parameters of the first generalised 
  lambda distribution, the second four numbers are the parameters of the second 
  generalised lambda distribution and the last value is the proportion of the 
  first generalised lambda distribution.}
  \item{value}{ The objective value of negative likelihood obtained. }
  \item{counts}{ A two-element integer vector giving the number of calls to 
  functions. Gradient is not used in this case. }
  \item{convergence}{ An integer code. \code{0} indicates successful 
  convergence. Error codes are:

          \code{1} indicates that the iteration limit 'maxit' had been
               reached.

          \code{10} indicates degeneracy of the Nelder-Mead simplex.}

      
  \item{message}{ A character string giving any additional information returned 
  by the optimizer, or \code{NULL}. }
}
\references{ Su (2007). Fitting Single and Mixture of Generalized Lambda 
Distributions to Data via Discretized and Maximum Likelihood Methods: 
GLDEX in R. Journal of Statistical Software: *21* 9.
}
\author{ Steve Su }
\note{ There is currently no guarantee of a global convergence. }
\seealso{ \code{\link{fun.bimodal.fit.ml}}, \code{\link{fun.bimodal.init}} }
\examples{
## Extract faithful[,2] into faithful2
# faithful2<-faithful[,2]

## Uses clara clustering method
# clara.faithful2<-clara(faithful2,2)$clustering

## Save into two different objects
# qqqq1.faithful2.cc<-faithful2[clara.faithful2==1]
# qqqq2.faithful2.cc<-faithful2[clara.faithful2==2]

## Find the initial values
# result.faithful2.init<-fun.bimodal.init(data1=qqqq1.faithful2.cc,
# data2=qqqq2.faithful2.cc, rs.leap1=3,fmkl.leap1=3,rs.init1 = c(-1.5, 1.5), 
# fmkl.init1 = c(-0.25, 1.5), rs.leap2=3,fmkl.leap2=3,rs.init2 = c(-1.5, 1.5), 
# fmkl.init2 = c(-0.25, 1.5))

## Find the final fits
# result.faithful2.rsrs<-fun.bimodal.fit.pml(data1=qqqq1.faithful2.cc,
# data2=qqqq2.faithful2.cc, result.faithful2.init[[2]],
# result.faithful2.init[[3]], result.faithful2.init[[1]],param1="rs",
# param2="rs",selc1="rs",selc2="rs")

## Output
# result.faithful2.rsrs
}
\keyword{smooth}