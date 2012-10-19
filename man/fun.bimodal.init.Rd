\name{fun.bimodal.init}
\alias{fun.bimodal.init}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Finds the initial values for optimisation in fitting the bimodal 
generalised lambda distribution. }
\description{
After classifying the data using \code{\link{fun.class.regime.bi}}, this 
function evaluates the temporary or initial solutions
by estimating each part of the bimodal distribution using the maximum likelihood 
estimation and starship method.
These initial solutions are then passed onto \code{\link{fun.bimodal.fit.ml}} or 
\code{\link{fun.bimodal.fit.pml}} to obtain the final fits.
}
\usage{
fun.bimodal.init(data1, data2, rs.leap1, fmkl.leap1, rs.init1, fmkl.init1, 
rs.leap2, fmkl.leap2, rs.init2, fmkl.init2,fun1="runif.sobol",
fun2="runif.sobol",no=10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data1}{ The first data obtained by the clustering algorithm. }
  \item{data2}{ The second data obtained by the clustering algorithm. }
  \item{rs.leap1}{ Scrambling (0,1,2,3) for the sobol sequence for the first RS 
  distribution fit. See scrambling/leap argument for \code{\link{runif.sobol}}, 
  \code{\link{runif.halton}} or \code{\link{QUnif}}. }
  \item{fmkl.leap1}{ Scrambling (0,1,2,3) for the sobol sequence for the first 
  FMKL distribution fit. See scrambling/leap argument for 
  \code{\link{runif.sobol}}, \code{\link{runif.halton}} 
  or \code{\link{QUnif}}. }
  \item{rs.init1}{ Inititial values (lambda3 and lambda4) for the first RS 
  generalised lambda distribution. \code{c(-1.5,1.5)} tends to work well. }
  \item{fmkl.init1}{ Inititial values (lambda3 and lambda4) for the first FMKL 
  generalised lambda distribution. \code{c(-0.25,1.5)} tends to work well }
  \item{rs.leap2}{ Scrambling (0,1,2,3) for the sobol sequence for the second RS 
  distribution fit. See scrambling/leap argument for \code{\link{runif.sobol}}, 
  \code{\link{runif.halton}} or \code{\link{QUnif}}. }
  \item{fmkl.leap2}{ Scrambling (0,1,2,3) for the sobol sequence for the second 
  FMKL distribution fit. See scrambling/leap argument for 
  \code{\link{runif.sobol}}, \code{\link{runif.halton}} or 
  \code{\link{QUnif}}. }
  \item{rs.init2}{ Inititial values (lambda3 and lambda4) for the second RS 
  generalised lambda distribution. \code{c(-1.5,1.5)} tends to work well. }
  \item{fmkl.init2}{ Inititial values (lambda3 and lambda4) for the second FMKL 
  generalised lambda distribution. \code{c(-0.25,1.5)} tends to work well }
  \item{fun1}{ A character string of either  \code{"runif.sobol"} (default), 
  \code{"runif.halton"} or \code{"QUnif"} for the first distribution fit. }
  \item{fun2}{ A character string of either  \code{"runif.sobol"} (default), 
  \code{"runif.halton"} or \code{"QUnif"} for the second distribution fit. }
  \item{no}{ Number of initial random values to find the best initial values 
  for optimisation. }
}
\details{
All three methods of fitting (RPRS, RMFMKL and STAR) will be given for each 
part of the bimodal distribution.
}
\value{
  \item{prop}{ Proportion of the number of observations in the first data in 
  relation to the entire data. }
  \item{first.fit}{ A matrix comprising the parameters of GLD obtained from 
  RPRS, RMFMKL and STAR for the first dataset. }
  \item{second.fit}{ A matrix comprising the parameters of GLD obtained from 
  RPRS, RMFMKL and STAR for the second dataset. }
}
\references{ Su (2007). Fitting Single and Mixture of Generalized Lambda 
Distributions to Data via Discretized and Maximum Likelihood Methods: 
GLDEX in R. Journal of Statistical Software: *21* 9.
}
\author{ Steve Su }
\note{ This is not designed to be called by the end user explicitly, 
the difficulties with RPRS parameterisation should be
noted by the users. 
}
\seealso{ \code{\link{fun.class.regime.bi}},\code{\link{fun.bimodal.fit.pml}},
\code{\link{fun.bimodal.fit.ml}} }
\examples{
## Split the first column of the faithful data into two using fun.class.regime.bi
# faithful1.mod<-fun.class.regime.bi(faithful[,1], 0.1, clara)

## Save the datasets
# qqqq1.faithful1.cc1<-faithful1.mod$data.a
# qqqq2.faithful1.cc1<-faithful1.mod$data.b

## Find the initial values for secondary optimisation.
# result.faithful1.init1<-fun.bimodal.init(data1=qqqq1.faithful1.cc1,
# data2=qqqq2.faithful1.cc1, rs.leap1=3,fmkl.leap1=3,rs.init1 = c(-1.5, 1.5), 
# fmkl.init1 = c(-0.25, 1.5), rs.leap2=3,fmkl.leap2=3,rs.init2 = c(-1.5, 1.5), 
# fmkl.init2 = c(-0.25, 1.5))

## These initial values are then passed onto fun,bimodal.fit.ml to obtain the 
## final fits.
}
\keyword{smooth}
