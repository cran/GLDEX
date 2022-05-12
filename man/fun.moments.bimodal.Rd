\name{fun.moments.bimodal}
\alias{fun.moments.bimodal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Finds the moments of fitted mixture of 
generalised lambda distribution by simulation. }
\description{
This functions compute the mean, variance, skewness and kurtosis of the fitted 
generalised lambda distribution mixtures using Monte
Carlo simulation.
}
\usage{
fun.moments.bimodal(result1, result2, prop1, prop2, len = 1000, 
no.test = 1000, param1, param2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{result1}{ A vector comprising four values for the first generalised 
  lambda distribution. }
  \item{result2}{ A vector comprising four values for the second generalised 
  lambda distribution. }
  \item{prop1}{ Proportion of the first generalised lambda distribution }
  \item{prop2}{ 1-prop1, this can be left unspecified. }
  \item{len}{ Length of object for each simulation run. }
  \item{no.test}{ Number of simulation run. }
  \item{param1}{ This can be \code{"rs"} or \code{"fmkl"}, specifying the type 
  of the first generalised lambda distribution. }
  \item{param2}{ This can be \code{"rs"} or \code{"fmkl"}, specifying the type 
  of the second generalised lambda distribution. }
}
\details{
There is also a theoretical computation of the moments in 
\code{\link{fun.theo.bi.mv.gld}}, it should be noted
that the theoretical moments may not exist. The length of object in \code{len} 
means how many observations should
be generated in each simulation run, with the number of simulation runs governed 
by \code{no.test}.
}
\value{
A matrix with four columns showing the mean, variance, skewness and kurtosis 
of the fitted generalised lambda distribution mixtures using Monte
Carlo simulation. Each row represents a simulation run.
}
\author{ Steve Su }
\seealso{ \code{\link{fun.theo.bi.mv.gld}}, \code{\link{fun.simu.bimodal}},
\code{\link{fun.rawmoments}} }
\examples{
\donttest{
# Fitting the first column of the Old Faithful Geyser data
 fit1<-fun.auto.bimodal.ml(faithful[,1],init1.sel="rmfmkl",init2.sel="rmfmkl",
 init1=c(-0.25,1.5),init2=c(-0.25,1.5),leap1=3,leap2=3)

# After fitting compute the monte carlo moments using fun.moments.bimodal
 fun.moments.bimodal(fit1$par[1:4],fit1$par[5:8],prop1=fit1$par[9],
 param1="fmkl",param2="fmkl")

# It is also possible to compare this with the moments of the original dataset:
 fun.moments(faithful[,1])
}
}
\keyword{univar}