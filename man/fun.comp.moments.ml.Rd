\name{fun.comp.moments.ml}
\alias{fun.comp.moments.ml}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Compare the moments of the data and the fitted univariate generalised 
lambda distribution. }
\description{
After fitting the distribution, it is often desirable to see whether the 
moments of the data matches with the fitted distribution. This function computes 
the theoretical and actual moments especially for \code{fun.data.fit.ml} 
function output.
}
\usage{
fun.comp.moments.ml(theo.obj, data, name = "ML")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theo.obj}{ Fitted distribution parameters, usually output from 
  \code{\link{fun.data.fit.ml}} }
  \item{data}{ Data set used }
  \item{name}{ Naming the method used in fitting the distribution, by default 
  this is "ML". }
}

\value{
 \item{r.mat}{ A matrix showing the mean, variance, skewness and kurtosis of the 
 fitted distribution in comparison to the data set. }
 \item{eval.mat}{ Absolute difference in each of the four moments from the data 
 under each of the distibutional fits. }
}
\author{ Steve Su }
\note{ Sometimes it is difficult to find RPRS type of fits to data set, so 
instead \code{\link{fun.comp.moments.ml.2}} is used to compare the theoretical 
moments of RMFMKL.ML and STAR methods
with respect to the dataset fitted.
}
\seealso{ \code{\link{fun.comp.moments.ml.2}} }
\examples{
\donttest{
# Generate random normally distributed observations.
 junk<-rnorm(1000,3,2)

# Fit the dataset using fun.data.ml
 fit<-fun.data.fit.ml(junk)

# Compare the resulting fits. It is usually the case the maximum likelihood 
# provides better estimation of the moments than the starship method.
 fun.comp.moments.ml(fit,junk)
}
}
\keyword{univar}