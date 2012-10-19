\name{fun.comp.moments.ml.2}
\alias{fun.comp.moments.ml.2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Compare the moments of the data and the fitted univariate generalised 
lambda distribution. Specialised funtion designed for RMFMKL.ML and STAR 
methods. }
\description{
After fitting the distribution, it is often desirable to see whether the 
moments of the data matches with the fitted distribution. This function computes 
the theoretical and actual moments for the FMKL GLD maximum likelihood 
estimation and starship method.
}
\usage{
fun.comp.moments.ml.2(theo.obj, data, name = "ML")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theo.obj}{ Fitted distribution parameters, there should be two sets, 
  both FMKL GLD. }
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
\note{ To compare all three fits under \code{\link{fun.data.fit.ml}} see 
\code{\link{fun.comp.moments.ml}} function.
}
\seealso{ \code{\link{fun.comp.moments.ml}} }
\examples{

## Generate random normally distributed observations.
# junk<-rnorm(1000,3,2)

## Fit the dataset using fun.data.ml
# fit<-cbind(fun.RMFMKL.ml(junk),starship(junk)$lambda)

## Compare the resulting fits. It is usually the case the maximum likelihood 
## provides better estimation of the moments than the starship method.
# fun.comp.moments.ml.2(fit,junk)
}
\keyword{univar}