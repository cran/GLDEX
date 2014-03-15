\name{fun.minmax.check.gld}
\alias{fun.minmax.check.gld}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Check whether the specified GLDs cover the minimum and the maximum values in a 
dataset
}
\description{
This function checks the lowest and highest quantiles of the specified GLDs 
against the 
specified dataset
}
\usage{
fun.minmax.check.gld(data, lambdas, param, lessequalmin = 1, 
greaterequalmax = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A vector of numerical dataset
}
  \item{lambdas}{
A matrix of four columns representing lambda 1 to lambda 4 of the GLD
}
  \item{param}{
Can be "rs", "fkml" or "fmkl"
}
  \item{lessequalmin}{
Can be 0 or 1
}
  \item{greaterequalmax}{
Can be 0 or 1
}
}
\details{

lessequalmin==1 means the lowest value of GLD <= minimum value of data
lessequalmin==0 means the lowest value of GLD < minimum value of data
greaterequalmin==1 means the highest value of GLD >= maximum value of data
greaterequalmin==0 means the highest value of GLD > maximum value of data

}

\value{
A vector of logical values indicating whether the specified data 
the specified GLDs cover the minimum and the maximum values in a dataset

}

\author{
Steve Su
}

\examples{
fun.minmax.check.gld(runif(100,.9,1),matrix(1:12,ncol=4),param="rs",0,0)
fun.minmax.check.gld(runif(100,.98,1),matrix(1:12,ncol=4),param="fkml",1,1)
}

\keyword{distribution}

