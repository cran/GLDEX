\name{fun.class.regime.bi}
\alias{fun.class.regime.bi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Classifies data into two groups using a clustering regime. }
\description{
This function is primarily designed to split a bimodal data vector into two 
groups to allow the fitting of mixture generalised lambda distributions. 
}
\usage{
fun.class.regime.bi(data, perc.cross, fun.cross)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ Data to be classified into two groups. }
  \item{perc.cross}{ Percentage of cross over from one data to the other, 
  usually set at 1\% }
  \item{fun.cross}{ Any clustering function such as \code{link{clara}}, 
  \code{\link{pam}}, \code{\link{fanny}} can be used here. Or a logical 
  vector indicating how data should be split. }
}
\details{
This function is part of the routine mixture fitting procedure provided in this 
package. The \code{perc.cross} argument or percentage of cross over is designed 
to allow the use of maximum likelihood estimation via EM algorithm for fitting 
bimodal data. When this is invoked, it will ensure both part of the data will 
contain both the minmum and maximum of the data set as well as a proportion (
specified in perc.cross argument) of observations from each other. If 1\% is 
required, then data.a will contains 1\% of the data.b and vice versa after the 
full data set has been classified into \code{data.a} and \code{data.b} by the 
\code{fun.cross} classification regime.
}
\value{
 
  \item{data.a}{ First group of data obtained by the classification algorithm. }
  \item{data.b}{ Second group of data obtained by the classification algorithm. }
}
\references{ 
Kaufman, L. and Rousseeuw, P. J. (1990). Finding Groups in Data: An Introduction 
to Cluster Analysis. Wiley, New York.

Su (2006) Maximum Log Likelihood Estimation using EM Algorithm and Partition 
Maximum Log Likelihood Estimation for Mixtures of Generalized Lambda 
Distributions. Working Paper. }
\author{ Steve Su }

\seealso{ \code{link{clara}}, \code{\link{pam}}, \code{\link{fanny}} }
\examples{
## Classify the faithful[,1] data into two categories with 10\% cross over mix.
# fun.class.regime.bi(faithful[,1],0.1,clara)

## Classify the faithful[,1] data into two categories with no mixing:
# fun.class.regime.bi(faithful[,1],0,clara)
}
\keyword{cluster}