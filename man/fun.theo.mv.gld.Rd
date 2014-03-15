\name{fun.theo.mv.gld}
\alias{fun.theo.mv.gld}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Find the theoretical first four moments of the generalised lambda 
distribution. }
\description{
  Computes the "mean","variance","skewness","kurtosis" statistics from a 
  given generalised lambda distribution.
}
\usage{
fun.theo.mv.gld(L1, L2, L3, L4, param, normalise="N")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{L1}{ Lambda 1. Or c(Lambda 1,Lambda 2,Lambda 3,Lambda 4). }
  \item{L2}{ Lambda 2. }
  \item{L3}{ Lambda 3. }
  \item{L4}{ Lambda 4. }
  \item{param}{ "rs" or "fmkl" or "fkml" }
  \item{normalise}{ "Y" if you want kurtosis to be calculated with reference
  to kurtosis = 0 under Normal distribution.
}
}

\value{
A vector listing the values of mean, variance, skewness and kurtosis.
}
\references{ Freimer, M., Mudholkar, G. S., Kollia, G. & Lin, C. T. (1988),  A
     study of the generalized tukey lambda family, Communications  in
     Statistics - Theory and Methods *17*, 3547-3567.

     Gilchrist, Warren G. (2000), Statistical Modelling with Quantile
     Functions, Chapman & Hall

     Karian, Z.A., Dudewicz, E.J., and McDonald, P. (1996), The
     extended  generalized lambda distribution system for fitting
     distributions to data:  history, completion of theory, tables,
     applications, the ``Final Word''  on Moment fits, Communications
     in Statistics - Simulation and Computation  *25*, 611-642.

     Karian,  Zaven A. and Dudewicz, Edward J. (2000), Fitting
     statistical distributions: the Generalized Lambda Distribution and
     Generalized Bootstrap methods, Chapman & Hall
 }
\author{ Steve Su }
\note{ Sometimes the theoretical moments may not exist, in those cases, 
\code{NA} is returned.
}
\seealso{ \code{\link{fun.comp.moments.ml}}, 
\code{\link{fun.comp.moments.ml.2}}, \code{\link{fun.lm.theo.gld}} }
\examples{
fun.theo.mv.gld(1, 2, 3, 4, "rs")
fun.theo.mv.gld(1, 2, 3, 4, "fmkl")
}
\keyword{univar} 