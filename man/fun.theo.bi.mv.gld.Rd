\name{fun.theo.bi.mv.gld}
\alias{fun.theo.bi.mv.gld}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculates the theoretical mean, variance, skewness and kurtosis for  
mixture of two generalised lambda distributions. }
\description{
This is the bimodal counterpart for \code{\link{fun.comp.moments.ml.2}} and 
\code{\link{fun.comp.moments.ml}}. }
\usage{
fun.theo.bi.mv.gld(L1, L2, L3, L4, param1, M1, M2, M3, M4, param2, p1,
normalise="N")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{L1}{ Location parameter of the first generalised lambda distribution. 
  Or all the parameters of mixture distribution in the form of
  c(L1,L2,L3,L4,M1,M2,M3,M4,p), you still must specify param1 and param2. }
  \item{L2}{ Scale parameter of the first generalised lambda distribution. }
  \item{L3}{ First shape parameter of the first generalised lambda 
  distribution. }
  \item{L4}{ Second shape parameter of the first generalised lambda 
  distribution. }
  \item{param1}{ \code{"rs"} or \code{"fmkl"} specifying the type of the first 
  generalised lambda distribution. }
  \item{M1}{ Location parameter of the second generalised lambda distribution }
  \item{M2}{ Scale parameter of the second generalised lambda distribution. }
  \item{M3}{ First shape parameter of the second generalised lambda 
  distribution. }
  \item{M4}{ Second shape parameter of the second generalised lambda 
  distribution. }
  \item{param2}{ \code{"rs"} or \code{"fmkl"} specifying the type of the 
  second generalised lambda distribution. }
  \item{p1}{ Proportion of the first generalisd lambda distribution. }
  \item{normalise}{ "Y" if you want kurtosis to be calculated with reference
  to kurtosis = 0 under Normal distribution. }
}

\value{
A vector showing the theoretical mean, variance, skewness and kurtosis for
mixture of two generalised lambda distributions. }
\references{  Freimer, M., Mudholkar, G. S., Kollia, G. & Lin, C. T. (1988), A
     study of the generalized tukey lambda family, Communications in
     Statistics - Theory and Methods *17*, 3547-3567.

     Karian, Zaven A. and Dudewicz, Edward J. (2000), Fitting
     statistical distributions: the Generalized Lambda Distribution and
     Generalized Bootstrap methods, Chapman & Hall

     Ramberg, J. S. & Schmeiser, B. W. (1974), An approximate method
     for generating asymmetric random variables, Communications of the
     ACM *17*, 78-82. }
\author{ Steve Su }
\note{ The theoretical moments may not always exist for generalised lambda 
distributions. }
\seealso{ \code{\link{fun.moments.bimodal}}, \code{\link{fun.simu.bimodal}}, 
\code{\link{fun.rawmoments}} }
\examples{
\donttest{
# Fits the Old Faithful geyser data (first column) using the maximum 
# likelihood.
 fit1<-fun.auto.bimodal.ml(faithful[,1],init1.sel="rmfmkl",init2.sel="rmfmkl",
 init1=c(-0.25,1.5),init2=c(-0.25,1.5),leap1=3,leap2=3)

# Find the theoretical moments of the fit
 fun.theo.bi.mv.gld(fit1$par[1],fit1$par[2],fit1$par[3],fit1$par[4],"fmkl",
 fit1$par[5],fit1$par[6],fit1$par[7],fit1$par[8],"fmkl",fit1$par[9])

# Compare this with the empirical moments from the data set.
 fun.moments(faithful[,1])
}}
\keyword{univar}