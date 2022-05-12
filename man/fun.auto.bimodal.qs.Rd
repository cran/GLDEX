\name{fun.auto.bimodal.qs}
\alias{fun.auto.bimodal.qs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fitting mixtures of generalied lambda distribtions to data using 
quantile matching method }
\description{
This function will fit mixture of generalised lambda distributions to dataset. 
It is restricted to two generalised lambda distributions. The method of fitting 
is quantile matching method. It is a two step optimization procedure, 
each unimodal part of the bimodal distribution is modelled using quantile 
matching method. The initial values obtained are then used to maximise the 
theoretical and empirical quantile match for the entire bimodal distribution. 
It fits mixture of the form p*(f1)+(1-p)*(f2) where f1 and f2 are pdfs of the 
generalised lambda distributions.
}
\usage{
fun.auto.bimodal.qs(data, per.of.mix = 0.01, clustering.m = clara, 
init1.sel = "rprs", init2.sel = "rprs", init1=c(-1.5, 1.5), init2=c(-1.5, 1.5), 
leap1=3, leap2=3, fun1 = "runif.sobol", fun2 = "runif.sobol", trial.n = 100, 
len = 1000, type = 7, no = 10000, maxit = 5000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ A numerical vector representing the dataset. }
  \item{per.of.mix}{ Level of mix between two parts of the distribution, 
  usually 1-2\% of cross mix is sufficient. }
  \item{clustering.m}{ Clustering method used in classifying the dataset into 
  two parts. Valid arguments include clara, fanny and pam from the cluster 
  library. Default is clara. Or a logical vector specifying how data should
  be split. }
  \item{init1.sel}{ This can be \code{"rprs"} or \code{"rmfmkl"}, 
  representing the choice (RS or FMKL) of the first distribution }
  \item{init2.sel}{ This can be \code{"rprs"} or \code{"rmfmkl"},  
  representing the choice (RS or FMKL) of the second distribution }
  \item{init1}{ Inititial values lambda3 and lambda4 for the first generalised 
  lambda distribution. }
  \item{init2}{ Inititial values lambda3 and lambda4 for the second generalised 
  lambda distribution. }
  \item{leap1}{ See scrambling argument in \code{\link{fun.gen.qrn}}. }
  \item{fun1}{ A character string of either  \code{"runif.sobol"} (default), \code{"runif.sobol.owen"}, 
  \code{"runif.halton"} or \code{"QUnif"}. }
  \item{leap2}{ See scrambling argument in \code{\link{fun.gen.qrn}}. }
  \item{fun2}{ A character string of either  \code{"runif.sobol"} (default), \code{"runif.sobol.owen"}, 
  \code{"runif.halton"} or \code{"QUnif"}. }
  \item{trial.n}{ Number of evenly spaced quantile ranging from 0 to 1 to be
  used in the checking phase, to find the best set of initial values for 
  optimisation, this is intended to be lower than \code{len} to speed up the 
  fitting algorithm. Default is 100. }
  \item{len}{ Number of evenly spaced quantile ranging from 0 to 1 to be used,
  default is 1000 }
  \item{type}{ Type of quantile to be used, default is 7, see \code{quantile} }
  \item{no}{ Number of initial random values to find the best initial values 
  for optimisation. }
  \item{maxit}{ Maximum number of iterations for numerical optimisation. Default
  is 5000. }
}
\details{
The initial values that work well for RPRS are \code{c(-1.5,1.5)} and for RMFMKL 
are \code{c(-0.25,1.5)}. For scrambling, if \code{1}, \code{2} or \code{3} the 
sequence is scrambled otherwise not. If \code{1}, Owen type type of scrambling 
is applied, if \code{2}, Faure-Tezuka type of scrambling, is applied, and if 
\code{3}, both Owen+Faure-Tezuka type of scrambling is applied. The \code{star} 
method uses the same initial values as \code{rmfmkl} since it uses the FMKL 
generalised lambda distribution. Nelder-Simplex algorithm is used in the 
numerical optimization. \code{rprs} stands for revised percentile method for 
RS generalised lambda distribution and "rmfmkl" stands for revised method of 
moment for FMKL generalised lambda distribution. These acronyms represents the 
initial optimization algorithm used to get a reasonable set of initial values 
for the subsequent optimization procedues.
}
\value{
\item{par}{ The best set of parameters found, the first four corresponds to the 
first distribution fit, the second four corresponds to the second distribution 
fit, the last value correspond to p for the first distribution fit. }
\item{value}{ The value of -ML for the paramters obtained. }
\item{counts}{ A two-element integer vector giving the number of calls to 
\code{fn} and \code{gr} respectively. This excludes those calls needed to 
compute the Hessian, if requested, and any calls to \code{fn} to compute a 
finite-difference approximation to the gradient. }
\item{convergence}{ \code{0} indicates successful convergence, \code{1} 
indicates the iteration limit \code{maxit} had been reached, \code{10} indicates 
degeneracy of the Nelder-Mead simplex. }
\item{message}{ A character string giving any additional information returned by 
the optimizer, or \code{NULL}. }
}
\references{ 
Bratley P. and Fox B.L. (1988) Algorithm 659: Implementing Sobol's quasi random 
sequence generator, ACM Transactions on Mathematical Software 14, 88-100.

Joe S. and Kuo F.Y. (1998) Remark on Algorithm 659: Implementing Sobol's quasi 
random Sequence Generator.

Nelder, J. A. and Mead, R. (1965) A simplex algorithm for function minimization. 
Computer Journal *7*, 308-313.

Su (2008). Fitting GLD to data via quantile matching method. (Book chapter to
appear)
}
\author{ Steve Su }
\note{ If the number of observations is small, 
\code{rprs} can sometimes fail as the percentiles may not exist for this data. 
Also, if the initial values do not span a valid generalised lambda distribution, 
try another set of initial values. }
\seealso{ \code{\link{fun.auto.bimodal.pml}}, \code{\link{fun.auto.bimodal.ml}},
\code{\link{fun.plot.fit.bm}}, 
\code{\link{fun.diag.ks.g.bimodal}} }
\examples{

\donttest{
# Fitting faithful data from the dataset library, with the clara clustering 
# regime. The first distribution is RS and the second distribution is fmkl. 
# The percentage of data mix is 1\%.

fun.auto.bimodal.qs(faithful[,1],per.of.mix=0.01,clustering.m=clara,
init1.sel="rprs",init2.sel="rmfmkl",init1=c(-1.5,1,5),init2=c(-0.25,1.5),
leap1=3,leap2=3)
}
}
\keyword{smooth}
