\name{fun.auto.bimodal.pml}
\alias{fun.auto.bimodal.pml}
\title{ Fitting mixture of generalied lambda distribtions to data using 
parition maximum likelihood estimation }
\description{
This function will fit mixture of generalised lambda distributions to dataset. 
It is restricted to two generalised lambda distributions. The method of fitting 
is parition maximum likelihood. It is a two step optimization procedure, each 
unimodal part of the bimodal distribution is modelled using the maximum 
likelihood method or the starship method (FMKL GLD only). These initial values 
the used to "maximise" the complete log likelihood for the entire bimodal 
distribution. It fits mixture of the form p*(f1)+(1-p)*(f2) where f1 and f2 are 
pdfs of the generalised lambda distributions.
}
\usage{
fun.auto.bimodal.pml(data, clustering.m = clara, init1.sel = "rprs", 
init2.sel = "rprs", init1=c(-1.5, 1.5), init2=c(-1.5, 1.5), leap1=3, leap2=3,
fun1="runif.sobol", fun2="runif.sobol",no=10000,max.it=5000, optim.further="Y")
}
\arguments{
  \item{data}{ A numerical vector representing the dataset. }
  \item{clustering.m}{ Clustering method used in classifying the dataset into 
  two parts. Valid arguments include clara, fanny and pam from the cluster 
  library. Default is clara. Or a logical vector specifying how data should
  be split. }
 \item{init1.sel}{ This can be \code{"rprs"}, \code{"rmfmkl"} or \code{"star"}, 
 the initial method used to fit the first distribution. }
  \item{init2.sel}{ This can be \code{"rprs"}, \code{"rmfmkl"} or \code{"star"}, 
  the initial method used to fit the second distribution. }
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
  \item{no}{ Number of initial random values to find the best initial values for 
  optimisation. }   
  \item{max.it}{ Maximum number of iterations for numerical optimisation. }
  \item{optim.further}{Whether to optimise the function further using full 
  maximum likelihood method, recommended setting is "Y"}
 
}
\details{
The initial values that work well for RPRS are \code{c(-1.5,1.5)} and for RMFMKL 
are \code{c(-0.25,1.5)}. For scrambling, if \code{1}, \code{2} or \code{3} the 
sequence is scrambled otherwise not. If \code{1}, Owen type type of scrambling 
is applied, if \code{2}, Faure-Tezuka type of scrambling, is applied, and if 
\code{3}, both Owen+Faure-Tezuka type of scrambling is applied. The \code{star} 
method uses the same initial values as \code{rmfmkl} since it uses the FMKL 
generalised lambda distribution. Nelder-Simplex algorithm is used in the 
numerical optimization. \code{rprs} stands for revised percentile method for RS 
generalised lambda distribution and "rmfmkl" stands for revised method of moment 
for FMKL generalised lambda distribution. These acronyms represents the initial 
optimization algorithm used to get a reasonable set of initial values for the 
subsequent optimization procedues. This function is an improvement from Su 
(2007) in Journal of Statistical Software.
}
\value{
\item{par}{ The best set of parameters found, the first four corresponds to the 
first distribution fit, the second four corresponds to the second distribution 
fit, the last value correspond to p for the first distribution fit. }
\item{value}{ The value of -PML for the paramters obtained. }
\item{counts}{ A two-element integer vector giving the number of calls to "fn" 
and "gr" respectively. This excludes those calls needed to compute the Hessian, 
if requested, and any calls to 'fn' to compute a finite-difference approximation 
to the gradient. }
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

Su (2007). Fitting Single and Mixture of Generalized Lambda Distributions to 
Data via Discretized and Maximum Likelihood Methods: GLDEX in R. Journal of 
Statistical Software: *21* 9.
}
\note{ If the number of observations is small, \code{rprs} can sometimes fail as 
the percentiles may not exist for this data. Also, if the initial values do not 
result in a valid generalised lambda distribution, try another set of initial 
values.
}
\seealso{ \code{\link{fun.auto.bimodal.ml}},\code{\link{fun.plot.fit.bm}},
\code{\link{fun.diag.ks.g.bimodal}} }
\examples{
\donttest{
# Fitting faithful data from the dataset library, with the clara clustering 
# regime. The first distribution is RS and the second distribution is fmkl. 

fun.auto.bimodal.pml(faithful[,1],clustering.m=clara,init1.sel="rprs",
init2.sel="rmfmkl",init1=c(-1.5,1,5),init2=c(-0.25,1.5),leap1=3,leap2=3)

}}
\keyword{smooth}
